# 
# Code to scan AHR ChIP peak fasta seqs, search for *all* potential DREs, and list them
# 
#Author  Kalyani Dhusia 3/22/2019
#email at kalyanidhusia.bhu@gmail.com
 
library(seqinr)
library(Biostrings)
library(TFBSTools)
library(seqLogo)
library(biomaRt)
library(JASPAR2018)
library(reshape2)
library(ggplot2)
library(nucleR)

# Get the Nrf2 PFM from JASPAR 2018 and convert to PWM
Nrf2_pfm <- getMatrixByName(JASPAR2018, "NFE2L2")
Nrf2_pwm <- toPWM(Nrf2_pfm)

Nrf2_pfm_mod <- Nrf2_pwm

Matrix(Nrf2_pfm_mod)

# Set threshold MS score for DREs (value from Dere et al, CRT 2011)
#MS_score_threshold = 0.85


##-----------
# Create an empty data frame to store long (19-bp) site, location and MS score
# for all putative DREs in sequence
longSites <- data.frame(seq = character(), 
                        chrom = character(),
                        start_coord = integer(),
                        end_coord = integer(),
                        stringsAsFactors=FALSE)


##----------
# Read in NRF2 peak seqs fasta file
peak_seqs <- read.fasta("nrf2_peaks.fa")
cat("--- No. of NRF2 peaks found = ", length(peak_seqs), "\n")


##----------
# Loop over all peak sequences read in
for (i in 1:length(peak_seqs)) {
  # Read in whole sequence
  mySeq <- peak_seqs[[i]]  
  
  peak_annot <- attr(mySeq, "Annot")  # Get the annotation of the peak
  peak_chrom1 <- strsplit(peak_annot, ":")[[1]][1]  # Get chromosome # of peak
  peak_chrom2 <- strsplit(peak_chrom1, ">")[[1]][2]
  
  peak_range <- strsplit(peak_annot, ":")[[1]][2]  # Get chromosomal range of peak
  peak_start <- as.numeric(strsplit(peak_range, "-")[[1]][1])  # Get start coord of peak
  
  
  # Convert sequence to DNA string
  mySeqDNAString <- DNAString(c2s(mySeq))
  
  # Search sequence for exct match to *modified* AhR PWM     
  siteSet <- searchSeq(Nrf2_pfm_mod, mySeqDNAString, 
                       seqname = "seq1", min.score = "70%", strand = "*")
  sites_in_myseq <- writeGFF3(siteSet)
  
  # If binding sites found, add all instances to dataframe of sites
  numSites_in_myseq <- dim(sites_in_myseq)[1]
  cat("* Found ", numSites_in_myseq, " putative DREs in peak # ", i, " \n")
  
  if (numSites_in_myseq > 0) {
    # initialize list of DRE start and end locations
    starts = numeric()
    ends = numeric()
    
    # Loop over DREs
    for (j in 1:numSites_in_myseq) {
      # Find start and end coords of site
      start <- sites_in_myseq$start[j]
      end <- sites_in_myseq$end[j]
    }
  
      
      # Add long_site info and MS score to longSites dataframe
      #cat("DRE no. = ", j, "\n")
      long_site_string <- toString(long_site)
      start_long_site = (start - 0)
      end_long_site = (end + 0) 
      
      # Add to list of starts and ends for current peak
      starts <- c(starts, start_long_site)
      ends <- c(ends, end_long_site)
      
      # Create temporary data frame of long sites
      temp_df <- data.frame(seq = long_site_string, 
                            chrom = peak_chrom2,
                            start_coord = start_long_site + peak_start, 
                            end_coord = end_long_site + peak_start,
                            peak_assoc_with = i)
      
      
      # Append to list of long sites
      longSites <- rbind(longSites, temp_df) 
      
    }  # end of loop over DREs
      

}  # end loop over peak seqs


##----------
# Save long-site seqs, scores and peak nos. to CSV file
write.csv(longSites2, file = "DHS-minus-NRF2_peaks_19bp.csv", row.names = FALSE)

# Savepeak info. to CSV file
write.csv(peaks_info, file = "NRF2peaks_features.csv", row.names = FALSE)

##----------
# Create bed file from long-site coords

# Create new data frame from longSites with only necessary fields
longSites2 <- data.frame(chr = longSites$chrom,
                         start = longSites$start_coord,
                         end = longSites$end_coord)

#head(longSites2)

# Write to bed file
f <- file("DHS-minus-NRF2_peaks_19bp.bed", "w") 
line = 'track name="unbound NRF2s"  description="NRF2s not under NRF2 peaks"  useScore=1'
writeLines(line, con = f)
close(f)
write.table(longSites2, file="DHS-minus-NRF2_peaks_19bp.bed",
            quote=F, sep="\t", row.names=F, col.names=F, append = T)