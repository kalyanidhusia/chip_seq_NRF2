# Code to scan NRF2 ChIP peak fasta seqs, search for *all* potential negative peaks, and list them
# 
#Author   Kalyani Dhusia 2/22/2019

 
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
longSites2 <- data.frame(seq = character(), 
                        chrom = character(),
                        start_coord = integer(),
                        end_coord = integer(),
                        stringsAsFactors=FALSE)

peaks_info <- data.frame(id = integer(),
                         num_dres = integer(),
                         density_dres = numeric(),
                         stringsAsFactors = FALSE)

##----------
# Get all false negative peak read counts from NRf2 peak bed file
system("bedtools subtract -a ENCFF126HBJ_nrf2.bed -b ENCFF767JPZ_dnase.bed > nrf2_negative.bed")


# Read in NRF2 peak seqs fasta file
peak_seqs <- read.fasta("nrf2_negative_sudin_peaks.fa")
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
  
  
  if (numSites_in_myseq > 0) {
    # initialize list of DRE start and end locations
    starts = numeric()
    ends = numeric()
    
    # Loop over DREs
    for (j in 1:numSites_in_myseq) {
      # Find start and end coords of site
      start <- sites_in_myseq$start[j]
      end <- sites_in_myseq$end[j]
      
      # Ignore sites too close to end of sequence string
      if ((start <= 1) || (end >= length(mySeqDNAString))) {
        next
      }
      
      if (sites_in_myseq$strand[j] == "+") {
        long_site <- mySeqDNAString[(start):(end)]
      } else {
        long_site <- reverseComplement(mySeqDNAString[(start):(end)])
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
      longSites2 <- rbind(longSites2, temp_df) 
      
    }  # end of loop over DREs
    
    # Calculate spread of DREs under peak
    spread_dres = max(ends) - min(starts)
    
    # Create temporary data frame of AHR peaks
    if (spread_dres != 0) {
      temp_df_peaks <- data.frame(id = i,
                                  num_dres = numSites_in_myseq,
                                  density_dres = numSites_in_myseq / spread_dres)
    } else {
      cat("**  Spread of DREs = 0!! \n")
    }
   
    
  } else {
    # No. of DREs in peak = 0
    temp_df_peaks <- data.frame(id = i,
                                num_dres = 0,
                                density_dres = 0)
  }
  
  # Append to list of peaks
  peaks_info <- rbind(peaks_info, temp_df_peaks) 
  
  
}  # end loop over peak seqs


##----------
# Save long-site seqs, scores and peak nos. to CSV file
write.csv(longSites2, file = "NRF2peak_including7_19bp_withPeakNum.csv", row.names = FALSE)

# Savepeak info. to CSV file
write.csv(peaks_info, file = "NRF2peaks_including7_features.csv", row.names = FALSE)

# Create histogram of no. of DREs among AHR peaks
A <- hist(peaks_info$num_dres, 
     main = "Distribution of NRF2 peaks",
     xlab = "Total No. of NRF2 binding sites", ylab = "No. of NRF2 peaks",
     col = "5.9",
     border = "blue")


##----------
# Create histogram of no. of DREs among NRF2 peaks with ggplot2
ggp <- ggplot(data = peaks_info, aes(x = num_dres)) + 
  geom_histogram(color="blue", fill="5.9", binwidth = 1) + 
  ggtitle("Distribution of NRF2 ChIP peaks") + 
  labs(x = "Total Number of NRF2 binding sites", y = "Number of NRF2 peak") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 14), axis.text=element_text(size = 12))
ggsave("ggp.png")

##----------
# System calls to Unix
# Remove header lines from NRf2 peaks bed file
system("tail -n +1 NRF2.bed > NRF2peaks_FDR1_noHeader.bed")

# Get peak read counts from NRf2 peak regions
system("bigWigAverageOverBed ENCFF973LBO.bigwig NRF2peaks_FDR1_noHeader.bed NRF2peaks_FDR1_readSums.tab")

peak_sums_df <- read.table("NRF2peaks_FDR1_readSums.tab", sep="\t")
colnames(peak_sums_df)[4] <- "read_sums"
colnames(peak_sums_df)[5] <- "read_avg"

# Add peak read sums and averages to peaks_info
peaks_info$peak_sum <- peak_sums_df$read_sums
peaks_info$peak_avg <- peak_sums_df$read_avg


##----------
# Plot scatter plot of peak sums

# Plot scatter plot of peak sums vs no. DREs
ggp2 <- ggplot(data = peaks_info, aes(x = num_dres, y = peak_sum)) + 
  geom_point(alpha = 1/3) + 
  labs(x = "Number of NRF2 peak", y = "Area under NRF2 binding sites") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggsave("ggp2.png")

# Plot scatter plot of peak sums vs density of DREs
ggp3 <- ggplot(data = peaks_info, aes(x = density_dres, y = peak_sum)) + 
  geom_point(alpha = 1/3) + 
  labs(x = "Density of NRF2 peak", y = "Area under NRF2 bidning site") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggsave("ggp3.png")

# --Filter out peaks with no AHREs before plotting histograms etc.
peaks_info_filt <- peaks_info[peaks_info$num_dres > 0, ]

# Plot scatter plot of peak sums vs no. DREs *after filtering out peaks with no DREs*
ggp4 <- ggplot(data = peaks_info_filt, aes(x = num_dres, y = peak_sum)) + 
  geom_point(colour = "darkblue", alpha = 1/3) + scale_y_log10() +
  labs(x = "Number of NRF2 peak", y = "Area under NRF2 bidning site") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggsave("ggp4.png")

# Plot scatter plot of peak sums vs density of DREs *after filtering out peaks with no DREs*
ggp5 <- ggplot(data = peaks_info_filt, aes(x = density_dres, y = peak_sum)) + 
  geom_point(colour = "green4", alpha = 1/3) + scale_y_log10(limits = c(10, 3.5e5)) +
  labs(x = "Density of NRF2 peak", y = "Area under NRF2 bidning site") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggsave("ggp5.png")

#--
# Plot *boxplot* of NRF2 peak areas grouped by no. DREs per peak
ggp4a <- ggplot(data = peaks_info_filt, aes(x = as.factor(num_dres), y = peak_sum)) + 
  geom_boxplot(fill='lightblue', color="black") + scale_y_log10() +
  labs(x = "Number of NRF2 peak", y = "Area under NRF2 bidning site") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggsave("ggp4a.png")

##----------
# Plot scatter plot of peak averages

# Plot scatter plot of peak averages vs no. DREs
ggp6 <- ggplot(data = peaks_info, aes(x = num_dres, y = peak_avg)) + 
  geom_point(alpha = 1/3) + 
  labs(x = "Number of average NRF2 peak", y = "Averaged Area under NRF2 bidning site") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggsave("ggp6.png")

# Plot scatter plot of peak averages vs no. DREs *after filtering out peaks with no DREs*
ggp7 <- ggplot(data = peaks_info_filt, aes(x = num_dres, y = peak_avg)) + 
  geom_point(colour = "red4", alpha = 1/3) + scale_y_log10() +
  labs(x = "Number of NRF2 peak", y = "Averaged Area under NRF2 bidning site") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggsave("ggp7.png")

# Plot scatter plot of peak averages vs density of DREs *after filtering out peaks with no DREs*
ggp8 <- ggplot(data = peaks_info_filt, aes(x = density_dres, y = peak_avg)) + 
  geom_point(colour = "6", alpha = 1/3) + scale_y_log10() +
  labs(x = "Density of NRF2 peak", y = "Averaged area under NRF2 binding site") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggsave("ggp8.png")

#--- Try this:
# Plot scatter plot of peak averages vs peak areas *after filtering out peaks with no DREs*
ggp9 <- ggplot(data = peaks_info_filt, aes(x = peak_sum, y = peak_avg)) + 
  geom_point(colour = "green4", alpha = 1/3) + scale_x_log10() +  scale_y_log10() +
  labs(x = "Area under NRF2 peak", y = "Averaged area under NRF2 binding site") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggsave("ggp9.png")


#------
# Get singleton AHREs (AHREs under peaks with only one AHRE)

# First, get peaks with only one AHRE
peaks_withSingleAHRE <- peaks_info[peaks_info$num_dres == 1, ]

# Next, get the AHREs under those peaks
singletonAHREs <- longSites2[longSites2$peak_assoc_with %in% peaks_withSingleAHRE$id, ]

# -- NOTE: 837 singleton_AHREs, but 844 peaks_withSingleAHRE => some AHREs shared by peaks?

#--- See if there is any correlation between AHRE MS_score and AHR peak area?

# Filter peaks_withSingleAHRE such that one peak assocaited with each singleton_AHRE
peaks_withSingleAHRE_filt <- peaks_withSingleAHRE[peaks_withSingleAHRE$id %in% 
                                                    singletonAHREs$peak_assoc_with, ]

# Add peak_sum column to singleton_AHRE
singletonAHREs$peak_sum <- peaks_withSingleAHRE_filt$peak_sum

# Draw scatter-plot of AHR peak area vs. AHRE MS score
# Plot scatter plot of peak sums vs no. DREs *after filtering out peaks with no DREs*
ggp10 <- ggplot(data = singletonAHREs, aes(x = peak_assoc_with, y = peak_sum)) + 
  geom_point(colour = "darkblue", alpha = 1/3) + scale_y_log10() +
  labs(x = "singelton NRF2 binding sites", y = "Area under NRF2 peaks") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggsave("ggp10.png")

# --> Little correlation between singleton_AHRE MS_score and AHR peak area

# Save singletonNRF2 dataframe to csv file
write.csv(singletonAHREs, file = "singletonNRF2_info.csv", row.names = FALSE)
