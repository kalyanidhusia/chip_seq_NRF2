#
# Code to group AHR ChIP peaks by "mode"
# 
#Author Kalyani Dhusia 03/25/2019
#email at kalyanidhusia.bhu@gmail.com


# System command to identify AHR and write out peaks not overlapping with DHSs
system("bedtools intersect -a NRF2.bed -b ENCFF535VLW.bed -v > NRF2peaks_not-in-DHSs.bed")

# Read in AHR peaks info
peaks_info <- read.csv("NRF2peaks_features.csv", header = TRUE)

# Read in AHR peaks bed file
peaks_bed <- read.table("NRF2.bed", sep="\t")

# Filter bed file by # of NRF2 peaks
peaks_noNRF2_bed <- peaks_bed[peaks_info$num_dres == 0, ]
peaks_oneNRF2_bed <- peaks_bed[peaks_info$num_dres == 1, ]
peaks_withNRF2_bed <- peaks_bed[peaks_info$num_dres > 0, ]

# Write different peak groups to beds file
write.table(peaks_noNRF2_bed, file="NRF2peaks_noNRF2s.bed",
            quote=F, sep="\t", row.names=F, col.names=F, append = T)

write.table(peaks_oneNRF2_bed, file="NRF2peaks_oneNRF2.bed",
            quote=F, sep="\t", row.names=F, col.names=F, append = T)

write.table(peaks_withNRF2_bed, file="NRF2peaks_withNRF2s.bed",
            quote=F, sep="\t", row.names=F, col.names=F, append = T)
