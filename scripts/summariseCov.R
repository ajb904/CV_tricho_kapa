#Take the full output from VarScan mpileup2cns (data generated for every position in the 
# assembly) and summarise coverage by contig.

library(plyr)
library(reshape2)
library(data.table)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)

infile <- args[1]
outfile <- args[2]
fastafile <- args[3]
                                                                              
cols <- c(1,2,5,6)
full.data <- fread(infile, header = T, sep = '\t', select = cols, data.table=F)

fullSummary <- ddply(full.data, 'Chrom', summarise,
                     medCov = median(Reads1+Reads2),
                     meanCov = mean(Reads1+Reads2),
                     sdCov = sd(Reads1+Reads2),
                     lqCov = quantile(Reads1+Reads2, 0.25), 
                     uqCov = quantile(Reads1+Reads2, 0.75),
                     minCov = min(Reads1+Reads2),
                     maxCov = max(Reads1+Reads2),
                     highCovPos = sum(Reads1+Reads2 > 250),
                     acceptCovPos = sum(Reads1+Reads2 >= 10 & Reads1+Reads2 <= 250))
# Positions with coverage < 8 are not included in VarScan file. Therefore, can't calculate low coverage
# positions just by adding those <10. Instead, get number of acceptably covered and high coverage positions
# and subtract from total length.

#The contig lengths given in the sequence IDs are not correct in all cases. Therefore, we need to calculate
#sequence lengths directly from the actual sequences

seqs <- readDNAStringSet(fastafile)
seqs.df <- data.frame(Chrom = names(seqs), ContigLength=width(seqs))

fullSummary <- merge(fullSummary, seqs.df, by='Chrom')
fullSummary[,'lowCovPos'] <- fullSummary$ContigLength - (fullSummary$highCovPos + fullSummary$acceptCovPos)

write.table(fullSummary, file = outfile, quote=F, row.names = F, col.names = T, sep='\t')
