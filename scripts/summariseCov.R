#Take the full output from VarScan mpileup2cns (data generated for every position in the 
# assembly) and summarise coverage by contig.

library(plyr)
library(reshape2)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

infile <- args[1]
outfile <- args[2]
                                                                              
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
                     lowCovPos = sum(Reads1+Reads2 < 10))

fullSummary <- cbind(colsplit(fullSummary$Chrom, '\\|size', c('Chrom','ContigLength')), fullSummary[,2:10])

write.table(fullSummary, file = outfile, quote=F, row.names = F, col.names = T, sep='\t')
