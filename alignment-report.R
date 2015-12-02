args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) path <- args[1] else stop("Supply path for alignment directory with --args")

rawfreqfile <- paste(path, "raw-freq.frq", sep = "/") # output of PLINK --freq before filtering
snplistfile <- paste(path, "snplist.snplist", sep = "/") # output of PLINK --write-snplist after filtering
liftbase <- paste(path, "chr", sep = "/") # add 1:22 to get imputation-alignment.py output...
freqfile <- paste(path, "alignment3-freq.frq", sep = "/") # output of PLINK --freq after alignment
reportmatchfile <- paste(path, "report-match.csv", sep = "/") # file to write report table
reportfreqfile <- paste(path, "report-freq.csv.gz", sep = "/") # file to write report table

allelesAB <- function(A1, A2, sep = "/") ifelse(A2 < A1, paste(A2, A1, sep = sep), paste(A1, A2, sep = sep))
freqB <- function(A1, A2, freq1) ifelse(A2 < A1, 1 - freq1, freq1)

rawfrq <- read.table(rawfreqfile, comment.char = "", header = TRUE, as.is = TRUE)
snplist <- scan(snplistfile, character(0))
##use nrow = 50 for efficiency, and assume all chr have the same headers
chrs <- which(paste("chr", 1:22, sep = "") %in% dir(path))
stopifnot(length(chrs) >= 1)
col.names <- sub("^#", "", unlist(strsplit(grep("^#SNP", read.table(paste(liftbase, chrs[1], sep = ""), comment = "", header = FALSE, sep="\n", as.is = TRUE, nrow = 50)[ , 1], value = TRUE), " ")))
lookup <- do.call("rbind", lapply(chrs, function(chr) read.table(paste(liftbase, chr, sep = ""), comment.char = "#", header = FALSE, col.names = col.names, as.is = TRUE)))
frq <- read.table(freqfile, comment.char = "", header = TRUE, as.is = TRUE)

rawfrq <- within(rawfrq, {passFilter <- SNP %in% snplist;
                    isAuto <- CHR %in% 1:22;
                    isACGT <- A1 %in% c("A", "C", "G", "T") & A2 %in% c("A", "C", "G", "T");
                    AB <- factor(paste("alleles", allelesAB(A1, A2)),
                                 paste("alleles", c("A/C", "A/G", "A/T", "C/G", "C/T", "G/T", "other")));
                    AB[is.na(AB)] <- "alleles other";
                    MAFbin <- factor(paste("MAF", ifelse(MAF == 0, "0",
                                                         ifelse(MAF <= 0.001, "(0, 0.001]",
                                                                ifelse(MAF <= 0.01, "(0.001, 0.01]",
                                                                       ifelse(MAF <= 0.05, "(0.01, 0.05]", "(0.05, 0.5]"))))),
                                     paste("MAF", c("0", "(0, 0.001]", "(0.001, 0.01]", "(0.01, 0.05]", "(0.05, 0.5]")));

                    lookupMatch <- match(SNP, lookup$SNP);
                    lookupStrand <- lookup$STRAND[lookupMatch]
                    lookupSign <- c(-1, +1)[match(lookupStrand, c("-", "+"))]
                    lookupResult <- c("allele-mismatch", "match-minus", "match-plus")[match(lookupStrand, c(".", "-", "+"))];
                    result <- factor(ifelse(!passFilter, "fail-filter",
                                            ifelse(!isAuto, "not-autosomal",
                                                   ifelse(!isACGT, "not-ACGT",
                                                          ifelse(is.na(lookupMatch), "no-match", lookupResult)))),
                         c("fail-filter", "not-autosomal", "not-ACGT", "no-match", "allele-mismatch", "match-minus", "match-plus"))})

resulttab <- with(rawfrq, rbind(table(result), table(MAFbin, result), table(AB, result)))
rownames(resulttab)[1] <- "all"

write.csv(resulttab, file = reportmatchfile,
          row.names = TRUE)

frq$allelesAB <- with(frq, allelesAB(A1, A2))
frq$freqB <- with(frq, freqB(A1, A2, MAF))
lookup$freqB <- with(lookup, freqB(ALT, REF, AC/AN))
lookup$SNP2 <- with(lookup, paste(CHROM, POS, sep = ":"))
lookup.match <- match(frq$SNP, lookup$SNP2)
frq$freqBref <- lookup$freqB[lookup.match] # should be same order
frq$OrigSNP <- lookup$SNP[lookup.match] # SNP name in raw data
  
write.csv(subset(frq, !is.na(freqB) & !is.na(freqBref), select = c("SNP", "OrigSNP", "allelesAB", "freqB", "freqBref")),
          file = gzfile(reportfreqfile),
          row.names = FALSE)


##pdftemp <- tempfile(fileext = ".pdf")
##pdf(pdftemp, width = 6, height = 6)
##with(frq, {plot(freqBref, freqB, type = "n",
##                xlim = c(0, 1), ylim = c(0, 1),
##                xlab = "Reference B allele frequency", ylab = "Target B allele frequency", las = 1);
##           text(freqBref, freqB, allelesAB(A1, A2), cex = 0.3)})
##dev.off()
#### reportplot does not exist
##system(paste("/usr/central/bin/convert -density 144x144", pdftemp, reportplot))
##system(paste("rm", pdftemp))


## post-imputation dataset will not be the only analysis dataset because target variants not present in reference are DROPPED
## hence exome content found in ESP but not 1KG will not make it...

#png(file = reportplot)
#with(frq, {plot(freqBref, freqB, type = "n",
#                xlim = c(0, 1), ylim = c(0, 1),
#                xlab = "Reference B allele frequency", ylab = "Target B allele frequency", las = 1);
#           text(freqBref, freqB, allelesAB(A1, A2), cex = 0.3)})
#dev.off()
#q()


