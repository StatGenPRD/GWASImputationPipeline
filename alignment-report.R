args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) path <- args[1] else stop("Supply path for alignment directory with --args")

# v3, assumes gzipped frequency and liftover files

cairo.okay <- getOption("bitmapType") == "cairo"

rawfreqfile <- file.path(path, "raw-freq.frq.gz") # output of PLINK --freq before filtering
snplistfile <- file.path(path, "snplist.snplist") # output of PLINK --write-snplist after filtering
liftbase <- file.path(path, "chr") # add 1:22 to get imputation-alignment.py output...
freqfile <- file.path(path, "alignment3-freq.frq.gz") # output of PLINK --freq after alignment
reportmatchfile <- file.path(path, "report-match.csv") # file to write report table
reportfreqfile <- file.path(path, "report-freq.csv.gz") # file to write report table
reportfreqpng <- file.path(path, "report-freq.png") # file to make png plot

allelesAB <- function(A1, A2, sep = "/") ifelse(A2 < A1, paste(A2, A1, sep = sep), paste(A1, A2, sep = sep))
freqB <- function(A1, A2, freq1) ifelse(A2 < A1, 1 - freq1, freq1)

rawfrq <- read.table(gzfile(rawfreqfile), comment.char = "", header = TRUE, as.is = TRUE)
snplist <- scan(snplistfile, character(0))
##use nrow = 50 for efficiency, and assume all chr have the same headers
chrs <- c(1:22, "X")[which(paste("chr", c(1:22, "X"), ".gz", sep = "") %in% dir(path))]
stopifnot(length(chrs) >= 1)
col.names <- sub("^#", "", unlist(strsplit(grep("^#SNP", read.table(gzfile(paste(liftbase, chrs[1], ".gz", sep = "")), comment = "", header = FALSE, sep="\n", as.is = TRUE, nrow = 50)[ , 1], value = TRUE), " ")))
##15Feb2013 changed to workaround for markers HumanOmni1-Quad_v1/hg18_b36_B with hash (#) symbol in names
#lookup <- do.call("rbind", lapply(chrs, function(chr) read.table(paste(liftbase, chr, sep = ""), comment.char = "#", header = FALSE, col.names = col.names, as.is = TRUE)))
lookup <- do.call("rbind", lapply(chrs, function(chr) read.table(gzfile(paste(liftbase, chr, ".gz", sep = "")), skip = 1, comment.char = "", header = FALSE, col.names = col.names, as.is = TRUE)))
frq <- read.table(gzfile(freqfile), comment.char = "", header = TRUE, as.is = TRUE)

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

if (cairo.okay) {
  res <- 300
  png(file = reportfreqpng, width = 6*res, height = 6*res, res = res)
} else {
  pdftemp <- tempfile(fileext = ".pdf")
  pdf(file = pdftemp, width = 6, height = 6)
}
par(mar = c(4, 4, 0, 0) + 0.1)
with(subset(frq, !is.na(freqB) & !is.na(freqBref)), {
  plot(freqBref, freqB, type = "n", xlim = c(0, 1), ylim = c(0, 1))
  text(freqBref, freqB, allelesAB, cex = 0.5)
})
dev.off()
if (cairo.okay) {
  # do nothing
} else {
  system(paste("/usr/central/bin/convert -density 150x150", pdftemp,
               "-crop 900x900+0+750", reportfreqpng))
  unlink(pdftemp)
}
