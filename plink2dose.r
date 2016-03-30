################################################################
#running command:R64-2.14.0 --vanilla --slave --args InfilePrefix OutputDir< path/plink2Dose.r 
#InfilePrefix: prefix(including full path) for plink *.raw and *.bim files
#OutputDir: create chr[i]chunkG.dose.gz, chr[i]chunkG.info.gz and chr[i]chunkG.done files under this directory
################################################################
options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {  
  prefix <- args[1]
  dir.out <- args[2]
  
} else {
  stop("Please input [infile prefix] [output directory]")
}

#prefix <- "I:/statgen/GXapp/GWASImputationPipelineUpdates/ChrX/HMaxiom_1KGP3/genotypes/no-imputation"
#dir.out <- "I:/statgen/GXapp/GWASImputationPipelineUpdates/ChrX/HMaxiom_1KGP3/imputed-20130502"
#setwd(dir.out)

#creates output dir if not exists
dir.create(dir.out)

#bim file: chr id cm bp a1 a2
#SNP = chr:pos
data.bim <- read.table(paste(prefix, "bim", sep = "."), sep = "\t",
                       comment = "", quote = "", header = F, as.is = TRUE)
names(data.bim) <- c("CHROM", "SNP", "CM", "POS", "Al1",  "Al2")

#raw file: FID IID PAT MAT SEX PHENOTYPE SNP_A1
data.raw <- read.table(paste(prefix, "raw", sep = "."), check.names = FALSE,
                       comment = "", quote = "", header = T, as.is = TRUE)
markers <- names(data.raw)[-c(1:6)] #chr:pos_A1
cols.info <- c("SNP", "Al1", "Al2", "Freq1", "MAF", "AvgCall", "Rsq", "Genotyped", "LooRsq", "EmpR", "EmpRsq", "Dose1", "Dose2")
data.info <-data.frame(MARKER = markers, stringsAsFactors =FALSE)
data.info <- cbind(data.info, read.table(text=markers, sep = "_", as.is = T, header = F, fill = T, col.names = c("SNP", "Al1")))
data.info<- cbind(data.info,  read.table(text=data.info$SNP, sep = ":",  as.is = T, header = F, fill = T, col.names = c("CHROM", "POS")))
data.info$Al2<- data.bim$Al2[match(data.info$SNP, data.bim$SNP)]
if(any(data.info$Al1==data.info$Al2)) 
  data.info$Al2[data.info$Al1==data.info$Al2]<- data.bim$Al1[match(data.info$SNP[data.info$Al1==data.info$Al2], data.bim$SNP)]
data.info$Freq1 <- unlist(lapply(data.raw[markers], mean, na.rm = T))/2
data.info$MAF<-data.info$Freq1
data.info$AvgCall<- unlist(lapply(data.raw[markers], function(x)length(which(!is.na(x)))))/nrow(data.raw)
data.info$Rsq<-1
data.info$Genotyped<-"Genotyped"
for(c in c("LooRsq",	"EmpR",	"EmpRsq",	"Dose1",	"Dose2")) {
  data.info[[c]]<-"-"
}
data.info<- data.info[order(data.info$CHROM, data.info$POS), ]



chrs<- sort(unique(data.info$CHROM))
for(chr in chrs){
  info<- data.info[data.info$CHROM == chr,]
  dose <- cbind(data.frame(SEX=data.raw$SEX,
                           ID = paste(data.raw$FID,"->", data.raw$IID, sep ="" ),
                           Dose = rep("DOSE", nrow(data.raw))), 
                data.frame(lapply(data.raw[match(info$MARKER, names(data.raw))], #fill missing data with mean
                                  function(x) {x[is.na(x)]<-mean(x, na.rm = T); return (x)})))
  names(dose)[-c(1:3)] <- info$MARKER
  if(chr != 23)  {
    out.prefix <- paste(dir.out, "/chr", chr, "chunkG", sep ="")
    write.table(info[cols.info], file = gzfile(paste(out.prefix, ".info.gz", sep = "")), 
                quote = F, row.names = F, sep = "\t")
    write.table(dose[-1], file = gzfile(paste(out.prefix, ".dose.gz", sep = "")),
                quote = F, row.names = F, col.names = F, sep = "\t")
    write("",file =paste(out.prefix, ".done", sep = "")) 
                     
  }else{
    #par: POS < 2699520 or POS > 154931043
    info.par <- info[info$POS <2699520 | info$POS >154931043, ]
    if(nrow(info.par)>0) {
      out.prefix <- paste(dir.out, "/chr", chr, "-parchunkG", sep ="")
      dose.par <- dose[c(2,3, match(info.par$MARKER, names(dose)))] 
      write.table(info.par[cols.info], file = gzfile(paste(out.prefix, ".info.gz", sep = "")), 
                  quote = F, row.names = F, sep = "\t")
      write.table(dose.par, file = gzfile(paste(out.prefix, ".dose.gz", sep = "")),
                  quote = F, row.names = F, col.names = F, sep = "\t")
      write("",file =paste(out.prefix, ".done", sep = ""))
    }
    
    #non-par: 2699520-154931043 seperate by gender
    info.nonpar <- info[info$POS >=2699520 | info$POS <= 154931043, ]
    if(nrow(info.nonpar)>0) {
      #1=male; 2=female
      for(i in 1:2) {
        out.prefix <- paste(dir.out, "/chr", chr,"-",c("M", "F")[i], "chunkG", sep ="")
        dose.nonpar <- dose[!is.na(dose$SEX)&dose$SEX==i, c(2,3, match(info.nonpar$MARKER, names(dose)))] 
        write.table(info.nonpar[cols.info], file = gzfile(paste(out.prefix, ".info.gz", sep = "")), 
                    quote = F, row.names = F, sep = "\t")
        write.table(dose.nonpar, file = gzfile(paste(out.prefix, ".dose.gz", sep = "")),
                    quote = F, row.names = F, col.names = F, sep = "\t")
        write("",file =paste(out.prefix, ".done", sep = ""))
      }
    }
  }
}










