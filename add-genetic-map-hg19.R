## R code to add genetic map positions to PLINK format .bim file, using
## Oxford/HapMap genetic maps for GRC Build 37/hg19 coordinates
##
## Supply name of .bim file as argument using --args
## e.g. R --vanilla --file=add-genetic-map-hg19.R --args plink.bim
##
## On completion, input .bim file will be saved with .bim.old extension
## and input .bim file will be overwritten with genetic map positions in Morgans in 3rd column
##
## The input .bim file may contain markers on multiple autosomes,
## but non-autosomal markers are not handled
##
## Toby.x.Johnson@gsk.com 04 Feb 2013

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) bimfile <- args[1] else stop("Supply input filename (e.g. plink.bim) with --args")

mappath <- "/GWD/appbase/projects/statgen/RD-MDD-GX_PUBLIC/HapMap/recombination/2011-01_phaseII_B37"
## Set to path to directory containing hg19 genetic maps

bim <- read.table(bimfile, comment = "", header = FALSE, as.is = TRUE,
                  col.names = c("CHR", "SNP", "MAP", "BP", "A1", "A2"))
options(scipen = 15) ## force fixed notation for integer BP positions
write.table(bim, file = paste(bimfile, "old", sep = "."),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

chrs <- unique(bim$CHR)
stopifnot(all(chrs %in% 1:23))
for (chr in chrs) {
  if(chr < 23) {
    map <- read.table(gzfile(file.path(mappath, paste("genetic_map_GRCh37_chr", chr, ".txt.gz", sep = ""))),
                    comment = "", header = TRUE, as.is = TRUE)
    names(map) <- sapply(names(map), function(nn) return(unlist(strsplit(nn, ".", fixed = TRUE))[1]))
    stopifnot(all(c("Chromosome", "Position", "Map") %in% names(map)))
    ## These maps start at non-zero physical position but zero genetic
    ## map position.  Thus, markers physically left of the first marker
    ## in the map will be interpolated as having zero genetic
    ## map positions, which is problematic because (i) HAPI-UR infers
    ## the existence of a genetic map if the genetic map position of the
    ## *first* marker (only) is non-zero, and (ii) HAPI-UR requires
    ## unique map positions for all markers.
    ##
    ## Hence we 'nudge' the whole map right, based on average cM/Mb:
    map$Map <- map$Map + map$Map[nrow(map)]/map$Position[nrow(map)]*map$Position[1]
  }
  if(chr == 23) {
    map.x <- NULL
    x.max <- 0
    for(i in c("X_par1", "X", "X_par2")){
      map <- read.table(gzfile(file.path(mappath, paste("genetic_map_GRCh37_chr", i, ".txt.gz", sep = ""))),
                        comment = "", header = TRUE, as.is = TRUE)
      names(map) <- sapply(names(map), function(nn) return(unlist(strsplit(nn, ".", fixed = TRUE))[1]))
      stopifnot(all(c("Chromosome", "Position", "Map") %in% names(map)))
      ## These maps start at non-zero physical position but zero genetic
      ## map position.  Thus, markers physically left of the first marker
      ## in the map will be interpolated as having zero genetic
      ## map positions, which is problematic because (i) HAPI-UR infers
      ## the existence of a genetic map if the genetic map position of the
      ## *first* marker (only) is non-zero, and (ii) HAPI-UR requires
      ## unique map positions for all markers.
      ##
      ## Hence we 'nudge' the whole map right, based on average cM/Mb:
      map$Map <- map$Map + map$Map[nrow(map)]/map$Position[nrow(map)]*map$Position[1] + x.max
      x.max<- map$Map[nrow(map)]
      map.x <- rbind(map.x, map)
    }
    map.x$Chromosome <-"chr23"
    map<- map.x
  }
  
  ## Interpolation, with 'safe' extremes, left at (0,0) and right based on average cM/Mb
  bim$MAP[bim$CHR == chr] <- 0.01*approx(c(0, map$Position, map$Position[nrow(map)]*1.5),
                                         c(0, map$Map, map$Map[nrow(map)]*1.5),
                                         bim$BP[bim$CHR == chr])$y
  ## Note 0.01* to convert cM to Morgans
}

## Rewrite markers with missing map positions in PLINK compatible way
bim.na <- is.na(bim$MAP)
bim$CHR[bim.na] <- 0
bim$MAP[bim.na] <- 0
bim$BP[bim.na] <- 0

options(scipen = 15) ## force fixed notation for integer BP positions
write.table(bim, file = bimfile,
            row.names = FALSE, col.names = FALSE, quote = FALSE, na = "0")

