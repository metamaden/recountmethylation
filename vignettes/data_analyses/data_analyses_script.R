#!/usr/bin/env R

# Script to generate env objects used in data_analyses.Rmd

savepath <- "data_analyses/env.RData"

suppressMessages(library(rhdf5))
suppressMessages(library(minfi))
suppressMessages(library(recountmethylation))
suppressMessages(library(knitr))

path <- "remethdb2.h5"
md <- data_mdpost(path)

# specify filter values
cn <- "tissue" # colname to query
labs <- c("blood", "brain") # labels of interest
regex.patt <- paste0("(^|;)", labs,"(;|$)") # pattern for strict label matches
# define and filter on sgroup variable
mdf$sgroup <- "NA"
sg.blood <- grepl(regex.patt[1], mdf[,cn]) & 
  !grepl(regex.patt[2], mdf[,cn])
sg.brain <- grepl(regex.patt[2], mdf[,cn]) & 
  !grepl(regex.patt[1], mdf[,cn])
mdf[sg.blood,]$sgroup <- labs[1]
mdf[sg.brain,]$sgroup <- labs[2]
mdf <- mdf[sg.blood|sg.brain,]
# select 50 samples per groups
set.seed(2)
ngsm.sg <- 50
gsm.blood <- which(mdf$sgroup == labs[1])
gsm.brain <- which(mdf$sgroup == labs[2])
gsm.indices <- c(sample(gsm.blood, ngsm.sg), sample(gsm.brain, ngsm.sg))
mdf <- mdf[gsm.indices,]

get_sst <- function(labs = c("blood", "brain")){
  # get summary stats by label
  sstat <- lapply(labs, function(x){
    mv <- mdf[mdf$sgroup==x,]
    v1 <- nrow(mv)
    numgsm.gse <- as.data.frame(table(mv$gseid))
    v2 <- round(mean(numgsm.gse[,2]), 3)
    v3 <- round(sd(numgsm.gse[,2]), 3)
    v4 <- nrow(numgsm.gse)
    # predage summaries
    pv <- as.numeric(mv$predage)
    v5 <- min(pv)
    v6 = max(pv)
    v7 <- round(mean(pv, na.rm = TRUE), 3)
    v8 <- round(sd(pv, na.rm = TRUE), 3)
    v9 <- length(pv[is.na(pv)])
    # predsex summaries
    v10 <- round(100*nrow(mv[mv$predsex=="F",])/nrow(mv), 3) # perc. F
    v11 <- nrow(mv[is.na(mv$predsex),])
    return(data.frame("ngsm" = v1, "meangsm.gse" = v2, "sdgsm.gse" = v3,  "numgse" = v4,
                      "min.predage" = v5, "max.predage" = v6, "mean.predage" = v7, 
                      "sd.predage" = v8, "numna.predage" = v9, "percfemale.predsex" = v10, 
                      "numna.predsex" = v11))
  })
  sst <- sapply(sstat, rbind) # bind label groups
  colnames(sst) <- labs
  return(sst)
}
sst1 <- get_sst()
save.image(file = savepath)