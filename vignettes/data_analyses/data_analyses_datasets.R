#!/usr/bin/env R

# Handle large data objects for the data_analyses vignette

set.seed(2020)

# h5 database file
fn <- "remethdb2.h5"

# metadata
path <- system.file("extdata", "metadata", package = "recountmethylation")
md <- get(load(list.files(path)[1]))

# get cg subset from manifest
ncg <- 5000
man.package <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
man <- get(data(Manifest, package = man.package))
cgrn <- rownames(man)
cgselect <- sample(cgrn, ncg)
msub <- man[cgselect,]
cgaddr <- unique(c(msub$AddressA, msub$AddressB))

#------------------------------------------
# analysis 1. age comparison, blood samples
#------------------------------------------
mdf <- md[grepl("(^|;)blood(;|$)", md$tissue),]
mdf$sgroup <- ifelse(as.numeric(mdf$predage) >= 50, "old", "young")
gsmv <- c(sample(mdf[mdf$sgroup == "young",]$gsm, 100),
          sample(mdf[mdf$sgroup == "old",]$gsm, 100))
mdf <- mdf[mdf$gsm %in% gsmv,]
rg1 <- getrg(gsmv = gsmv, cgv = cgaddr, all.cg = FALSE)
save(rg1, file = "rgraw_analysis1.rda")

# validat a couple of samples

#---------------------------------------
# analysis 2. brain versus blood samples
#---------------------------------------
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
ngsm.sg <- 100
gsm.blood <- which(mdf$sgroup == labs[1])
gsm.brain <- which(mdf$sgroup == labs[2])
gsm.indices <- c(sample(gsm.blood, ngsm.sg), sample(gsm.brain, ngsm.sg))
mdf <- mdf[gsm.indices,]
rg2 <- getrg(gsmv = mdf$gsm)
save(rg2, file = "rgraw_analysis2.rda")