#!/usr/bin/env R

library(recountmethylation)

# path to h5 db
dbn = "remethdb_test.h5"
h5ls(dbn) # summarize db objects

# get full gsm metadata
mdp = data.mdpost(dbn.path = dbn)
dim(mdp)

# check which samples available in h5 db
dsn = "redsignal"
rs.gsm = rhdf5::h5read(dbn, paste(dsn, "rownames", sep = ".")) # redsignal rownames (gsm ids)
rs.gsm = gsub("\\..*", "", rs.gsm)
mdf = mdp[mdp$gsm %in% rs.gsm,]
unique(unlist(strsplit(mdf$disease, ";"))) # available disease terms
unique(unlist(strsplit(mdf$tissue, ";"))) # available tissue terms

# get sample id query by tissue term
termi = "blood"
var.query = "tissue"
which.index = 1:5
which.gsm = which(grepl(paste0(".*", termi, ".*"), 
                        mdf[,var.query]))[which.index]
gsmvi = mdf$gsm[which.gsm]

# get cg address query by annotation filt
# library(anno.name)
#data(Manifest); man = Manifest
#data(Locations); loc = Locations

anno.name = "IlluminaHumanMethylation450kanno.ilmn12.hg19"
man = eval(parse(text = paste(anno.name, "Manifest", sep = "::")))
loc = eval(parse(text = paste(anno.name, "Locations", sep = "::")))
identical(rownames(loc), rownames(man))

chrname = "chr9"
cgfilt = rownames(loc[grepl(paste0("^", chrname, "$"), loc$chr),])
cgvi = unique(c(man[cgfilt,]$AddressA, man[cgfilt,]$AddressB))

# make df list from h5 db query
ldat.rgi = getrg(gsmv = gsmvi, dbn = dbn, dat.type = "df")
se.rgi = getrg(gsmv = gsmvi, dbn = dbn, dat.type = "se")

# basic preprocessing
lpre = list(se = list(),
            mv = list(), umv = list(),
            bv = list())
lpre[["se"]][["raw"]] = se.raw = preprocessRaw(se.rgi)
lpre[["se"]][["nb"]] = se.noob = preprocessNoob(se.rgi)

lpre[["mv"]][["raw"]] = getMeth(se.raw)
lpre[["mv"]][["nb"]] = getMeth(se.noob)

lpre[["umv"]][["raw"]] = getUnmeth(se.raw)
lpre[["umv"]][["nb"]] = getUnmeth(se.noob)

lpre[["bv"]][["raw"]] = getBeta(se.raw)
lpre[["bv"]][["nb"]] = getBeta(se.noob)

bv.raw.means = apply(lpre$bv$raw, 1, mean)
bv.noob.means = apply(lpre$bv$nb, 1, mean)
bv.dif = bv.noob.means - bv.raw.means
summary(bv.dif)

hist(bv.dif, main = "Noob - Raw\n(mean Beta-value differences)")



