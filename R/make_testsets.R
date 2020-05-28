library(recountmethylation)
library(HDF5Array)

#-----------------
# gr h5se test set
#-----------------
# make new test set
# note: scrape metadata from originating object
dn <- "remethdb-h5se_gr_0-0-1_1590090412"
dv <- unlist(strsplit(dn, "_"))
class <- paste0(dv[2], "-test") # "gr-test"
version <- dv[3] # "0-0-1"
ts <- dv[4] # "1590090412
dntest <- paste("remethdb-h5se", class, version, ts, sep = "_")

gr <- loadHDF5SummarizedExperiment(dn)
gsmv <- colnames(gr)[1:2]
cgv <- getAnnotation(gr)
cgv <- rownames(cgv[cgv$chr == "chr22",])
grtest <- gr[cgv,gsmv]
dim(grtest) # [1] 8552    2
saveHDF5SummarizedExperiment(grtest, dir = dntest)

# validate
rg.geo <- gds_idat2rg(gsub("\\..*", "", gsmv))
gr.geo <- preprocessNoob(rg.geo)
bv.geo <- getBeta(gr.geo)
bv.h5se <- as.matrix(getBeta(grtest))
# format colnames
colnames(bv.geo) <- gsub("_.*", "", colnames(bv.geo))
colnames(bv.h5se) <- gsub("\\..*", "", colnames(bv.h5se))
# subset rownames
bv.geo <- bv.geo[rownames(bv.geo) %in% rownames(bv.h5se),]
match.rows <- order(match(rownames(bv.geo), rownames(bv.h5se)))
match.cols <- order(match(colnames(bv.geo), colnames(bv.h5se)))
bv.geo <- bv.geo[match.rows, match.cols]
# check rows/cols identical
identical(rownames(bv.geo), rownames(bv.h5se))
identical(colnames(bv.geo), colnames(bv.h5se))
# check objects
identical(bv.geo, bv.h5se)

# copy test file to package repo
# syspath <- system.file("extdata", package = "recountmethylation")
newdn <- dntest
syspath <- paste("/Users", "maden", "Documents", "GitHub", 
                 "recountmethylation", "inst", "extdata", sep = "/")
dn.syspath <- paste(syspath, newdn, sep = "/")
dir.create(dn.syspath)
lfcp <- list.files(dntest)
for(f in lfcp){
  f.oldpath <- paste(dntest, f, sep = "/")
  f.syspath <- paste(dn.syspath, f, sep = "/")
  file.copy(f.oldpath, f.syspath)
}

#---------------
# h5 test object
#---------------
# make new test set
# note: scrape metadata from originating object
require(rhdf5)
require(minfi)

dn <- "remethdb-h5_rg_1590090412_0-0-1.h5"
dv <- unlist(strsplit(dn, "_"))
class <- paste0(dv[2], "-test") # "gr-test"
version <- "0-0-1"
ts <- "1590090412"
dntest <- paste(dv[1], class, version, ts, sep = "_")
ext <- ".h5"
dntest <- paste0(dntest, ext)

# gr <- loadHDF5SummarizedExperiment(dn)
gsmv <- c("GSM1038309", "GSM1038308")
rg <- getrg(gsmv, dbn = dn)
cgv <- getAnnotation(rg)
cgv <- cgv[cgv$chr == "chr22",]
cg.addrv <- unique(c(cgv$AddressA, cgv$AddressB))
rgtest <- rg[rownames(rg) %in% cg.addrv, ]
dim(rgtest) # [1] 11162     2

# make the test h5 set
{
  h5dbn <- dntest
  h5createFile(h5dbn)
  ngsm <- ncol(rgtest)
  ncg <- nrow(rgtest)
  gsmv <- colnames(rgtest)
  cgv <- rownames(rgtest)
  # get data
  red <- as.matrix(getRed(rgtest))
  grn <- as.matrix(getGreen(rgtest))
  # make new datasets
  dsv <- c("redsignal", "greensignal")
  for(d in dsv){
    if(d == "redsignal"){
      dat <- t(getRed(rgtest))
    } else{
      dat <- t(getGreen(rgtest))
    }
    h5createDataset(h5dbn, d, dims = list(ngsm, ncg), 
                    maxdims = c(H5Sunlimited(), H5Sunlimited()), 
                    storage.mode = "double", level = 5, chunk = c(2, 500))
    h5createDataset(h5dbn, paste(d, "rownames", sep = "."), dims = list(ngsm), 
                    maxdims = c(H5Sunlimited()), storage.mode = "character", 
                    size = 256, level = 5, chunk = c(2))
    h5createDataset(h5dbn, paste(d, "colnames", sep = "."), dims = list(ncg), 
                    maxdims = c(H5Sunlimited()), storage.mode = "character", 
                    size = 256, level = 5, chunk = c(500))
    h5write(gsmv, file = h5dbn, name = paste(d, "rownames", sep = "."), 
            index = list(1:ngsm))
    h5write(cgv, file = h5dbn, name = paste(d, "colnames", sep = "."), 
            index = list(1:ncg))
    h5write(dat, file = h5dbn, name = d, 
            index = list(1:ngsm, 1:ncg))
  }
  # add metadata
  md <- as.matrix(pData(rg))
  h5createDataset(h5dbn, "mdpost", dims = c(nrow(md), ncol(md)),
                  maxdims = c(H5Sunlimited(), H5Sunlimited()),  storage.mode = "character", 
                  level = 5, chunk = c(2, 5), size = 256)
  h5createDataset(h5dbn, paste("mdpost", "colnames", sep = "."), dims = ncol(md),
                  maxdims = c(H5Sunlimited()), storage.mode = "character",
                  level = 5, chunk = c(5), size = 256)
  h5closeAll()
  h5ls(h5dbn)
}

# get data
{
  # geo
  rg.geo <- gds_idat2rg(gsmv)
  # h5 test
  rg.h5 <- getrg(gsmv, cgv = cg.addrv, dbn = h5dbn, metadata = FALSE, verbose = TRUE)
  
}

# match sets
{
  rows.filt <- rownames(rg.geo) %in% rownames(rg.h5)
  rgf.geo <- rg.geo[rows.filt,]
  colnames(rgf.geo) <- gsub("\\_.*", "", colnames(rgf.geo))
  match.rows <- order(match(rownames(rgf.geo), rownames(rg.h5)))
  match.cols <- order(match(colnames(rgf.geo), colnames(rg.h5)))
  rgf.geo <- rgf.geo[match.rows, match.cols]
  # get matrices
  lh5 <- list("red" = as.matrix(getRed(rg.h5)), "grn" = as.matrix(getGreen(rg.h5)))
  lgeo <- list("red" = as.matrix(getRed(rgf.geo)), "grn" = as.matrix(getGreen(rgf.geo)))
  for(i in 1:2){class(lh5[[i]]) <- class(lgeo[[i]]) <- "integer"}
  # check matrices
  identical(lgeo[["red"]], lh5[["red"]])
  identical(lgeo[["grn"]], lh5[["grn"]])
}

# copy test file to package repo
# syspath <- system.file("extdata", package = "recountmethylation")
newdn <- dntest
syspath <- paste("/Users", "maden", "Documents", "GitHub", 
                 "recountmethylation", "inst", "extdata", sep = "/")
dn.syspath <- paste(syspath, newdn, sep = "/")
dir.create(dn.syspath)
lfcp <- list.files(dntest)
f.oldpath <- newdn
f.syspath <- paste(dn.syspath, dntest, sep = "/")
file.copy(f.oldpath, f.syspath)