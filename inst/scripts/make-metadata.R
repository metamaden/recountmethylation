require(recountmethylation)

# get datasets info
dn <- RCurl::getURL("https://recount.bio/data/", 
                    .opts = list(ssl.verifypeer = FALSE))
sm <- servermatrix(dn)
smfn <- sm[,1]
tsnum <- as.numeric(gsub(".*\\_|\\..*", "", smfn), smfn)
tsfilt <- which(tsnum == max(tsnum, na.rm = T))
smf <- sm[tsfilt,]
smf1 <- smf[,1]

# get field info
linfo <- list()
linfo[["h5se_gr"]] <- list(fn = "remethdb-h5se_gr-test_0-0-1_1590090412",
                           info = "h5se, GenomicRatioSet, noob-normalized, version 0.0.1")
linfo[["h5se_gm"]] <- list(fn = "remethdb-h5se_gm_0-0-1_1590090412",
                           info = "h5se, MethylSet, raw/unnormalized, version 0.0.1")
linfo[["h5se_rg"]] <- list(fn = "remethdb_h5se-rg_0-0-1_1590090412",
                           info = "h5se, RGChannelSet, raw/unnormalized, version 0.0.1")
linfo[["h5se_gr-test"]] <- list(fn = "remethdb-h5se_gr-test_0-0-1_159009041",
                           info = "h5se, GenomicRatioSet, noob-normalized, version 0.0.1")
linfo[["h5_rg"]] <- list(fn = "remethdb-h5_rg_0-0-1_1590090412.h5",
                         info = "h5, red and green signal, raw/unnormalized, version 0.0.1")
linfo[["h5_rg-test"]] <- list(fn = "remethdb-h5_rg-test_0-0-1_1590090412.h5",
                         info = "h5, red and green signal, raw/unnormalized, version 0.0.1")
ndb <- 6
version <- "0.0.1"
biocv <- BiocManager::version()
gbld <- "hg19"
stype <- "tar.gz"
surl <- "http://www.ncbi.nlm.nih.gov/geo/query/"
svers <- "Mar 31 2019"
spp <- "Homo sapiens"
taxid <- "9606"
coordb <- "NA"
dprov <- "GEO/GDS"
maint <- maintainer("recountmethylation")
rdclass <- c(rep("HDF5-SummarizedExperiment", 4), rep("HDF5Database", 2))
dsclass <- c(rep("HDF5-SummarizedExperiment", 4), rep("H5File", 2))
rdpath <- paste("https://recount.bio", "data", 
                as.character(unlist(lapply(linfo, function(x){x$fn}))), 
                sep = "/")
tagsv <- rep("", ndb)

meta <- data.frame(Title = as.character(unlist(lapply(linfo, function(x){x$fn}))),
                   Description = as.character(unlist(lapply(linfo, function(x){x$info}))),
                   BiocVersion = rep(biocv, ndb),
                   Genome = rep(gbld, ndb),
                   SourceType = rep(stype, ndb),
                   SourceURL = rep(surl, ndb),
                   SourceVersion = rep(svers, ndb),
                   Species = rep(spp, ndb),
                   TaxonomyId = rep(taxid, ndb),
                   Coordinate_1_based = rep(coordb, ndb),
                   DataProvider = rep(dprov, ndb),
                   Maintainer = rep(maint, ndb), 
                   RDataClass = rdclass,
                   DispatchClass = dsclass,
                   RDataPath = rdpath,
                   Tags = tagsv,
                   Version = rep(version, ndb),
                   stringsAsFactors = FALSE)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
