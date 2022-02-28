require(recountmethylation)

# get datasets info
dn <- RCurl::getURL("https://recount.bio/data/", 
                    .opts = list(ssl.verifypeer = FALSE))
sm <- servermatrix(dn)
smfn <- sm[,1]
tsnum <- as.numeric(gsub(".*\\_|\\..*", "", smfn), smfn)
tsfilt <- which(tsnum == max(tsnum, na.rm = TRUE))
smf <- sm[tsfilt,]
smf1 <- smf[,1]

# get field info
linfo <- list()

# version 0.0.1
linfo[["h5se_gr"]] <- list(fn = "remethdb-h5se_gr_0-0-1_1590090412",
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

# version 0.0.2
linfo[["h5se-gr_hm450k"]] <- list(fn = "remethdb_h5se-gr_hm450k_0-0-2_1607018051",
                           info = "h5se, HM450K platform samples, GenomicRatioSet, 
                           noob-normalized, version 0.0.2")
linfo[["h5se-gm_hm450k"]] <- list(fn = "remethdb_h5se-gm_hm450k_0-0-2_1607018051",
                           info = "h5se, HM450K platform samples, MethylSet, 
                           raw/unnormalized, version 0.0.2")
linfo[["h5se-rg_hm450k"]] <- list(fn = "remethdb_h5se-rg_hm450k_0-0-2_1607018051",
                           info = "h5se, HM450K platform samples, RGChannelSet, 
                           raw/unnormalized, version 0.0.2")
linfo[["h5-rg_hm450k"]] <- list(fn = "remethdb_h5-rg_hm450k_0-0-2_1607018051.h5",
                         info = "h5, HM450K platform samples, red and green signal, 
                         raw/unnormalized, version 0.0.2")
linfo[["h5se-gr_epic"]] <- list(fn = "remethdb_h5se-gr_epic_0-0-2_1607018051",
                           info = "h5se, EPIC/HM850K platform samples, 
                           GenomicRatioSet, noob-normalized, version 0.0.2")
linfo[["h5se-gm_epic"]] <- list(fn = "remethdb_h5se-gm_epic_0-0-2_1589820348",
                           info = "h5se, EPIC/HM850K platform samples, 
                           MethylSet, raw/unnormalized, version 0.0.2")
linfo[["h5se-rg_epic"]] <- list(fn = "remethdb_h5se-rg_epic_0-0-2_1589820348",
                           info = "h5se, EPIC/HM850K platform samples, 
                           RGChannelSet, raw/unnormalized, version 0.0.2")
linfo[["h5-rg_epic"]] <- list(fn = "remethdb_h5-rg_epic_0-0-2_1589820348.h5 ",
                         info = "h5, EPIC/HM850K platform samples, 
                         red and green signal, raw/unnormalized, version 0.0.2")


ndb <- length(linfo);tagsv <- rep("", ndb);biocv <- BiocManager::version()
dsclass<-rep("FilePath",ndb);locprefix<-rep("https://recount.bio/data/",ndb)
version<-"0.0.2";gbld<-"hg19";stype<-"HDF5";coordb<-"NA";dprov<-"GEO/GDS"
surl<-"http://www.ncbi.nlm.nih.gov/geo/query/";spp<-"Homo sapiens"
taxid<-"9606";svers <- "Nov 7 2020";maint <- maintainer("recountmethylation")
rdclass <- ifelse(grepl("h5se", names(linfo)), "HDF5-SummarizedExperiment",
                  "HDF5Database")
rdpath <- as.character(unlist(lapply(linfo, function(x){x$fn})))

meta <- data.frame(Title = as.character(unlist(lapply(linfo, function(x){x$fn}))),
                   Description = as.character(unlist(lapply(linfo, function(x){x$info}))),
                   BiocVersion = rep(biocv, ndb),
                   Genome = rep(gbld, ndb),
                   SourceType = rep(stype, ndb),
                   SourceUrl = rep(surl, ndb),
                   SourceVersion = rep(svers, ndb),
                   Species = rep(spp, ndb),
                   TaxonomyId = rep(taxid, ndb),
                   Coordinate_1_based = rep(coordb, ndb),
                   DataProvider = rep(dprov, ndb),
                   Maintainer = rep(maint, ndb), 
                   RDataClass = rdclass,
                   DispatchClass = dsclass,
                   Location_Prefix = locprefix,
                   RDataPath = rdpath,
                   Tags = tagsv,
                   Version = rep(version, ndb),
                   stringsAsFactors = FALSE)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)

write.csv(meta, file="metadata.csv", row.names=FALSE)
