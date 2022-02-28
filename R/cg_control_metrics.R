#!/usr/bin/env R

# Functions for quality control of Illumina's BeadArray platforms, including 
# HM450K and EPIC.

#' bactrl
#'
#' Get the BeadArray control metrics from HM450K platform red/grn signals.
#'
#' @param cdf Control probe annotations (optional, NULL, data.frame, cols = properties, rows = probes). 
#' @param baset Either the most informative BeadArray metrics (a.k.a. "reduced" set, 
#' 5 metrics, recommended), or the full set of 17 metrics ("all").
#' @param rg An RGChannelSet object (optional, NULL).
#' @param rs Red signal data (optional, NULL, data.frame, columns are probes, rows are samples, 
#' column names are addresses, rownames are samples/GSM IDs).
#' @param gs Green signal data (optional, NULL, data.frame, columns are probes, rows are samples, 
#' column names are addresses, rownames are samples/GSM IDs).
#' @param baseline Baseline measure for signals (integer, 3000).
#' @param biotin.baseline Baseline to use for biotin controls (integer, 1).
#' @details This function calculates the BeadArray quality metrics based on Illumina's 
#' documentation and previous work (see references). Based on previous work, this function
#' can calculate the full set of 17 metrics (e.g. baset = 'all'), or the reduced set of
#' 5 metrics (e.g. baset = 'reduced'), where the latter is recommended for most purposes. 
#' For additional details, consult the BeadArray metrics vignette in this package.
#' @returns Matrix of BeadArray signal values
#' @references 
#' 
#' 1. S. K. Maden, et. al., "Human methylome variation across Infinium 450K data on the Gene 
#' Expression Omnibus," NAR Genomics and Bioinformatics, Volume 3, Issue 2, June 2021
#' 
#' 2. Illumina, “Illumina Genome Studio Methylation Module v1.8,” Nov. 2010.
#' 
#' 3. Illumina, “BeadArray Controls Reporter Software Guide,” Oct. 2015.
#' 
#' 4. J. A. Heiss and A. C. Just, “Identifying mislabeled and contaminated DNA methylation 
#' microarray data: an extended quality control toolset with examples from GEO,” Clinical 
#' Epigenetics, vol. 10, June 2018.
#''
#' @examples 
#' dir <- system.file("extdata", "bactrl", package = "recountmethylation")
#' rgf <- get(load(file.path(dir, "rgf-cgctrl-test_hm450k-minfidata.rda")))
#' mba <- bactrl(rg = rgf)
#' @seealso bathresh
#' @export
bactrl <- function(cdf = NULL, baset = "reduced", 
                   rg = NULL, rs = NULL, gs = NULL, 
                   baseline = 3000, biotin.baseline = 1){
  message("Getting control probes info...")
  if(is.null(cdf)){
    message("Loading control probes info...")
    dir <- system.file("extdata", "bactrl", package = "recountmethylation")
    cdf <- get(load(file.path(dir, "dfaddr_control-probes.rda")))
  }
  message("Prepping signal matrices...")
  if(is.null(rg)){
    if(is.null(rs)|is.null(gs)){
      stop("Error: provide an RGChannelSet or signal matrices.")}
  } else{rs <- minfi::getRed(rg); gs <- minfi::getGreen(rg)}
  rs <- rs[rownames(rs) %in% cdf$Address,]
  gs <- gs[rownames(gs) %in% cdf$Address,]
  if(nrow(rs) == 0 | nrow(gs) == 0){
    stop(paste0("Error: one or both signal matrices don't contain ",
                "any control probe addresses. Are the rows labeled",
                " with valid probe addresses?"))}
  rm <- data.frame(sample_id = colnames(rs)) # return matrix
  addr.bkg <- recountmethylation:::ba.background(cdf)
  if(length(addr.bkg) == 1 & is.na(addr.bkg[1])){
    stop("Error: couldn't get addresses for background signal probes.")}
  message("Calculating BeadArray metrics...")
  batot <- ifelse(baset == "reduced", 5, 17)
  if(baset == "reduced"){
    message("Calculating reduced metric set...")} else{
    message("Calculating the full BeadArray metric set...")
  }
  rm <- recountmethylation:::ba.biotinstaining.red(rs = rs, 
    biotin.baseline = biotin.baseline, rm = rm, cdf = cdf, mtot = batot) 
  rm <- recountmethylation:::ba.biotinstaining.grn(gs = gs, rm = rm, cdf = cdf, 
    biotin.baseline = biotin.baseline, mtot = batot) 
  rm <- recountmethylation:::ba.nonpolymorphic.red(rs = rs, rm = rm, cdf = cdf, 
    mtot = batot) 
  rm <- recountmethylation:::ba.nonpolymorphic.grn(gs = gs, rm = rm, cdf = cdf, 
    mtot = batot) 
  rm <- recountmethylation:::ba.bisulfiteconv1.red(rs = rs, rm = rm, cdf = cdf, 
    mtot = batot) 
  if(baset == "reduced"){
    message("Returning reduced metric set...")
    return(rm)
  }
  rm <- recountmethylation:::ba.restoration(gs = gs, rm = rm, cdf = cdf, 
    addr.bkg = addr.bkg) 
  rm <- recountmethylation:::ba.specificity1.red(rs = rs, rm = rm, cdf = cdf)
  rm <- recountmethylation:::ba.specificity1.grn(gs = gs, rm = rm, cdf = cdf)
  rm <- recountmethylation:::ba.specificity2(rs = rs, gs = gs, rm = rm, 
    cdf = cdf)
  rm <- recountmethylation:::ba.extension.red(rs = rs, rm = rm, cdf = cdf)
  rm <- recountmethylation:::ba.extension.grn(gs = gs, rm = rm, cdf = cdf)
  rm <- recountmethylation:::ba.hybridization.hi.vs.med(gs = gs, rm = rm, 
    cdf = cdf)
  rm <- recountmethylation:::ba.hybridization.med.vs.low(gs = gs, rm = rm, 
    cdf = cdf)
  rm <- recountmethylation:::ba.targetremoval1(gs = gs, rm = rm, cdf = cdf, 
    baseline = baseline)
  rm <- recountmethylation:::ba.targetremoval2(gs = gs, rm = rm, cdf = cdf,
    baseline = baseline)
  rm <- recountmethylation:::ba.bisulfiteconv1.grn(gs = gs, rm = rm, cdf = cdf)
  rm <- recountmethylation:::ba.bisulfiteconv2(rs = rs, gs = gs, rm = rm, 
    cdf = cdf)
  message("Returning the full BeadArray metric set...")
  return(rm)
}

#' Get BeadArray control outcomes from a matrix of metric signals
#'
#' Get Illumina's prescribed minimum quality thresholds for BeadArray metrics. 
#'
#' @returns Data frame of minimum BeadArray quality thresholds.
#' @references 
#' 1. Illumina, “Illumina Genome Studio Methylation Module v1.8,” Nov. 2010.
#'
#' 2. Illumina, “BeadArray Controls Reporter Software Guide,” Oct. 2015.
#'
#' @examples 
#' dfthresh <- bathresh()
#' @seealso bactrl
#' @export
bathresh <- function(){
  data.frame(restoration.grn = 0, biotin.stain.red = 5, biotin.stain.grn = 5,
                   specificityI.red = 1, specificityI.grn = 1, specificityII = 1,
                   extension.red = 5, extension.grn = 5, hyb.hi.med = 1, hyb.med.low = 1,
                   target.removal.1 = 1, target.removal.2 = 1, bisulfite.conv.I.red = 1,
                   bisulfite.conv.I.grn = 1, bisulfite.conv.II = 1, nonpolymorphic.red = 5, 
                   nonpolymorphic.grn = 5, stringsAsFactors = F)
}

#' get_qcsignal
#'
#' Get the medians of the log2-transformed M and U signals. This function 
#' uses the DelayedMatrixStats implementations of colMedians for rapid
#' calculations on DelayedArray-formatted matrices.
#' 
#' @param se Valid SummarizedExperiment object, such as a MethylSet or similar object
#' for which getMeth() and getUnmeth() methods are defined (optional).
#' @param mm Matrix of methylated/M signals (optional, not required if se provided).
#' @param mu Matrix of unmethylated/U signals (optional, not required if se provided).
#' @param sample_idv Vector of sample IDs to label rows in the returned data frame 
#' (optional, uses mm colnames instead if not provided).
#' @returns Data frame of signal summaries.
#' @details Calculates the log2 of median signal for methylated/M and 
#' unmethylated/U signals separately.
#' @examples
#' library(minfiData)
#' data(MsetEx)
#' se <- MsetEx
#' class(se)
#' # [1] "MethylSet"
#' # attr(,"package")
#' # [1] "minfi"
#' ms <- get_qcsignal(se)
#' @seealso bactrl
#' @export
get_qcsignal <- function(se = NULL, mm = NULL, mu = NULL, sample_idv = NULL){
  if(is.null(se)){
    if(is.null(mm)|is.null(mu)){
      stop("Error: Provide either a valid se object, or both mm and mu")}
  } else{
    if(!is(se, "MethylSet")){stop("Error: se must be a MethylSet.")}
    mm <- minfi::getMeth(se); mu <- minfi::getUnmeth(se)}
  lsignal <- list(meth = mm, unmeth = mu)
  message("Calculating log2 median signals...")
  ms <- as.data.frame(do.call(cbind, lapply(seq(2), function(ii){
    log2(DelayedMatrixStats::colMedians(lsignal[[ii]]))
  })), stringsAsFactors = F)
  colnames(ms) <- c("log2_meth", "log2_unmeth")
  message("Labeling rows with sample ids...")
  if(is.null(sample_idv)){
    rownames(ms) <- colnames(mm)
    } else{
      if(length(sample_idv)==nrow(ms)){
        rownames(ms) <- sample_idv
        } else{
        message("Warning: `sample_idv`` length not equal to nrow",
                " in ms. Using mm colnames instead.")
      }
    }
  return(ms)
}

#' get_crossreactive_cpgs
#' 
#' Get cross-reactive CpG probe IDs for Illumina BeadArray platforms.
#' 
#' @param probeset Specify the set of probes to filter ("all", "hm450k", "epic",
#' "chen", "pidlsey", "illumina").
#' @returns Vector of cross-reactive CpG probe IDs.
#' @details Prior work showed significant cross-reactivity at subsets of CpG probes on 
#' Illumina's BeadArray platforms, including HM450K and EPIC. This was primarily due to 
#' the probe sequence, as the targeted 50-bp sequence can be either too short or too 
#' degenerate to bind a particular DNA region with high specificity. This can cause 
#' cross-reaction with off-target DNA locations, including at entirely different 
#' chromosomes than the target sequence. Consult the individual publication sources for 
#' details about the identification and consequences of cross-reactive CpG probes.
#' 
#' You can retrieve a cross-reactive probe set in a variety of ways. For instance, 
#' declare the publication source with either "chen" (for Chen et al 2013), "pidsley" 
#' (for Pidsley et al 2016), or "illumina" (for official Illumina documentation), or 
#' declare the platform category as either "all" (both HM450K and EPIC), "hm450k", or 
#' "epic." 
#' 
#' @references 
#' 1. Yi-an Chen, Mathieu Lemire, Sanaa Choufani, Darci T. Butcher, Daria Grafodatskaya, 
#' Brent W. Zanke, Steven Gallinger, Thomas J. Hudson & Rosanna Weksberg (2013) 
#' Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium
#' HumanMethylation450 microarray, Epigenetics, 8:2, 203-209, DOI: 10.4161/epi.23470
#'
#' 2. Pidsley, R., Zotenko, E., Peters, T.J. et al. Critical evaluation of the Illumina
#' MethylationEPIC BeadChip microarray for whole-genome DNA methylation profiling.
#' Genome Biol 17, 208 (2016). https://doi.org/10.1186/s13059-016-1066-1
#'
#' @examples
#' length(get_crossreactive_cpgs("all"))      # 46324
#' length(get_crossreactive_cpgs("hm450k"))   # 30540
#' length(get_crossreactive_cpgs("epic"))     # 43410
#' length(get_crossreactive_cpgs("chen"))     # 29233
#' length(get_crossreactive_cpgs("pidsley"))  # 43254
#' length(get_crossreactive_cpgs("illumina")) # 1031
#'
#' @seealso bactrl, get_qcsignal
#' @export
get_crossreactive_cpgs <- function(probeset = "all"){
  valid.setv <- c("all", "hm450k", "epic", "chen", "pidsley", "illumina")
  if(probeset %in% valid.setv){
    data.path <- system.file("extdata", "crossreactive_cpgprobes",
                             package = "recountmethylation")
    cgidv.fpath <- file.path(data.path, 
                             paste0(probeset, "_crx_cgv.rda"))
    cgidv <- get(load(cgidv.fpath))
  } else{
    stop("Error: Didn't recognize provided `probeset` argument. Is it a ",
         "valid option (e.g. either ", paste0("'", valid.setv, "'", 
                                              collapse = ", "), ")?")
  }
  return(cgidv)
}
