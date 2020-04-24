head(list.files("/pfs/ccle_rnaseq_kallisto_0_46_1"))

library(PharmacoGx)
library(genefu)
library(data.table)
library(reshape2)
library(Biobase)

options(stringsAsFactors = FALSE)


    matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
      sapply(ids, function(x) {
        myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
        if(length(myx) > 1){
          stop("Something went wrong in curating ids, we have multiple matches")
        }
        if(length(myx) == 0){return(NA_character_)}
        return(tbl[myx, returnColumn])
      })
    }
    
myDirPrefix <- "/pfs/"


args = commandArgs(trailingOnly=TRUE)
rnaseq_select <- args
print(rnaseq_select)
rnaseq_results <- list()
ORCESTRA_ID = tail(rnaseq_select, n=1)

	  
tools <- grep(pattern = 'Kallisto|Salmon', x = rnaseq_select)
tools <- rnaseq_select[tools]
tools <- gsub("-", "_", tools)
transcriptome <- grep(pattern = 'Gencode|Ensembl', x = rnaseq_select)
transcriptome <- rnaseq_select[transcriptome]
tool_path = expand.grid(a = tools,b = transcriptome)
tool_path = paste0(tool_path$a, "_",tool_path$b)
	  
print(tool_path)	  
	  
	  
	  
	  
	  
	  
	  
    badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
    ## drug information
    message("Read drug information")
    druginfo <- read.csv("/pfs/get_CCLE/CCLE_NP24.2009_profiling_2012.02.20.csv", stringsAsFactors=FALSE)
    druginfo[druginfo == "" | druginfo == " "] <- NA
    
    ##########################################################################
    ## manual drug name curation
    druginfo[druginfo[ , "Compound..code.or.generic.name."] == "Panobinostat\xa0\xa0", "Compound..code.or.generic.name."] <- "Panobinostat"
    #druginfo[druginfo[ , "Compound..code.or.generic.name."] == "PF-2341066", "Compound..code.or.generic.name."] <- "Crizotinib"
    druginfo[druginfo[ , "Compound..code.or.generic.name."] == "PD-0332991 ", "Compound"] <- "PD-0332991"
    
    ##########################################################################
    
    
    druginfo <- data.frame(druginfo, "drugid"=paste("drugid", toupper(gsub(pattern=badchars, replacement="", x=toupper(druginfo[ , "Compound..code.or.generic.name."]))), sep="_"), check.names=FALSE)
    rownames(druginfo) <- druginfo[ , "drugid"]
    
    ## profiles for the drugs
    message("Read drug sensitivity measurements")
    ## drug sensitivity data from the addendum in Nature
    
    
    #      drugpheno <- gdata::read.xls(xls=file.path(path.drug, "ccle_drug_pheno_file.xls"), sheet=12)
    drugpheno <- read.csv("/pfs/get_CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv")
    drugpheno[drugpheno == "" | drugpheno == " "] <- NA
    drugpheno[drugpheno[ , "Compound"]=="ZD-6474", "Compound"] <- "Vandetanib"
    drugpheno[drugpheno[ , "Compound"]=="PF2341066", "Compound"] <- "PF-2341066"
    
    
    
    drugpheno <- data.frame(drugpheno, "drugid"=paste("drugid", as.character(drugpheno[ , "Compound"]), sep="_"), "cellid"=as.character(drugpheno[ , "Primary.Cell.Line.Name"]), check.names=FALSE)
    
    ## drug screening concentrations
    ll <- sapply(strsplit(drugpheno[ , "Doses..uM."], ","), length)
    drugconc <- lapply(strsplit(drugpheno[ , "Doses..uM."], ","), function(x) {
      xx <- as.numeric(x)
      xx2 <- c(0.0025, 0.0080, 0.0250, 0.0800, 0.2500, 0.8000, 2.5300, 8.0000)
      if(any(!is.element(xx, xx2))) { stop("Unexpected drug screening concentrations!") }
      xx3 <- rep(NA, length(xx2))
      names(xx3) <- xx2
      xx3[match(xx, names(xx3))] <- xx
      return(xx3)
    })
    drugconc <- do.call(rbind, drugconc)
    drugconc <- data.frame("cellid"=as.character(drugpheno[ , "cellid"]), "drugid"=as.character(drugpheno[ , "drugid"]), "nbr.conc.tested"=ll, drugconc, check.names=FALSE)
    dimnames(drugconc) <- list(paste(drugconc[ , "drugid"], drugconc[ , "cellid"], sep="_"), c("cellid", "drugid", "nbr.conc.tested", sprintf("Dose%i.uM", 1:(ncol(drugconc) - 3))))
    
    ## combine all drugs
    dix <- sort(unique(c(as.character(druginfo[ , "drugid"]), as.character(drugpheno[ , "drugid"]))))
    ## update druginfo
    druginfo2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(druginfo), dimnames=list(dix, colnames(druginfo))), check.names=FALSE)
    newlev <- sapply(druginfo, levels)
    newlev$drugid <- dix
    druginfo2 <- setcolclass.df(df=druginfo2, colclass=sapply(druginfo, class), factor.levels=newlev)
    druginfo2[match(as.character(druginfo[ , "drugid"]), dix), colnames(druginfo)] <- druginfo
    druginfo2[ , "drugid"] <- newlev$drugid
    druginfo <- druginfo2
    celln <- unique(as.character(drugpheno[ , "cellid"]))
    drugpheno.all <- NULL
    ## IC50
    drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
    for(i in 1:length(celln)) {
      tt <- as.numeric(drugpheno[drugpheno[ , "cellid"] == celln[i], "IC50..uM."])
      names(tt) <- as.character(drugpheno[drugpheno[ , "cellid"] == celln[i] , "drugid"])
      drugphenot[celln[i],names(tt)] <- tt
    }
    colnames(drugphenot) <- paste(colnames(drugphenot), "IC50", sep="_")
    if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- cbind(drugpheno.all, drugphenot, check.names=FALSE) }
    ## EC50..uM.
    drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
    for(i in 1:length(celln)) {
      tt <- as.numeric(drugpheno[drugpheno[ , "cellid"] == celln[i] , "EC50..uM."])
      names(tt) <- as.character(drugpheno[drugpheno[ , "cellid"] == celln[i] , "drugid"])
      drugphenot[celln[i],names(tt)] <- tt
    }
    colnames(drugphenot) <- paste(colnames(drugphenot), "EC50", sep="_")
    if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- cbind(drugpheno.all, drugphenot, check.names=FALSE) }
    ## Amax
    drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
    for(i in 1:length(celln)) {
      tt <- as.numeric(drugpheno[drugpheno[ , "cellid"] == celln[i] , "Amax"])
      names(tt) <- as.character(drugpheno[drugpheno[ , "cellid"] == celln[i] , "drugid"])
      drugphenot[celln[i],names(tt)] <- tt
    }
    colnames(drugphenot) <- paste(colnames(drugphenot), "Amax", sep="_")
    if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- cbind(drugpheno.all, drugphenot, check.names=FALSE) }
    ## ActArea
    drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
    for(i in 1:length(celln)) {
      tt <- as.numeric(drugpheno[drugpheno[ , "cellid"] == celln[i] , "ActArea"])
      names(tt) <- as.character(drugpheno[drugpheno[ , "cellid"] == celln[i] , "drugid"])
      drugphenot[celln[i],names(tt)] <- tt
    }
    colnames(drugphenot) <- paste(colnames(drugphenot), "ActivityArea", sep="_")
    if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- cbind(drugpheno.all, drugphenot, check.names=FALSE) }
    ## Dose
    drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
    for(i in 1:length(celln)) {
      tt <- drugpheno[drugpheno[ , "cellid"] == celln[i] , "Doses..uM."]
      names(tt) <- as.character(drugpheno[drugpheno[ , "cellid"] == celln[i] , "drugid"])
      drugphenot[celln[i],names(tt)] <- tt
    }
    colnames(drugphenot) <- paste(colnames(drugphenot), "Doses", sep="_")
    if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- cbind(drugpheno.all, drugphenot, check.names=FALSE) }
    ## Activity.Data..median
    drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
    for(i in 1:length(celln)) {
      tt <- drugpheno[drugpheno[ , "cellid"] == celln[i] , "Activity.Data..median."]
      names(tt) <- as.character(drugpheno[drugpheno[ , "cellid"] == celln[i] , "drugid"])
      drugphenot[celln[i],names(tt)] <- tt
    }
    colnames(drugphenot) <- paste(colnames(drugphenot), "ActivityData", sep="_")
    if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- cbind(drugpheno.all, drugphenot, check.names=FALSE) }
    ## Activity.SD
    drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
    for(i in 1:length(celln)) {
      tt <- drugpheno[drugpheno[ , "cellid"] == celln[i] , "Activity.SD"]
      names(tt) <- as.character(drugpheno[drugpheno[ , "cellid"] == celln[i] , "drugid"])
      drugphenot[celln[i],names(tt)] <- tt
    }
    colnames(drugphenot) <- paste(colnames(drugphenot), "Activity.SD", sep="_")
    if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- cbind(drugpheno.all, drugphenot, check.names=FALSE) }
    ## FitType
    drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
    for(i in 1:length(celln)) {
      tt <- drugpheno[drugpheno[ , "cellid"] == celln[i] , "FitType"]
      names(tt) <- as.character(drugpheno[drugpheno[ , "cellid"] == celln[i] , "drugid"])
      drugphenot[celln[i],names(tt)] <- tt
    }
    colnames(drugphenot) <- paste(colnames(drugphenot), "FitType", sep="_")
    if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- data.frame(drugpheno.all, drugphenot, check.names=FALSE) }
    ## save the new spreadsheet
    drugpheno <- data.frame("cellid"=rownames(drugpheno.all), drugpheno.all, check.names=FALSE)
    
    ## info about each experiment
    message("Read sample information")
    sampleinfo <- read.csv("/pfs/get_CCLE/CCLE_sample_info_file_2012-10-18.txt", sep="\t")
    sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
    sampleinfo <- data.frame("cellid"=as.character(sampleinfo[ , "Cell.line.primary.name"]), sampleinfo, check.names=FALSE)
    ## remove duplicated cell line hybridization
    ## only the most recent experiment (as determine by hyridization date or CEL file timestamp) will be kept for each cell line
    # sampleinfo <- sampleinfo[!duplicated(sampleinfo[ , "cellid"]), ,drop=FALSE]
    #     sampleinfo[ , "cellid"] <- as.character(sampleinfo[ , "cellid"])
    rownames(sampleinfo) <- as.character(sampleinfo[ , "cellid"])
    
    
    
    iix <- which(!is.element(rownames(drugpheno), rownames(sampleinfo)))
    celline.ccle <- data.frame(matrix(NA, nrow=nrow(sampleinfo) + length(iix), ncol=ncol(sampleinfo), dimnames=list(c(rownames(sampleinfo), rownames(drugpheno)[iix]), colnames(sampleinfo))), check.names=FALSE)
    celline.ccle[rownames(sampleinfo), colnames(sampleinfo)] <- sampleinfo
    celline.ccle[rownames(drugpheno)[iix], "Cell.line.primary.name"] <- rownames(drugpheno)[iix]
    celline.ccle[rownames(drugpheno)[iix], "cellid"] <- rownames(drugpheno)[iix]
    ## add url based on CCLE name
    uurl <- paste("http://www.broadinstitute.org/ccle/cell%20lines/", celline.ccle[ , "CCLE.name"], sep="")
    uurl[is.na(celline.ccle[ , "CCLE.name"])] <- NA
    celline.ccle <- data.frame("cellid"=celline.ccle[ , "cellid"], "link"=uurl, celline.ccle[ , !is.element(colnames(celline.ccle), "cellid")], check.names=FALSE)
    
    cellnall <- sort(unique(c(as.character(sampleinfo[ , "cellid"]), as.character(drugpheno[ , "cellid"]))))
    drugnall <- sort(unique(c(rownames(druginfo), as.character(drugconc[ , "drugid"]), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))
    ## update drugpheno
    ## IC50 in microM
    iix <- grep("_IC50$", colnames(drugpheno))
    iixn <- gsub("_IC50$", "", colnames(drugpheno)[iix])
    dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
    dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
    drugpheno.ic50 <- dd
    
    
    iix <- grep("_Amax$", colnames(drugpheno))
    iixn <- gsub("_Amax$", "", colnames(drugpheno)[iix])
    dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
    dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
    drugpheno.amax <- dd
    
    iix <- grep("_ActivityArea$", colnames(drugpheno))
    iixn <- gsub("_ActivityArea$", "", colnames(drugpheno)[iix])
    auct <- drugpheno[ , iix]
    colnames(auct) <- iixn
    ## division by the number of concentrations tested
    myx <- sapply(strsplit(colnames(drugpheno), "_"), function(x) { return(x[length(x)] == "Doses") })
    ndoses <- drugpheno[ ,myx,drop=FALSE]
    ndoses <- apply(ndoses, 2, function(x) {
      return(sapply(x, function(x) {
        return(sapply(strsplit(as.character(x), ","), function(x) { if(is.na(x[1])) { return(NA) } else { return(length(x)) } }))
      }))
    })
    colnames(ndoses) <- gsub("_Doses", "", colnames(ndoses))
    ndoses <- ndoses[rownames(auct), colnames(auct)]
    
    ## compute AUC
    iix <- grep("_ActivityData$", colnames(drugpheno))
    names(iix) <- gsub("_ActivityData$", "", colnames(drugpheno)[iix])
    iix2 <- grep("_Doses$", colnames(drugpheno))
    names(iix2) <- gsub("_Doses$", "", colnames(drugpheno)[iix2])
    myx <- intersect(names(iix), names(iix2))
    iix <- iix[myx]
    iix2 <- iix2[myx]
    xx <- array(NA, dim=c(nrow(drugpheno), length(iix), 2), dimnames=list(rownames(drugpheno), names(iix), c("ActArea", "Doses")))
    xx[ , , "ActArea"] <- apply(drugpheno[ , iix, drop=FALSE], 2, as.character)
    xx[ , , "Doses"] <- apply(drugpheno[ , iix2, drop=FALSE], 2, as.character)
    auct <- apply(X=xx, MARGIN=c(1,2), FUN=function(x) {
      aa <- - as.numeric(unlist(strsplit(x[1], ","))) / 100
      dd <- as.numeric(unlist(strsplit(x[2], ","))) / 100
      ## average
      xx <- sum(aa) / sum(!is.na(aa))
      ## compute AUC using trapezoidal estimation
      # xx <- MESS::auc(x=1:length(aa), y=aa, type="linear") / 8
      ## compute AUC using Simpson's rule
      # xx <- MESS::auc(x=dd, y=aa, type="spline")
      return(xx)
    })
    dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
    dd[rownames(auct), colnames(auct)] <- data.matrix(auct)
    dd[dd < 0] <- 0
    dd[dd > 1] <- 1
    drugpheno.auc <- dd
    
    nms <- sapply(colnames(drugpheno.ic50), paste, rownames(drugpheno.ic50), sep="_")
    ic50 <- c(drugpheno.ic50)
    names(ic50) <- nms
    auc <- c(drugpheno.auc)
    names(auc) <- nms
    amax <- c(drugpheno.amax)
    names(amax) <- nms
    
    profiles <- cbind("ic50_published"=ic50, "aac_published"=auc, "amax_published"=amax)
    drugpheno.auc.ccle <-drugpheno.auc
    drugpheno.ic50.ccle <- drugpheno.ic50
    drugpheno.amax.ccle <- drugpheno.amax
    
    
    
    profiles <- profiles[complete.cases(profiles),]
    
    load("/pfs/ccleProfilesAssemble/profiles.RData")
    load("/pfs/ccleRawSensitivity/drug_norm_post.RData")

    recomputed <- res
    
    drugconc <- drugconc[rownames(raw.sensitivity),]
    duration <- rep(x=72, length=nrow(drugconc))
    sensitivityInfo <- cbind(drugconc, "duration_h"=duration)
    profiles <- profiles[rownames(sensitivityInfo),]
    
    print("Compute AMAX")
    #Amax <- NULL
    #for (exp in names(raw.sensitivity)){
      #Amax <- c(Amax, computeAmax(conc=raw.sensitivity[exp, , "Dose"], viability=raw.sensitivity[exp, , "Viability"]))
    #}
    #names(Amax) <- names(raw.sensitivity)

    
    profiles <- cbind(profiles, recomputed[rownames(profiles),])
    print("colnames of profiles")
    print(colnames(profiles))
    sens.profiles <-  data.frame("ic50_published" = as.numeric(profiles[,"ic50_published"]),
                        "aac_published" = as.numeric(profiles[,"aac_published"]),
                        "amax_published" = as.numeric(profiles[,"amax_published"]),
                        "aac_recomputed" = as.numeric(profiles[,"AAC"]), 
                        "ic50_recomputed"=as.numeric(profiles[,"IC50"]), 
                        "HS"=as.numeric(profiles[,"HS"]), 
                        "E_inf"=as.numeric(profiles[,"E_inf"]), 
                        "EC50"=as.numeric(profiles[,"EC50"]))

    #profiles <- cbind(profiles, "amax_recomputed"= Amax)
    sens.profiles$aac_recomputed <- sens.profiles$aac_recomputed/100
    rownames(sens.profiles) <- rownames(profiles)
    profiles <- sens.profiles
    print(head(profiles))
	  
    print("Profiles done")
    ### Temporary solution while we wait for the release of PharmacoDb!!!
    
    druginfo <- cbind(druginfo, "drug.name"= druginfo[ , "Compound..code.or.generic.name."])
    
    cell_all <- read.csv("/pfs/downAnnotations/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
    
    curationCell <- cell_all[which(!is.na(cell_all[ , "CCLE.cellid"]) | !is.na(cell_all[, "CCLE_rnaseq.cellid"])),]
    curationTissue <- cell_all[which(!is.na(cell_all[ , "CCLE.cellid"]) | !is.na(cell_all[, "CCLE_rnaseq.cellid"])),]
    
    curationCell <- curationCell[ , c("unique.cellid", "CCLE.cellid", "CCLE_rnaseq.cellid")]
    curationTissue <- curationTissue[ , c("unique.tissueid", "CCLE.tissueid", "CCLE_rnaseq.tissueid")]
    
    rownames(curationTissue) <- curationCell[ , "unique.cellid"]
    rownames(curationCell) <- curationCell[ , "unique.cellid"]
    
    print("Curation Cell Done")
    drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))
    curationDrug <- drug_all[which(!is.na(drug_all[ , "CCLE.drugid"])),]
    curationDrug <- curationDrug[ , c("unique.drugid", "CCLE.drugid")]
    rownames(curationDrug) <- curationDrug[ , "unique.drugid"]
    
    print("Curation Drug Done")
    curationDrug[which(curationDrug[ , "CCLE.drugid"]=="PD-0332991 "), "CCLE.drugid"] <- "PD-0332991"
    
    #change cell slot rownames to unique ids
    celline.ccle <- celline.ccle[celline.ccle$cellid %in% curationCell$CCLE.cellid, ]
    rownames(celline.ccle) <- rownames(curationCell)[match(celline.ccle[ , "cellid"], curationCell[ , "CCLE.cellid"])]
    celline.ccle[ , "cellid"] <- rownames(celline.ccle)
    #integrate cell slot with all the extra rnaseq cells
    celline.ccle <- cbind(celline.ccle, "CCLE_rnaseq.cellid" = curationCell[rownames(celline.ccle), "CCLE_rnaseq.cellid"])
    tt <- as.data.frame(matrix(NA, ncol=ncol(celline.ccle), nrow=length(which(!(rownames(curationCell) %in% rownames(celline.ccle))))), stringsAsFactors=FALSE)
    colnames(tt) <- colnames(celline.ccle)
    tt[ , "CCLE_rnaseq.cellid"] <- rownames(curationCell)[which(!rownames(curationCell) %in% rownames(celline.ccle) )]
    rownames(tt) <- rownames(curationCell)[match(tt[ , "CCLE_rnaseq.cellid"], curationCell[ , "unique.cellid"])]
    celline.ccle <- rbind(celline.ccle,tt)
    
    ##update sensitivity info cellid to be unique.cellid
    sensitivityInfo[ , "cellid"] <- as.character(matchToIDTable(ids=sensitivityInfo[ , "cellid"], tbl=curationCell, column = "CCLE.cellid", returnColumn = "unique.cellid"))
    ###!!!!!!!!!!!!!!!!!!!!!!!
    ##check for PD-0332991 next time
    ###!!!!!!!!!!!!!!!!!!!!!!!
    ##update drugid to be unique dugid
    rownames(curationDrug) <- curationDrug[ , "unique.drugid"]
    druginfo$drug.name[which(druginfo$drug.name =="PD-0332991 ")] <- "PD-0332991"
    sensitivityInfo[ , "drugid"] <- curationDrug[match(gsub("drugid_", "",sensitivityInfo[ , "drugid"]), curationDrug[ , "CCLE.drugid"]), "unique.drugid"]
    druginfo <- druginfo[druginfo$drug.name %in% curationDrug$CCLE.drugid, ]
    rownames(druginfo) <- curationDrug$unique.drugid[match(druginfo$drug.name, curationDrug$CCLE.drugid)]
    
    
    ##update tissueid to be unique tissueid
    celline.ccle[ , "tissueid"] <- curationTissue[rownames(celline.ccle), "unique.tissueid"]

    summarizeRnaSeq <- function (dir, 
                             features_annotation,
                             samples_annotation,
			      method) {
  library(Biobase)
  library(readr)
  library(tximport)
  
  load(features_annotation)
    
  tx2gene <- as.data.frame(cbind("transcript"=tx2gene$transcripts, "gene"=tx2gene$genes))
  
  files <- list.files(dir, recursive = TRUE, full.names = T)
  if(method=="kallisto"){
  resFiles <- grep("abundance.h5", files)
  }else{
  resFiles <- grep("quant.sf", files)
  }
  resFiles <- files[resFiles]
  length(resFiles)
  names(resFiles) <- basename(dirname(resFiles))
  
  if(features_annotation == "/pfs/downAnnotations/Ensembl.v99.annotation.RData"){
  txi <- tximport(resFiles, type=method, tx2gene=tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)
  } else{
  txi <- tximport(resFiles, type=method, tx2gene=tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = FALSE)	  
  }
	  
  head(txi$counts[,1:5])
  dim(txi$counts)
	  
  xx <- txi$abundance
  gene.exp <- Biobase::ExpressionSet(log2(xx + 0.001))
  fData(gene.exp) <- features_gene[featureNames(gene.exp),]
  pData(gene.exp) <- samples_annotation[sampleNames(gene.exp),]
  annotation(gene.exp) <- "rnaseq"
  
  xx <- txi$counts
  gene.count <- Biobase::ExpressionSet(log2(xx + 1))
  fData(gene.count) <- features_gene[featureNames(gene.count),]
  pData(gene.count) <- samples_annotation[sampleNames(gene.count),]
  annotation(gene.count) <- "rnaseq"
  
  txii <- tximport(resFiles, type=method, txOut=T)
  
  if(features_annotation == "/pfs/downAnnotations/Ensembl.v99.annotation.RData"){
  #remove non-coding transcripts in ensembl 	  
  rownames(txii$abundance) <-  gsub("\\..*","",rownames(txii$abundance))
  txii$abundance[which(!rownames(txii$abundance)  %in% features_transcript$transcript_id)]
  missing_transcript <- rownames(txii$abundance)[which(!rownames(txii$abundance)  %in% features_transcript$transcript_id)]
  txii$abundance <- txii$abundance [-which(rownames(txii$abundance) %in% missing_transcript),]
  }
  	  
  xx <- txii$abundance
  transcript.exp <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 0.001))
  if(features_annotation == "/pfs/downAnnotations/Gencode.v33.annotation.RData" || features_annotation == "/pfs/downAnnotations/Gencode.v33lift37.annotation.RData"){
  featureNames(transcript.exp) <- gsub("\\|.*","",featureNames(transcript.exp))
  fData(transcript.exp) <- features_transcript[featureNames(transcript.exp),]
  }else{
  fData(transcript.exp) <- features_transcript[featureNames(transcript.exp),]
  }
  pData(transcript.exp) <- samples_annotation[sampleNames(transcript.exp),]
  annotation(transcript.exp) <- "isoforms"
  
	  
  if(features_annotation == "/pfs/downAnnotations/Ensembl.v99.annotation.RData"){
  #remove non-coding transcripts in ensembl
  rownames(txii$counts) <-  gsub("\\..*","",rownames(txii$counts))
  txii$counts <- txii$counts [-which(rownames(txii$counts) %in% missing_transcript),]	  
  }	  
  xx <- txii$counts
  transcript.count <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 1))
  if(features_annotation == "/pfs/downAnnotations/Gencode.v33.annotation.RData" || features_annotation == "/pfs/downAnnotations/Gencode.v33lift37.annotation.RData"){
  featureNames(transcript.count) <- gsub("\\|.*","",featureNames(transcript.count))
  fData(transcript.count) <- features_transcript[featureNames(transcript.count),]
  }else{
  fData(transcript.count) <- features_transcript[featureNames(transcript.count),]
  }
  pData(transcript.count) <- samples_annotation[sampleNames(transcript.count),]
  annotation(transcript.count) <- "isoforms"
  
  pData(gene.exp)[ ,"batchid"] <- NA
  pData(gene.count)[ ,"batchid"] <- NA	  
  pData(transcript.exp)[ ,"batchid"] <- NA
  pData(transcript.count)[ ,"batchid"] <- NA
	    
  return(list("rnaseq"=gene.exp, 
              "rnaseq.counts"=gene.count, 
              "isoforms"=transcript.exp, 
              "isoforms.counts"=transcript.count))
}

    

    rnaseq.sampleinfo <- read.csv(file="/pfs/downAnnotations/ccle_rnaseq_meta.csv", stringsAsFactors=FALSE)
    rownames(rnaseq.sampleinfo) <- rnaseq.sampleinfo$Run
    rnaseq.sampleinfo[ , "cellid"] <- as.character(matchToIDTable(ids=rnaseq.sampleinfo[ , "Cell_Line"], tbl=curationCell, column = "CCLE.cellid", returnColumn = "unique.cellid"))
   
    for (r in 1:length(tool_path)){
  print(tool_path[r])
  if (length(grep(pattern = 'Kallisto', x = tool_path[r])) > 0){
    tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
    tdir = paste0("ccle_rnaseq_",gsub(".","_",tolower(tool), fixed = T), "/",  tool, "/", tool, "/")  
    rnatool="kallisto"	  
  } else {
    tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
    tdir = paste0("ccle_rnaseq_",gsub(".","_",tolower(tool), fixed = T), "/",  tool, "/", tool, "/")
    rnatool="salmon"	  
  }
  
  
  if (length(grep(pattern = 'lift37', x = tool_path[r])) > 0){
    annot = "/pfs/downAnnotations/Gencode.v33lift37.annotation.RData"
  } else if (length(grep(pattern = 'v33', x = tool_path[r])) > 0){
    annot = "/pfs/downAnnotations/Gencode.v33.annotation.RData"
  } else {
    annot = "/pfs/downAnnotations/Ensembl.v99.annotation.RData"
  }
  print(annot)
  
  print(tdir)
  print(file.path(paste0(myDirPrefix, tdir, tool_path[r])))
  rnaseq <- summarizeRnaSeq(dir=file.path(paste0(myDirPrefix, tdir, tool_path[r])),
                            features_annotation=annot,
                            samples_annotation=rnaseq.sampleinfo,
			    method = rnatool)
  rnaseq_results <- c(rnaseq_results,c(
    rnaseq <- setNames(rnaseq,  paste0(tool,".", names(rnaseq)))
  )
  )
}
    
print("finished rnaseq summarization")

rnaseq_cellid_all <- pData(rnaseq_results[[1]])[,"cellid"]
    
    #microarray (rna)
    load("/pfs/processCCLE_Brain-Array/ccle_ge_brainarray_rma.RData")
    rna <- eset
    rna$cellid <- matchToIDTable(ids=rna$Cell.line.primary.name, tbl=curationCell, column = "CCLE.cellid", returnColumn = "unique.cellid")
 
    
    #cnv
    ttt <- readRDS("/pfs/download_ccle_molec/CCLE_molecular/CCLE_molecular/2015/CNV/CCLE_CN.gene.RDS")
    ccle.eset <- ttt
    ccle.eset <- ccle.eset[ , which(!is.na(pData(ccle.eset)[, "Cell line primary name"]))]
    tt <- rownames(pData(ccle.eset))
    pData(ccle.eset) <- as.data.frame(apply(pData(ccle.eset), MARGIN=2, as.character), stringsAsFactors=FALSE)
    rownames(pData(ccle.eset)) <- tt
    pData(ccle.eset)[ , "cellid"] <- as.character(matchToIDTable(ids=pData(ccle.eset)[, "Cell line primary name"], tbl=cell_all, column = "CCLE.cellid", returnColumn = "unique.cellid"))
    pData(ccle.eset)[,"batchid"] <- NA
    annot <- read.csv("/pfs/downAnnotations/annot_ensembl_all_genes.csv", stringsAsFactors=FALSE, check.names=FALSE, header=TRUE, row.names=1)
    tt <- annot[match(rownames(fData(ccle.eset)), annot$gene_name), c("gene_id", "EntrezGene.ID", "gene_name", "gene_biotype")]
    rownames(tt) <- rownames(fData(ccle.eset))
    colnames(tt) <- c("EnsemblGeneId", "EntrezGeneId", "Symbol", "GeneBioType")
    fData(ccle.eset) <- tt
    annotation(ccle.eset) <- "cnv"
    tt <- exprs(ccle.eset)  
    temp <- rownames(tt)
    tt <- apply(tt, MARGIN=2, as.numeric)
    rownames(tt) <- temp
    exprs(ccle.eset) <- tt
    data <- exprs(ccle.eset)
    tt <- apply(data, MARGIN = 1, function(x){table(is.na(x))["FALSE"]})
    features <- setdiff(rownames(exprs(ccle.eset)), names(tt)[which(is.na(tt))])
    ccle.eset <- ccle.eset[features, ]
    
    #mutation
    
    ## from hybrid capture
  mut <- read.csv("/pfs/download_ccle_molec/CCLE_molecular/CCLE_molecular/2015/Mutation/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf", sep="\t")
  mut <- mut[ , c("Hugo_Symbol", "Tumor_Sample_Barcode", "Protein_Change"), drop=FALSE]
  mut[!is.na(mut) & mut == ""] <- NA
  mut[is.na(mut[ , "Protein_Change"]) | mut[ , "Protein_Change"] == "", "Protein_Change"] <- "wt"
  mut[!is.na(mut[ , "Protein_Change"]) & (mut[ , "Protein_Change"] == "p.?" | mut[ , "Protein_Change"] == "p.0?"), "Protein_Change"] <- NA
  mut <- mut[complete.cases(mut), , drop=FALSE]
  myx <- !duplicated(paste(mut[ , c("Tumor_Sample_Barcode")], mut[ , c("Hugo_Symbol")], mut[ , c("Protein_Change")], sep="///"))
  mut <- mut[myx, , drop=FALSE]
  ## from oncomap
  mut2 <- read.csv("/pfs/download_ccle_molec/CCLE_molecular/CCLE_molecular/2015/Mutation/CCLE_Oncomap3_2012-04-09.maf", sep="\t")
  mut2 <- mut2[ , c("Hugo_Symbol", "Tumor_Sample_Barcode", "Protein_Change"), drop=FALSE]
  mut2[!is.na(mut2) & mut2 == ""] <- NA
  # mut2[is.na(mut2[ , "Protein_Change"]) | mut2[ , "Protein_Change"] == "", "Protein_Change"] <- "wt"
  mut2[!is.na(mut2[ , "Protein_Change"]) & (mut2[ , "Protein_Change"] == "p.?" | mut2[ , "Protein_Change"] == "p.0?"), "Protein_Change"] <- NA  
  mut2 <- mut2[complete.cases(mut2), , drop=FALSE]
  myx <- !duplicated(paste(mut2[ , c("Tumor_Sample_Barcode")], mut2[ , c("Hugo_Symbol")], mut2[ , c("Protein_Change")], sep="///"))
  mut2 <- mut2[myx, , drop=FALSE]
  ## merge mutations discovered by hybrid capture and oncomap
  mutation <- rbind(mut, mut2)
  ucell <- sort(unique(mutation[ , "Tumor_Sample_Barcode"]))
  ugene <- sort(unique(mutation[ , "Hugo_Symbol"]))
  
  dd <- matrix("wt", nrow=length(ucell), ncol=length(ugene), dimnames=list(ucell, ugene))
  ## maximum number of mutations per gene
  # mutmax <- max(table(paste(mutation[ , c("Tumor_Sample_Barcode")], mutation[ , c("Hugo_Symbol")], sep="::")))
  mm <- 1:nrow(mutation)
  ff <- TRUE
  while(length(mm) > 1) {
    myx <- !duplicated(paste(mutation[mm, c("Tumor_Sample_Barcode")], mutation[mm, c("Hugo_Symbol")], sep="///"))
    if(ff) {
      dd[as.matrix(mutation[mm[myx], c("Tumor_Sample_Barcode", "Hugo_Symbol")])] <- mutation[mm[myx], "Protein_Change"]
      ff <- FALSE
    } else {
      dd[as.matrix(mutation[mm[myx], c("Tumor_Sample_Barcode", "Hugo_Symbol")])] <- paste(dd[as.matrix(mutation[mm[myx], c("Tumor_Sample_Barcode", "Hugo_Symbol")])], mutation[mm[myx], "Protein_Change"], sep="///")
    }
    mm <- mm[!myx]
  }
  ## check for inconsistencies (wt + mutations)
  iix <- grep("///", dd)
  for(iii in iix) {
    x <- sort(unique(unlist(strsplit(dd[iii], split="///"))))
    if(length(x) > 1) { x <- x[!is.element(x, "wt")] }
    dd[iii] <- paste(x, collapse="///")
  }
  nn <- sampleinfo[match(rownames(dd), sampleinfo[ , "CCLE.name"]), "cellid"]
  ## remove if we do not have cell line identifier
  dd <- dd[!is.na(nn), , drop=FALSE]
  #nn <- as.character(matchToIDTable(ids=nn, tbl=curationCell, column = "CCLE.cellid", returnColumn = "unique.cellid"))
  rownames(dd) <- nn[!is.na(nn)]
  mutation <- dd
  
  save(drugpheno, file="/pfs/out/drug_pheno.RData")
    
  cellnall <- sort(unique(c(as.character(sampleinfo[ , "cellid"]), as.character(drugpheno[ , "cellid"]), rownames(mutation))))
  dd <- matrix(NA, ncol=ncol(mutation), nrow=length(cellnall), dimnames=list(cellnall, colnames(mutation)))
  dd[rownames(mutation), colnames(mutation)] <- mutation
  mutation <- dd
  #rownames(mutation) <- as.character(matchToIDTable(ids=rownames(mutation), tbl=curationCell, column = "CCLE.cellid", returnColumn = "unique.cellid"))
  #mutation <- mutation[unique(rownames(mutation)),]
  
  MutationEset <- ExpressionSet(t(mutation))
  geneMap <- read.csv("/pfs/downAnnotations/annot_ensembl_all_genes.csv", as.is=TRUE)
  geneInfoM <- geneMap[na.omit(match(rownames(MutationEset),geneMap[ , "gene_name"]) ), c('gene_id', 'gene_biotype', 'gene_name', 'EntrezGene.ID')]
  rownames(geneInfoM) <- geneInfoM[ , "gene_name"]
  geneInfoM <- geneInfoM[rownames(MutationEset),]
  rownames(geneInfoM) <- rownames(MutationEset)
  colnames(geneInfoM) <- c("EnsemblGeneId", "GeneBioType", "Symbol", "EntrezGeneId")
  fData(MutationEset) <- geneInfoM
  tttt <- data.frame(row.names=colnames(MutationEset), colnames(MutationEset))
  colnames(tttt) <- 'cellid'
  pData(MutationEset) <- tttt
  annotation(MutationEset) <- "mutation"
  pData(MutationEset)[ , "batchid"] <- NA
  pData(MutationEset)[ , "cellid"] <- as.character(matchToIDTable(ids=pData(MutationEset)[ , "cellid"], tbl=curationCell, column="CCLE.cellid", returnColumn = "unique.cellid"))
    
  #add missing celllines and drugs to cell/drug info
    
cellnall <- unionList(rownames(celline.ccle), 
			ccle.eset$cellid, 
		        rna$cellid, 
		        rnaseq_cellid_all,
		        MutationEset$cellid)
    
newcells <- setdiff(cellnall, rownames(celline.ccle))
newRows <- matrix(NA_character_, nrow=length(newcells), ncol=ncol(celline.ccle))
# newRows <- cell.info[newcells,]

rownames(newRows) <- newcells
colnames(newRows) <- colnames(celline.ccle)
newRows[,"unique.cellid"] <- newcells

celline.ccle <- rbind(celline.ccle, newRows)

cellsPresent <- sort(unionList(sensitivityInfo$cellid, 
				rna$cellid, 
				MutationEset$cellid,
				ccle.eset$cellid,
		    		rnaseq_cellid_all))
    
celline.ccle <- celline.ccle[cellsPresent,]

celline.ccle$tissueid <- curationTissue[rownames(celline.ccle), "unique.tissueid"]
celline.ccle$cellid <- rownames(celline.ccle)

curationTissue <- curationTissue[rownames(celline.ccle),]
curationCell <- curationCell[rownames(celline.ccle),]

#add cellosaurus disease type to cell-info
disease <- cell_all$Cellosaurus.Disease.Type[match(celline.ccle$cellid, cell_all$unique.cellid)]
celline.ccle$Cellosaurus.Disease.Type <- disease
	  

#add cellosaurus assession to cell-info
assession <- cell_all$Cellosaurus.Accession.id[match(celline.ccle$cellid, cell_all$unique.cellid)]
celline.ccle$Cellosaurus.Accession.id <- assession
		 
#add pharmacodb id to cell-info
pdb <- cell_all$PharmacoDB.id[match(celline.ccle$cellid, cell_all$unique.cellid)]
celline.ccle$PharmacoDB.id <- pdb

#add study tissue id to cell_info
study_tissue <- cell_all$unique.tissueid.fromstudies[match(celline.ccle$cellid, cell_all$unique.cellid)]
celline.ccle$unique.tissueid.fromstudies <- study_tissue
		 
#add study cell-line type to cell_info
cell_type <- cell_all$CellLine.Type[match(celline.ccle$cellid, cell_all$unique.cellid)]
celline.ccle$CellLine.Type <- cell_type
		 
#add metastatic info to cell_info		 
metastatic <- cell_all$Metastatic[match(celline.ccle$cellid, cell_all$unique.cellid)]
celline.ccle$Metastatic <- metastatic


drugsPresent <- sort(unique(sensitivityInfo$drugid))

druginfo <- druginfo[drugsPresent,]
	  
	  
drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))
drug_all <- drug_all[which(!is.na(drug_all[ , "CCLE.drugid"])),]
drug_all <- drug_all[ , c("unique.drugid", "CCLE.drugid","smiles","inchikey","cid","FDA")]
rownames(drug_all) <- drug_all[ , "unique.drugid"]

drug_all <- drug_all[rownames(druginfo),]
druginfo[,c("smiles","inchikey","cid","FDA")] <- drug_all[,c("smiles","inchikey","cid","FDA")]
    
curationDrug <- curationDrug[rownames(druginfo),]
druginfo$drugid <- rownames(druginfo)

standardizeRawDataConcRange <- function(sens.info, sens.raw){
	unq.drugs <- unique(sens.info$drugid)

	conc.m <- data.table(melt(sens.raw[,,1], as.is=TRUE))
	conc.m[,drugid := sens.info$drugid[match(Var1, rownames(sens.info))]]
	conc.ranges <- conc.m[,.(l = min(value, na.rm=T), r = max(value, na.rm=T)), c("drugid", "Var1")]
	conc.ranges[,Var1 := NULL]
	conc.ranges <- conc.ranges[,unique(.SD), drugid]	
	# conc.ranges[,N := .N, drugid]
	conc.ranges.disj <- conc.ranges[, {sq <- sort(unique(c(l,r))); 
				                       l = sq[seq(1,length(sq)-1)];
				                       r = sq[seq(2,length(sq))];
				                       .(l=l,r=r)}, drugid]
    ## Function below returns all consecutive ranges of ints between 1 and N
    returnConsInts <- function(N) {
        stopifnot(N>0)
        unlist(sapply(seq(1,N), function(ii) return(sapply(seq(ii, N), function(jj) return(seq(ii,jj))))), recursive=FALSE)
    }
    rangeNoHoles <- function(indicies, lr.tbl){
        if(length(indicies) == 1) return(TRUE)
        sq <- seq(indicies[1], indicies[length(indicies)]-1)
        all(lr.tbl[["l"]][sq+1] <= lr.tbl[["r"]][sq])
    }
    per.drug.range.indicies <- sapply(conc.ranges.disj[,.N,drugid][,N], returnConsInts)

    names(per.drug.range.indicies) <- conc.ranges.disj[,unique(drugid)] ## checked this: conc.ranges.disj[,.N,drugid][,drugid] == conc.ranges.disj[,unique(drugid)]
    

    # Check if there are any holes in the chosen range combination
    per.drug.range.indicies <- sapply(names(per.drug.range.indicies), function(drug){

        lr.tbl <- conc.ranges.disj[drugid == drug]
        per.drug.range.indicies[[drug]][sapply(per.drug.range.indicies[[drug]], rangeNoHoles, lr.tbl = lr.tbl)]

        })
    per.drug.range.indicies.2 <- sapply(names(per.drug.range.indicies), function(drug){

        lr.tbl <- conc.ranges.disj[drugid == drug]
        res <- t(sapply(per.drug.range.indicies[[drug]], function(x) return(c(lr.tbl[x[1],l], lr.tbl[x[length(x)],r]))))
        colnames(res) <- c("l", "r")
        res <- data.frame(res)
        res <- cbind(drugid = drug, res)
        }, simplify=FALSE)
    per.drug.range.indicies.dt <- rbindlist(per.drug.range.indicies.2)
    
    conc.ranges <- conc.m[,.(l = min(value, na.rm=T), r = max(value, na.rm=T)), c("drugid", "Var1")]
    setkey(conc.m, Var1)
    conc.m <- na.omit(conc.m)
    setkey(conc.m, drugid, Var1, value)
    setkey(conc.ranges, drugid, l, r)
    # tic()
    ## NOTE:: Data.table used for maximum speed. Probably possible to do this more intelligently by 
    ## NOTE:: being aware of which conditions overlap, but its fast enough right now as it is.
    chosen.drug.ranges <- lapply(unq.drugs, function(drug){
        num.points.in.range <- apply(per.drug.range.indicies.dt[drugid==drug, .(l,r)], 1, function(rng){
            conc.m[drugid==drug][conc.ranges[drugid==drug][l<=rng["l"]][r>=rng["r"],Var1], on="Var1"][value >= rng["l"]][value <= rng["r"],.N]
            # conc.m[drugid==drug][, Var1]
            })
        max.ranges <- per.drug.range.indicies.dt[drugid==drug][which(num.points.in.range==max(num.points.in.range))]
        max.ranges[which.max(log10(r) - log10(l)), ]
    })
    # toc()
    names(chosen.drug.ranges) <- sapply(chosen.drug.ranges, `[[`, "drugid")
    removed.experiments <- unlist(lapply(unq.drugs, function(drug){
        rng <- unlist(chosen.drug.ranges[[drug]][,.(l,r)])
        exp.out.range <- conc.ranges[drugid==drug][l>rng["l"] | r<rng["r"],Var1]
        return(exp.out.range)
        }))

    sens.raw[removed.experiments,,] <- NA_real_
    conc.ranges.kept <- conc.ranges[!Var1 %in% removed.experiments]

    for(drug in unq.drugs){
        rng <- unlist(chosen.drug.ranges[[drug]][,.(l,r)])
        myx <- conc.ranges.kept[drugid==drug,Var1]
        doses <- sens.raw[myx, ,"Dose"]
        which.remove <- (doses < rng["l"] | doses > rng["r"])
        sens.raw[myx, ,"Dose"][which(which.remove,arr.ind=TRUE)] <- NA_real_
        sens.raw[myx, ,"Viability"][which(which.remove,arr.ind=TRUE)] <- NA_real_

        ## Annotate sens info with chosen range
        sens.info[sens.info$drugid==drug,"chosen.min.range"] <- rng["l"]
        sens.info[sens.info$drugid==drug,"chosen.max.range"] <- rng["r"]
    }
    sens.info$rm.by.conc.range <- FALSE
    sens.info[removed.experiments,"rm.by.conc.range"] <- TRUE

    return(list("sens.info" = sens.info, sens.raw = sens.raw))
}
		 

    
    
  z <- list()
  
  z <- c(z,c(
    "rna"=rna,
    rnaseq_results,
    "mutation"=MutationEset,
    "cnv"=ccle.eset)
  )
    
sensitivityInfo <- as.data.frame(sensitivityInfo)
standardize <- standardizeRawDataConcRange(sens.info = sensitivityInfo, sens.raw = raw.sensitivity)
       
    
CCLE <- PharmacoSet(molecularProfiles=z, name="CCLE", cell=celline.ccle, drug=druginfo, sensitivityInfo=standardize$sens.info, sensitivityRaw=standardize$sens.raw, sensitivityProfiles=profiles, sensitivityN=NULL,  curationCell=curationCell, curationDrug=curationDrug, curationTissue=curationTissue, datasetType="sensitivity")

filterNoisyCurves2 <- function(pSet, epsilon=25 , positive.cutoff.percent=.80, mean.viablity=200, nthread=1) {
acceptable <- mclapply(rownames(sensitivityInfo(pSet)), function(xp) {
  #for(xp in rownames(sensitivityInfo(pSet))){
  drug.responses <- as.data.frame(apply(pSet@sensitivity$raw[xp , ,], 2, as.numeric), stringsAsFactors=FALSE)
  if (!all(is.na(drug.responses))){
    
  
  drug.responses <- drug.responses[complete.cases(drug.responses), ]
  doses.no <- nrow(drug.responses)
  drug.responses[,"delta"] <- .computeDelta(drug.responses$Viability)
  
  delta.sum <- sum(drug.responses$delta, na.rm = TRUE)
  
  max.cum.sum <- .computeCumSumDelta(drug.responses$Viability)
  
  if ((table(drug.responses$delta < epsilon)["TRUE"] >= (doses.no * positive.cutoff.percent)) &
      (delta.sum < epsilon) &
      (max.cum.sum < (2 * epsilon)) &
      (mean(drug.responses$Viability) < mean.viablity)) {
    return (xp)
  }
  }
  
}, mc.cores=nthread)
acceptable <- unlist(acceptable)
noisy <- setdiff(rownames(sensitivityInfo(pSet)), acceptable)
return(list("noisy"=noisy, "ok"=acceptable))
}

.computeDelta <- function(xx ,trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc)
  {
    return(c(pmin(100, xx[2:length(xx)]) - pmin(100, xx[1:length(xx)-1]), 0))
  }else{
    return(c(xx[2:length(xx)] - xx[1:length(xx)-1]), 0)
  }
}

#' @importFrom utils combn
.computeCumSumDelta <- function(xx, trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc) {
    xx <- pmin(xx, 100)
  }
  tt <- t(combn(1:length(xx), 2 , simplify = TRUE))
  tt <- tt[which(((tt[,2] - tt[,1]) >= 2) == TRUE),]
  cum.sum <- unlist(lapply(1:nrow(tt), function(x){xx[tt[x,2]]-xx[tt[x,1]]}))
  return(max(cum.sum))
}
		 
noisy_out <- filterNoisyCurves2(CCLE)
print("filter done")
CCLE@sensitivity$profiles[noisy_out$noisy, ] <- NA		 
saveRDS(CCLE, file="/pfs/out/CCLE.rds")		 
		 

