## ----pre, echo=FALSE, warning=FALSE, include=FALSE-----------------------
library(knitr)
is_check <- ("CheckExEnv" %in% search()) || any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv()))
knitr::opts_chunk$set(echo = TRUE, eval=!is_check, cache=TRUE)
#rmarkdown::render("vignettes/TcGSA_userguide.Rmd")

## ----GS_import, message=FALSE, cache=TRUE--------------------------------
temp <- tempfile()
utils::download.file("http://doi.org/10.1371/journal.pcbi.1004310.s007", destfile = temp, mode = "wb")
load(unz(temp, "ReproducibleRFiles/GMTs_PLOScb.RData", open = "r"))
unlink(temp)
rm(temp)

## ----GEOquery, include=FALSE, message=FALSE, cache=FALSE-----------------
if (!requireNamespace("GEOquery", quietly = TRUE)) {
	if (!requireNamespace("BiocManager", quietly = TRUE)){
    	install.packages("BiocManager")
	}
	BiocManager::install("GEOquery")
}

## ----import_dalia, message=FALSE, cache=TRUE, results='hide'-------------
GEOquery::getGEOSuppFiles('GSE46734', filter_regex="(*NonParamCombat*)|(*DESIGN*)")

## ----design_dalia, cache=TRUE--------------------------------------------
design <- read.delim(gzfile("GSE46734/GSE46734_DALIA1longitudinalTranscriptome_DESIGN_anonym.txt.gz"))
design_preATI <- design[-which(design$TimePoint<0 | design$TimePoint==16 | design$TimePoint>22), ]
head(design_preATI,5)

## ---- include=FALSE------------------------------------------------------
stopifnot(nrow(design_preATI)==90 & ncol(design_preATI)==6)

## ----expr_dalia, cache=TRUE----------------------------------------------
expr_preATI <- read.delim(gzfile("GSE46734/GSE46734_DALIA1longitudinalTranscriptome_PALO01_PreATI_NEQC_NonParamCombat.txt.gz"))
rownames(expr_preATI) <- expr_preATI$PROBE_ID
expr_preATI <- expr_preATI[,as.character(design_preATI$Sample_name)]

expr_preATI[1:4,1:4]

## ---- include=FALSE------------------------------------------------------
stopifnot(nrow(expr_preATI)==32978 & ncol(expr_preATI)==90)

## ----dim_expr_DALIA------------------------------------------------------
identical(ncol(expr_preATI), nrow(design_preATI))

## ----LR_ST, message=FALSE, warning=FALSE, cache=TRUE---------------------
tcgsa_result <- TcGSA::TcGSA.LR(expr = expr_preATI, 
								gmt = gmt_modulesV2, 
								design = design_preATI, 
								subject_name = "Patient_ID", 
								time_name = "TimePoint")

## ----tcgsa_result, cache=TRUE, echo=FALSE--------------------------------
tcgsa_result

## ----summary_dalia, cache=TRUE-------------------------------------------
summary(tcgsa_result)

## ----signifLRT_ST, cache=TRUE--------------------------------------------
head(TcGSA::signifLRT.TcGSA(tcgsa_result)$mixedLRTadjRes)

## ----multtest_ST, cache=TRUE---------------------------------------------
head(TcGSA::multtest.TcGSA(tcgsa_result))

## ----plot1GS_ST, message=FALSE, warning=FALSE, fig.keep='all', cache=TRUE----
TcGSA::plot1GS(expr = expr_preATI, 
			   #plot1GS(expr = tcgsa_result$Estimations,
			   gmt = gmt_modulesV2, 
			   Subject_ID = design_preATI$Patient_ID, 
			   TimePoint = design_preATI$TimePoint,
			   clustering = FALSE, 
			   time_unit = "W", 
			   geneset.name = "M3.2", 
			   title="",
			   margins=0.4,
			   lab.cex=0.37,
			   axis.cex=0.37,
			   line.size=0.45,
			   gg.add=list(ggplot2::theme(legend.position="none"),
			   			ggplot2::ylim(-1.26,1.26)
			   ))

## ----import_ober, message=FALSE, warning=FALSE, cache=TRUE---------------
utils::download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30101/soft/GSE30101_family.soft.gz", destfile = "GSE30101_family.soft.gz", mode = "wb", cacheOK = FALSE)
gse.soft <- GEOquery::getGEO(filename="GSE30101_family.soft.gz")

## ----expr_ober, cache=TRUE-----------------------------------------------
probesIDs <- GEOquery::Table(GEOquery::GSMList(gse.soft)[[1]])$ID
data.matrix <- do.call('cbind', lapply(GEOquery::GSMList(gse.soft),
									   function(x) {
									   	tab <- GEOquery::Table(x)
									    mymatch <- match(probesIDs,tab$ID_REF)
									    return(tab$VALUE[mymatch])
									   }))
rownames(data.matrix) <- probesIDs
expr.All.ChaussVac <- apply(X = data.matrix, MARGIN = 2, FUN = as.numeric)
rownames(expr.All.ChaussVac) <- probesIDs

## ----design_ober, cache=TRUE---------------------------------------------
design_list <- lapply(GEOquery::GSMList(gse.soft), 
					  function(x){GEOquery::Meta(x)$characteristics_ch1})
design <- data.frame(row.names = names(design_list))
design$sample_ID <- names(design_list)
s_id <- unlist(lapply(design_list, function(x){gsub("subject id: ", "", x[grep("subject id: ", x)])}))
design$Subject_ID <- as.character(paste("P", s_id[design$sample_ID], sep=""))

time <- unlist(lapply(design_list, function(x){gsub("day: ", "", x[grep("day: ", x)])}))
time[which(time %in% c("-7", "0.5", "1", "7", "10", "14",  "21", "28"))] <-
	paste("D", time[which(time %in% c("-7", "0.5", "1", "7", "10", "14", "21", "28"))], sep="")
time[which(time %in% c("-168", "1.5", "6", "9", "12", "15", "24", "36", "48"))] <-
	paste("H", time[which(time %in% c("-168", "1.5", "6", "9", "12", "15", "24", "36", "48"))], sep="")
design$Time <- as.character(time[design$sample_ID])

vac <- unlist(lapply(design_list, function(x){
	gsub("vaccine: ", "", x[grep("vaccine: ", x)])
}))
vac <- as.factor(vac)
levels(vac) <- c("influenza", "influenza", "influenza", "influenza", "saline", 
				 "pneumo", "pneumo", "pneumo", "saline", "saline")
design$Vaccine <- as.character(vac[design$sample_ID])

sampSet <- unlist(lapply(design_list, function(x){
	gsub("sample set: ", "", x[grep("sample set: ", x)])
}))
design$sampSet <- as.character(sampSet[design$sample_ID])

design$Time[which(design$sampSet=="Training_Set_Vein" & design$Time %in% c("0", "3"))] <-
	paste("D", design$Time[which(design$sampSet=="Training_Set_Vein" & design$Time %in% c("0", "3"))], sep="")
design$Time[which(design$sampSet=="Training_Set_Finger" & design$Time %in% c("0", "3"))] <-
	paste("H", design$Time[which(design$sampSet=="Training_Set_Finger" & design$Time %in% c("0", "3"))], sep="")
design$Time[which(design$sampSet=="Test_Set_Vein" & design$Time %in% c("0", "3"))] <-
	paste("D", design$Time[which(design$sampSet=="Test_Set_Vein" & design$Time %in% c("0", "3"))], sep="")
design$Time[which(design$sampSet=="Test_Set_Finger" & design$Time %in% c("0", "3"))] <-
	paste("D", design$Time[which(design$sampSet=="Test_Set_Finger" & design$Time %in% c("0", "3"))], sep="")
design$Time[which(design$sampSet=="Validation_Vein" & design$Time %in% c("0", "3"))] <-
	paste("D", design$Time[which(design$sampSet=="Validation_Vein" & design$Time %in% c("0", "3"))], sep="")

design$Day <- gsub("D", "", design$Time)
design$Day[grep("H", design$Day)] <- as.numeric(gsub("H", "", design$Day[grep("H", design$Day)]))/24
design$Day <- as.numeric(design$Day)

design.All.ChaussVac <- design



# Avg Baseline -----
design.All.ChaussVac.trainSetVein <- design.All.ChaussVac[which(design.All.ChaussVac$sampSet=="Training_Set_Vein"),]
samplesSaline2rmv <- design.All.ChaussVac.trainSetVein[162:214,"sample_ID"]
design.All.ChaussVac.noDup <- design.All.ChaussVac.trainSetVein[-which(design.All.ChaussVac.trainSetVein$sample_ID%in%samplesSaline2rmv),]

design.All.ChaussVac.AvgBl <- design.All.ChaussVac.noDup[which(design.All.ChaussVac.noDup$Day!=0),]
design.All.ChaussVac.AvgBl[which(design.All.ChaussVac.AvgBl$Day==-7),"Day"] <- 0
design.All.ChaussVac.AvgBl[which(design.All.ChaussVac.AvgBl$Time=="D-7"),"Time"] <- "D0"

expr.All.ChaussVac.AvgBl <- expr.All.ChaussVac[, design.All.ChaussVac.AvgBl$sample_ID]
for(p in unique(design.All.ChaussVac.AvgBl$Subject_ID)){
	if(length(which(design.All.ChaussVac.noDup$Subject_ID==p & (design.All.ChaussVac.noDup$Day==0 | design.All.ChaussVac.noDup$Day==-7)))>1){
		expr.All.ChaussVac.AvgBl[, which(design.All.ChaussVac.AvgBl$Subject_ID==p & design.All.ChaussVac.AvgBl$Day==0)] <-
			apply(X=cbind(expr.All.ChaussVac[, design.All.ChaussVac.noDup[which(design.All.ChaussVac.noDup$Subject_ID==p & design.All.ChaussVac.noDup$Day==0), "sample_ID"]],
						  expr.All.ChaussVac[, design.All.ChaussVac.noDup[which(design.All.ChaussVac.noDup$Subject_ID==p & design.All.ChaussVac.noDup$Day==-7), "sample_ID"]]),
				  MARGIN=1, FUN=mean, na.rm=TRUE)
	}
}
rownames(expr.All.ChaussVac.AvgBl) <- probesIDs

if(!all.equal(as.character(design.All.ChaussVac.AvgBl$sample_ID), colnames(expr.All.ChaussVac.AvgBl))){stop("\n\n\nWARNING: EXPRESSION FILE ORDER NOT MATCHING DESIGN FILE\n\n\n")}
design.All.ChaussVac.AvgBl$Subject_ID <- as.factor(design.All.ChaussVac.AvgBl$Subject_ID)


design.PNEUMOvsSALINE.ChaussVac.AvgBl <- design.All.ChaussVac.AvgBl[which(design.All.ChaussVac.AvgBl$Vaccine!="influenza"), ]
design.PNEUMOvsSALINE.ChaussVac.AvgBl$Vaccine <- as.factor(as.character(design.PNEUMOvsSALINE.ChaussVac.AvgBl$Vaccine))

expr.PNEUMOvsSALINE.ChaussVac.AvgBl <- expr.All.ChaussVac.AvgBl[,design.PNEUMOvsSALINE.ChaussVac.AvgBl$sample_ID]

## ----LR_MT, cache=TRUE, message=FALSE, warning=FALSE---------------------
tcgsa_result_MT <- TcGSA::TcGSA.LR(expr = expr.PNEUMOvsSALINE.ChaussVac.AvgBl, 
								   gmt = gmt_modulesV2,
								   design = design.PNEUMOvsSALINE.ChaussVac.AvgBl, 
								   subject_name = "Subject_ID", 
								   time_name = "Day", 
								   group_name = "Vaccine")
summary(tcgsa_result_MT)
head(TcGSA::signifLRT.TcGSA(tcgsa_result_MT)$mixedLRTadjRes)

## ----clust_MT, message=FALSE, warning=FALSE, cache=TRUE------------------
clust <- TcGSA::clustTrend(tcgs = tcgsa_result_MT, 
						   expr=tcgsa_result_MT$Estimations,
						   Subject_ID=design_gen$Patient_ID,
						   TimePoint=design_gen$TimePoint,
						   baseline = 0,
						   group_of_interest="pneumo")
clust

## ----heatmap_MT, message=FALSE, cache=TRUE, results='asis'---------------
plot(x=tcgsa_result_MT, expr=tcgsa_result_MT$Estimations,
	 Subject_ID=design_gen$Patient_ID,
	 TimePoint=design_gen$TimePoint,
	 group_of_interest="pneumo",
	 clust_trends=clust,
	 legend.breaks=seq(from=-2,to=2, by=0.01), time_unit="D",
	 subtitle="Pneumo vs Saline", cex.label.row=0.5, cex.label.col=1, cex.main=0.7,
	 heatmap.width=0.2, dendrogram.size=0.3, margins=c(2,3),
	 heatKey.size=0.8)

## ----dl_GEOquery, warning=FALSE, message=FALSE, eval=FALSE---------------
#  if (!requireNamespace("GEOquery", quietly = TRUE)) {
#  	if (!requireNamespace("BiocManager", quietly = TRUE)){
#  		install.packages("BiocManager")
#  	}
#  	BiocManager::install("GEOquery")
#  }

