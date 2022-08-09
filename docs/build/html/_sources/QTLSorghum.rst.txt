================================================================================
A PBGL Snakemake Workflow and Quatitative Trait Locus on Bulk Segregant Analysis 
================================================================================

=====================
:Author: Michael Hall
:Date: 08/09/2022
====================



Software Prerequisites in an Environment .YAML Form
===================================================

.. code:: shell

	name: QTLseqr

	channels:
  	- default
  	- bioconda
  	- conda-forge
  	- r

	dependencies:
  	- r-modeest >= 2.3.2
  	- r-ggplot2 >= 2.2.0
  	- r-gtools
  	- r-dplyr
  	- r-readr
  	- r-tidyr
  	- r-Rcpp
  	- r-locfit
  	- r-knitr
  	- r-rmarkdown
  	- r-kableExtra
  	- r-devtools
  	- r-data.table
  	- r-vcfR
  	- r-optparse
  	- r-base >= 4.1.0
  
Install manually in RStudio
----------------------------

.. code:: shell

	#install the dependencies
	install.packages(c("data.table", 
                   "modeest",
                   "locfit",
                   "dplyr", 
                   "tidyr", 
                   "vcfR", 
                   "ggplot2"), dependencies=TRUE)

	# install devtools
	install.packages("devtools", dependencies=TRUE)

	# use devtools to install QTLseqr from github
	library("devtools")
	#devtools::install_github("bmansfeld/QTLseqr")
	#devtools::install_github("warthmann/QTLseqr")
	devtools::install_github("pbgl/QTLseqr")



	#load dependencies
	library("devtools")
	library("data.table")
	library("dplyr")
	library("tidyr")
	library("vcfR")
	library("ggplot2")

	#load the package
	library("QTLseqr")


Data Analysis
-------------

.. code:: shell
   
	QTLseqr::importFromVCF(file = "freebayes_D2.filtered.vcf",highBulk = "D2_F2_tt",lowBulk = "D2_F2_TT",
	chromList = c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"),
	,filter=FALSE, outfile=TRUE)
	#Set High bulk and Low bulk sample names and parser generated file name
	#The file name is generated from the QTLParser_1_MH function in line 119

	HighBulk <- "D2_F2_tt"
	LowBulk <- "D2_F2_TT"
	file <- "Hall.csv"

#Choose which chromosomes/contigs will be included in the analysis,

	Chroms <- c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")

	df <-
   QTLseqr::importFromTable(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = chromList
    sep = ","
   ) 
  
	#Filter SNPs based on some criteria
	df_filt <-
    QTLseqr::filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.20,
    minTotalDepth = 100,
    maxTotalDepth = 400,
    minSampleDepth = 40,
             minGQ = 99,
    verbose = TRUE
   )
  
	#Run G' analysis
	df_filt<-QTLseqr::runGprimeAnalysis(
    SNPset = df_filt,
    windowSize = 5000000,
    outlierFilter = "deltaSNP",
    filterThreshold = 0.1)
  
	#Run QTLseq analysis
	df_filt <- QTLseqr::runQTLseqAnalysis(
    SNPset = df_filt,
    windowSize = 5000000,
    popStruc = "F2",
    bulkSize = c(45, 38),
    replications = 10000,
    intervals = c(95, 99)
	)

	#Plot
	QTLseqr::plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
	QTLseqr::plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals  = TRUE)

	#export summary CSV
	QTLseqr::getQTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")
