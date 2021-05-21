## --------------------------------------------------------------------------
##
## This file is part of the miRNA-QC-and-Diagnosis software package.
##
## Version 1.1.2 - May 2021
##
##
## The miRNA-QC-and-Diagnosis package is free software; you can use it,
## redistribute it, and/or modify it under the terms of the GNU General
## Public License version 3 as published by the Free Software Foundation.
## The full text of the license can be found in the file LICENSE.txt at the top
## level of the package distribution.
##
## Authors:
##	Michele Castelluzzo (1), Alessio Perinelli (2), Simone Detassis (3),
##	Michela A. Denti (3) and Leonardo Ricci (1,2)
##	(1) Department of Physics, University of Trento, 38123 Trento, Italy
##	(2) CIMeC, Center for Mind/Brain Sciences, University of Trento,
##		38068 Rovereto, Italy
##	(3) Department of Cellular, Computational and Integrative Biology
##		(CIBIO), University of Trento, 38123 Trento, Italy
##
##	michele.castelluzzo@unitn.it
##	alessio.perinelli@unitn.it
##	michela.denti@unitn.it
##	leonardo.ricci@unitn.it
##	https://github.com/LeonardoRicci/
##	https://nse.physics.unitn.it/
##
##
## If you use the miRNA-QC-and-Diagnosis package for your analyses, please cite:
##
##	L. Ricci, V. Del Vescovo, C. Cantaloni, M. Grasso, M. Barbareschi and
##	M. A. Denti, Statistical analysis of a Bayesian classifier based on the
##	expression of miRNAs, BMC Bioinformatics 16:287 (2015).
##	DOI: 10.1186/s12859-015-0715-9
##
##
## --------------------------------------------------------------------------

#' Pre-processing of datasets.
#'
#' This function carries out the pre-processing required by the other functions of the miRNA-QC-and-Diagnosis package.
#'
#' @param inputDataset Dataset (data frame) to be pre-processed. The data frame must contain the columns 'Subject', 'miRNA', 'Value' and possibly 'Class'. Any other column is ignored, and any missing column forbids execution. Please note that using the character '-' within the dataset causes undefined behaviour (even if data were correctly loaded by 'read.table').
#' @param multipletSize Size of the multiplets to be considered. Any multiplet of different size is ignored.
#'
#' @return A pre-processed data frame, containing the columns 'Subject', 'miRNA', 'Mean', 'StdDev', 'SampleSize', and possibly 'Class'.
#'
#' Please refer to the user manual installed in "/path-to-library/MiRNAQCD/doc/manual.pdf" for detailed function documentation. The path "/path-to-library" can be shown from R by calling ".libPaths()"
#'
#' @examples
#' requiredFile = paste(system.file(package="MiRNAQCD"), "/extdata/test_dataset_alpha.dat", sep='')
#' myDataFrame <- read.table(file=requiredFile, header=TRUE)
#' myPreprocessedDataFrame <- miRNA_expressionPreprocessing(myDataFrame, 3)
#'

#' @export
miRNA_expressionPreprocessing <- function(inputDataset, multipletSize) {

	if (!(("miRNA" %in% colnames(inputDataset)) & ("Subject" %in% colnames(inputDataset)) & ("Value" %in% colnames(inputDataset))))  {
		stop("ERROR: unsuitable dataset format. Dataset must contain at least columns 'miRNA', 'Subject', 'Value'.\n")
	}

	if (length(inputDataset[1,]) > 4) {
		warning("WARNING: more than 4 dataset columns. Columns other than 'miRNA', 'Subject', 'Value', 'Class' will be ignored.\n")
	}

	if ("Class" %in% colnames(inputDataset)) {
		selectedDataset <- subset(inputDataset, select=c("Subject", "miRNA", "Value", "Class"))
		columnClassExists <- TRUE
	} else {
		selectedDataset <- subset(inputDataset, select=c("Subject", "miRNA", "Value"))
		columnClassExists <- FALSE
	}

	selectedDataset$Value <- as.numeric(as.vector(selectedDataset$Value))
	selectedDataset <- selectedDataset[is.numeric(selectedDataset$Value) & !is.na(selectedDataset$Value), ]

	setOfSubjects <- unique(as.character(selectedDataset$Subject))

	for (subject in setOfSubjects) {
		singleSubjectFrame <- subset(selectedDataset, selectedDataset$Subject == subject)
		if (length(unique(singleSubjectFrame$Class)) > 1) {
			stop("ERROR: multiple classes for some of the subjects. Each subject's entry must report the same class.\n")
		}
	}

	setOfFeatures <- unique(selectedDataset$miRNA)
	setOfFeatures <- setOfFeatures[order(setOfFeatures)]

	for (feature in setOfFeatures) {
		singleFeatureFrame <- subset(selectedDataset, selectedDataset$miRNA == feature)

		columnSubjectFeatures <- cbind(sort(unique(as.character(singleFeatureFrame$Subject))), rep(feature, length(unique(singleFeatureFrame$Subject))))
		columnValueMean <- as.vector(tapply(singleFeatureFrame$Value, as.character(singleFeatureFrame$Subject), mean))
		columnValueStd <- as.vector(tapply(singleFeatureFrame$Value, as.character(singleFeatureFrame$Subject), stats::sd))
		columnSampleSize <- as.vector(tapply(singleFeatureFrame$Value, as.character(singleFeatureFrame$Subject), length))

		if (columnClassExists) {
			columnClass <- as.vector(tapply(as.character(singleFeatureFrame$Class), as.character(singleFeatureFrame$Subject), unique))
			tempRows <- cbind(columnSubjectFeatures, columnValueMean, columnValueStd, columnSampleSize, columnClass)
		} else {
			tempRows <- cbind(columnSubjectFeatures, columnValueMean, columnValueStd, columnSampleSize)
		}

		if (!exists("tempFrame"))
			tempFrame <- tempRows
		else
			tempFrame <- rbind (tempFrame, tempRows)
	}
	processedDataset <- data.frame(tempFrame)

	if (columnClassExists) {
		names(processedDataset) <- c("Subject", "miRNA", "Mean", "StdDev", "SampleSize", "Class")
	} else {
		names(processedDataset) <- c("Subject", "miRNA", "Mean", "StdDev", "SampleSize")
	}

	processedDataset$Mean <- as.numeric(as.character(processedDataset$Mean))
	processedDataset$StdDev <- as.numeric(as.character(processedDataset$StdDev))
	processedDataset$SampleSize <- as.numeric(as.character(processedDataset$SampleSize))

	processedDataset$Mean <- round(processedDataset$Mean,4)
	processedDataset$StdDev <- round(processedDataset$StdDev,4)
	processedDataset$SampleSize <- round(processedDataset$SampleSize,0)

	processedDataset <- processedDataset[processedDataset$SampleSize == multipletSize, ]
	processedDataset <- data.frame(processedDataset[with(processedDataset, order(Subject)),])

	return(processedDataset)

}
