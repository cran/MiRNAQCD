## --------------------------------------------------------------------------
##
## This file is part of the miRNA-QC-and-Diagnosis software package.
##
## Version 1.0 - June 2020
##
##
## The miRNA-QC-and-Diagnosis package is free software; you can use it,
## redistribute it, and/or modify it under the terms of the GNU General
## Public License version 3 as published by the Free Software Foundation.
## The full text of the license can be found in the file LICENSE.txt at the top
## level of the package distribution.
##
## Authors:
##	Michele Castelluzzo (1), Alessio Perinelli (1), Simone Detassis (2),
##	Michela A. Denti (2) and Leonardo Ricci (1,3)
##	(1) Department of Physics, University of Trento, 38123 Trento, Italy
##	(2) Department of Cellular, Computational and Integrative Biology
##		(CIBIO), University of Trento, 38123 Trento, Italy
##	(3) CIMeC, Center for Mind/Brain Sciences, University of Trento,
##		38068 Rovereto, Italy
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

#' Classification of a dataset (diagnosis).
#'
#' This function classifies the entries of the input dataset as either target or versus by using the chosen classifier and given the corresponding disgnostic threshold value.
#'
#' @param inputDataset Dataset (data frame) to be classified. The data frame must comply with the output format of the quality control functions (miRNA_expressionPreprocessing and miRNA_removeOutliers), thus containing the columns 'Subject', 'miRNA', 'Mean', 'StdDev', 'SampleSize'. Any other column is ignored, and any missing column forbids execution.
#' @param inputMiRNAList List of miRNAs to be used by the classifier. The chosen miRNAs must be present in the 'miRNA' column of the inputDataset.
#' @param coeffList List of coefficients for the classifier. The number of coefficients must be the same as the number of used miRNAs and listed in the same order.
#' @param inputThreshold Diagnostic threshold data frame for the classifier. The data frame must comply with the output format of the classifier setup function (miRNA_classifierSetup), thus containing the columns 'Threshold', 'DeltaThreshold', 'ChiUp', 'DChiUp', 'ChiDown', 'DChiDown'. Any other column is ignored.
#' @param saveOutputFile Boolean option setting whether results are written to file (TRUE) or not (FALSE). Default is FALSE.
#' @param outputFileBasename Name of the output file where the diagnosis results are to be stored. If not assigned, a filename is automatically generated.
#' @param sep Field separator character for the output file; the default is tabulation.
#' @param plotFormat String specifying the format of generated graphic files (plots): can either be "pdf" (default) or "png".
#'
#' @return A data frame containing the columns 'Subject', 'Diagnosis' and 'Score'.
#'
#' @examples
#' requiredDataFile = paste(system.file(package="MiRNAQCD"),
#'			"/extdata/dataset_beta_clean.dat", sep='')
#' myDataFrame <- read.table(file=requiredDataFile, header=TRUE)
#' requiredThresholdFile = paste(system.file(package="MiRNAQCD"),
#'			"/extdata/dataset_alpha_threshold.txt", sep='')
#' thresholdDataFrame <- read.table(file=requiredThresholdFile, header=TRUE)
#' mirnaToUse <- c("FX", "FZ")
#' coefficientsToUse <- c(1.0, -1.0)
#'
#' ## Classification
#' classifiedDataset <- miRNA_diagnosis(myDataFrame, mirnaToUse, coefficientsToUse,
#'					thresholdDataFrame)

#' @export
miRNA_diagnosis <- function(inputDataset, inputMiRNAList, coeffList, inputThreshold, saveOutputFile=FALSE, outputFileBasename="", sep='\t', plotFormat="pdf") {

	## Input validation and preprocessing

	if (!(("Subject" %in% colnames(inputDataset)) & ("miRNA" %in% colnames(inputDataset)) & ("Mean" %in% colnames(inputDataset)) & ("StdDev" %in% colnames(inputDataset)) & ("SampleSize" %in% colnames(inputDataset))))  {
		stop("ERROR: unsuitable dataset format. Dataset must contain columns 'Subject', 'miRNA', 'Mean', 'StdDev', 'SampleSize'.\n")
	}

	if (length(inputDataset[1,]) > 6) {
		warning("WARNING: more than 6 dataset columns. Columns other than 'Subject', 'miRNA', 'Mean', 'StdDev', 'SampleSize', 'Class' will be ignored.\n")
	}

	if ("Class" %in% colnames(inputDataset)) {
		selectedDataset <- subset(inputDataset, select=c("Subject", "miRNA", "Mean", "StdDev", "SampleSize", "Class"))
	} else {
		selectedDataset <- subset(inputDataset, select=c("Subject", "miRNA", "Mean", "StdDev", "SampleSize"))
	}

	if (sep == "") {
		sep <- " "
	}

	# Input validation: check Features
	availableFeatures <- unique(selectedDataset$miRNA)
	listOfFeature <- unique(inputMiRNAList)
	if (length(listOfFeature) != length(inputMiRNAList)){
		warning("WARNING: feature list presents some duplicates which will be ignored.\n")
	}
	inputMiRNAList <- unique(inputMiRNAList)
	listOfFeature <- intersect(listOfFeature, availableFeatures)
	if (length(listOfFeature) != length(inputMiRNAList)){
		stop("ERROR: some selected features are not present in the dataset.\n")
	}
	if (length(listOfFeature) != length(coeffList)) {
		stop("ERROR: unsuitable function arguments. Size of feature list and coefficient list do not match.\n")
	}

	# Input validation: check Subjects
	listOfSubjects <- unique(selectedDataset$Subject)
	for (subject in listOfSubjects) {
		subjectFrame <- selectedDataset[selectedDataset$Subject == subject,]
		availableSubjectFeatures <- unique(subjectFrame$miRNA)
		if (length(listOfFeature) != length(intersect(availableSubjectFeatures, listOfFeature))) {
			if (!exists("subjectsToRemove")) {
				subjectsToRemove <- subject
			} else {
				subjectsToRemove <- rbind(subjectsToRemove, subject)
			}
		}
	}
	if (exists("subjectsToRemove")) {
		subjectsToRemove <- unique(subjectsToRemove)
		correctSubjects <- setdiff(listOfSubjects, subjectsToRemove)
	} else {
		correctSubjects <- listOfSubjects
	}
	if (length(correctSubjects) == 0) {
		stop("ERROR: no available subjects for this classifier.\n")
	}
	availableDataset <- selectedDataset[selectedDataset$Subject %in% correctSubjects, ]
	availableDataset <- availableDataset[availableDataset$miRNA %in% listOfFeature, ]

	# Bookkeeping
	for (feature in listOfFeature) {
		columnSubjectMean <- availableDataset[availableDataset$miRNA == feature,]
		columnSubjectMean <- subset(columnSubjectMean, select=c("Subject", "Mean"))
		if (!exists("featureFrame")) {
			featureFrame <- columnSubjectMean
		} else {
			featureFrame <- merge(featureFrame, columnSubjectMean, by = "Subject")
		}
	}
	names(featureFrame) <- c("Subject", listOfFeature)

	## Actual Diagnosis

	threshold <- as.numeric(as.vector(inputThreshold$Threshold))

	# Compute score as linear combination of features
	dataFrameTemp <- featureFrame[2:(1+length(listOfFeature))]
	for (feature in listOfFeature) {
		dataFrameTemp[,feature] <- dataFrameTemp[,feature] * as.numeric(coeffList[which(listOfFeature == feature)])
	}
	classifierDataFrame <- cbind(featureFrame, Score=rowSums(dataFrameTemp))
	classifierDataFrame <- data.frame(classifierDataFrame[with(classifierDataFrame, order(Score)),])

	# Create new Diagnosis columns and writes "target" or "versus" at each row, then join all together
	targetFrame <- classifierDataFrame[classifierDataFrame$Score > threshold, ]
	targetFrame <- cbind(targetFrame, Diagnosis=rep("target", length(targetFrame[,1])))

	versusFrame <- classifierDataFrame[classifierDataFrame$Score <= threshold, ]
	versusFrame <- cbind(versusFrame, Diagnosis=rep("versus", length(versusFrame[,1])))

	outputDataFrame <- rbind(targetFrame, versusFrame)
	outputDataFrame <- subset(outputDataFrame, select=c("Subject", "Diagnosis", "Score", listOfFeature))
	outputDataFrame <- outputDataFrame[with(outputDataFrame, order(Score, decreasing=FALSE)),]

	processedFrameToWrite <- subset(outputDataFrame, select=c("Subject", "Diagnosis", "Score"))

	## Output

	if (saveOutputFile) {

		# Write classified dataset to file
		classifierLabel <- ""
		for (i in 1:length(coeffList)) {
			classifierLabel <- paste(classifierLabel, "_", coeffList[i], "_", listOfFeature[i], sep="")
		}
		if (outputFileBasename == "") {
			outputFileBasename <- paste("diagnosis", classifierLabel, sep="")
			plotFileName <- paste("diagnosis", classifierLabel, sep="")
		} else {
			plotFileName <- outputFileBasename
		}
		outputFile <- paste(outputFileBasename, ".dat", sep="")

		if (file.exists(outputFile) & file.access(outputFile, mode=2)) {
			stop("ERROR: cannot write ", outputFile, ". Check write permission.\n")
		}

		utils::write.table(format(processedFrameToWrite, drop0trailing=FALSE), file=outputFile, sep=sep, row.names=FALSE, quote=FALSE)
		message("LOG:\tDiagnosis written to ", outputFile, " successfully.\n", sep="")

		# Plot Scores
		thresholdPlot <- miRNA_plotThresholds(outputDataFrame, inputThreshold, plotFileName, plotFormat=plotFormat)

		message("LOG:\tScore classification plot saved to ", plotFileName, "_score.", plotFormat, " successfully.\n", sep="")
	}

	return(processedFrameToWrite)
}
