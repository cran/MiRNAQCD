## --------------------------------------------------------------------------
##
## This file is part of the miRNA-QC-and-Diagnosis software package.
##
## Version 1.1.1 - February 2021
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

#' Classification of a dataset (diagnosis).
#'
#' This function classifies the entries of the input dataset as either target or versus by using the chosen classifier and given the corresponding disgnostic threshold value.
#'
#' This function can also run in 'Performance analysis mode' to evaluate the performance of a classifier by running it on an already-classified dataset.
#' In order to carry out performance analysis, inputDataset has to contain a 'Class' column. Moreover, a list of Target classes has to be provided to the function via the inputTargetList argument.
#'
#' @param inputDataset Dataset (data frame) to be classified. The data frame must comply with the output format of the quality control functions (miRNA_expressionPreprocessing and miRNA_removeOutliers), thus containing the columns 'Subject', 'miRNA', 'Mean', 'StdDev', 'SampleSize'. Any other column is ignored, and any missing column forbids execution. If the 'Performance analysis mode' is selected (see inputTargetList), the dataset has to contain the 'Class' column as well.
#' @param inputMiRNAList List of miRNAs to be used by the classifier. The chosen miRNAs must be present in the 'miRNA' column of the inputDataset.
#' @param coeffList List of coefficients for the classifier. The number of coefficients must be the same as the number of used miRNAs and listed in the same order.
#' @param inputThreshold Diagnostic threshold data frame for the classifier. The data frame must comply with the output format of the classifier setup function (miRNA_classifierSetup), thus containing the columns 'Threshold', 'DeltaThreshold', 'ChiUp', 'DChiUp', 'ChiDown', 'DChiDown'. Any other column is ignored.
#' @param inputTargetList List of classes to use as target. Providing this argument corresponds to selecting the 'Performance analysis mode'. Consequently, inputDataset is expected to contain the 'Class' column as well. The chosen target must correspond to at least one of the classes present in the 'Class' column of the inputDataset.
#' @param inputVersusList List of classes to use as versus in 'Performance analysis mode'. If the argument is left empty, all classes present in the 'Class' column of the inputDataset, minus the Target classes, are used as Versus.
#' @param saveOutputFile Boolean option setting whether results are written to file (TRUE) or not (FALSE). Default is FALSE.
#' @param outputFileBasename Name of the output file where the diagnosis results are to be stored. If not assigned, a filename is automatically generated.
#' @param sep Field separator character for the output file; the default is tabulation.
#' @param plotFormat String specifying the format of generated graphic files (plots): can either be "pdf" (default) or "png".
#' @param scorePlotAscending Boolean option to set the direction in which samples are ordered: TRUE corresponds to samples ordered by ascending score, FALSE corresponds to samples ordered by descending score. Default is TRUE. This argument is meaningful only if saveOutputFile is set to TRUE.
#' @param scorePlotParameters String specifying the y-axis parameters of the score plot. If empty, the axis is configured by assessing suitable parameters from the data.  This argument is meaningful only if saveOutputFile is set to TRUE. The string has to comply with the format "yl_yu_yt", where: yl is the lower y limit; yu is the upper y limit; yt is the interval between tics along the axis.
#' @param histogramParameters (Used in 'Performance analysis mode' only). String specifying the parameters used to build the histogram. If empty, the histogram is built by assessing suitable parameters from the data. This parameter is meaningful only if saveOutputFile is set to TRUE. The string has to comply with the format "xl_xu_bw", where: xl is the lower boundary of the leftmost bin; xu is the upper boundary of the rightmost bin; bw is the bin width.
#' @param colorComplementFlag Boolean option to switch between the default palette (FALSE) and its inverted version (TRUE). Default is FALSE, corresponding to target samples reported in blue and versus samples in red. This argument is meaningful only if saveOutputFile is set to TRUE.
#'
#' @return A data frame containing the columns 'Subject', 'Diagnosis' and 'Score'.
#'
#' Please refer to the user manual installed in "/path-to-library/MiRNAQCD/doc/manual.pdf" for detailed function documentation. The path "/path-to-library" can be shown from R by calling ".libPaths()"
#'
#' @examples
#' requiredDataFile = paste(system.file(package="MiRNAQCD"),
#'			"/extdata/test_dataset_beta_clean.dat", sep='')
#' myDataFrame <- read.table(file=requiredDataFile, header=TRUE)
#' requiredThresholdFile = paste(system.file(package="MiRNAQCD"),
#'			"/extdata/test_dataset_alpha_threshold.txt", sep='')
#' thresholdDataFrame <- read.table(file=requiredThresholdFile, header=TRUE)
#' mirnaToUse <- c("FX", "FZ")
#' coefficientsToUse <- c(1.0, -1.0)
#'
#' ## Classification
#' classifiedDataset <- miRNA_diagnosis(myDataFrame, mirnaToUse, coefficientsToUse,
#'					thresholdDataFrame)

#' @export
miRNA_diagnosis <- function(inputDataset, inputMiRNAList, coeffList, inputThreshold, inputTargetList=character(), inputVersusList=character(), saveOutputFile=FALSE, outputFileBasename="", sep='\t', plotFormat="pdf", scorePlotParameters=character(), scorePlotAscending=TRUE, colorComplementFlag=FALSE, histogramParameters=character()) {

	## Input validation and preprocessing

	if (!(("Subject" %in% colnames(inputDataset)) & ("miRNA" %in% colnames(inputDataset)) & ("Mean" %in% colnames(inputDataset)) & ("StdDev" %in% colnames(inputDataset)) & ("SampleSize" %in% colnames(inputDataset))))  {
		stop("ERROR: unsuitable dataset format. Dataset must contain columns 'Subject', 'miRNA', 'Mean', 'StdDev', 'SampleSize'.\n")
	}

	if (length(inputDataset[1,]) > 6) {
		warning("WARNING: more than 6 dataset columns. Columns other than 'Subject', 'miRNA', 'Mean', 'StdDev', 'SampleSize', 'Class' will be ignored.\n")
	}

	availableClasses <- unique(inputDataset$Class)

	classIsAvailable <- FALSE
	if ("Class" %in% colnames(inputDataset)) {
		selectedDataset <- subset(inputDataset, select=c("Subject", "miRNA", "Mean", "StdDev", "SampleSize", "Class"))
		classIsAvailable <- TRUE
	} else {
		selectedDataset <- subset(inputDataset, select=c("Subject", "miRNA", "Mean", "StdDev", "SampleSize"))
		classIsAvailable <- FALSE
	}

	# Decide if posterior analysis can be performed, check Target/Versus lists
	performPosteriorAnalysis <- FALSE
	if (length(inputTargetList) == 0) {
		performPosteriorAnalysis <- FALSE
	} else {
		if (classIsAvailable) {
			performPosteriorAnalysis <- TRUE
			message("LOG:\tmiRNA_diagnosis() will perform posterior (performance) analysis.\n")
			listOfTargets <- unique(inputTargetList)
			if (length(listOfTargets) != length(inputTargetList)){
				warning("WARNING: target list presents some duplicates which will be ignored.\n")
			}
			inputTargetList <- unique(inputTargetList)
			listOfTargets <- intersect(listOfTargets, availableClasses)
			if (length(listOfTargets) != length(inputTargetList)){
				warning("WARNING: some of the target classes are not present in the dataset and will be ignored.\n")
			}

			if (length(inputVersusList) == 0) {
				inputVersusList <- setdiff(availableClasses, listOfTargets)
				inputVersusList <- unique(inputVersusList)
			}
			listOfVersus <- unique(inputVersusList)
			if (length(listOfVersus) != length(inputVersusList)){
				warning("WARNING: versus list presents some duplicates which will be ignored.\n")
			}
			inputVersusList <- unique(inputVersusList)
			listOfVersus <- intersect(listOfVersus, availableClasses)
			if (length(listOfVersus) != length(inputVersusList)){
				warning("WARNING: some of the versus classes are not present in the dataset and will be ignored.\n")
			}

			if (length(listOfTargets)==0 | length(listOfVersus)==0) {
				stop("ERROR: unsuitable function arguments. The requested target and/or versus classes are empty, or not present in the dataset.\n")
			} else if (length(intersect(listOfTargets, listOfVersus))) {
				stop("ERROR: conflicting function arguments; target set and versus set have non-empty intersection.\n")
			}
		} else {
			warning("WARNING: posterior analysis is requested (inputTargetList is non-empty), but no 'Class' column is available. Posterior (performance) analysis will NOT be carried out.\n")
			performPosteriorAnalysis <- FALSE
		}
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

	# Bookkeeping, if Class is present
	if (classIsAvailable) {
		columnSubjectClass <- subset(availableDataset, select=c("Subject", "Class"))
	}

	# Prepare frame to perform diagnosis
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
	if (scorePlotAscending) {
		classifierDataFrame <- data.frame(classifierDataFrame[with(classifierDataFrame, order(Score)),])
	} else {
		classifierDataFrame <- data.frame(classifierDataFrame[with(classifierDataFrame, order(-1.0*Score)),])
	}

	# Create new Diagnosis columns and writes "target" or "versus" at each row, then join all together
	targetFrame <- classifierDataFrame[classifierDataFrame$Score > threshold, ]
	targetFrame <- cbind(targetFrame, Diagnosis=rep("target", length(targetFrame[,1])))

	versusFrame <- classifierDataFrame[classifierDataFrame$Score <= threshold, ]
	versusFrame <- cbind(versusFrame, Diagnosis=rep("versus", length(versusFrame[,1])))

	outputDataFrame <- rbind(targetFrame, versusFrame)
	outputDataFrame <- subset(outputDataFrame, select=c("Subject", "Diagnosis", "Score", listOfFeature))
	if (scorePlotAscending) {
		outputDataFrame <- unique(outputDataFrame[with(outputDataFrame, order(Score)),])
	} else {
		outputDataFrame <- unique(outputDataFrame[with(outputDataFrame, order(-1.0*Score)),])
	}

	processedFrameToWrite <- subset(outputDataFrame, select=c("Subject", "Diagnosis", "Score"))
	if (classIsAvailable) {
		processedFrameToWrite <- unique(merge(processedFrameToWrite, columnSubjectClass, by = "Subject"))
	}

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
		if (!performPosteriorAnalysis) {
			thresholdPlot <- miRNA_plotThresholds(outputDataFrame, inputThreshold, plotFileName, plotFormat=plotFormat, scorePlotParameters=scorePlotParameters, colorComplementFlag=colorComplementFlag)
			message("LOG:\tScore classification plot saved to ", plotFileName, "_score.", plotFormat, " successfully.\n", sep="")
		}
	}

	## Posterior (performance) analysis - if available and requested
	if (performPosteriorAnalysis) {

		if (saveOutputFile) {
			performanceFileName=paste(outputFileBasename, "_performance.txt", sep="")
		} else {
			performanceFileName=tempfile()
		}

		dataFrameTemp_T <- unique(processedFrameToWrite[processedFrameToWrite$Class %in% listOfTargets, ])
		names(dataFrameTemp_T) <- c("Subject", "ClassifierDiagnosis", "Score", "Class")
		nr_target_true <- NROW(dataFrameTemp_T[dataFrameTemp_T$ClassifierDiagnosis == "target", ])
		nr_target_false <- NROW(dataFrameTemp_T[dataFrameTemp_T$ClassifierDiagnosis == "versus", ])
		completeDataFrame <- cbind(dataFrameTemp_T, Diagnosis=rep("target", length(dataFrameTemp_T[,1])))

		dataFrameTemp_V <- unique(processedFrameToWrite[processedFrameToWrite$Class %in% listOfVersus, ])
		names(dataFrameTemp_V) <- c("Subject", "ClassifierDiagnosis", "Score", "Class")
		nr_versus_true <- NROW(dataFrameTemp_V[dataFrameTemp_V$ClassifierDiagnosis == "versus", ])
		nr_versus_false <- NROW(dataFrameTemp_V[dataFrameTemp_V$ClassifierDiagnosis == "target", ])
		completeDataFrame <- rbind(completeDataFrame, cbind(dataFrameTemp_V, Diagnosis=rep("versus", length(dataFrameTemp_V[,1]))))

		names(completeDataFrame) <- c("Subject", "ClassifierDiagnosis", "Score", "Class", "Diagnosis")

		performance_sensitivity <- nr_target_true / (nr_target_true + nr_target_false)
		performance_specificity <- nr_versus_true / (nr_versus_true + nr_versus_false)
		performance_accuracy <- (nr_target_true + nr_versus_true) / (nr_target_true + nr_versus_true + nr_target_false + nr_versus_false)
		performance_F1score <- 2.0*nr_target_true / (2.0*nr_target_true + nr_target_false + nr_versus_false)

		# ROC and corresponding area under curve
		pred <- completeDataFrame$Score
		obs <- completeDataFrame$Diagnosis
		pROC_data <- suppressMessages(pROC::roc(obs, pred))
		ciAuc <- suppressMessages(pROC::ci.auc(pROC_data))
		performance_AUC <- ciAuc[2];
		xT <- mean(completeDataFrame$Score[completeDataFrame$Diagnosis == "target"])
		sT <- stats::sd(completeDataFrame$Score[completeDataFrame$Diagnosis == "target"])
		xV <- mean(completeDataFrame$Score[completeDataFrame$Diagnosis == "versus"])
		sV <- stats::sd(completeDataFrame$Score[completeDataFrame$Diagnosis == "versus"])
		d <- xT - xV
		dd <- sqrt((sT*sT + sV*sV)/2.)
		dee <- d/dd

		cat("#### Classifier performance\n\n", append=FALSE, file=performanceFileName)
		cat("Confusion matrix:\tPredicted\tPredicted\n", append=TRUE, file=performanceFileName)
		cat("\t\t\tTarget\t\tVersus\n", append=TRUE, file=performanceFileName)
		cat(paste("Actual Target\t\t", nr_target_true, "\t\t", nr_target_false, "\n", sep=""), append=TRUE, file=performanceFileName)
		cat(paste("Actual Versus\t\t", nr_versus_false, "\t\t", nr_versus_true, "\n", sep=""), append=TRUE, file=performanceFileName)

		cat(paste("\n\nAccuracy:\t", sprintf("%.3f", performance_accuracy), "\n", sep=""), append=TRUE, file=performanceFileName)
		cat(paste("Specificity:\t", sprintf("%.3f", performance_specificity), "\n", sep=""), append=TRUE, file=performanceFileName)
		cat(paste("Sensitivity:\t", sprintf("%.3f", performance_sensitivity), "\n", sep=""), append=TRUE, file=performanceFileName)
		cat(paste("F1-score:\t", sprintf("%.3f", performance_F1score), "\n", sep=""), append=TRUE, file=performanceFileName)
		cat(paste("dPrime:\t\t", sprintf("%.3f", performance_F1score), "\n", sep=""), append=TRUE, file=performanceFileName)
		cat(paste("ROC AUC:\t", sprintf("%.3f", performance_AUC), "\n", sep=""), append=TRUE, file=performanceFileName)

		if (saveOutputFile) {
			if (scorePlotAscending) {
				completeDataFrame <- data.frame(completeDataFrame[with(completeDataFrame, order(Score)),])
			} else {
				completeDataFrame <- data.frame(completeDataFrame[with(completeDataFrame, order(-1.0*Score)),])
			}
			message("LOG:\tclassifier performance metrics written to ", performanceFileName, " successfully.\n", sep="")
			thresholdPlot <- miRNA_plotThresholds(completeDataFrame, inputThreshold, plotFileName, plotFormat=plotFormat, scorePlotParameters=scorePlotParameters, colorComplementFlag=colorComplementFlag)
			message("LOG:\tScore classification plot saved to ", plotFileName, "_score.", plotFormat, " successfully.\n", sep="")
			histPlot <- miRNA_plotHistograms(completeDataFrame, inputThreshold, plotFileName, plotFormat=plotFormat, histogramParameters=histogramParameters, colorComplementFlag=colorComplementFlag)
			message("LOG:\tHistogram plot saved to ", plotFileName, "_histogram.", plotFormat, " successfully.\n", sep="")
			rocPlot <- miRNA_plotROC(completeDataFrame, plotFileName, plotFormat=plotFormat)
			message("LOG:\tROC plot saved to ", plotFileName, "_ROC.", plotFormat, " successfully.\n", sep="")
		} else {
			cat(readLines(performanceFileName), sep="\n")
		}
	}

	return(processedFrameToWrite)
}
