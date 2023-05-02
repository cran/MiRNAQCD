## --------------------------------------------------------------------------
##
## This file is part of the miRNA-QC-and-Diagnosis software package.
##
## Version 1.1.3 - April 2023
##
##
## The miRNA-QC-and-Diagnosis package is free software; you can use it,
## redistribute it, and/or modify it under the terms of the GNU General
## Public License version 3 as published by the Free Software Foundation.
## The full text of the license can be found in the file LICENSE.txt at the top
## level of the package distribution.
##
## Authors:
##	Michele Castelluzzo (1), Alessio Perinelli (1), Simone Detassis (3),
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

#' Assessment of quality threshold values.
#'
#' This function assesses a set of quality threshold values (standard deviations), one for each miRNA, out of a dataset.
#'
#' @param inputDataset Dataset (data frame) to be used for the assessment. The data frame must comply with the output format of the preprocessing function (miRNA_expressionPreprocessing), thus containing the columns 'Subject', 'miRNA', 'Mean', 'StdDev', 'SampleSize' and possibly 'Class'. Any other column is ignored, and any missing column forbids execution.
#' @param significanceLevel Significance level to be used for the assessment (must be greater than zero and less than one). Default is 0.05 (i.e. 5 percent).
#' @param saveOutputFile Boolean option setting whether results are written to file (TRUE) or not (FALSE). Default is FALSE.
#' @param outputFileName Name of the output file where the quality threshold values are to be stored. If not assigned, a filename is automatically generated.
#' @param sep Field separator character for the output files; the default is tabulation.
#'
#' @return A data frame of quality threshold values, containing the columns 'miRNA' and 'QualityThreshold'.
#'
#' Please refer to the user manual installed in "/path-to-library/MiRNAQCD/doc/manual.pdf" for detailed function documentation. The path "/path-to-library" can be shown from R by calling ".libPaths()"
#'
#' @examples
#' requiredFile = paste(system.file(package="MiRNAQCD"),
#'			"/extdata/test_dataset_alpha_prep.dat", sep='')
#' myDataFrame <- read.table(file=requiredFile, header=TRUE)
#' qt <- miRNA_assessQualityThreshold(myDataFrame, significanceLevel=0.05)

#' @export
miRNA_assessQualityThreshold <- function(inputDataset, significanceLevel=0.05, saveOutputFile=FALSE, outputFileName="", sep='\t') {

	## Input validation

	if (!(("Subject" %in% colnames(inputDataset)) & ("miRNA" %in% colnames(inputDataset)) & ("Mean" %in% colnames(inputDataset)) & ("StdDev" %in% colnames(inputDataset)) & ("SampleSize" %in% colnames(inputDataset)))) {
		stop("ERROR: unsuitable dataset format. Dataset must contain columns 'Subject', 'miRNA', 'Mean', 'StdDev', 'SampleSize'.\n")
	}

	if (length(inputDataset[1,]) > 7) {
		warning("WARNING: more than 7 dataset columns. Columns other than 'Subject', 'miRNA', 'Mean', 'StdDev', 'SampleSize', 'Class' will be ignored.\n")
	}

	if ("Class" %in% colnames(inputDataset)) {
		selectedDataset <- subset(inputDataset, select=c("Subject", "miRNA", "Mean", "StdDev", "SampleSize", "Class"))
		columnClassExists <- TRUE
	} else {
		selectedDataset <- subset(inputDataset, select=c("Subject", "miRNA", "Mean", "StdDev", "SampleSize"))
		columnClassExists <- FALSE
	}

	if (sep == "") {
		sep <- " "
	}

	multipletSize <- unique(selectedDataset$SampleSize)
	if (length(multipletSize) > 1) {
		stop("ERROR: multiple sample sizes detected. Dataset must contain a unique value of SampleSize.\n")
	}
	multipletSize <- as.numeric(as.vector(multipletSize))

	## Define custom ks-test function

	ksTestPValue <- function(StdDev, varianceVector){
		res <- stats::ks.test(varianceVector/(StdDev*StdDev)*(multipletSize-1), stats::pchisq, df=multipletSize-1, alternative = "t")
		return(res[[2]])
	}

	## Assessment of critical sigma

	setOfFeatures <- unique(selectedDataset$miRNA)
	setOfFeatures <- setOfFeatures[order(setOfFeatures)]
	for (feature in setOfFeatures) {
		tempFrame <- subset(selectedDataset, selectedDataset$miRNA == feature)
		varianceColumn <- tempFrame$StdDev*tempFrame$StdDev
		varianceColumn <- data.frame(varianceColumn)
		names(varianceColumn) <- "Variance"
		tempFrame <- cbind(tempFrame, varianceColumn)

		# Find the value of StdDev in interval that gives maximum value of ksTestPValue

		ks_v <- as.numeric(stats::optimize(ksTestPValue, varianceVector=tempFrame$Variance, interval=c(0, 0.5), maximum=TRUE))


		sigmaLimit <- ks_v[1] * sqrt(-2*log(significanceLevel))

		tempRow <- cbind(feature, round(sigmaLimit,2))
		if (!exists("criticalValuesTemp")) criticalValuesTemp <- tempRow
		else criticalValuesTemp <- rbind (criticalValuesTemp, tempRow)
	}

	qualityThresholdFrame <- data.frame(criticalValuesTemp)
	names(qualityThresholdFrame) <- c("miRNA", "QualityThreshold")

	## Output
	if (saveOutputFile) {
		if (outputFileName == "") {
			outputFileName <- paste("quality_threshold_m_", multipletSize, "_significance_", significanceLevel, ".dat", sep="")
		}
		if (file.exists(outputFileName) & file.access(outputFileName, mode=2)) {
			stop("ERROR: cannot write ", outputFileName, ". Check write permission.\n")
		}
		cat(c("# multiplet size =", multipletSize, "\n"), file=outputFileName)
		cat(c("# significance level =", significanceLevel, "\n"), file=outputFileName, append=TRUE)


		utils::write.table(format(qualityThresholdFrame, drop0trailing=FALSE), file=outputFileName, append=TRUE, sep=sep, col.names=TRUE, row.names=FALSE, quote=FALSE)


		message("LOG:\tQuality threshold values written to ", outputFileName, " successfully.\n", sep="")
	}

	return(qualityThresholdFrame)
}
