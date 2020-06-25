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

#' Load quality threshold values.
#'
#' This function loads from file a data frame containing the quality threshold values for a set of miRNAs.
#'
#' @param inputFileName Name of the file to be loaded. The file has to contain at least the columns 'miRNA', 'QualityThreshold' (not necessarily in this order).
#' @param sep Field separator character; the default is any white space (one or more spaces or tabulations).
#'
#' @return A data frame containing the columns 'miRNA' and 'QualityThreshold'.
#'
#' @examples
#' requiredFile = paste(system.file(package="MiRNAQCD"), "/extdata/dataset_alpha_qt.dat", sep='')
#' qtDataFrame <- miRNA_loadQualityThreshold(requiredFile)

#' @export
miRNA_loadQualityThreshold <- function(inputFileName, sep="") {

	if (!file.exists(inputFileName)) {
		stop("ERROR: cannot read ", inputFileName, ". No such file or directory.\n")
	} else if (file.access(inputFileName, mode=4)) {
		stop("ERROR: cannot read ", inputFileName, ". Check read permission.\n")
	}

	if (sep == "") {
		readFrame <- utils::read.table(file=inputFileName, fileEncoding="UTF-8", header=T)
	} else {
		readFrame <- utils::read.table(file=inputFileName, sep=sep, fileEncoding="UTF-8", header=T)
	}

	if (!(("miRNA" %in% colnames(readFrame)) & ("QualityThreshold" %in% colnames(readFrame)))) {
		stop("ERROR: unsuitable dataset format. Dataset must contain columns 'miRNA', 'QualityThreshold'.\n")
	}

	if (length(readFrame[1,]) > 2) {
		warning("WARNING: more than 2 dataset columns. Columns other than 'miRNA', 'QualityThreshold' will be ignored.\n")
	}

	outputFrame <- subset(readFrame, select=c("miRNA", "QualityThreshold"))

	return(outputFrame)
}
