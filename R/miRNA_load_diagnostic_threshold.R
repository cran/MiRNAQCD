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

#' Load diagnostic threshold values.
#'
#' This function loads from file a data frame containing the diagnostic threshold values of a trained classifier.
#'
#' @param inputFileName Name of the file to be loaded. The file has to contain at least the columns 'Threshold', 'DeltaThreshold', 'ChiUp', 'DChiUp', 'ChiDown', 'DChiDown' (not necessarily in this order).
#' @param sep Field separator character; the default is any white space (one or more spaces or tabulations).
#'
#' @return A data frame containing all the columns present in the file.
#'
#' Please refer to the user manual installed in "/path-to-library/MiRNAQCD/doc/manual.pdf" for detailed function documentation. The path "/path-to-library" can be shown from R by calling ".libPaths()"
#'
#' @examples
#' requiredFile = paste(system.file(package="MiRNAQCD"),
#'		"/extdata/test_dataset_alpha_threshold.txt", sep='')
#' threshold <- miRNA_loadDiagnosticThreshold(requiredFile)

#' @export
miRNA_loadDiagnosticThreshold <- function(inputFileName, sep="") {

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

	if (!(("Threshold" %in% colnames(readFrame)) & ("DeltaThreshold" %in% colnames(readFrame)) & ("ChiUp" %in% colnames(readFrame)) & ("DChiUp" %in% colnames(readFrame)) & ("ChiDown" %in% colnames(readFrame)) & ("DChiDown" %in% colnames(readFrame)) )) {
		stop("ERROR: unsuitable dataset format. Dataset must contain columns 'Threshold', 'DeltaThreshold', 'ChiUp', 'DChiUp', 'ChiDown', 'DChiDown'.\n")
	}

	return(readFrame)
}
