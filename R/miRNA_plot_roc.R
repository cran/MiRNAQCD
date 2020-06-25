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

#' Plot of classifier ROC curve.
#'
#' Generates a plot of the ROC curve out of the input dataset.
#'
#' @param inputDataset Dataset (data frame) to be used for the plot.
#' @param outputFileLabel Label to be used to build the name of the output file.
#' @param plotFormat String to set the format of the output file. Can either be 'pdf' (default) or 'png'.
#'
#' @return A ggplot object containing the plot.
#'
#' This function is not exported to the package NAMESPACE, but it is called by other functions of the same package.

miRNA_plotROC <- function(inputDataset, outputFileLabel, plotFormat="pdf") {

	title <- "ROC curve"
	size <- 24

	pred <- inputDataset$Score

	obs <- inputDataset$Diagnosis

	pROC_data <- suppressMessages(pROC::roc(obs, pred))

	pROC_data <- data.frame(pROC_data[3:2])
	pROC_data[,1] <- 1 - pROC_data[,1]
	pROC_data <- pROC_data[length(pROC_data[,1]):1,]
	names(pROC_data) <- c("specificities", "sensitivities")

	pROC_data$group <- rep(factor(1), times=c(nrow(pROC_data)))

	plotObject <- ggplot2::ggplot(data=pROC_data, ggplot2::aes(x=pROC_data$specificities, y=pROC_data$sensitivities, colour=pROC_data$group)) +
		ggplot2::geom_line(size = 2, alpha = 0.9) +
		ggplot2::labs(x = "1 - Specificity", y = "Sensitivity") +
		ggplot2::scale_colour_manual(name= "ROC curve", values=c("green","red"), labels="") +
		ggplot2::theme(axis.text=ggplot2::element_text(size=14,face="bold",color=1), axis.title=ggplot2::element_text(size=20,face="bold")) +
		ggplot2::theme(legend.position = "none", legend.justification = c(0, 1), legend.background = ggplot2::element_rect(colour = NA, fill = "white"), legend.title=ggplot2::element_text(size=20), legend.text=ggplot2::element_text(size=20))

	switch(plotFormat,
		png = suppressMessages(ggplot2::ggsave(paste(sep="", outputFileLabel, "_ROC.png"), device="png")),		# png case
		suppressMessages(ggplot2::ggsave(paste(sep="", outputFileLabel, "_ROC.pdf"), device=grDevices::cairo_pdf))	# default pdf
	)

	return(plotObject)
}
