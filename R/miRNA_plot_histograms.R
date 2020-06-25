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

#' Plot of classifier distribution histograms (Target and Versus).
#'
#' Generates a plot of the histograms of the score values for the target and versus sets.
#'
#' @param inputDataset Dataset (data frame) to be used for the plot.
#' @param thresholdFrame Diagnostic threshold values (data frame) to be used for the plot. If omitted, no threshold is drawn on the histogram.
#' @param outputFileLabel Label to be used to build the name of the output file.
#' @param plotFormat String to set the format of the output file. Can either be 'pdf' (default) or 'png'.
#'
#' @return A ggplot object containing the plot.
#'
#' This function is not exported to the package NAMESPACE, but it is called by other functions of the same package.
#' @importFrom stats density

miRNA_plotHistograms <- function(inputDataset, thresholdFrame=character(), outputFileLabel, plotFormat="pdf") {

	title <- "Histograms of Target and Versus"
	size <- 24

	dnorm.count <- function(x, mean = 0, sd = 1, log = FALSE, n = 1, binwidth = 1) {
		n * stats::dnorm(x = x, mean = mean, sd = sd, log = log)
	}

	dataTarget <- subset(inputDataset, inputDataset$Diagnosis == 'target')
	dataVersus <- subset(inputDataset, inputDataset$Diagnosis == 'versus')

	if (length(thresholdFrame) > 0) {
		thresholdFlag <- TRUE
		chi <- as.numeric(as.vector(thresholdFrame$Threshold))
		dchi <- as.numeric(as.vector(thresholdFrame$DeltaThreshold))
	} else {
		thresholdFlag <- FALSE
	}

	nT <- length(inputDataset$Score[inputDataset$Diagnosis == "target"])
	xT <- mean(inputDataset$Score[inputDataset$Diagnosis == "target"])
	sT <- stats::sd(inputDataset$Score[inputDataset$Diagnosis == "target"])

	nV <- length(inputDataset$Score[inputDataset$Diagnosis == "versus"])
	xV <- mean(inputDataset$Score[inputDataset$Diagnosis == "versus"])
	sV <- stats::sd(inputDataset$Score[inputDataset$Diagnosis == "versus"])

	if (nT == 0 && nV != 0) {
		maxPlot <- round(max(inputDataset$Score[inputDataset$Diagnosis == "versus"])+0.5,0)
	} else {
		maxPlot <- round(max(inputDataset$Score[inputDataset$Diagnosis == "target"], inputDataset$Score[inputDataset$Diagnosis == "versus"])+0.5,0)
	}

	if (nV == 0 && nT != 0) {
		minPlot <- round(min(inputDataset$Score[inputDataset$Diagnosis == "target"])-0.5,0)
	} else {
		minPlot <- round(min(inputDataset$Score[inputDataset$Diagnosis == "target"], inputDataset$Score[inputDataset$Diagnosis == "versus"])-0.5,0)
	}

	bwT = (max(inputDataset$Score[inputDataset$Diagnosis == "target"]) - min(inputDataset$Score[inputDataset$Diagnosis == "target"]))/round(1.5 + log2(nT));
	bwV = (max(inputDataset$Score[inputDataset$Diagnosis == "versus"]) - min(inputDataset$Score[inputDataset$Diagnosis == "versus"]))/round(1.5 + log2(nV));


	if (nT == 0) {
		plotObject <- ggplot2::ggplot(data=inputDataset, ggplot2::aes_string(x = "Score")) +
			ggplot2::theme(legend.position="none") +
			ggplot2::ggtitle(title) +
			ggplot2::theme(plot.title = ggplot2::element_text(lineheight=.8, size=size, face="bold", hjust = 0.25 + (30-size)*0.1, vjust = -5)) +
			ggplot2::xlab("x") +
			ggplot2::theme(axis.text=ggplot2::element_text(size=14,face="bold",color=1), axis.title=ggplot2::element_text(size=20,face="bold")) +
			ggplot2::scale_x_continuous(breaks=round(seq(minPlot,maxPlot, length.out=10), 1)) +
			ggplot2::scale_y_continuous()+
			ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(density)), data=dataVersus, fill = "red", alpha = 0.5, binwidth=bwV) +
			ggplot2::stat_function(fun = dnorm.count, args = c(mean = xV, sd = sV, n = 1.0, binwidth=bwV), colour = "red")
		if (thresholdFlag) {
		plotObject <- plotObject + ggplot2::geom_vline(xintercept=chi, linetype="solid", size=1, colour="green") +
			ggplot2::geom_vline(xintercept=c(chi-dchi,chi+dchi), linetype=4, size=1, colour="green")
		}
	} else if (nV == 0) {
		plotObject <- ggplot2::ggplot(data=inputDataset, ggplot2::aes_string(x = "Score")) +
			ggplot2::theme(legend.position="none") +
			ggplot2::ggtitle(title) +
			ggplot2::theme(plot.title = ggplot2::element_text(lineheight=.8, size=size, face="bold", hjust = 0.25 + (30-size)*0.1, vjust = -5)) +
			ggplot2::xlab("x") +
			ggplot2::theme(axis.text=ggplot2::element_text(size=14,face="bold",color=1), axis.title=ggplot2::element_text(size=20,face="bold")) +
			ggplot2::scale_x_continuous(breaks=round(seq(minPlot,maxPlot, length.out=10), 1)) +
			ggplot2::scale_y_continuous()+
			ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(density)), data=dataTarget, fill = "blue", alpha = 0.5, binwidth=bwT) +
			ggplot2::stat_function(fun = dnorm.count, args = c(mean = xT, sd = sT, n = 1.0, binwidth=bwT), colour = "blue")
		if (thresholdFlag) {
		plotObject <- plotObject + ggplot2::geom_vline(xintercept=chi, linetype="solid", size=1, colour="green") +
			ggplot2::geom_vline(xintercept=c(chi-dchi,chi+dchi), linetype=4, size=1, colour="green")
		}
	} else {
		plotObject <- ggplot2::ggplot(data=inputDataset, ggplot2::aes_string(x = "Score")) +
			ggplot2::theme(legend.position="none") +
			ggplot2::ggtitle(title) +
			ggplot2::theme(plot.title = ggplot2::element_text(lineheight=.8, size=size, face="bold", hjust = 0.25 + (30-size)*0.1, vjust = -5)) +
			ggplot2::xlab("x") +
			ggplot2::theme(axis.text=ggplot2::element_text(size=14,face="bold",color=1), axis.title=ggplot2::element_text(size=20,face="bold")) +
			ggplot2::scale_x_continuous(breaks=round(seq(minPlot,maxPlot, length.out=10), 1)) +
			ggplot2::scale_y_continuous()+
			ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(density)), data=dataTarget, fill="blue", alpha=0.5, binwidth=bwT, show.legend = FALSE) +
			ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(density)), data=dataVersus, fill="red", alpha=0.5, binwidth=bwV, show.legend = FALSE) +
			ggplot2::stat_function(fun = dnorm.count, args = c(mean = xT, sd = sT, n = 1.0, binwidth=bwT), colour = "blue") +
			ggplot2::stat_function(fun = dnorm.count, args = c(mean = xV, sd = sV, n = 1.0, binwidth=bwV), colour = "red")
		if (thresholdFlag) {
		plotObject <- plotObject + ggplot2::geom_vline(xintercept=chi, linetype="solid", size=1, colour="green") +
			ggplot2::geom_vline(xintercept=c(chi-dchi,chi+dchi), linetype=4, size=1, colour="green")
		}
	}
	switch(plotFormat,
		png = suppressMessages(ggplot2::ggsave(paste(sep="", outputFileLabel, "_histogram.png"), device="png")),		# png case
		suppressMessages(ggplot2::ggsave(paste(sep="", outputFileLabel, "_histogram.pdf"), device=grDevices::cairo_pdf))	# default pdf
	)


	return(plotObject)
}
