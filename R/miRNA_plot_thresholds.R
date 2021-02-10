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

#' Plot of scores and thresholds of a Bayes classifier.
#'
#' Generates a plot of the classifier scores of a dataset, as well as the corresponding classifier thresholds.
#'
#' @param inputDataset Dataset (data frame) to be used for the plot.
#' @param thresholdsFrame Diagnostic threshold values (data frame) to be used for the plot.
#' @param outputFileLabel Label to be used to build the name of the output file.
#' @param plotFormat String to set the format of the output file. Can either be 'pdf' (default) or 'png'.
#' @param scorePlotParameters String specifying the parameters of the score plot y-axis. If empty, the axis is configured by assessing suitable parameters from the data. The string has to comply with the format "yl_yu_yt", where: yl is the lower y limit; yu is the upper y limit; yt is the interval between tics along the axis.
#' @param colorComplementFlag Boolean option to switch between the default palette (FALSE) and its inverted version (TRUE). Default is FALSE, corresponding to target samples reported in blue and versus samples in red.
#'
#' @return A ggplot object containing the plot.
#'
#' This function is not exported to the package NAMESPACE, but it is called by other functions of the same package.

miRNA_plotThresholds <- function(inputDataset, thresholdsFrame, outputFileLabel, plotFormat="pdf", scorePlotParameters=character(), colorComplementFlag=FALSE) {

	title <- "Diagnosis results"
	size <- 24

	dnorm.count <- function(x, mean = 0, sd = 1, log = FALSE, n = 1, binwidth = 1) {
		n * binwidth * stats::dnorm(x = x, mean = mean, sd = sd, log = log)
	}

	dataTarget <- subset(inputDataset, inputDataset$Diagnosis == 'target')
	dataVersus <- subset(inputDataset, inputDataset$Diagnosis == 'versus')

	if (length(scorePlotParameters) > 0) {
		parseParameters <- strsplit(scorePlotParameters, "_", fixed = TRUE)
		minPlot <- as.double(parseParameters[[1]][1])
		maxPlot <- as.double(parseParameters[[1]][2])
		ytics <- as.double(parseParameters[[1]][3])
	} else {
		maxPlot <- round(max(inputDataset$Score)+0.5,0)
		minPlot <- round(min(inputDataset$Score)-0.5,0)
		ytics <- (maxPlot - minPlot) / 6.0
	}

	if (dim(inputDataset)[1] <= 40) {
		xtics <- 5
	} else if (dim(inputDataset)[1] <= 80) {
		xtics <- 10
	} else if (dim(inputDataset)[1] <= 160) {
		xtics <- 20
	} else {
		xtics <- 40
	}

	chi <- as.numeric(as.vector(thresholdsFrame$Threshold))
	dchi <- as.numeric(as.vector(thresholdsFrame$DeltaThreshold))
	chiUp <- as.numeric(as.vector(thresholdsFrame$ChiUp))
	dchiUp <- as.numeric(as.vector(thresholdsFrame$DChiUp))
	chiDown <- as.numeric(as.vector(thresholdsFrame$ChiDown))
	dchiDown <- as.numeric(as.vector(thresholdsFrame$DChiDown))


	versusColor <- "red"
	targetColor <- "blue"
	rectanglePalette <- c("green", "lightgreen", "yellow", "orange")
	if (colorComplementFlag) {
		versusColor <- "blue"
		targetColor <- "red"
		rectanglePalette <- rev(rectanglePalette)
	}

	rects <- data.frame(ystart = c(-Inf,chiDown,chi,chiUp), yend = c(chiDown,chi,chiUp,Inf), Diagnosis = c("4", "3", "2", "1"))
	plotObject <- ggplot2::ggplot() +
		ggplot2::theme(legend.title=ggplot2::element_text(size=14)) +
		ggplot2::theme(legend.text=ggplot2::element_text(size=14)) +
		ggplot2::geom_rect(data = rects, ggplot2::aes(xmin=-Inf, xmax=Inf, ymin=rects$ystart, ymax=rects$yend, fill=rects$Diagnosis), alpha = 0.4) +
		ggplot2::scale_fill_manual(name= "Score", values=rectanglePalette, labels=c("target", "", "", "versus")) +
		ggplot2::ggtitle(title) +
		ggplot2::theme(plot.title = ggplot2::element_text(lineheight=.8, size=size, face="bold", hjust = 0.25, vjust = -5)) +
		ggplot2::xlab("sample #") + ggplot2::ylab("y") +
		ggplot2::theme(axis.text=ggplot2::element_text(size=14,face="bold",color=1), axis.title=ggplot2::element_text(size=20,face="bold")) +
		ggplot2::scale_x_continuous(breaks=seq(0,dim(inputDataset)[1],by=xtics)) +
		ggplot2::scale_y_continuous(breaks=seq(minPlot,maxPlot,ytics), limits=c(minPlot, maxPlot)) +
		ggplot2::geom_point(data=inputDataset, ggplot2::aes(x=1:nrow(inputDataset), y=inputDataset$Score, color=inputDataset$Diagnosis), size=3) +
		ggplot2::scale_colour_manual(name="Diagnosis", values=c(targetColor, versusColor)) +
		ggplot2::geom_hline(yintercept=max(chi), linetype=1, color="black", size = 1.2) +
		ggplot2::geom_hline(yintercept=max(chi+dchi), linetype=4, color="black", size = 0.6) +
		ggplot2::geom_hline(yintercept=max(chi-dchi), linetype=4, color="black", size = 0.6) +
		ggplot2::geom_hline(yintercept=max(chiUp), linetype=1, color=targetColor, size = 1.0) +
		ggplot2::geom_hline(yintercept=max(chiUp+dchiUp), linetype=4, color=targetColor, size = 0.6) +
		ggplot2::geom_hline(yintercept=max(chiUp-dchiUp), linetype=4, color=targetColor, size = 0.6) +
		ggplot2::geom_hline(yintercept=max(chiDown), linetype=1, color=versusColor, size = 1.0) +
		ggplot2::geom_hline(yintercept=max(chiDown+dchiDown), linetype=4, color=versusColor, size = 0.6) +
		ggplot2::geom_hline(yintercept=max(chiDown-dchiDown), linetype=4, color=versusColor, size = 0.6) +

	switch(plotFormat,
		png = {if (capabilities("png")) {
			suppressMessages(ggplot2::ggsave(paste(sep="", outputFileLabel, "_scores.png"), device="png"))		# png case
		} else {
			warning("WARNING: could not save a png image for the graph; check with 'capabilities(\"png\")' whether you R setup has the capability of saving png pictures.\n")
		}},
		{if (capabilities("cairo")) {
			suppressMessages(ggplot2::ggsave(paste(sep="", outputFileLabel, "_scores.pdf"), device=grDevices::cairo_pdf))	# default pdf
		} else {
			tryCatch(
				{suppressMessages(ggplot2::ggsave(paste(sep="", outputFileLabel, "_scores.pdf"), device="pdf"))},
				error=function(cond) {
					message("WARNING: Could not save a pdf image; no pdf-capable grDevice is available.\n")
				},
				warning=function(cond) {
					message("WARNING: Output pdf file might have issues.\n")
				}
			)
		}}
	)


	return(plotObject)
}
