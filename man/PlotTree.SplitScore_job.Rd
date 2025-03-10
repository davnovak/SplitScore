% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_Visualise.R
\name{PlotTree.SplitScore_job}
\alias{PlotTree.SplitScore_job}
\title{Visualise SplitScore metaclustering using tree topology}
\usage{
\method{PlotTree}{SplitScore_job}(ssj, no_plot = FALSE, palette = NULL)
}
\arguments{
\item{ssj}{\code{SplitScore_job} object}

\item{palette}{optional character vector of colours for distinguishing between metaclusters on the dendrogram plot (with at least as many colours as there were metaclusters generated). Default value is \code{NULL}}
}
\description{
Plots a hierarchical-clustering dendrogram, visualises split scores at each node and shows division into metaclusters.
Metaclustering results are shown by distinctive colouring of subtrees representing metaclusters created by pruning the dendrogram.
The input is an \code{SplitScore_job} object to which both \code{AgglomerativeClustering} and \code{Prune} were applied before.
}
