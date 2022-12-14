% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TIGER.R
\name{el2regulon}
\alias{el2regulon}
\title{Convert a bipartite edgelist to regulon}
\usage{
el2regulon(el)
}
\arguments{
\item{el}{An edge list dataframe with three columns. First column is TF name,
second column is gene name, and third column is edge weight.}
}
\value{
A VIPER required regulon object
}
\description{
Convert a bipartite edgelist to regulon
}
