% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/im_mat.R
\name{im_mat}
\alias{im_mat}
\title{Plot a heatmap of a matrix in red and blue}
\usage{
im_mat(x, high = "red", xaxt = "n", yaxt = "n", sort = FALSE,
  scale = FALSE, ballance = FALSE, zlim = NA, ...)
}
\arguments{
\item{x}{matrix}

\item{high}{string of either 'red' or 'blue' used to show higher values}

\item{xaxt}{string indicating how to display the x axis. Suppress x axis with 'n'}

\item{yaxt}{string indicating how to display the y axis. Suppress y axis with 'n'}

\item{sort}{logical indicating whether the columns of the matrix should 
be sorted in decreasing order of their means}

\item{scale}{logical indicating whether the matrix should be z scaled to have 
columns with norm zero and standard deviation one.}

\item{ballance}{logical indicating whether to expand the range so it stays centered at zero}

\item{zlim}{numeric bounds on the max and min range for colors.}

\item{...}{other graphical parameters passed to image}
}
\value{
none
}
\description{
Displays a heatmap of a matrix using red and blue colors. Options to scale
and sort as well as any other graphical parameters with ...
}
\examples{
im_mat(matrix(1:12, nrow=3, ncol=4), sort=FALSE, scale=FALSE)
im_mat(matrix(1:12, nrow=3, ncol=4), sort=TRUE, scale=FALSE)
im_mat(matrix(1:12, nrow=3, ncol=4), sort=FALSE, scale=TRUE)
im_mat(matrix(1:12, nrow=3, ncol=4), sort=TRUE, scale=TRUE)
}

