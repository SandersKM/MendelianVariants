\name{plotScoreDistribution}
\alias{plotScoreDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot distribution of number of variant scores passed
}
\description{
Plots the number of variants that pass each number of scores. Overlayed with gene.
}
\usage{
plotScoreDistribution(df, gene.symbol, save.plot = FALSE, output.dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
 Data frame from sorted scored variant file.
}
  \item{gene.symbol}{
 HGNC gene symbol of the gene that the variants are within.
}
  \item{save.plot}{
 TRUE/FALSE
}
  \item{output.dir}{
 Needed if save.plot == TRUE. The directory that the image should be saved to.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A distribution of the variant scores is plotted and optionally saved to a directory.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

data("MUC1_variants")
cleanVariants <- removeUnwantedVariants(MUC1_variants)
annotatedVariants <- addVariantAnnotations(cleanVariants, "MUC1")
scoreSortedVariants <- sortScoredVariants(annotatedVariants)
plotScoreDistribution(scoreSortedVariants, "MUC1")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~hplot }% use one of  RShowDoc("KEYWORDS")
