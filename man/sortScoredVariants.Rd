\name{sortScoredVariants}
\alias{sortScoredVariants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sorts variants by the number of score cutoffs passed
}
\description{
Takes annotated variants with scores, compares the scores of each variant to score cutoffs, and organizes
variants based on how variant scores are passed.
}
\usage{
sortScoredVariants(df, cadd.min = 20, fitCon.min = 0.4, phastCon.min = 0.55, mcap.min = 0.025, AFgnomAD.max = 0.01, sift.max = 0.05, polyphen.min = 0.909)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
 Annoted variant dataframe.
}
  \item{cadd.min}{
 Minimum CADD score accepted
}
  \item{fitCon.min}{
Minimum fitCon score accepted
}
  \item{phastCon.min}{
Minimum phastCon score accepted
}
  \item{mcap.min}{
Minimum M-CAP score accepted
}
  \item{AFgnomAD.max}{
Maximum gnomAD Allele Frequency accepted
}
  \item{sift.max}{
Maximum SIFT score accepted

}
  \item{polyphen.min}{
Miniumum PolyPhen accepted

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame of variants sorted by the number of variant scores passed.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kate Sanders
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~manip }% use one of  RShowDoc("KEYWORDS")
