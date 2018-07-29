\name{addVariantAnnotations}
\alias{addVariantAnnotations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotate variant file from gnomAD.
}
\description{
Add readable, descriptive annotations and scores to a gnomAD variant file.
}
\usage{
addVariantAnnotations(df, gene.of.interest.symbol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
  gnomAD variant data frame
}
  \item{gene.of.interest.symbol}{
 HGNC Symbol
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame with annotated variants.
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

data("MUC1_variants")
cleanVariants <- removeUnwantedVariants(MUC1_variants)
annotatedVariants <- addVariantAnnotations(cleanVariants, "MUC1")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~database }% use one of  RShowDoc("KEYWORDS")
