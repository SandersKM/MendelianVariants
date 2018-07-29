\name{removeUnwantedVariants}
\alias{removeUnwantedVariants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove variants from gnomAD CSV that don't pass filters
}
\description{
Takes in a CSV of variants from gnomAD, as well as optional filters for annotations and flags. Returns
a data frame with the variants that don't contain the filters and flags
}
\usage{
removeUnwantedVariants(df, unwanted_annotations = c("splice region", "synonymous", "3' UTR", "5' UTR", "downstream gene","intron", "upstream gene", "non coding transcript exon"),
unwanted_flags = c("SEGDUP", "LC LoF"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
  Data frame of the gnomAD CSV
}
  \item{unwanted_annotations}{
  Vector containing gnomAD variant annotations that should be filtered out.
  }
  \item{unwanted_flags}{
  Vector containing gnomAD flags that should be filtered out.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame containing variants that passed the flag and variant filters
}
\references{
http://gnomad.broadinstitute.org
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~manip }% use one of  RShowDoc("KEYWORDS")
