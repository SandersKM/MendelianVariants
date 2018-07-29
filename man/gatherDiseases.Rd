\name{gatherDiseases}
\alias{gatherDiseases}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gather diseases from MedGen.
}
\description{
gatherDiseases takes a phenotype keyword, queries MedGen, and returns all unique diseases
that have related genes.
}
\usage{
gatherDiseases(keyword)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{keyword}{
    Phenotype of interest. To avoid overflow, use "Kidney Disease" rather than "Kidney"
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Data frame with all diseases from MedGen with associated genes
\references{
 https://www.ncbi.nlm.nih.gov/medgen/
}
\author{
Kat Sanders
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

gatherDiseases("Kidney Disease")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~database }% use one of  RShowDoc("KEYWORDS")

