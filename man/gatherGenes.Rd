\name{gatherGenes}
\alias{gatherGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotate disease-related genes
}
\description{
Takes the disease file and annotates each gene using NCBI Gene and Ensembl.
}
\usage{
gatherGenes(disease_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{disease_file}{
The disease data frame returned from gatherDiseases(keyword)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame where each row is a gene. Includes information about phenotype, location, and exons.
}
\references{
https://www.ncbi.nlm.nih.gov/gene/
http://grch37.ensembl.org/biomart/martview/a751c8ace82f5c63091516237797a2fe
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

## The function is currently defined as

diseases <- gatherDiseases("Kidney Disease")
genes <- gatherGenes(diseases)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~database }% use one of  RShowDoc("KEYWORDS")
