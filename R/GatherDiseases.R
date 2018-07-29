library(rentrez)
library(XML)

gatherDiseases <- function(keyword) {
  disease.search.count <- entrez_search(db = "medgen", term = keyword, retmax = 0)$count
  disease.search <- entrez_search(db = "medgen", term = keyword, retmax = disease.search.count)
  disease_file <- data.frame(id = disease.search$ids, stringsAsFactors = FALSE)
  disease_file$summary <- entrez_summary(db = "medgen", id = disease_file$id)
  disease_file$conceptid <- extract_from_esummary(disease_file$summary, "conceptid")
  disease_file$title <- extract_from_esummary(disease_file$summary, "title")
  disease_file$definition <- sapply(extract_from_esummary(disease_file$summary, "definition"), "[", 1)
  disease_file$definition <- unlist(lapply(disease_file$definition, toString))
  disease_file <- disease_file[, !(names(disease_file) %in% c("summary"))]
  disease_genes <- entrez_link(dbfrom = "medgen", id = disease_file$id, db = "gene")$links[1]$medgen_gene_diseases
  disease_genes_symbols <- extract_from_esummary(entrez_summary(db = "gene", id = disease_genes), c("uid", "name"))
  disease_gene_by_id <- entrez_link(dbfrom = "medgen", id = disease_file$id, db = "gene", by_id = TRUE)
  disease_file$genes_symbol <- character(dim(disease_file)[1])
  disease_file$genes <- character(dim(disease_file)[1])
  # put together genes and gene symbol information for each disease
  for (i in 1:dim(disease_file)[1]) {
    rownum <- match(xmlValue(disease_gene_by_id[i][[1]]$file[2]$IdList[1]$Id[1]$text), disease_file$id)
    disease_file$genes[rownum] <- paste0(disease_gene_by_id[i][[1]]$links$medgen_gene_diseases, collapse = "; ")
    disease_file$genes_symbol[rownum] <- paste0(disease_genes_symbols["name", disease_gene_by_id[i][[1]]$links$medgen_gene_diseases], 
                                                collapse = "; ")
  }
  # take out diseases without associated genes
  disease_file <- disease_file[!disease_file$genes == "", ]
  # Get rid of Summary Column
  disease_file <- disease_file[, !(names(disease_file) %in% c("summary"))]
  disease_file$inheritance <- sapply(disease_file$definition, function(n) {
    def <- strsplit(tolower(n), split = " ")
    if (length(def) > 0) {
      if ("autosomal" %in% def[[1]]) {
        if ("dominant" %in% def[[1]]) {
          return("AD")
        }
        if ("recessive" %in% def[[1]]) {
          return("AR")
        }
        return("A")
      }
      if (("x" %in% def[[1]]) || ("linked" %in% def[[1]]) || ("x-linked" %in% def[[1]])) {
        return("X")
      }
    }
    return("")
  })
  disease_file$inheritance <- unlist(disease_file$inheritance)
  return(disease_file)
}



