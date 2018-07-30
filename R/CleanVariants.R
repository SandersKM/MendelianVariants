removeUnwantedVariants <- function(df, unwanted_annotations = c("splice region", "synonymous", "3' UTR", "5' UTR", "downstream gene",
                                                                "intron", "upstream gene", "non coding transcript exon"),
                                   unwanted_flags = c("SEGDUP", "LC LoF")) {
  df$temp <- TRUE
  for (n in 1:dim(df)[1]) {
    if (df$Annotation[n] %in% unwanted_annotations) {
      df$temp[n] <- FALSE
    }
    if(length(unwanted_flags) > 0 ){
      if(df$Flags[n] %in% unwanted_flags | df$Flags[n] == "SEGDUP LC LoF"){
        df$temp[n] <- FALSE
      }
    }
  }

  df <- df[!df$temp == FALSE, !(names(df) %in% c("clean", "Flags", "temp"))]
  return(df)
}

addVariantAnnotations <- function(df, gene.of.interest.symbol){
  try({
    api.response <- fromJSON(paste("http://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/", gene.of.interest.symbol,
                                   "?content-type=application/json;expand=1", sep = ""))

    gene.of.interest.ensembl_gene_id <- api.response$id
    gene.of.interest.start <- api.response$start
    gene.of.interest.end <- api.response$end
    gene.of.interest.strand <- api.response$strand
    gene.of.interest.ch <- api.response$seq_region_name
    transcripts <- as.data.frame(api.response)
    canonical.rownum <- which(transcripts$Transcript.is_canonical == 1)
    gene.of.interest.cannonical_transcript <- transcripts$Transcript.id[canonical.rownum]
    gene.of.interest.exon <- cbind(transcripts$Transcript.Exon[canonical.rownum][[1]]$start,
                                   transcripts$Transcript.Exon[canonical.rownum][[1]]$end)
  }, silent = TRUE)

  # find the cannonical exon (if any) each variant falls into
  df$exon <- sapply(df$Position, function(position){
    closest_exon <- 0
    exon_dist <- abs(position - gene.of.interest.exon[1,1])
    symb <- "+"
    for(i in 1:length(gene.of.interest.exon[,1])){
      if((position >= gene.of.interest.exon[i, 1]) && (position <= gene.of.interest.exon[i,2])){
        return(i)
      }
      if(abs(position - gene.of.interest.exon[i, 1]) < exon_dist){
        closest_exon <- i
        exon_dist <- abs(position - gene.of.interest.exon[i,2])
        symb <- "+"
      }
      if(abs(position - gene.of.interest.exon[i,2]) < exon_dist){
        closest_exon <- i
        exon_dist <- abs(position - gene.of.interest.exon[i,2])
        symb <- "-"
      }
    }
    return(paste(closest_exon, symb, exon_dist, sep = ""))
  })
  df$ancestors <- sapply(1:dim(df)[1], get_ancestors, df = df)
  drops <- c("Allele.Count.African", "Allele.Number.African", "Homozygote.Count.African",
             "Allele.Count.Ashkenazi.Jewish", "Allele.Number.Ashkenazi.Jewish", "Homozygote.Count.Ashkenazi.Jewish",
             "Allele.Count.East.Asian", "Allele.Number.East.Asian", "Homozygote.Count.East.Asian",
             "Allele.Count.European..Finnish.", "Allele.Number.European..Finnish.", "Homozygote.Count.European..Finnish.",
             "Allele.Count.European..Non.Finnish.", "Allele.Number.European..Non.Finnish.", "Homozygote.Count.European..Non.Finnish.",
             "Allele.Count.Latino", "Allele.Number.Latino", "Homozygote.Count.Latino",
             "Allele.Count.South.Asian", "Allele.Number.South.Asian", "Homozygote.Count.South.Asian",
             "Allele.Count.Other", "Allele.Number.Other", "Homozygote.Count.Other")
  df <- df[ , !(names(df) %in% drops)]

  # function to get gnomAD website for specific variant
  df$gnomAD.website <- sapply(1:dim(df)[1], function(n){
    base_url <- "http://gnomad.broadinstitute.org/variant/"
    variantID <- paste(df$Chrom[n], df$Position[n], df$Reference[n], df$Alternate[n], sep = "-")
    return(paste(base_url,variantID,sep = ""))
  })


  df <- df[order(df$Position),]
  # used to get the row number of the scores
  get_position_offset <- function(n){
    offset <- df$Position[n] - gene.of.interest.start + 1
    return(offset)
  }
  df$distance.from.start <- sapply(1:dim(df)[1], get_position_offset)

  tryCatch({library(phastCons100way.UCSC.hg19)}, 
           error = function(e){
             source("https://bioconductor.org/biocLite.R")
             biocLite("phastCons100way.UCSC.hg19")
             library(phastCons100way.UCSC.hg19)
           })
  tryCatch({library(fitCons.UCSC.hg19)}, 
           error = function(e){
             source("https://bioconductor.org/biocLite.R")
             biocLite("fitCons.UCSC.hg19")
             library(fitCons.UCSC.hg19)
           })
 
  gr <- GRanges(seqnames=gene.of.interest.ch,
                IRanges(start=gene.of.interest.start:gene.of.interest.end, width=1))

  # phastCons100way.UCSC.hg19 - phastCon scores are derived from the alignment of the human genome (hg19)
  # and 99 other vertabrate species
  phastCon.scores <- scores(phastCons100way.UCSC.hg19, gr)
  df$phastCon.score <- phastCon.scores[df$distance.from.start]$scores

  # fitCons.UCSC.hg19 - fitCons scores measure the fitness consequences of function annotation for the
  # human genome (hg19)
  fitCon.scores <- scores(fitCons.UCSC.hg19, gr)
  df$fitCon.score <- fitCon.scores[df$distance.from.start]$scores

  # This is to get the position of the scores for GScores where there are multiple allele options
  alleles <- c("A", "C", "G", "T")
  df$temp <- sapply(1:dim(df)[1], function(n){
    which(alleles[which(alleles != df$Reference[n])] == df$Alternate[n])[1]
  })

  # cadd.v1.3.hg19 - fitCons scores measure the fitness consequences of function annotation for the
  # human genome (hg19)
  # These scores are rounded to provide faster lookup
  if (!exists("cadd")){
    cadd <- getGScores("cadd.v1.3.hg19")
  }

  cadd.scores <- scores(cadd, gr)
  df$cadd.score <- sapply(1:dim(df)[1], function(n){
    if(is.na(df$temp[n])){return(NULL)}
    if(df$temp[n] == 1){
      return(cadd.scores[df$distance.from.start[n]]$scores1[1])
    }
    if(df$temp[n] == 2){
      return(cadd.scores[df$distance.from.start[n]]$scores2[1])
    }
    if(df$temp[n] == 3){
      return(cadd.scores[df$distance.from.start[n]]$scores3[1])
    }
  })
  df$cadd.score <- lapply(df$cadd.score , toString)
  df$cadd.score <- unlist(df$cadd.score)
  df$cadd.score <- lapply(df$cadd.score, as.integer)

  # get mcap scores
  if (!exists("mcap")){
    mcap <- getGScores("mcap.v1.0.hg19")
  }

  mcap.scores <- scores(mcap, gr)
  df$mcap.score <- sapply(1:dim(df)[1], function(n){
    if(is.na(df$temp[n])){return(NULL)}
    if(df$temp[n] == 1){
      return(mcap.scores[df$distance.from.start[n]]$scores1[1])
    }
    if(df$temp[n] == 2){
      return(mcap.scores[df$distance.from.start[n]]$scores2[1])
    }
    if(df$temp[n] == 3){
      return(mcap.scores[df$distance.from.start[n]]$scores3[1])
    }
  })
  df$mcap.score <- lapply(df$mcap.score , toString)
  df$mcap.score <- unlist(df$mcap.score)
  df$mcap.score <- lapply(df$mcap.score, as.numeric)


  # Fetch sift and polyphen scores from Ensembl

  df$polyphen.score <- character(dim(df)[1])
  df$sift.score <- character(dim(df)[1])
  df$polyphen.prediction <- character(dim(df)[1])
  df$sift.prediction <- character(dim(df)[1])
  df$consequences.all <- character(dim(df)[1])
  df$impact.all <- character(dim(df)[1])

  server <- "http://grch37.rest.ensembl.org"
  ext <- "/vep/human/hgvs"
  i <- 1
  while(i < dim(df)[1]){
    j <- i + 290 # Ensembl takes at most 300 requests at a time.
    if(j > dim(df)[1]){
      j = dim(df)[1]
    }
    rest.api.response <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"),
                              body = paste('{ "hgvs_notations" : [', paste0(gene.of.interest.symbol, ":",  df$Transcript.Consequence[i:j],
                                                                            collapse = "\",\""), ' ] }', sep = "\""))

    rest.api.info <- fromJSON(toJSON(content(rest.api.response)))
    for(k in 1:length(rest.api.info$transcript_consequences)){
      rest.info <- rest.api.info$transcript_consequences[[k]]
      n <- which(df$Transcript.Consequence == strsplit(rest.api.info$input[[k]],split = ":")[[1]][2])
      info <- rest.info[rest.info$gene_symbol == gene.of.interest.symbol &
                          rest.info$transcript_id == gene.of.interest.cannonical_transcript,]

      if( "polyphen_prediction" %in% colnames(info)){
        df$polyphen.prediction[n] <- info$polyphen_prediction
      }
      if( "sift_prediction" %in% colnames(info)){
        df$sift.prediction[n] <- info$sift_prediction
      }
      if( "polyphen_score" %in% colnames(info)){
        df$polyphen.score[n] <- info$polyphen_score
      }
      if( "sift_score" %in% colnames(info)){
        df$sift.score[n] <- info$sift_score
      }
      if("consequence_terms" %in% colnames(info)){
        df$consequences.all[n] <- info$consequence_terms
      }
      if("impact" %in% colnames(info)){
        df$impact.all[n] <- info$impact
      }
    }
    i <- j + 1
  }

  df$cadd.score <- unlist(df$cadd.score)
  df$mcap.score <- unlist(df$mcap.score)
  df$polyphen.score <- as.numeric(unlist(lapply(df$polyphen.score, toString)))
  df$sift.score <- as.numeric(unlist(lapply(df$sift.score, toString)))
  df$impact.all <- unlist(lapply(df$impact.all, toString))
  df$consequences.all <- unlist(lapply(df$consequences.all, toString))
  df$sift.prediction <- unlist(lapply(df$sift.prediction, toString))
  df$polyphen.prediction <- unlist(lapply(df$polyphen.prediction, toString))

  return(df)
}

# functions to get the ancestry of each variant in a nice format
get_ancestors <- function(df, n){
  ancestors = ""
  African <- df$Allele.Count.African[n] - df$Homozygote.Count.African[n]
  Ashkenazi_Jewish <- df$Allele.Count.Ashkenazi.Jewish[n] - df$Homozygote.Count.Ashkenazi.Jewish[n]
  East_Asian <- df$Allele.Count.East.Asian[n] - df$Homozygote.Count.East.Asian[n]
  European_Finnish <- df$Allele.Count.European..Finnish.[n] - df$Homozygote.Count.European..Finnish.[n]
  European <- df$Allele.Count.European..Non.Finnish.[n] - df$Homozygote.Count.European..Non.Finnish.[n]
  Latino <- df$Allele.Count.Latino[n] - df$Homozygote.Count.Latino[n]
  South_Asian <- df$Allele.Count.South.Asian[n] - df$Homozygote.Count.South.Asian[n]
  Other <- df$Allele.Count.Other[n] - df$Homozygote.Count.Other[n]
  ancestors <- paste(write_ancestors("African", African),write_ancestors("Ashkenazi Jewish", Ashkenazi_Jewish),
                     write_ancestors("East Asian", East_Asian), write_ancestors("Latino", Latino),
                     write_ancestors("European Finnish", European_Finnish),write_ancestors("European", European),
                     write_ancestors("South Asian", South_Asian),write_ancestors("Other", Other), sep = "")
  return(substr(ancestors, 1, nchar(ancestors) - 2))
}

write_ancestors <- function(name, number){
  ancestor <- ""
  if(length(number) > 0 && number != 0){
    ancestor <- paste(name," (", number, "); ", sep = "")
  }
  return(ancestor)
}

sortScoredVariants <- function(df, cadd.min = 20, fitCon.min = .4, phastCon.min = .55,
                               mcap.min = .025, AFgnomAD.max = .01, sift.max = .05, polyphen.min = .909){
  for(n in 1:dim(df)[1]){
    pass <- 0
    fail <- 0
    na <- 0
    # SIFT
    if(is.na(df$sift.score[n])){
      na <- na + 1
    }
    else if(df$sift.score[n] <= sift.max){
      pass <- pass + 1
    }
    else{
      fail <- fail + 1
    }
    # PolyPhen
    if(is.na(df$polyphen.score[n])){
      na <- na + 1
    }
    else if(df$polyphen.score[n] >= polyphen.min){
      pass <- pass + 1
    }
    else{
      fail <- fail + 1
    }
    # gnomAD AF
    if(is.na(df$Allele.Frequency[n])){
      na <- na + 1
    }
    else if(df$Allele.Frequency[n] <= AFgnomAD.max){
      pass <- pass + 1
    }
    else{
      fail <- fail + 1
    }
    # M-CAP
    if(is.na(df$mcap.score[n])){
      na <- na + 1
    }
    else if(df$mcap.score[n] >= mcap.min){
      pass <- pass + 1
    }
    else{
      fail <- fail + 1
    }
    # phastCons
    if(is.na(df$phastCon.score[n])){
      na <- na + 1
    }
    else if(df$phastCon.score[n] >= phastCon.min){
      pass <- pass + 1
    }
    else{
      fail <- fail + 1
    }
    # fitCons
    if(is.na(df$fitCon.score[n])){
      na <- na + 1
    }
    else if(df$fitCon.score[n] >= fitCon.min){
      pass <- pass + 1
    }
    else{
      fail <- fail + 1
    }
    # CADD
    if(is.na(df$cadd.score[n])){
      na <- na + 1
    }
    else if(df$cadd.score[n] >= cadd.min){
      pass <- pass + 1
    }
    else{
      fail <- fail + 1
    }
    df$num.pass[n] <- pass
    df$num.na[n] <- na
    df$num.fail[n] <- fail
  }
  df <- df[ , !(names(df) %in% c("rest.api.url", "rest.api.info", "source", "distance.from.start", "temp",
                                                   "rest.api.transcript.consequences","Consequence", "Filters...exomes", "Filters...genomes"))]
return(df[order(-df$num.pass),])
}

plotScoreDistribution <- function(df, gene.symbol, save.plot = FALSE, output.dir){
  if(!exists("ensembl")){
    ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  }
  geneInfo <- getBM(attributes = c("start_position", "end_position"), filters = "hgnc_symbol",
                    values = gene.symbol, mart = ensembl)
  gene.image <- makeGene(id = gene.symbol, type = "hgnc_symbol", biomart = ensembl)
  genomeAxis <- makeGenomeAxis(add53 = TRUE, add35=TRUE)
  expres <- makeGenericArray(intensity = as.matrix(df$num.pass), probeStart = as.numeric(
    df$Position), dp = DisplayPars(type = "dot", lwd = 2, pch = "o")) # shows number of scores passed
  if(save.plot){
    jpeg(file = paste(output.dir, gene.symbol, "_variants", ".jpeg", sep = ""), width = 1000,
         height = 900)
    gdPlot(list(Exons=gene.image, "Number of Scores Passed"= expres, "BP" = genomeAxis),
           minBase = geneInfo$start_position, maxBase =geneInfo$end_position, labelCex = 1.5) # plots all 3 images
    dev.off()
  }
  gdPlot(list(Exons=gene.image, "Number of Scores Passed"= expres, "BP" = genomeAxis),
         minBase = geneInfo$start_position, maxBase =geneInfo$end_position, labelCex = 1.5) # plots all 3 images

}

