#!/usr/bin/env Rscript
# Map and pathway enrichment analysis.

library("jsonlite")
library("optparse")

# parse options
option_list <- list(
  make_option(c("-o", "--output_json"),
              type = "character",
              default = "wikipathways.json",
              help = "An output JSON of fetched data [default= %default]",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Massive SPARQL query:
# https://www.wikipathways.org/index.php/Help:WikiPathways_Sparql_queries#List_
# of_WikiPathways_for_HGNC_symbols
# with the Homo sapiens condition added

sparql_url <- paste0("http://sparql.wikipathways.org/sparql?&query=select+",
                     "distinct+%3FpathwayRes+str%28%3Fwpid%29+as+%3Fpathway+",
                     "str%28%3Ftitle%29+as+%3FpathwayTitle+",
                     "fn%3Asubstring%28%3FhgncId%2C36%29+as+%3FHGNC+",
                     "where+%7B%0D%0A",
                     "++%3Fgene+a+wp%3AGeneProduct+%3B%0D%0A",
                     "++++dcterms%3Aidentifier+%3Fid+%3B%0D%0A",
                     "++++dcterms%3AisPartOf+%3FpathwayRes+%3B%0D%0A",
                     "++++wp%3AbdbHgncSymbol+%3FhgncId+.%0D%0A",
                     "++%3FpathwayRes+a+wp%3APathway+%3B%0D%0A",
                     "++++wp%3AorganismName+",
                     "%22Homo%20sapiens%22%5E%5Exsd%3Astring+%3B%0D%0A",
                     "++++dcterms%3Aidentifier+%3Fwpid+%3B%0D%0A",
                     "++++dc%3Atitle+%3Ftitle+.%0D%0A%7D",
                     "&format=text%2Fcsv&timeout=0&debug=on")

# Retrieve the table whole opening/closing the connection gracefully
con <- url(sparql_url)
res <- try(read.table(con, sep = ",", header = T, stringsAsFactors = F),
           silent = T)
if (class(res) == "try-error") {
  message(paste0("Cannot read from ", sparql_url))
  message(res)
  close(con)
}

# Compile a similar JSON with the WikiPathways pathways

# Melt the table by the unique pathway name
wp_list <- lapply(unique(res$pathwayTitle),
                  function(x) list(name = x,
                                   hgncs = unique(res[res$pathwayTitle == x,
                                                      "HGNC"])))

message("Writing WikiPathways to JSON...")
cat(toJSON(wp_list), file = opt$output_json)

message("Done.")
