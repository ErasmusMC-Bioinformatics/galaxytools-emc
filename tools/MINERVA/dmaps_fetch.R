#!/usr/bin/env Rscript
### Map and pathway enrichment analysis.

library('jsonlite')
library('httr')
library('optparse')


# parse options
option_list = list(
  make_option(c("-i", "--input_config"), type="character", default=NULL, 
              help="A comma seperated file listing the resources to fetch", metavar="character"),
  make_option(c("-o", "--output_json"), type="character", default="dmaps_out.json", 
              help="An output JSON of fetched data [default= %default]", metavar="character")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input_config)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Check parameter values

if (!file.exists(opt$input_config)){
  stop(paste('File', opt$input_config, 'does not exist.n'), call.=FALSE)
}

### A convenience function to fetch data from MINERVA
ask_GET <- function(furl, fask) {
  resp <- httr::GET(url = paste0(furl, fask),
                    httr::add_headers('Content-Type' = "application/x-www-form-urlencoded"))
  if(httr::status_code(resp) == 200) {
    return(httr::content(resp, as = "text"))
  }
  return(NULL)
}

### In case we want to pass parameters as an argument
args = commandArgs(trailingOnly=TRUE)

### Config listing the resources, which we want to fetch
config <- opt$input_config
suppressWarnings(config <- read.table(config, sep = ",", stringsAsFactors = F, header = T))

### for a given model id (for all maps/submaps in the project)
### Fetch HGNC symbols for proteins, RNAs and genes.
get_hgncs <- function(fmodelid, all_models) {
  ### MINERVA call fetching elements with relevant parameters
  fmodel <- fromJSON(ask_GET(paste0(mnv_base,"models/",fmodelid,"/"), "bioEntities/elements/?columns=id,name,type,references,bounds"), flatten = F)
  
  ### Model name used in naming gene sets
  model_name <- all_models[all_models$idObject == fmodelid,"name"]
  
  ### If no annotations, break
  if(all(is.null(unlist(fmodel$references)))) { return(NULL) }
  
  ### List all HGNCs in maps and submaps
  ret <- list(
    list(name = model_name, 
         hgncs = unlist(sapply(fmodel$references[fmodel$type %in% c("Protein", "RNA", "Gene")], function(x) x[x$type == "HGNC_SYMBOL", "resource"]))))
  
  ### List all HGNCs in pathways (smaller areas within a single diagram)
  if(length(unique(fmodel[fmodel$type == "Pathway","name"])) > 0) {
    ### Calculating the inclusion of the elements per bioentity set
    for(pname in unique(me[fmodel$type == "Pathway","name"])) {
      res <- apply(me[fmodel$name == pname,"bounds"], 1, 
                   function(x) fmodel$bounds$x >= x["x"] & fmodel$bounds$x <= (x["x"] +  x["width"]) & 
                     fmodel$bounds$y >= x["y"] & fmodel$bounds$y <= (x["y"] +  x["height"]))
      
      hgncs <- unlist(sapply(fmodel$references[res & fmodel$type %in% c("Protein", "RNA", "Gene")], 
                             function(x) x[x$type == "HGNC_SYMBOL", "resource"]))
      ret <- c(ret, list(name = paste0(model_name,":",pname), hgncs = hgncs))
    }
  }
  return(ret)
}

message("Fetching contents of disease maps in config")

all_maps <- list()

for(map in config[config$type == "map","resource"]) {
  message(paste0("Querying ", map))
  cfg <- ask_GET(map, "configuration/")
  ### For erroneous response, skip to next map
  if(is.null(cfg)) { next }
  cfg <- fromJSON(cfg)
  project_id <- cfg$options[cfg$options$type == "DEFAULT_MAP","value"]
  
  mnv_base <- paste0(map,"projects/",project_id,"/")
  
  ### Ask for models
  message(paste0("Asking for models: ", mnv_base, "models/"))
  models <- ask_GET(mnv_base, "models/")
  ### For erroneous response, skip to next map
  if(is.null(models)) { next }
  models <- fromJSON(models, flatten = F)
  
  message(paste0("Fetching HGNCs for individual models..."))
  
  model_hgncs <- lapply(models$idObject, get_hgncs, models)
  
  #Remove null entries
  model_hgncs <- model_hgncs[!sapply(model_hgncs, is.null)]
  
  all_maps <- c(all_maps, list(map_name = map, map_elements = model_hgncs))
  
  message("Done.")
}

message("Writing disease map(s) to JSON...")
cat(toJSON(all_maps), file = opt$output_json)

