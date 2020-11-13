suppressPackageStartupMessages({
    library(RCurl)
    library(RColorBrewer)
    library(rWikiPathways)
    library(RCy3)
    library(optparse)
})

# inputs
option_list <- list(
    make_option(c("-p", "--wikiPathway"),
                type = "character",
                default = 'WP528',
                help = "wikiPathway Identtifier",
                metavar = "character"),
    make_option(c("-g", "--genes"),
                type = "character",
                help = "Path to differentially expressed genes file",
                metavar = "character"),
    make_option(c("-l", "--header"),
                type = "logical",
                help = "File has header line?",
                metavar = "logical")
    )

opt_parser <- OptionParser(usage = "%prog [options] file",
                           option_list = option_list)
opt <- parse_args(opt_parser)

wp_id <- opt$wikiPathway
headers <- opt$header
input_data <- opt$genes


#fetch and modify data
dat <- read.table(input_data, header = headers, sep="\t", fill=TRUE)
if(headers)
{
  colnames(dat)[1] <- "geneid"
  colnames(dat)[2] <- "fc"
} else
{
  colnames(dat) <- c("geneid", "fc", "p-value", "adj p-value")
}

dat <- data.frame(lapply(dat, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))
min.fc = min(dat["fc"],na.rm=TRUE)
max.fc = max(dat["fc"],na.rm=TRUE)
abs.fc = max(abs(min.fc),abs(max.fc))
data.values = c(-abs.fc,0,abs.fc)

#open in browser
rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
fill_hex <- function(x) {
  out <- vector(mode="character", length=length(x))
  for (i in seq_along(x))
  {
    if(is.na(x[[i]]))
    {
      out[[i]] <- "#ffffff"
    }
    else
    {
      rel_value = 70+(185*( 1 - (abs(x[[i]])/abs.fc)))
      if(x[[i]] < 0)
      {
        if(abs(min.fc) == 0)
        {
          out[[i]] <- "#ffffff"
        }
        else{
          out[[i]] <- rgb2hex(rel_value, rel_value, 255)
        }
      }
      else
      {
        if(abs(min.fc) == 0)
        {
          out[[i]] <- "#ffffff"
        }
        else{
          out[[i]] <- rgb2hex(255, rel_value, rel_value)
        }
      }
    }


  }
  return(out)
}
dat$hex <- fill_hex(dat[["fc"]])
url_params <- paste(with(dat, paste(substr(hex,2,nchar(hex)), geneid, sep="=HGNC_")), collapse="&")
url <- paste("https://pathway-viewer.toolforge.org/?id=", wp_id, "&", url_params, sep="")
download.file(url, destfile='pathwayview.html', method="libcurl")

