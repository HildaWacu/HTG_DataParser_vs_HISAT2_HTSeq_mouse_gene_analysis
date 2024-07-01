#Install curl if missing
if(!require(curl, quietly = TRUE)){
  install.packages("curl")
  library(curl)
}

# Parameters
#   accession: ArrayExpress/BioStudies accession of dataset to download.
#   path: File path to download folder (default is work directory).
getBio <- function(accession, path = "."){
  if(!dir.exists(path)) stop("Path '", path, "' does not exist")
  zipfile <- paste0(path, "/", accession, ".zip")
  curl_download(
    url = paste0("https://www.ebi.ac.uk/biostudies/files/", accession),
    destfile = zipfile)
  unzip(zipfile = zipfile, exdir = path)
  file.remove(zipfile)
  return(list.files(path = path, full.names = TRUE))
}
