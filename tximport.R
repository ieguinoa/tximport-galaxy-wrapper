# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("getopt")
#library("tools")
options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "quiet", "q", 0, "logical",
  "help", "h", 0, "logical",
  "base_dir", "w", 1, "character",
  "out_file", "o", 1, "character",
  "countsFiles", "n", 1, "character",
  "countsFromAbundance", "r", 1, "character",
  "format", "v", 1, "character",
  "gff_file", "H", 0, "logical",
  "tx2gene", "f", 0, "character",
  "geneIdCol", "l", 0, "character",
  "txIdCol" , "p", 1, "character",
  "abundanceColt", "i", 0, "logical",
  "countsCol", "y", 1, "character",
  "lengthCol", "x", 1, "character"),
  byrow=TRUE, ncol=4)

opt <- getopt(spec)

countsFiles

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}
	



library(tximport)
args = commandArgs(trailingOnly=TRUE)
#length(args)
#if (length(args) < 3) {
#  stop("At least 3 arguments must be supplied: gene2tx-Table, output-File, (samples)xN", call.=FALSE)
#}




### if the input is a gff/gtf file first need to create the tx2gene table
if (!is.null(gff_file)) {
    suppressPackageStartupMessages({
        library("GenomicFeatures")
    })
    txdb <- makeTxDbFromGFF(gff_file)
    k <- keys(txdb, keytype = "TXNAME")
    tx2gene <- select(txdb, keys=k, columns="GENEID", keytype="TXNAME")
    # Remove 'transcript:' from transcript IDs (when gffFile is a GFF3 from Ensembl and the transcript does not have a Name)
    tx2gene$TXNAME <- sub('^transcript:', '', tx2gene$TXNAME)

} else if (!is.null(tx2gene)){
    tx2gene_table <- read.table(tx2gene,header=FALSE)

# parse sample list
samples <- character()
for (v in 3:length(countsFiles))
  samples <- c(samples, args[v])

names(samples) <- paste0("sample", 1:(length(args)-2))
names(tx2gene) <- c('tx','gene')
out <- tximport(samples, type="salmon", tx2gene=tx2gene)

# write count as table
write.table(out$counts, file=out_file, row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")


