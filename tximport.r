# params are:  
# 1 = gene2tx_table
# 2 = out_file
# 3..n = sample files 

library(tximport)
args = commandArgs(trailingOnly=TRUE)
#length(args)
if (length(args) < 3) {
  stop("At least 3 arguments must be supplied: gene2tx-Table, output-File, (samples)xN", call.=FALSE)
}

# tx2gene links transcript IDs to gene IDs for summarization
tx2gene <- read.table(file.path(args[1]),header=FALSE)

# parse sample list
samples <- character()
for (v in 3:length(args))
  samples <- c(samples, args[v])

names(samples) <- paste0("sample", 1:(length(args)-2))
names(tx2gene) <- c('tx','gene')
out <- tximport(samples, type="salmon", tx2gene=tx2gene)

# write count as table
write.table(out$counts, file=args[2], row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")


