
# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check that input and output files are supplied
if (length(args) != 2) {
	stop("Please supply input and output file paths.", call.=TRUE)
}

print(paste0("Input path: ", args[1]))
print(paste0("Output path: ", args[2]))


# Load in Biostrings library
library(Biostrings)

# Assign a number to each unique sequence
seq <- readDNAStringSet(args[1])
seq <- unique(seq)
seq <- setNames(seq, paste0(seq_along(seq), "_", names(seq)))

# Create a new fasta file with numbered sequences
writeXStringSet(seq, args[2])
write.table(col.names = FALSE, row.names = FALSE, data.frame(names(seq)), paste0(args[2],".names.txt"))