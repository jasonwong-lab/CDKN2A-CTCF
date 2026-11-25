args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript CDKN2B_parse_bam_readcount_ChIP.r INPUT.readcount OUTPUT.txt")
}

lines <- readLines(args[1], warn = FALSE)
if (length(lines) == 0) {
  stop("No read counts in ", args[1])
}

parse_line <- function(line) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(fields) < 5 || !grepl("^[0-9]+$", fields[2]) ||
      !grepl("^[0-9]+$", fields[4])) {
    stop("Invalid bam-readcount row in ", args[1], ": ", line)
  }
  bases <- strsplit(fields[-(1:4)], ":", fixed = TRUE)
  counts <- vapply(c("A", "C", "G", "T"), function(base) {
    hits <- Filter(function(x) length(x) >= 2 && x[1] == base, bases)
    if (length(hits) != 1 || !grepl("^[0-9]+$", hits[[1]][2])) {
      stop("Missing or invalid ", base, " count in ", args[1], ": ", line)
    }
    as.numeric(hits[[1]][2])
  }, numeric(1))

  data.frame(
    chr = fields[1], pos = as.numeric(fields[2]), ref_base = fields[3],
    Total_count = as.numeric(fields[4]),
    A = counts["A"], C = counts["C"], G = counts["G"], T = counts["T"],
    row.names = NULL
  )
}

result <- do.call(rbind, lapply(lines, parse_line))
write.table(result, file = args[2], sep = "\t", row.names = FALSE, quote = FALSE)
