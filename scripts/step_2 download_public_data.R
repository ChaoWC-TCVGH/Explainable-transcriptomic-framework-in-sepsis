library(readr)
.script_dir <- "D:/Gene_analysis/test"
source(file.path(.script_dir, "step_1 define and dictionary.R"))

options(timeout = max(7200, getOption("timeout")))
man <- read_csv(PUBLIC_MANIFEST_CSV)

fetch_recount3 <- function(project, dest) {
  assert(requireNamespace("recount3", quietly = TRUE), "package recount3 is required for this step")
  projects <- recount3::available_projects(organism = "human")
  info <- projects[projects$project == project, , drop = FALSE]
  assert(nrow(info) >= 1L, "recount3 project not found: ", project)
  rse <- recount3::create_rse(info[1, , drop = FALSE], type = "gene")
  saveRDS(rse, dest)
}

check <- list()
for (i in seq_len(nrow(man))) {
  dest <- file.path(RAW_DIR, man$file[i])
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  comparable <- isTRUE(as.logical(man$sha256_comparable[i]))
  have <- file.exists(dest) && (!comparable || sha256_file(dest) == man$sha256[i])
  action <- "present"
  if (!have) {
    action <- "downloaded"
    cat(sprintf("   downloading %s\n", man$file[i]))
    if (grepl("^https://", man$source[i])) {
      utils::download.file(man$source[i], dest, mode = "wb", quiet = TRUE)
    } else {
      fetch_recount3(sub("^recount3_(.*)_rse\\.rds$", "\\1", basename(man$file[i])), dest)
    }
  }
  sha <- sha256_file(dest)
  status <- if (!comparable) "not comparable (serialised R object)" else
            if (sha == man$sha256[i]) "MATCH" else "MISMATCH"
  cat(sprintf("   %-62s %-10s %s\n", man$file[i], action, status))
  check[[i]] <- data.frame(dataset = man$dataset[i], file = man$file[i], action = action,
                           bytes = file.info(dest)$size, expected_bytes = man$bytes[i],
                           sha256 = sha, expected_sha256 = man$sha256[i], status = status,
                           stringsAsFactors = FALSE)
}
check <- do.call(rbind, check)
write_csv(check, file.path(TABLE_DIR, "public_data_check.csv"))
bad <- check$file[check$status == "MISMATCH"]
assert(!length(bad), "downloaded file differs from the file the authors analysed: ",
       paste(bad, collapse = ", "),
       " -- GEO may have revised it; results can differ. See README, section 'If a hash differs'.")

