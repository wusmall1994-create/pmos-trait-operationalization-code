root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION"))) stop("Run from repository root.")

all_files <- list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE)
relative <- substring(normalizePath(all_files, winslash = "/", mustWork = FALSE), nchar(root) + 2)
excluded <- grepl("^(\\.git|data/public|data/restricted|outputs)(/|$)", relative)
files <- all_files[!excluded & file.info(all_files)$isdir %in% FALSE]
relative_files <- relative[!excluded & file.info(all_files)$isdir %in% FALSE]

forbidden_extensions <- c("doc", "docx", "pdf", "png", "jpg", "jpeg", "tif", "tiff",
                          "xpt", "sav", "dta", "sas7bdat", "rds", "rdata", "xlsx")
extensions <- tolower(tools::file_ext(relative_files))
bad_extensions <- relative_files[extensions %in% forbidden_extensions]
if (length(bad_extensions)) stop("Forbidden release file types: ", paste(bad_extensions, collapse = ", "))

read_text <- function(path) {
  tryCatch(paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n"),
           error = function(e) "")
}
content <- setNames(vapply(files, read_text, character(1)), relative_files)
patterns <- c(
  manuscript_heading = "(^|\\n)#{1,3} (Abstract|Introduction|Discussion|Conclusion)(\\n|$)",
  github_token = "gh[pousr]_[A-Za-z0-9]{20,}",
  zenodo_token = "(?i)zenodo[_-]?(access[_-]?)?token\\s*[=:]\\s*[A-Za-z0-9._-]{12,}",
  private_key = "BEGIN (RSA |EC |OPENSSH )?PRIVATE KEY",
  email = "[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}",
  telephone = "\\+?86[ -]?1[3-9][0-9][ -]?[0-9]{4}[ -]?[0-9]{4}",
  windows_absolute_path = "[A-Za-z]:[\\\\/]Users[\\\\/]",
  local_codex_path = "[A-Za-z]:[\\\\/]Codex[\\\\/]"
)
findings <- unlist(lapply(names(patterns), function(label) {
  hit <- names(content)[vapply(content, grepl, logical(1), pattern = patterns[[label]], perl = TRUE)]
  if (length(hit)) paste(label, hit, sep = ": ") else character()
}))
if (length(findings)) stop("Release audit failed:\n", paste(findings, collapse = "\n"))

cat("Release audit passed:", length(files), "code/metadata files checked.\n")
