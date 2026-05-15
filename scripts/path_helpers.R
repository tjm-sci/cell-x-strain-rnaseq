# Resolve project-local paths from the repository root while preserving absolute
# paths supplied on the command line.
resolve_project_path <- function(path, must_work = FALSE) {
  expanded_path <- path.expand(path)
  is_absolute <- grepl("^/", expanded_path) || grepl("^[A-Za-z]:[/\\\\]", expanded_path)

  if (is_absolute) {
    return(normalizePath(expanded_path, winslash = "/", mustWork = must_work))
  }

  normalizePath(here::here(expanded_path), winslash = "/", mustWork = must_work)
}
