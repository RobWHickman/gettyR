#' Create directory structure for plots from a cell sorting
#' @param session The recording session ID (date and session) for the directory to live in
#' @param cells The names of the sorted cell clusters (can be set to NULL to ignore)
#' @param dir The base directory where the session are found. Defaults to Dan's lab laptop
#'
#' @author Robert Hickman
#' @export create_spike_dirs

create_spike_dirs <- function(session, cells, dir = "C:/Users/DHill/Desktop/robert_cell_sorting/sessions") {
  save_dir_root <- file.path(dir, session, "figures")

  #create figs directory
  if(!dir.exists(file.path(save_dir_root))) {
    dir.create(save_dir_root)
  }

  #create dir for unsorted getty response figures
  if(!dir.exists(file.path(save_dir_root, "unsorted_responses"))) {
    dir.create(file.path(save_dir_root, "unsorted_responses"))
  }

  #create separate dir for each cell cluster
  if(!is.null(cells)) {
    for(cell in cells) {
      cell_dir <- file.path(save_dir_root, cell)
      if(!dir.exists(cell_dir)) {
        dir.create(cell_dir)
      }
    }
  }

  message("directories created successfully")
}
