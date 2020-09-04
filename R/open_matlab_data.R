open_mat <- function(file) {
  x <- R.matlab::readMat(file)

  #get trial data
  data <- x$savefile
  trial_data_row <- which(rownames(data) == "trial")
  all_trial_data <- data[[trial_data_row]]
  #hack to find out number of trials
  n_trials <- length(all_trial_data) / 12

  for (i in seq(n_trials)) {
    single_trial_data <- all_trial_data[, , i]

    duration <-
      as.numeric(single_trial_data[[which(names(single_trial_data) == "duration")]])
    situation <-
      as.numeric(single_trial_data[[which(names(single_trial_data) == "situation")]])

    addvals <- data.frame(names = name_params("addvals"),
                          vals = single_trial_data[[which(names(single_trial_data) == "addvals")]])
    trial_addvals <- addvals[!grepl("^prv_", addvals$names),]
    prv_trial_addvals <- addvals[grepl("^prv_", addvals$names),]

    spike_matrix <- t(single_trial_data[[which(names(single_trial_data) == "neuron")]])
    if(length(spike_matrix) == 0 ) { spike_matrix <- matrix(NA)}

    spikes <-
      data.frame(trial_spike_time_ms = spike_matrix,
                 spike = 1)

    bits <- data.frame(
      names = rep(name_params("bits"), each = 3),
      index = rep(c("initial", "upat", "downat"), length(name_params("bits"))),
      vals = as.matrix(single_trial_data[[which(names(single_trial_data) == "bit")]])
    )

    null_indexes <- !(sapply(bits$vals, length))
    bits$vals[null_indexes] <- NA
    expanded_bits <- bits[rep(row.names(bits), lengths(bits$vals)), seq(ncol(bits))]

    expanded_bits$vals <- unlist(bits$vals)

    trial_df <- dplyr::tibble(
      trial = i,
      duration,
      situation,
      spikes = list(spikes),
      trial_addvals = list(trial_addvals),
      prv_trial_addvals = list(prv_trial_addvals),
      bits = list(expanded_bits)
    )

    if (i == 1) {
      session_df <- trial_df
    } else {
      session_df <- rbind(session_df, trial_df)
    }
  }
  return(session_df)
}

open_sorted_spikes <- function(file) {
  x <- R.matlab::readMat(file)

  data <- x[[names(x)]]

  times <- data[[which(rownames(data) == "times")]]
  codes <- data[[which(rownames(data) == "codes")]][,1]

  coded_spikes <- data.frame(time = times*1000,
                             cell = codes)
  return(coded_spikes)
}

open_cell_trace <- function(dir, sorted_file = "sorted_spikes.mat") {
  message(paste("looking for", sorted_file, "in", dir))

  x <- R.matlab::readMat(file.path(dir, session, sorted_file))

  data <- x[[names(x)]]
  traces <- data[[which(rownames(data) == "values")]]
  codes <- data[[which(rownames(data) == "codes")]][,1]
  df <- data.frame(traces) %>%
    mutate(cell = paste0("cell", codes)) %>%
    pivot_longer(cols = starts_with("X")) %>%
    mutate(time = as.numeric(gsub("^X", "", name))) %>%
    select(-name)

  return(df)
}

open_spike_data <- function(dir, session, sorted_file = "sorted_spikes.mat") {
  message(paste("looking for .mat session file and", sorted_file, "in", dir))

  data <- gettyR::open_mat(file.path(dir, session, dir(file.path(dir, session)) %>% .[grepl("[0-9]\\.mat", .)])) %>%
    mutate(cumulative_duration_ms = cumsum(duration),
           trial_start_time_ms = lag(cumulative_duration_ms, default = 0),
           trial_vals = lead(prv_trial_addvals)) %>%
    select(-cumulative_duration_ms, -prv_trial_addvals)
  #change the names from previous trial to the trial it is on now
  data$trial_vals <- lapply(data$trial_vals, function(x) {
    x$names <- gsub("^prv_", "", x$names)
    return(x)
  })

  sorted_spikes <- gettyR::open_sorted_spikes(file.path(dir, session, sorted_file)) %>%
    mutate(cell = paste0("cell", cell),
           trial = map_dbl(time, function(s) last(which(data$trial_start_time_ms < s)))) %>%
    rename(spike_time_ms = time)

  data <- left_join(data, sorted_spikes, by = "trial") %>%
    mutate(trial_spike_time_ms = spike_time_ms - trial_start_time_ms) %>%
    tidyr::nest(sorted_spikes = c(spike_time_ms, trial_spike_time_ms, cell))

  return(data)
}

