#' Plot the trace data of sorted spike trace data
#' @param trace_data The trace data exported from getty via Spike2. A df of 3 variables an n rows
#' @param hz The frequency of the getty sampling of neuron data. Should be set at 22kHz
#' @param ci The confidence interval to plot on the error of the neuron shape. Defaults to 1.96
#' @param specific_cell A specific cell to plot or to highlight all cells. Defaults to NULL
#'
#' @author Robert Hickman
#' @export plot_cell_trace

plot_cell_trace <- function(trace_data, cell, hz = 22000, ci = 1.96, specific_cell = NULL) {
  if(!is.null(specific_cell)) trace_data <- dplyr::filter(trace_data, cell == specific_cell)

  trace <- trace_data %>%
    group_by(cell, time) %>%
    mutate(time = time * 1000 / hz) %>%
    summarise(mu = mean(value),
              sd = sd(value),
              n = n()) %>%
    mutate(cell_title = paste0(cell, ": ", n, " spikes"))

  trace_plot <- ggplot(trace, aes(x = time, y = mu)) +
    geom_ribbon(aes(ymin = mu - ci*sd, ymax = mu + ci*sd), alpha = 0.5, fill = "grey80") +
    geom_line(aes(colour = cell), size = 2) +
    scale_colour_discrete(guide = FALSE) +
    labs(
      title = "cell trace",
      x = "time (/ms)",
      y = "voltage change"
    ) +
    theme_minimal() +
    facet_wrap(~cell_title)

  return(trace_plot)
}

#' Plot the ISI of sorted spike trace data
#' @param spikes The sorted spike nested data column
#' @param binwidth The binwidth of the histogram produced
#' @param specific_cell A specific cell to plot or to highlight all cells. Defaults to NULL
#'
#' @author Robert Hickman
#' @export plot_isi_histogram

plot_isi_histogram <- function(spikes, bindwith = 30, specific_cell = NULL) {
  isi_data <- do.call(rbind, spikes) %>%
    arrange(spike_time_ms) %>%
    group_by(cell) %>%
    mutate(isi = lead(spike_time_ms) - spike_time_ms) %>%
    filter(!is.na(cell) & !is.na(isi) & isi < 1000)

  if(!is.null(specific_cell)) isi_data <- dplyr::filter(isi_data, cell == specific_cell)

  isi_histogram <- ggplot(isi_data, aes(x = isi)) +
    geom_histogram(aes(fill = cell), alpha = 0.8) +
    scale_fill_discrete(guide = FALSE) +
    labs(
      title = "interspike interval histogram",
      x = "isi (/ms)",
      y = "count"
    ) +
    theme_minimal() +
  facet_wrap(~cell, scales = "free_y")

  return(isi_histogram)
}

#' Plot the autocorrelelogram of sorted spike trace data
#' @param spikes The sorted spike nested data column
#'
#' @author Robert Hickman
#' @export plot_autocorr

#' Plot the PCA of sorted spike trace data
#' @param trace_data The trace data exported from getty via Spike2. A df of 3 variables an n rows
#' @param specific_cell A specific cell to plot or to highlight all cells. Defaults to NULL
#' @param xaxis Which principal component to plot on the x axis. PC1-5
#' @param yaxis Which principal component to plot on the y axis. PC1-5
#'
#' @author Robert Hickman
#' @export plot_pca

plot_pca <- function(trace_data, specific_cell = NULL, xaxis = "PC1", yaxis = "PC2") {
  if(!xaxis %in% paste0("PC", 1:5) | !yaxis %in% paste0("PC", 1:5)) {
    errorCondition("axes must be PC1-5")
  }

  pca_data <- trace_data %>%
    dplyr::mutate(time = paste0("t", time)) %>%
    tidyr::pivot_wider(id_cols = c(cell, spike_no), names_from = time, values_from = value) %>%
    recipe(~., data = .) %>%
    update_role(cell, spike_no, new_role = "id") %>%
    recipes::step_center(all_predictors()) %>%
    recipes::step_scale(all_predictors()) %>%
    recipes::step_pca(all_predictors()) %>%
    recipes::prep() %>%
    recipes::juice() %>%
    dplyr::rename(x_dim = !!xaxis, ydim = !!yaxis)

  pca_plot <- ggplot(pca_data, aes(x = x_dim, y = ydim)) +
    geom_point(data = dplyr::select(pca_data, -cell), colour = "grey80", alpha = 0.65) +
    theme_minimal() +
    labs(
      title = "pca of cell trace",
      x = xaxis,
      y = yaxis
    )

  if(is.null(specific_cell)) {
    pca_plot <- pca_plot +
      geom_point(aes(colour = cell)) +
      facet_wrap(~cell)
  } else {
    pca_plot <- pca_plot +
      geom_point(data = dplyr::filter(pca_data, cell == specific_cell), colour = "dodgerblue")
  }
}


#' Plot the unsorted responses from getty
#' @param data A df of data on spikes. Should include the situation, bits and spikes
#'
#' @author Robert Hickman
#' @export plot_getty_responses

plot_getty_responses <- function(data) {
  warning("function is deprecated and will be removed")
  getty_spike_responses <- data %>%
    select(trial, trial_start_time_ms, situation, bits, spikes) %>%
    unnest(bits) %>%
    filter(index == "upat") %>%
    unnest(spikes) %>%
    select(-spike)

  free_reward_getty <- getty_spike_responses %>%
    filter(situation == 4 & names == "budget_tap") %>%
    mutate(post_event_time_ms = trial_spike_time_ms - vals) %>%
    filter(post_event_time_ms < 2000 & post_event_time_ms > -1000) %>%
    ggplot() +
    geom_histogram(aes(x = post_event_time_ms), binwidth = 20) +
    geom_vline(xintercept = 0, colour = "red", linetype = "dashed", size = 2) +
    labs(title = "getty spike responses to Free Reward",
         x = "post event time (/ms)",
         y = "spike count") +
    theme_minimal()

  fix_cross_getty <- getty_spike_responses %>%
    filter(situation == 2 & names == "fixation_cross") %>%
    mutate(post_event_time_ms = trial_spike_time_ms - vals) %>%
    filter(post_event_time_ms < 2000 & post_event_time_ms > -1000) %>%
    ggplot() +
    geom_histogram(aes(x = post_event_time_ms), binwidth = 20) +
    geom_vline(xintercept = 0, colour = "red", linetype = "dashed", size = 2) +
    labs(title = "getty spike responses to Fixation Cross",
         x = "post event time (/ms)",
         y = "spike count") +
    theme_minimal()

  win_lose_getty <- getty_spike_responses %>%
    filter(situation == 2 & names == "win_lose") %>%
    mutate(post_event_time_ms = trial_spike_time_ms - vals) %>%
    filter(post_event_time_ms < 2000 & post_event_time_ms > -1000) %>%
    ggplot() +
    geom_histogram(aes(x = post_event_time_ms), binwidth = 20) +
    geom_vline(xintercept = 0, colour = "red", linetype = "dashed", size = 2) +
    labs(title = "getty spike responses to Winning Bids",
         x = "post event time (/ms)",
         y = "spike count") +
    theme_minimal()

  plot <- gridExtra::grid.arrange(
    free_reward_getty,
    fix_cross_getty,
    win_lose_getty,
    nrow = 3
  )
  return(plot)
}

#' Plot the raster and firing rate of the unsorted getty data
#' @param data A df of data on spikes. Should include the situation, bits and spikes
#'
#' @author Robert Hickman
#' @export plot_getty_responses2

plot_getty_responses <- function(data, trial_situations, trial_bits, back_window = -1000, front_window = 2000, binwidth = 20, raster_width = 5) {
  situation_trials <- data %>%
    dplyr::filter(situation == trial_situations) %>%
    dplyr::select(trial, trial_start_time_ms, situation, bits, trial_vals, spikes) %>%
    tidyr::unnest(bits) %>%
    dplyr::filter(index == "upat" & names %in% c(trial_bits, "error")) %>%
    dplyr::mutate(names = case_when(
      names == "error" ~ "error",
      TRUE ~ "interest_bit_time"
    )) %>%
    tidyr::pivot_wider(names_from = "names", values_from = "vals") %>%
    dplyr::filter(is.na(error)) %>%
    tidyr::unnest(spikes) %>%
    dplyr::select(-spike, -error) %>%
    dplyr::mutate(post_event_time_ms = trial_spike_time_ms - interest_bit_time) %>%
    dplyr::filter(post_event_time_ms < front_window & post_event_time_ms > back_window)

  if(trial_bits == "win_lose") {
    unique_trials <- situation_trials %>%
      dplyr::filter(!duplicated(trial))
    trial_data <- map2_df(
      unique_trials$trial,
      unique_trials$trial_vals,
      function(i, v) {
        if(length(v) > 1) {
          x <- v %>%
            dplyr::filter(names %in% c("monkey_bid", "computer_bid")) %>%
            tidyr::pivot_wider(names_from = names, values_from = vals) %>%
            mutate(trial = i)
          return(x)
        } else {
          return(data.frame(computer_bid = NA, monkey_bid = NA, trial = i))
        }
      }) %>%
      dplyr::mutate(win = case_when(
        monkey_bid > computer_bid ~ 1,
        TRUE ~ 0
      ))
    situation_trials <- left_join(situation_trials, trial_data, by = "trial") %>%
      dplyr::filter(win == 1)
  }

  raster_plot <- situation_trials %>%
    ggplot() +
    geom_vline(xintercept = 0, colour = "red", linetype = "dotted") +
    geom_rect(aes(xmin = post_event_time_ms, xmax = post_event_time_ms + raster_width, ymin = trial, ymax = trial + 1)) +
    scale_x_continuous(breaks = seq(back_window, front_window, 500)) +
    labs(
      title = paste("unsorted nba responses to", trial_bits),
      y = "trial count") +
    theme_minimal()

  firing_rate <- situation_trials %>%
    dplyr::mutate(bin_time = cut(post_event_time_ms, seq(back_window, front_window, binwidth))) %>%
    dplyr::mutate(bin_time = as.numeric(gsub("(^\\()(.*)(,.*$)", "\\2", as.character(bin_time))) + binwidth/2) %>%
    dplyr::group_by(trial, bin_time) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete(trial, bin_time, fill = list(n = 0)) %>%
    dplyr::group_by(bin_time) %>%
    dplyr::summarise(mean_spikes = mean(n),
                     sem = sd(n) / sqrt(n())) %>%
    dplyr::mutate_at(c("mean_spikes", "sem"), ~. * (1000 / binwidth))

  fr_plot <- firing_rate %>%
    ggplot(aes(x = bin_time, y = mean_spikes)) +
    geom_vline(xintercept = 0, colour = "red", linetype = "dotted") +
    geom_ribbon(aes(ymin = mean_spikes - sem, ymax = mean_spikes + sem), fill = "grey80", colour = "black", linetype = "dotted") +
    geom_line(size = 2) +
    scale_x_continuous(breaks = seq(back_window, front_window, 500), "time (ms)") +
    labs(y = "firing rate (Hz)") +
    theme_minimal()

  plot <- raster_plot / fr_plot

  return(plot)

}
