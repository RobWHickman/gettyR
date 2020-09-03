#' Plot the trace data from getty
#' @param trace_data The trace data exported from getty via Spike2. A df of 3 variables an n rows
#' @param hz The frequency of the getty sampling of neuron data. Should be set at 22kHz
#' @param ci The confidence interval to plot on the error of the neuron shape. Defaults to 1.96
#'
#' @author Robert Hickman
#' @export plot_cell_trace

plot_cell_trace <- function(trace_data, hz = 22000, ci = 1.96) {
  trace <- trace_data %>%
    group_by(cell, time) %>%
    mutate(time = time * 1000 / hz) %>%
    summarise(m = mean(value),
              sd = sd(value),
              n = n()) %>%
    mutate(cell = paste0(cell, ": ", n, " spikes"))

  trace_plot <- ggplot(trace, aes(x = time, y = m)) +
    geom_ribbon(aes(ymin = m - ci*sd, ymax = m + ci*sd), alpha = 0.5) +
    geom_line(aes(colour = cell)) +
    scale_colour_discrete(guide = FALSE) +
    labs(title = "cell traces",
         x = "time (/ms)",
         y = "voltage change") +
    theme_minimal() +
    facet_wrap(~cell)

  return(trace_plot)
}

#' Plot the ISI from getty
#' @param spikes The sorted spike nested data column
#'
#' @author Robert Hickman
#' @export plot_isi_histogram

plot_isi_histogram <- function(spikes, bindwith = 30) {
  isi_data <- do.call(rbind, spikes) %>%
    arrange(spike_time_ms) %>%
    group_by(cell) %>%
    mutate(isi = lead(spike_time_ms) - spike_time_ms) %>%
    filter(!is.na(cell) & !is.na(isi) & isi < 1000)

  isi_histogram <- ggplot(isi_data, aes(x = isi)) +
    geom_histogram(aes(fill = cell), alpha = 0.8) +
    scale_fill_discrete(guide = FALSE) +
    labs(title = "interspike interval histogram",
         x = "ISI (/ms)",
         y = "count") +
    theme_minimal() +
    facet_wrap(~cell, scales = "free_y")

  return(isi_histogram)
}


#' Plot the unsorted responses from getty
#' @param data A df of data on spikes. Should include the situation, bits and spikes
#'
#' @author Robert Hickman
#' @export plot_getty_responses

plot_getty_responses <- function(data) {
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

plot_getty_responses <- function(data, situation, bit) {

}
