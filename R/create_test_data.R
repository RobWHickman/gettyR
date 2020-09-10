#' Create fake spike data for testing
#'
#' @example
#' test_spikes <- purrr::imap_dfr(create_test_data(), function(x, y) data.frame(spike_time = x, trial = y))
#'
#' @author Robert Hickman
#' @export create_test_data

create_test_data <- function(t1 = 1.2, t1str = 3, t2 = 3, t2str = 1.5, t3 = 4.5, t3str = 2, ttotal = 8, cells = 3, baseline_cell_rate = 1, reaction_time = 0.25, trials = 50) {
  trial_spikes <- c()
  for(t in seq(trials)) {
    all_ms <- ttotal * 1000

    spikes <- c()
    for(c in seq(cells)) {
      baseline_spike_rate <- (baseline_cell_rate + runif(min = -1, max = 2, 1)) * ttotal
      spikes <- append(spikes, sample(all_ms, ceiling(baseline_spike_rate)))

      epoch1_range <- (t1 * 1000):((t1 + reaction_time) * 1000)
      epoch1_spikes <- sample(epoch1_range, round(baseline_spike_rate * (t1str + runif(1, min = -1, max = 1))) * reaction_time)
      spikes <- append(spikes, epoch1_spikes)

      epoch2_range <- (t2 * 1000):((t2 + reaction_time) * 1000)
      epoch2_spikes <- sample(epoch1_range, round(baseline_spike_rate * (t2str + runif(1, min = -1, max = 1))) * reaction_time)
      spikes <- append(spikes, epoch1_spikes)

      epoch3_range <- (t3 * 1000):((t3 + reaction_time) * 1000)
      epoch3_spikes <- sample(epoch1_range, round(baseline_spike_rate * (t3str + runif(1, min = -1, max = 1))) * reaction_time)
      spikes <- append(spikes, epoch1_spikes)

    }
    trial_spikes[[t]] <- spikes
  }
  return(trial_spikes)
}

