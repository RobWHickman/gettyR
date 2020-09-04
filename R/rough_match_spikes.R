#' Match sorted and unsorted spike times to remove collection lag
#' @param trial The trial number of a trial in a dataset
#' @param spikes The unsorted getty spike times
#' @param sorted_spikes The sorted Spike2 spike times
#'
#' @author Robert Hickman
#' @export rough_match_spikes

rough_match_spikes <- function(trial, spikes, sorted_spikes) {
  warning("function is deprecated with the insertion of the rad trial data and should not be used")
  unsorted_spike_times <- spikes$trial_spike_time_ms
  sorted_spike_times <- sorted_spikes$trial_spike_time_ms

  if(is.na(sorted_spike_times)) {
    return(data.frame(trial, lag = 0, Freq = 0))
  } else {

    vec <- c()
    for(i in seq(length(sorted_spike_times))) {
      v <- sorted_spike_times[i] - unsorted_spike_times
      vec <- c(vec, v)
    }

    df <- data.frame(table(round(vec))) %>%
      arrange(-Freq) %>%
      .[1,] %>%
      mutate(trial) %>%
      mutate(lag = as.numeric(as.character(Var1))) %>%
      select(trial, lag, Freq)

    return(df)
  }
}
