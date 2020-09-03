#function to get the names of the parameters in the getty recording
name_params <- function(parameters) {
  if(parameters == "addvals") {
    #the addval names
    names <- c(
      "number_of_addvals",
      "trial_no",
      "trial_no2",
      "getty_setting_trial_duration",
      "getty_setting_trial_duration2",
      "situation",
      "getty_task",
      "getty_subtask",
      "fractal_value",
      "budget_magnitude",
      "trial_start_bid",
      "trial_computer_bid",
      "prv_budget_liquid",
      "prv_budget_magnitude",
      "prv_computer_bid",
      "prv_correct",
      "prv_free_reward",
      "prv_monkey_bid",
      "prv_reward_chance",
      "prv_reward_liquid",
      "prv_reward_magnitude"
    )
    #the bit names
  } else if(parameters == "bits") {
    names <- c(
      "trial_onset",
      "fixation_cross",
      "fractal_display",
      "bid_start",
      "bid_stable",
      "win_lose",
      "reward_epoch_end",
      "budget_epoch_end",
      "free_reward",
      "reward_tap",
      "budget_tap",
      "third_tap",
      "trial_end",
      "error"
    )
  }
  return(names)
}
