adjust_test_positivity_data = function(test_positivity_data, model){
  #' function to adjust test positivity data with a gam model of smoothed test positivity rates
  #' 
  #' @param test_positivity_data dataframe of test positivity data aggregated at a geographic unit (e.g., states)
  #' @param model gam model to smooth test positivity rates
  #' @return adjusted_test_positivity_data dataframe of observed and adjusted test positivity data
  test_positivity_data[["smoothed_positivity_rate"]] = predict(model, newdata = test_positivity_data, type = "response")
  adjusted_test_positivity_data = test_positivity_data %>%
    mutate(adjusted_test_positivity_rate = positive_increase * sqrt(smoothed_positivity_rate)) %>%
    select(date, state_name, positive_increase, adjusted_test_positivity_rate) %>%
    gather(variable, value, -date, -state_name) %>%
    mutate(variable = if_else(variable == "adjusted_test_positivity_rate",
                              true = "Adjusted Test Positivity Rate",
                              false = "Observed Test Positivity Rate")) %>%
    group_by(state_name, variable) %>%
    arrange(date) %>%
    ungroup() %>%
    mutate(value = value / value[n()] * 100)

  return(adjusted_test_positivity_data)
}  

#' adjustment themes and label

adjustment_theme = theme(panel.grid.major.x = element_blank(),
                         panel.grid.minor.x = element_blank(),
                         panel.grid.minor.y = element_blank(),
                         legend.position = "bottom")

adjustment_labels = labs(x = "", 
                         colour = "", 
                         y = "New Daily Confirmed Cases",
                         subtitle = "",
                         caption = "rolling seven-day average, sourced from covidtracking.com")
