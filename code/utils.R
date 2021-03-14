generate_smooth_data = function(state_df, model){
  state_df[["pos_rate_smoothed"]] = predict(model, newdata = state_df, type = "response")
  state_df = state_df %>%
  mutate(adj_pos = positive_increase * sqrt(pos_rate_smoothed)) %>%
  select(date, state_name, positive_increase, adj_pos) %>%
  gather(variable, value, -date, -state_name) %>%
  mutate(variable = if_else(variable == "adj_pos",
                            true = "Adjusted for test positivity rate",
                            false = "Original")) %>%
  group_by(state_name, variable) %>%
  arrange(date) %>%
  ungroup() %>%
  mutate(value = value / value[n()] * 100)
  return(state_df)
}  

#-----------------Common themes and labels-------

the_theme = theme(axis.text.y = element_blank(),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.minor.y = element_blank())

the_labs = labs(x = "", 
                 colour = "", 
                 y = "New daily confirmed cases",
                 caption = "Source: data from covidtracking.com, positivity adjustment by http://freerangestats.info")
