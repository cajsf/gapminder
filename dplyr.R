library(readr)

library(dplyr)

library(ggplot2)

gapminder_df <- readr::read_csv('data/gapminder.csv')

gapminder_df |> glimpse()

gapminder_df |>
  select(-lifeExp)