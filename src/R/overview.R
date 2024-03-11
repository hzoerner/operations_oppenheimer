library(tidyverse)

julia_dir <- "operations_oppenheimer/"

snf_s <- read.csv(file = paste0(julia_dir,"snf_stored.csv"))

nc_s <- read.csv(file= paste0(julia_dir,"nc_stored.csv"))

snf_t <- read.csv(file = paste0(julia_dir,"snf_shipped.csv"))

snf_t_1 <- snf_t %>%
  filter(SNF > 0)

nc_t <- read.csv(file = paste0(julia_dir,"nc_shipped.csv"))

nc_t_1 <- nc_t %>%
  filter(NC > 0)


nc_s %>%
  group_by(year) %>%
  summarise(total = sum(NC))

snf_s %>%
  group_by(year) %>%
  summarise(total = sum(SNF))

ggplot(snf_s, aes(x = year, y = SNF, fill = node)) +
  
  geom_col()

ggplot(nc_s, aes(x = year, y = NC, fill = node)) +
  
  geom_col()
