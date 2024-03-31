library(tidyverse)
library(ggtheme)

DIR <- "C:/Users/ACER/Desktop/Uni/OR-INF/3-29/"
EXCEL_PATH <- "C:/Users/ACER/Desktop/Uni/OR-INF/operations_oppenheimer/data/ExtendedNuclearData.xlsx"

cap_util <- data.frame(cisf = character(), total = numeric(), utilization = numeric(), n = character())

for (n in c(1, 3, 5)){
  snf <- read.csv(file = paste0(DIR, "factor5_", n, "cisf/snf_stored.csv")) %>%
    rename("total" = "SNF")
  nc <- read.csv(file = paste0(DIR, "factor5_", n, "cisf/nc_stored.csv")) %>%
    rename("total" = "NC")
  
  cisf <- read.csv(file = paste0(DIR, "factor5_", n, "cisf/cisf_build.csv")) %>%
    filter(year == 2099 & build == 1)
  
  total = bind_rows(snf, nc) %>%
    filter(node %in% cisf$cisf) %>%
    rename("cisf" = "node") %>%
    group_by(cisf, year) %>%
    summarise(total = sum(total)) %>%
    arrange(year) %>%
    ungroup() %>%
    mutate(utilization = round(total / 1500, 2), 
           n = as.character(n))
  
  at_the_end <- total %>%
    filter(year == 2099) %>%
    select(-year)
  
  cap_util = bind_rows(cap_util, at_the_end)
  rm(at_the_end, total, cisf, snf, nc)
}

util_plot <- ggplot(cap_util, aes(cisf, utilization, fill = n)) +
  
  geom_col(color = "black", position = position_dodge(width = .6), width = .5) +
  
 # geom_text(aes(label = utilization), position = position_dodge(width = .6)) +
  
  labs(y = "Capacity utilization", fill = "Total number of\nCISF in scenario:", x = "CISF") 

ggsave(plot = util_plot, filename = "cap_util_plot.png")
