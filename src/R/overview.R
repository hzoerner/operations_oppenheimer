library(tidyverse)

julia_dir <- "operations_oppenheimer/"

snf_s <- read.csv(file = paste0(julia_dir,"snf_stored.csv"))

nc_s <- read.csv(file= paste0(julia_dir,"nc_stored.csv"))

snf_t <- read.csv(file = paste0(julia_dir,"snf_shipped.csv"))

nc_end <- read.csv(file = paste0(julia_dir,"end_transport.csv"))

snf_t_1 <- snf_t %>%
  filter(SNF > 0)

nc_t <- read.csv(file = paste0(julia_dir,"nc_shipped.csv"))

nc_t_1 <- nc_t %>%
  filter(NC > 0)

nc_s %>%
  group_by(year, node) %>%
  summarise(total = sum(NC))

snf_over_time <- snf_s %>%
  filter(node %in% c("Reaktor 1", "Reaktor 2"))

ggplot(snf_s, aes(x = year, y = SNF, fill = node)) +
  
  geom_col()

ggplot(nc_s, aes(x = year, y = NC, fill = node)) +
  
  geom_col() +
  
  labs(y = "Accumulated amount of casks", fill = "Node") +
  
  theme_bw()

cisf_build <- read.csv(paste0(julia_dir, "cisf_build.csv"))

cisf_build %>%
  filter(build != 0)

library(xlsx)

options(scipen = 100, digits = 4)

cisf <- c("Finkenkrug", "Gorleben", "Bitterfeld-Wolfen")

#### transport costs

transport_nc <- read.xlsx(file = paste0(julia_dir, "data/NuclearData.xlsx"), sheetName = "Transport") %>%
  filter(is_possible != 0) %>%
  unite(col = "key", from, to) %>%
  right_join(mutate(nc_t_1, key = paste(from, to, sep = "_"))) %>%
  rename("transport" = "NC")

transport_snf <- read.xlsx(file = paste0(julia_dir, "data/NuclearData.xlsx"), sheetName = "Transport") %>%
  filter(is_possible != 0) %>%
  unite(col = "key", from, to) %>%
  right_join(mutate(snf_t_1, key = paste(from, to, sep = "_"))) %>%
  rename("transport" = "SNF")

transport <- bind_rows(transport_nc, transport_snf) %>%
  select(-key) %>%
  mutate(total_costs = costs * transport)

transport_summary <-transport %>%
  group_by(year) %>%
  summarise(transport = sum(total_costs))

rm(transport, transport_nc, transport_snf)

general <- read.xlsx(file = paste0(julia_dir, "data/NuclearData.xlsx"), sheetName = "General")

#### store costs

store <- snf_s %>%
  rename(units = SNF) %>%
  bind_rows(rename(nc_s, "units" = "NC")) %>%
  group_by(node, year) %>%
  summarise(units = sum(units)) %>%
  ungroup() %>%
  mutate(is_cisf = ifelse(node %in% cisf, T, F),
         costs = ifelse(is_cisf, units * general[1,2], units * general[4,2]))

store_summary <- store %>%
  group_by(year) %>%
  summarise(store = sum(costs))

rm(store)

#### cisf costs
hot_cells <- c("Hot Cell 1", "Hot Cell 2")

refill_costs <- snf_t %>%
  filter(to %in% hot_cells) %>%
  filter(SNF != 0) %>%
  select(-from) %>%
  mutate(costs = SNF * general[4,2])

refill_summary <- refill_costs %>%
  group_by(year) %>%
  summarise(refill = sum(costs))
  

costs <- data.frame(year = c(2030:2099)) %>%
  left_join(transport_summary, by = "year") %>%
  left_join(store_summary, by = "year") %>%
  left_join(refill_summary, by = "year") %>%
  pivot_longer(cols = c("store", "transport", "refill"), names_to = "type", values_to = "costs_per_year") %>%
  replace_na(list(costs_per_year = 0))

rm(store_summary, transport_summary, refill_summary)

ggplot(costs, aes(year, costs_per_year, fill = type)) +
  
  geom_col() +
  
  labs(x = "Year", y = "Costs per year in Euro", fill = "Type:") +
  
  theme_bw()

