library(tidyverse)
library(DemoKin)
library(furrr)

# LAC
all_LAC_loc <- c(
  32, 68, 76, 152, 170, 188, 192, 214, 218, 222,
  320, 332, 340, 484, 558, 591, 600, 604, 858, 862,
  44, 388, 662, 533, 52, 531, 308, 312, 474, 630,
  670, 850, 28, 740, 780, 660, 535, 92, 136, 212,
  238, 500, 652, 659, 663, 534, 796, 84, 254, 328
)
lt_female <- data.table::fread("D:/data+/WPP2024_Life_Table_Complete_Medium_Female_1950-2023.csv") %>% 
  filter(LocID %in% all_LAC_loc) %>% 
  bind_rows(
    data.table::fread("D:/data+/WPP2024_Life_Table_Complete_Medium_Female_2024-2100.csv") %>% 
      filter(LocID %in% all_LAC_loc))
lt_male <-   data.table::fread("D:/data+/WPP2024_Life_Table_Complete_Medium_Male_1950-2023.csv") %>% 
  filter(LocID %in% all_LAC_loc) %>% 
  bind_rows(
    data.table::fread("D:/data+/WPP2024_Life_Table_Complete_Medium_Male_2024-2100.csv") %>% 
      filter(LocID %in% all_LAC_loc))
asfr <- data.table::fread("D:/data+/WPP2024_Fertility_by_Age1.csv") %>% 
  filter(LocID %in% all_LAC_loc, Variant == "Medium")

# plot demographics
demo <- bind_rows(
  lt_female %>% 
    summarise(`e(0)` = ex[AgeGrpStart == 0],
              `q(0,5)` = 1-lx[AgeGrpStart == 5]/lx[AgeGrpStart == 0], .by = c(Time, LocID, Location)) %>% 
    pivot_longer(`e(0)`:`q(0,5)`, names_to = "Indicator", values_to = "Value"),
  asfr %>% 
    summarise(TFR = sum(ASFR/1000),
              MAC = sum(AgeGrpStart * ASFR/1000)/TFR, .by = c(Time, LocID, Location))  %>% 
    pivot_longer(TFR:MAC, names_to = "Indicator", values_to = "Value")
)
# save demographic data for later use
save(lt_female |> filter(LocID %in% my_locations), 
        lt_male |> filter(LocID %in% my_locations), 
        asfr |> filter(LocID %in% my_locations), 
        file = "data/demo_data_LA.rda")

# build kin network --------------------------------------------------------

load("data/demo_data_LA.rda")
countries <- unique(lt_female$Location)

# create trees
future::plan(future::multisession, workers = max(1, future::availableCores() - 1))
kin_net_country_years <- future_map_dfr(countries, function(this_country){
      pf <- lt_female %>% filter(Location == this_country) %>% 
        select(AgeGrpStart, Time, px) %>% 
        pivot_wider(names_from = Time, values_from = px) %>% 
        select(-AgeGrpStart) %>% as.matrix()
      pm <- lt_male %>% filter(Location == this_country) %>% 
        select(AgeGrpStart, Time, px) %>% 
        pivot_wider(names_from = Time, values_from = px) %>% 
        select(-AgeGrpStart) %>% as.matrix()
      ff <- asfr %>% filter(Location == this_country) %>% 
        select(AgeGrpStart, Time, ASFR) %>% 
        mutate(ASFR = ASFR/1000) %>% 
        right_join(expand.grid(AgeGrpStart = 0:100, Time = 1950:2100)) %>% 
        replace(is.na(.), 0) %>% 
        arrange(Time, AgeGrpStart) %>% 
        pivot_wider(names_from = Time, values_from = ASFR) %>% 
        select(-AgeGrpStart) %>% as.matrix()
      kin_net <- kin2sex(pf, pm, ff, ff, birth_female = .48, 
                                       time_invariant = F, sex_focal = "f",
                                       output_cohort = cohorts)$kin_full %>% 
        mutate(kin = case_when(kin %in% c("os","ys") ~ "s",
                               kin %in% c("oa","ya") ~ "a",
                               kin %in% c("coa","cya") ~ "c",
                               kin %in% c("nos","nys") ~ "n",
                               T ~ kin)) %>%
        summarise(living = sum(living), dead = sum(dead), 
                  .by = c(cohort, year, age_focal, kin, age_kin, sex_kin)) %>% 
        mutate(country = this_country, .before = 1)
      return(kin_net)
    }, .options = furrr_options(seed = TRUE)) %>% 
  left_join(demokin_codes %>% rename(kin = DemoKin))
future::plan(future::sequential)
save(kin_net_country_years, file = "data/kin_data_LA.rda")
