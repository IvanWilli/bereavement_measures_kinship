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
  filter(LocID %in% all_LAC_loc) 
  # %>% 
  # bind_rows(
  #   data.table::fread("D:/data+/WPP2024_Life_Table_Complete_Medium_Female_2024-2100.csv") %>% 
  #     filter(LocID %in% all_LAC_loc))
lt_male <-   data.table::fread("D:/data+/WPP2024_Life_Table_Complete_Medium_Male_1950-2023.csv") %>% 
  filter(LocID %in% all_LAC_loc) 
  # %>% 
  # bind_rows(
  #   data.table::fread("D:/data+/WPP2024_Life_Table_Complete_Medium_Male_2024-2100.csv") %>% 
  #     filter(LocID %in% all_LAC_loc))
asfr <- data.table::fread("D:/data+/WPP2024_Fertility_by_Age1.csv") %>% 
  filter(LocID %in% all_LAC_loc, Variant == "Medium", Time %in% unique(lt_male$Time))

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
save(lt_female, lt_male, asfr, demo, file = "data/demo_data_LA.rda")

# build kin network --------------------------------------------------------

load("data/demo_data_LA.rda")
countries <- unique(lt_female$Location)
cohorts <- 1950
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
        right_join(expand.grid(AgeGrpStart = 0:100, Time = 1950:2023)) %>% 
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
    }, .options = furrr_options(seed = TRUE), .progress = TRUE) 
future::plan(future::sequential)
kin_net_country_years <- kin_net_country_years %>% left_join(demokin_codes %>% rename(kin = DemoKin))
save(kin_net_country_years, file = "data/kin_data_LA.rda")

# indicators -------------------------------------------------------------

countries_above1M <- c(
  "Argentina",
  "Bolivia (Plurinational State of)",
  "Brazil",
  "Chile",
  "Colombia",
  "Costa Rica",
  "Cuba",
  "Dominican Republic",
  "Ecuador",
  "El Salvador",
  "Guatemala",
  "Haiti",
  "Honduras",
  "Jamaica",
  "Mexico",
  "Nicaragua",
  "Panama",
  "Paraguay",
  "Peru",
  "Puerto Rico",
  "Trinidad and Tobago",
  "Uruguay",
  "Venezuela (Bolivarian Republic of)"
)

load("data/kin_data_LA.rda")
kin_net_country_years <- kin_net_country_years |>
  filter(country %in% countries_above1M)

# survival cohort
surv_df <- lt_female %>% 
  bind_rows(lt_male) %>% 
  filter(Location %in% countries_above1M) %>%
  select(country = Location, age = AgeGrpStart, year = Time, px, sex = Sex) %>% 
  mutate(cohort = year - age,
         sex = ifelse(sex == "Female", "f", "m"))%>% 
  arrange(country, desc(cohort), age)

# e0
ex <- surv_df %>% 
  mutate(S = cumprod(px),
         ex = rev(cumsum(rev(S))) / S, .by = c(country, cohort, sex))

cohort_example <- 1950
# joint life expectancy
e0x <- map_df(
  expand.grid(cohort = cohort_example, sex = c("f", "m"), country = countries) %>% 
    split(list(.$cohort, .$sex, .$country)), function(csc){
      # csc <- expand.grid(cohort = cohorts, sex = c("f", "m"), country = countries) %>% slice(1)
      focal_cohort <- csc$cohort
      sex_focal <- "f"
      sex_kin   <- csc$sex
      this_country <- csc$country
      # tabla de supervivencia de Focal desde edad 0
      focal_surv <- surv_df %>%
        filter(cohort == focal_cohort, sex == sex_focal, country == this_country) %>%
        arrange(age) %>%
        mutate(
          t  = row_number(),
          Sf = cumprod(px)
        ) %>%
        select(t, age_focal = age, Sf)
      # edades iniciales posibles del kin en el momento en que Focal nace
      res_e0x <- tibble(x = 0:max(surv_df$age, na.rm = TRUE)) %>%
        mutate(
          kin_cohort = focal_cohort - x,
          e0x = map_dbl(x, \(x0) {
            kin_surv <- surv_df %>%
              filter(cohort == focal_cohort - x0,
                     sex == sex_kin, country == this_country,
                     age >= x0) %>%
              arrange(age) %>%
              mutate(
                t  = row_number(),
                Sk = cumprod(px)
              ) %>%
              select(t, age_kin = age, Sk)
            both <- focal_surv %>%
              inner_join(kin_surv, by = "t")
            if (nrow(both) == 0) return(NA_real_)
            sum(both$Sf * both$Sk, na.rm = TRUE)
          })
        )
      res_e0x %>% select(kin_age = x, e0x) %>% 
        mutate(country = this_country, sex_kin = sex_kin, focal_cohort = focal_cohort)
    }
)

### Burden
plot_data <- kin_net_country_years %>% 
  filter(cohort == cohort_example) %>% 
  left_join(ex %>% 
              filter(age == 0, sex == "f", cohort == cohort_example) %>% 
              select(cohort, e0 = ex, country)) %>% 
  filter(age_focal < e0) %>% 
  summarise(l = sum(living), d = sum(dead), e0 = unique(e0),
    .by = c(kin, age_focal, country)
  ) %>% 
  mutate(
    D = cumsum(d), L = lag(D, default = 0) + l,
    .by = c(kin, country)
  ) %>%
  filter(kin %in% c("d", "m", "s", "gm"))

labels <- plot_data %>%
  summarise(B = sum(D) / sum(L), x = mean(age_focal),
    .by = c(kin, country)
  ) %>% 
  left_join(demokin_codes %>% 
              rename(kin = DemoKin)) %>%
  mutate(label = paste0(round(B*100,0), "%"))

#------------------------------------------------------

table_S <- kin_net_country_years |> 
  arrange(country, cohort, age_focal, kin) |> 
 filter(cohort == cohort_example) %>% 
      # focal surviving until kin_age
      left_join(ex %>% 
                  filter(country %in% countries, sex =="f") %>% 
                  select(country, age_focal = age, cohort, S, ex)) %>%
      # joint life expectancy for each kin_age
      left_join(e0x %>% 
                  filter(country %in% countries) %>% 
                  select(country, age_kin = kin_age, sex_kin, e0x, cohort = focal_cohort)) %>% 
      mutate(S_x = dead * S * e0x) %>%
      summarise(
                # mean_age = sum(dead * age_focal)/sum(dead),
                # D_mean_age = sum(dead[age_focal < trunc(mean_age)]),
                D_e0 = sum(dead[age_focal < trunc(ex[age_focal == 0])], na.rm = T),
                # mean_age_kin = sum(dead[age_focal == trunc(mean_age)] * age_kin[age_focal == trunc(mean_age)])/sum(dead[age_focal == trunc(mean_age)]),
                # e0x = min(e0x[age_kin == trunc(mean_age_kin) & sex_kin == "f"]), 
                S = sum(S_x, na.rm = T), 
                S_mean = S/D_e0, .by = c(Labels_2sex, country)) |> 
  filter(Labels_2sex  %in% c("Children", "Siblings", "Parents", "Grandparents")) |> 
  select(Country = country, Kin = Labels_2sex, D_e0, S, S_mean)

# ---------------------------------------------------------

O_data <- kin_net_country_years %>% 
  filter(cohort == cohort_example) %>% 
  left_join(ex %>% 
              filter(age == 0, sex == "f", cohort == cohort_example) %>% 
              select(cohort, e0 = ex, country)) %>% 
  filter(age_focal < e0)
  m_gm <- O_data  %>%
  filter(kin %in% c("m","gm"), cohort == cohort_example) %>%
  summarise(`d(x)` = sum(dead), .by = c(country, kin, age_focal))
  m_gm_o <- m_gm  %>%
  summarise(overl = min(`d(x)`),
            total = sum(`d(x)`), .by = c(country, age_focal))
  # overlapping s and m
  s_m <- O_data  %>%
  filter(kin %in% c("s","m"), cohort == cohort_example) %>%
  summarise(`d(x)` = sum(dead), .by = c(country, kin, age_focal))
  s_m_o <- s_m  %>%
  summarise(overl = min(`d(x)`),
            total = sum(`d(x)`), .by = c(country, age_focal))
# overlapping d and m
d_m <- O_data  %>%
  filter(kin %in% c("d","m"), cohort == cohort_example) %>%
  summarise(`d(x)` = sum(dead), .by = c(country, kin, age_focal))
d_m_o <- d_m  %>%
  summarise(overl = min(`d(x)`),
            total = sum(`d(x)`), .by = c(country, age_focal))
overl <- bind_rows(
  m_gm_o %>% mutate(Related = "Parents - Grandparent"),
  s_m_o %>% mutate(Related = "Siblings - Parents"),
  d_m_o %>% mutate(Related = "Children - Parents"),
)
ddistr <- bind_rows(
  m_gm %>% mutate(Related = "Parents - Grandparent"),
  s_m %>% mutate(Related = "Siblings - Parents"),
  d_m %>% mutate(Related = "Children - Parents")) %>% 
  left_join(demokin_codes %>% 
              rename(kin = DemoKin))
O_table <- left_join(
  ddistr %>% 
    summarise(D = sum(`d(x)`), 
              mad = sum(`d(x)`* age_focal)/D, .by = c(Related, kin, country)),
  overl %>% 
    summarise(dOverl = sum(overl),
              dTotal = sum(total), 
              Overl = dOverl*2/dTotal, .by = c(Related, country)))

# maps: B, shared time, and overlap --------------------------------------
# Run after plot_data, table_S, and O_table exist in the environment.
if (!requireNamespace("sf", quietly = TRUE) ||
    !requireNamespace("rnaturalearth", quietly = TRUE) ||
    !requireNamespace("patchwork", quietly = TRUE)) {
  stop("Install packages `sf`, `rnaturalearth`, and `patchwork` to draw the maps.")
}

map_name_recode <- c(
  "Bolivia (Plurinational State of)" = "Bolivia",
  "Falkland Islands (Malvinas)" = "Falkland Islands",
  "Saint Martin (French part)" = "Saint Martin",
  "Venezuela (Bolivarian Republic of)" = "Venezuela"
)
map_name_recode <- c(
  map_name_recode,
  setNames(
    paste0("Saint-Barth", intToUtf8(0xe9), "lemy"),
    paste0("Saint Barth", intToUtf8(0xe9), "lemy")
  )
)

add_map_name <- function(dat, country_col = "country") {
  dat %>%
    mutate(map_name = recode(.data[[country_col]], !!!map_name_recode,
                             .default = .data[[country_col]]))
}

map_shapes <- rnaturalearth::ne_countries(
  scale = "medium", type = "map_units", returnclass = "sf"
) %>%
  select(admin, name, name_long, brk_name, sovereignt, geometry)

join_map_shapes <- function(dat) {
  wanted <- unique(dat$map_name)
  
  out <- map_shapes %>%
    mutate(map_name = case_when(
      admin %in% wanted ~ admin,
      name %in% wanted ~ name,
      name_long %in% wanted ~ name_long,
      brk_name %in% wanted ~ brk_name,
      sovereignt %in% wanted ~ sovereignt,
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(map_name)) %>%
    select(map_name, geometry) %>%
    left_join(dat, by = "map_name")
  
  missing_shapes <- setdiff(wanted, unique(out$map_name))
  if (length(missing_shapes) > 0) {
    warning("No map shape found for: ", paste(missing_shapes, collapse = ", "))
  }
  
  out
}

map_gradient <- c("#dadad9ff", "#828283ff", "#141414ff")

make_map <- function(dat, legend_title, labels) {
  ggplot(join_map_shapes(dat)) +
    geom_sf(aes(fill = value), color = "white", linewidth = 0.15) +
    coord_sf(xlim = c(-120, -30), ylim = c(-60, 35), expand = FALSE) +
    scale_fill_gradientn(
      colors = map_gradient, na.value = "grey90",
      labels = labels, name = legend_title
    ) +
    theme_void() +
    theme(legend.position = "bottom")
}

make_map_list <- function(dat, legend_title, labels) {
  dat %>%
    split(.$facet) %>%
    map(~make_map(.x, legend_title, labels))
}

make_map_faceted <- function(map_list, ncol = 2) {
  patchwork::wrap_plots(map_list, ncol = ncol)
}

kin_facet_levels <- c("Children", "Parents", "Siblings", "Grandparents")

B_map_data <- labels %>%
  filter(kin %in% c("d", "m", "s", "gm")) %>%
  transmute(country, facet = Labels_2sex, value = B) %>%
  mutate(facet = factor(facet, levels = kin_facet_levels)) %>%
  add_map_name()

S_map_data <- table_S %>%
  transmute(country = Country, facet = Kin, value = S_mean) %>%
  mutate(facet = factor(facet, levels = kin_facet_levels)) %>%
  add_map_name()

O_map_data <- O_table %>%
  distinct(country, Related, Overl) %>%
  transmute(
    country,
    facet = recode(str_squish(Related),
                   "Parents - Grandparent" = "Parents - Grandparents"),
    value = Overl
  ) %>%
  mutate(facet = factor(
    facet,
    levels = c("Children - Parents", "Siblings - Parents", "Parents - Grandparents")
  )) %>%
  add_map_name()

B_maps <- make_map_list(B_map_data, "B", scales::label_percent(accuracy = 1))
S_maps <- make_map_list(S_map_data, "S_mean", scales::label_number(accuracy = 0.1))
O_maps <- make_map_list(O_map_data, "Overl", scales::label_percent(accuracy = 1))

B_maps_faceted <- make_map_faceted(B_maps)
S_maps_faceted <- make_map_faceted(S_maps)
O_maps_faceted <- make_map_faceted(O_maps)

B_maps
S_maps
O_maps
B_maps_faceted
S_maps_faceted
O_maps_faceted


# country hierarchies -----------------------------------------------------

make_country_hierarchy <- function(dat, legend_title, labels,
                                   scale_value = 1,
                                   facet_levels = NULL) {
  highlight_cols <- c("Argentina" = "#0072B2", "Guatemala" = "#E69F00")
  
  if (is.null(facet_levels)) {
    facet_levels <- if (is.factor(dat$facet)) {
      levels(droplevels(dat$facet))
    } else {
      unique(dat$facet)
    }
  }
  
  dat <- dat %>%
    mutate(facet = factor(as.character(facet), levels = facet_levels)) %>%
    filter(!is.na(facet))
  
  first_facet <- levels(dat$facet)[1]
  
  country_order <- dat %>%
    filter(facet == first_facet) %>%
    summarise(value = first(value), .by = country) %>%
    arrange(desc(value)) %>%
    pull(country)
  
  plot_dat <- dat %>%
    mutate(
      value_plot = value * scale_value,
      segment_col = coalesce(highlight_cols[as.character(country)], "grey80"),
      country = factor(country, levels = rev(country_order))
    ) %>%
    filter(!is.na(country))
  
  ggplot(plot_dat, aes(value_plot, country, color = value_plot)) +
    geom_segment(aes(x = 0, xend = value_plot, yend = country,
                     color = I(segment_col)),
                 linewidth = 0.5, show.legend = FALSE) +
    geom_point(size = 2.1) +
    facet_wrap(vars(facet), scales = "free_x", drop = TRUE) +
    scale_x_continuous(labels = labels) +
    scale_color_gradientn(
      colors = map_gradient,
      labels = labels,
      name = legend_title
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 18) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

B_country_hierarchy <- make_country_hierarchy(
  B_map_data, "B (per 100)",
  scales::label_number(accuracy = 1),
  scale_value = 100
)
S_country_hierarchy <- make_country_hierarchy(
  S_map_data, "S_mean",
  scales::label_number(accuracy = 0.1)
)
O_country_hierarchy <- make_country_hierarchy(
  O_map_data, "Overl (per 100)",
  scales::label_number(accuracy = 1),
  scale_value = 100,
  facet_levels = c("Children - Parents", "Siblings - Parents", "Parents - Grandparents")
)

B_country_hierarchy
S_country_hierarchy
O_country_hierarchy

ggsave(
  "plots/B_country_hierarchy.pdf",
  B_country_hierarchy,
  width = 11, height = 12
)
ggsave(
  "plots/S_country_hierarchy.pdf",
  S_country_hierarchy,
  width = 11, height = 12
)
ggsave(
  "plots/O_country_hierarchy.pdf",
  O_country_hierarchy,
  width = 11, height = 10
)
