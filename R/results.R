# compute kin bereavement measures for Argentina and Guatemala

library(tidyverse)
library(DemoTools)
library(geomtextpath)
library(kableExtra)
library(DemoKin)
library(DDSQLtools)
source("R/funs.R")
options(tibble.print_min=30)

# get data and plot demographics --------------------------------------------------

# LA locations
get_locations()
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
my_locations <- c(32, 320, 904)
plot_demo <- demo |> 
  filter(Time < 2024) |> 
  ggplot(aes(Time, Value, group = Location)) +
  geom_line(data = . %>% filter(!LocID %in% my_locations), col = "grey80") +
  geom_line(data = . %>% filter(LocID %in% my_locations), aes(col = Location)) +
  facet_wrap(~Indicator, scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Year", y = "")
ggsave(plot = plot_demo, 
      filename = "plots/plot_demo.pdf", 
      width = 8, height = 6)

# save demographic data for later use
save(lt_female |> filter(LocID %in% my_locations), 
        lt_male |> filter(LocID %in% my_locations), 
        asfr |> filter(LocID %in% my_locations), file = "data/demo_data.rda")

# build kinship network --------------------------------------------------------

load("data/demo_data.rda")

# params
countries <- c("Argentina", "Guatemala")
cohorts <- seq(1950, 2010, 30)

# create trees
kin_net_country_years <- map_df(countries, function(this_country){
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
    }) %>% 
  left_join(demokin_codes %>% rename(kin = DemoKin))
save(kin_net_country_years, file = "data/kin_data.rda")

# kin death analysis ----------------------------------------------------

load("data/kin_data.rda")

# summary by age
kin_death <- kin_net_country_years %>% 
  filter(kin %in% c("d", "s", "m", "gm"), cohort %in% c(1950)) %>%
  summarise(l = sum(living), d = sum(dead), .by = c(country, year, kin, age_focal)) %>% 
  mutate(D = cumsum(d), .by = c(country, year, kin)) %>% 
  mutate(L = lag(D)+l) %>% 
  mutate(kin = case_when(kin == "d" ~ "Children", 
                         kin == "s" ~ "Siblings", 
                         kin == "m" ~ "Parents", 
                         kin == "gm" ~ "Granparents", T ~ "Other"),
         year = as.character(year)) 

# cohort
cohort_example <- 1950

# e0 cohort
e0_df <- lt_female %>% 
  filter(Location %in% countries) %>%
  select(country = Location, age = AgeGrpStart, year = Time, px, sex = Sex) %>% 
  mutate(cohort = year - age,
         sex = ifelse(sex == "Female", "f", "m"))%>% 
  arrange(country, desc(cohort), age) %>% 
  mutate(S = cumprod(px),
         ex = rev(cumsum(rev(S))) / S, .by = c(country, cohort, sex)) |> 
  filter(cohort == cohort_example, sex == "f", age == 0)

# plot
plot_death_kin <- kin_death %>% 
  filter(age_focal < 80) %>%
  ggplot(aes(age_focal, D, col = kin, label = kin)) +
  geom_textline() +
  geom_vline(data = e0_df, aes(xintercept = ex), linetype = 2) +
  theme_bw() +
  guides(color = "none") +
  scale_y_continuous(name = "Deaths") +
  scale_x_continuous(name = "Age of Focal",
                     breaks = seq(0,75,10), labels =  seq(0,75,10), expand = c(0,0)) +
  facet_grid(cols = vars(country), scales = "free") + 
  theme(panel.grid.minor.x = element_blank())
ggsave(plot = plot_death_kin, filename = "plots/plot_death_kin.pdf", width = 8, height = 6)

# calculate indicators -------------------------------------------------------------

load("data/kin_data.rda")
load("data/demo_data.rda")
countries <- c("Argentina", "Guatemala")
cohort_example <- 1950

# survival cohort
surv_df <- lt_female %>% bind_rows(lt_male) %>% 
  select(country = Location, age = AgeGrpStart, year = Time, px, sex = Sex) %>% 
  mutate(cohort = year - age,
         sex = ifelse(sex == "Female", "f", "m"))%>% 
  arrange(country, desc(cohort), age)

# e0
ex <- surv_df %>% 
  mutate(S = cumprod(px),
         ex = rev(cumsum(rev(S))) / S, .by = c(country, cohort, sex))

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

# Burden

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

plot_B <- ggplot(plot_data %>% 
         left_join(demokin_codes %>% rename(kin = DemoKin)), 
       aes(x = age_focal)) +
  geom_area(aes(y = D, fill = "Accumulated deaths"), linewidth = 0.8) +
  geom_line(aes(y = L, color = "Kin ever met"), linewidth = 0.8) +
  geom_vline(data = plot_data %>% distinct(kin, country, e0), 
             aes(xintercept = e0), col = "grey", linetype = 2) +
  geom_text(
    data = labels,
    aes(x = x, y = Inf, label = label),
    vjust = 1.2,
    size = 3.5) +
  facet_grid(cols = vars(country), rows = vars(Labels_2sex))+
  labs(x = "Focal Age", y  ="") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_manual(
    name = NULL,
    values = c("Accumulated deaths" = "grey80")
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Kin ever met" = "black")
  )
ggsave(filename = "plots/plot_I.pdf",
       plot = plot_B, 
       width = 8, height = 6)

# Shared Time

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

# table
library(kableExtra)
table_S |>
  mutate(
    across(c(D_e0, S, S_mean), ~round(.x, 1))
  ) |>
  kbl(
    format = "latex",
    booktabs = TRUE,
    longtable = FALSE,
    align = c("l", "l", "r", "r", "r"),
    col.names = c("Country", "Kin", "$D(e(0))$", "$S(0)$", "$overline{S}(0)$")
  ) |>
  kable_styling(
    latex_options = c("hold_position"),
    font_size = 10
  ) |>
  collapse_rows(columns = 1, valign = "top", latex_hline = "major")

# Overlapping

# overlapping m and ggm
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
 
# table
O_table <- left_join(
  ddistr %>% 
    summarise(D = sum(`d(x)`), 
              mad = sum(`d(x)`* age_focal)/D, .by = c(Related, kin, country)),
  overl %>% 
    summarise(dOverl = sum(overl),
              dTotal = sum(total), 
              Overl = dOverl*2/dTotal, .by = c(Related, country)))

# plot
O_labels <- O_table %>% 
  distinct(Related, country, Overl) %>% 
  mutate(
    age_focal = mean(c(0, 77)),
    label_y = Inf,
    label = paste0(round(Overl*100, 0), "%")
  )
O_plot <- ggplot() +
  geom_area(
    data = overl,
    aes(age_focal, 2*overl, fill = "Area")
  ) +
  geom_textline(
    data = overl,
    aes(age_focal, total),
    linetype = 2, label = "T", size = 3
  ) +
  geom_textline(
    data = ddistr, 
    aes(age_focal, `d(x)`, label = substr(Labels_2sex, 1, 1), col = kin),
    size = 3
  ) +
  geom_text(
    data = O_labels,
    aes(age_focal, label_y, label = label),
    vjust = 1.2,
    hjust = 0.5,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(name = "Deaths", expand = c(0, 0)) +
  scale_x_continuous(
    name = "Age of Focal",
    breaks = seq(0, 77, 10),
    labels = seq(0, 77, 10),
    limits = c(0, 77),
    expand = c(0, 0)
  ) +
  geom_vline(data = plot_data %>% distinct(kin, country, e0), 
             aes(xintercept = e0), col = "grey", linetype = 2) +
  scale_fill_manual(
    name = NULL,
    values = c("Area" = "grey90")
  ) +
  guides(color = "none", fill = "none") +
  theme_bw() +
  facet_grid(rows = vars(Related), cols = vars(country))

ggsave(filename = "plots/O_plot.pdf",
       plot = O_plot, 
       width = 8, height = 6)


# end --------------------------------------------------------------------




















# ### all indicators
# bereavement_indicators <- map_df(
#   expand.grid(cohort = cohorts, country = countries) %>% 
#     split(list(.$cohort, .$country)), function(cc){
      
#     # cc <- expand.grid(cohort = cohorts, country = countries) %>% slice(4)
#     this_country <- cc$country
#     this_cohort <- cc$cohort
    
#     # árbol
#     this_kin_net <- kin_net_country_years %>% 
#       filter(country == this_country, cohort == this_cohort)
    
#     # focal
#     e0_focal <- ex %>% filter(country == this_country, cohort == this_cohort, 
#                               age == 0, sex == "f") %>% pull(ex)
    
#     ### I
#     I_table <- this_kin_net %>% 
#       filter(age_focal < e0_focal) %>% 
#       summarise(l = sum(living),
#                 d = sum(dead), .by = c(kin, age_focal)) %>% 
#       mutate(D = cumsum(d), 
#              L = lag(D, default = 0) + l, .by = kin) 
    
#     # indicador
#     I_index <- I_table %>% 
#       summarise(ll = sum(l), LL = sum(L), I = 1-ll/LL, .by = kin)
    
#     # I plot
#     I_plot <- I_table %>% 
#       select(age_focal, l, L, kin) %>% 
#       pivot_longer(l:L, names_to = "Living", values_to = "n") %>% 
#       ggplot(aes(age_focal, n, col = Living)) + geom_line() + facet_wrap(~kin)
    
#     ### S
#     S_index <- this_kin_net %>%
#       # focal surviving until kin_age
#       left_join(ex %>% 
#                   filter(country == this_country, sex =="f") %>% 
#                   select(age_focal = age, cohort, S)) %>%
#       # joint life expectancy for each kin_age
#       left_join(e0x %>% 
#                   filter(country == this_country) %>% 
#                   select(age_kin = kin_age, sex_kin, e0x,
#                          cohort = focal_cohort)) %>% 
#       mutate(S_x = dead * S * e0x) %>%
#       summarise(S = sum(S_x, na.rm = T), .by = kin)
    
#     ### O overlapping/unnatural: ggm>gm, gm>m, m>s, m>d, a>d, d>gd, gd>ggd, a>c, c>d
#     O_index <- this_kin_net  %>% 
#       filter(age_focal < e0_focal) %>% 
#       summarise(d = sum(dead), .by = c(kin, age_focal)) %>% 
#       pivot_wider(names_from = kin, values_from = d) %>% 
#       rowwise() %>% 
#       mutate(
#         gm_ggm_overlap = min(ggm, gm), gm_ggm_total = sum(ggm, gm),
#         m_gm_overlap = min(gm, m),     m_gm_total = sum(gm, m),
#         a_gm_overlap = min(gm, a),     a_gm_total = sum(gm, a),
#         s_m_overlap = min(s, m),       s_m_total = sum(s, m),
#         d_m_overlap = min(d, m),       d_m_total = sum(d, m),
#         gd_d_overlap = min(d, gd),     gd_d_total = sum(d, gd),
#         ggd_gd_overlap = min(gd, ggd), ggd_gd_total = sum(gd, ggd),
#         c_a_overlap = min(a, c),       c_a_total = sum(c, a),
#         n_s_overlap = min(n, s),       n_s_total = sum(n, s)) %>% 
#       ungroup %>% 
#       summarise(
#         ggm = 0,
#         gm = sum(gm_ggm_overlap)/sum(gm_ggm_total),
#         m = sum(m_gm_overlap)/sum(m_gm_total),
#         a = sum(a_gm_overlap)/sum(a_gm_total),
#         s = sum(s_m_overlap)/sum(s_m_total),
#         d = sum(d_m_overlap)/sum(d_m_total),
#         gd = sum(gd_d_overlap)/sum(gd_d_total),
#         ggd = sum(ggd_gd_overlap)/sum(ggd_gd_total),
#         c = sum(c_a_overlap)/sum(c_a_total),
#         n = sum(n_s_overlap)/sum(n_s_total)
#       ) %>% 
#       pivot_longer(ggm:n, names_to = "kin", values_to = "O")
    
#     ## indicators
#     I_index %>% 
#       left_join(S_index, by = "kin") %>% 
#       left_join(O_index, by = "kin") %>% 
#       mutate(country = this_country,
#              cohort = this_cohort, .before = 1)
# })

# # analysis
# bereavement_indicators %>% 
#   filter(kin %in% c("d", "s", "m", "gm")) %>% 
#   left_join(demokin_codes %>% select(kin = DemoKin, name_kin = Labels_2sex)) %>% 
#   select(country, cohort, name_kin, I, S, O) %>% 
#   pivot_longer(I:O, names_to = "Indicator", values_to = "Value") %>%
#   ggplot(aes(cohort, Value, col = country)) +
#   geom_line() + geom_point() +
#   facet_grid(rows = vars(Indicator), 
#              cols = vars(name_kin), scales = "free_y") +
#   theme_bw()

# ######################################################################








# # build kin: focal is female
# country_years <- expand.grid(country = c("Argentina", "Guatemala"), years = seq(1950, 2100, 5))
# kin_net_country_years <- map_df(1:nrow(country_years) , function(i){
#   # i = 4
#   print(i)
#   qf <- lt_female %>% filter(Location == country_years$country[i], 
#                              Time == as.character(country_years$year[i])) %>% pull(qx)
#   pf <- 1 - qf
#   qm <- lt_male %>% filter(Location == country_years$country[i], 
#                            Time == as.character(country_years$year[i])) %>% pull(qx)
#   pm <- 1 - qm
#   f <- asfr %>% filter(Location == country_years$country[i], 
#                        Time == as.character(country_years$year[i])) %>% pull(ASFR)/1000
#   f <- c(rep(0,15), f, rep(0,51))
  
#   # Markov Chain and ex
#   ages <- length(pf)
#   age <- 0:(ages-1)
#   Uf <- Um <- matrix(0, ages, ages)
#   Uf[row(Uf)-1==col(Uf)] <- pf[-ages]
#   Um[row(Um)-1==col(Um)] <- pm[-ages]
#   ex_m <- data.frame(age = age, ex = solve(diag(1,ages)-Um) %>% colSums())
#   ex_f <- data.frame(age = age, ex = solve(diag(1,ages)-Uf) %>% colSums())
  
#   # Markov Chain joint-survival
#   exy0 <- exy_t0(pf, pm, sex_x = "f", sex_y = c("f","m"))
  
#   # kin network - time invariant
#   px_kin <- bind_rows(
#     data.frame(sex_kin = rep("m", ages), age_kin = age, px_kin = pm),
#     data.frame(sex_kin = rep("f", ages), age_kin = age, px_kin = pf))
#   ex_kin <- bind_rows(
#     data.frame(sex_kin = rep("m", ages), age_kin = age, ex_kin = ex_m$ex),
#     data.frame(sex_kin = rep("f", ages), age_kin = age, ex_kin = ex_f$ex))
#   kin_out <- kin2sex(pf, pm, f, f, birth_female = .5, sex_focal = "f")$kin_full %>%
#     mutate(kin = case_when(kin %in% c("os","ys") ~ "s",
#                            kin %in% c("oa","ya") ~ "a",
#                            kin %in% c("coa","cya") ~ "c",
#                            kin %in% c("nos","nys") ~ "n",
#                            T ~ kin)) %>%
#     summarise(living = sum(living), dead = sum(dead), .by = c(kin, age_kin, sex_kin, age_focal)) %>%
#     left_join(data.frame(age_focal = age, px_focal = pf), by = "age_focal") %>%
#     left_join(data.frame(age_focal = age, ex_focal = ex_f$ex), by = "age_focal") %>%
#     left_join(px_kin, by = c("age_kin", "sex_kin")) %>%
#     left_join(ex_kin, by = c("age_kin", "sex_kin")) %>% 
#     left_join(exy0 %>% select(age_focal = x0, age_kin = y0, sex_kin = sex_y, exy = exy0), 
#               by = c("age_kin", "sex_kin", "age_focal")) %>%
#     mutate(country = country_years$country[i], year = country_years$years[i], .before = 1)
  
#   # return
#   return(kin_out)
# })
# # save(kin_net_country_years, file = "data/kin_net_country_years.rda")
# load("data/kin_net_country_years.rda")

# beareavement patterns ------------------------------------------------------

# calculation
# load("data/kin_net_country_years.rda")
# x_focal <- 0
# country_year <- kin_net_country_years %>% distinct(country, year)
# country_year_indicators <- map_df(1:nrow(country_year),function(i){
#   # i = 60
  
#   # arbol
#   this_kin_net <- kin_net_country_years %>% 
#     filter(country == country_year$country[i], year == country_year$year[i])
  
#   # focal
#   e0_focal <- this_kin_net %>% filter(age_focal == x_focal) %>% distinct(ex_focal) %>% pull() %>% round()
  
#   sss <- this_kin_net %>% filter(age_focal <= e0_focal, kin == "d") %>% 
#     summarise(l = sum(living),
#               d = sum(dead), .by = age_focal) %>% 
#     mutate(D = cumsum(d), 
#            L = lag(D, default = 0) + l)
  
#   sss %>% 
#     select(age_focal, l, L) %>% 
#     pivot_longer(l:L, names_to = "Living", values_to = "n") %>% 
#     ggplot(aes(age_focal, n, col = Living)) + geom_line()
  
#   sss %>% summarise(ll = sum(l), LL = sum(L), I = 1-ll/LL)
  
#   # burden and timming
#   D_e0 <- this_kin_net %>% 
#     filter(age_focal <= e0_focal) %>%
#     summarise(l = sum(living), d = sum(dead), .by = c(kin, age_focal)) %>% 
#     # porcentaje de muerte hasta los 50 años (completó nuevos kins, no harbá más, excepto ggd)
#     summarise(L = sum(l[age_focal < e0_focal]),
#               D = sum(d[age_focal < e0_focal]),
              
              
#               mean_age = sum((d[age_focal < e0_focal]) * (age_focal[age_focal < e0_focal]))/sum(d[age_focal < e0_focal]),
#               mean_age_trunc = trunc(mean_age),
#               L = l[age_focal == e0_focal], 
#               D = sum(d[age_focal < mean_age_trunc]),
#               Total = L + D,
#               pD = D/Total,
#               aD = 1 - mean_age/e0_focal,
#               Index = pD * aD,
#               .by = c(kin)) %>% 
#     arrange(kin) %>% 
#     select(kin, D, L, Total, pD, , aD, Index, mean_age) %>% 
#     mutate(e0 = e0_focal, .before = 1)
  
#   # S time loss share
#   lx_focal <- this_kin_net %>% distinct(age_focal, px_focal)
#   lx_focal$cum_px_focal <- cumprod(lx_focal$px_focal)
#   lx_focal$lx_focal <-c(1, lx_focal$cum_px_focal[-nrow(lx_focal)]) 
#   S <- this_kin_net %>%
#     filter(age_focal >= x_focal) %>%
#     left_join(lx_focal %>% select(age_focal, lx_focal), by = "age_focal") %>% 
#     arrange(kin, age_focal, x_focal) %>%
#     select(kin, sex_kin, age_focal, age_kin, dead, lx_focal, exy) %>% 
#     mutate(px_y = lx_focal/min(lx_focal[age_focal==x_focal]),
#            S_x = dead * px_y * exy) %>%
#     group_by(kin) %>%
#     summarise(S = sum(S_x, na.rm = T))
  
#   # O overlapping/unnatural: ggm>gm, gm>m, m>s, m>d, a>d, d>gd, gd>ggd, a>c, c>d
#   O <- this_kin_net  %>% 
#     filter(age_focal>=x_focal, age_focal<e0_focal) %>% 
#     summarise(d = sum(dead), .by = c(kin, age_focal)) %>% 
#     pivot_wider(names_from = kin, values_from = d) %>% 
#     rowwise() %>% 
#     mutate(
#       gm_ggm_overlap = min(ggm, gm), gm_ggm_total = sum(ggm, gm),
#       m_gm_overlap = min(gm, m),     m_gm_total = sum(gm, m),
#       a_gm_overlap = min(gm, a),     a_gm_total = sum(gm, a),
#       s_m_overlap = min(s, m),       s_m_total = sum(s, m),
#       d_m_overlap = min(d, m),       d_m_total = sum(d, m),
#       gd_d_overlap = min(d, gd),     gd_d_total = sum(d, gd),
#       ggd_gd_overlap = min(gd, ggd), ggd_gd_total = sum(gd, ggd),
#       c_a_overlap = min(a, c),       c_a_total = sum(c, a),
#       n_s_overlap = min(n, s),       n_s_total = sum(n, s)) %>% 
#     ungroup %>% 
#     summarise(
#       ggm = 0,
#       gm = sum(gm_ggm_overlap)/sum(gm_ggm_total),
#       m = sum(m_gm_overlap)/sum(m_gm_total),
#       a = sum(a_gm_overlap)/sum(a_gm_total),
#       s = sum(s_m_overlap)/sum(s_m_total),
#       d = sum(d_m_overlap)/sum(d_m_total),
#       gd = sum(gd_d_overlap)/sum(gd_d_total),
#       ggd = sum(ggd_gd_overlap)/sum(ggd_gd_total),
#       c = sum(c_a_overlap)/sum(c_a_total),
#       n = sum(n_s_overlap)/sum(n_s_total)
#     ) %>% 
#     pivot_longer(ggm:n, names_to = "kin", values_to = "O")
  
#   ## indicators
#   D_e0 %>% 
#     left_join(S, by = "kin") %>% 
#     left_join(O, by = "kin") %>% 
#     mutate(country = country_year$country[i],
#            year = country_year$year[i], .before = 1)
# })
# save(country_year_indicators, file = "data/country_year_indicators.rda")

# # analysis
# country_year_indicators %>% 
#   filter(kin %in% c("d", "s", "m", "gm")) %>% 
#   left_join(demokin_codes %>% select(kin = DemoKin, name_kin = Labels_2sex)) %>% 
#   select(country, year, name_kin, pD, aD, Index, mean_age, S, O) %>% 
#   pivot_longer(pD:O, names_to = "Indicator", values_to = "Value") %>%
#   filter(Indicator %in% c("pD", "aD", "Index")) %>% 
#   ggplot(aes(year, Value, col = country)) +
#   geom_line() +
#   facet_grid(rows = vars(Indicator), 
#              cols = vars(name_kin), scales = "free_y") +
#   theme_bw()

# country_year_indicators %>% 
#   filter(kin %in% c("d", "s", "m", "gm")) %>% 
#   left_join(demokin_codes %>% select(kin = DemoKin, name_kin = Labels_2sex)) %>% 
#   select(country, year, name_kin, pD, mean_age) %>% 
#   ggplot(aes(mean_age, pD, col = name_kin, linetype = country)) +
#   geom_line() +
#   theme_bw()

# # indicators
# kin_net_country_years %>% 
#   filter(kin %in% c("d", "s", "m", "gm"), 
#          age_focal == 50, year < 2020) %>%
#   summarise(living = sum(living), .by = c(country, year, kin)) %>% 
#   left_join(demo %>%
#               filter(Indicator %in% c("MAC", "TFR")) %>% 
#               pivot_wider(names_from = Indicator, values_from = Value) %>% 
#               rename(year = Time, country = Location)) %>% 
#   ggplot(aes(TFR, MAC)) +
#   geom_point(aes(col = country, size = living)) +
#   theme_bw() +
#   facet_wrap(~kin)
  

# joint life expectancy ---------------------------------------------------

# https://chatgpt.com/share/69849a75-2db8-8001-9d22-d9fff9b7ebab

# exy_t0 <- function(pf, pm, sex_x = "f", sex_y = c("f","m")) {
#   ages <- length(pf)
#   age  <- 0:(ages-1)
  
#   # pick survival vectors once
#   pX <- if (sex_x == "m") pm else pf
  
#   grid <- expand.grid(
#     sex_x = sex_x,
#     sex_y = sex_y,
#     x0 = age,
#     y0 = age,
#     KEEP.OUT.ATTRS = FALSE,
#     stringsAsFactors = FALSE
#   )
  
#   # compute e_xy(t=0) per row (very light work, no matrices)
#   grid$exy0 <- mapply(function(x0, y0, sy) {
#     pY <- if (sy == "m") pm else pf
    
#     # available steps until max age reached
#     L <- ages - max(x0, y0)
#     if (L <= 0) return(1)
    
#     # joint 1-year survivals along the diagonal
#     pxy <- pX[(x0+1):(x0+L)] * pY[(y0+1):(y0+L)]
    
#     if (L == 1) return(1)
#     1 + sum(cumprod(pxy[1:(L-1)]))
#   }, grid$x0, grid$y0, grid$sex_y)
  
#   as_tibble(grid)
# }
