# compute kin bereavement measures for Argentina

library(tidyverse)
library(DemoTools)
library(geomtextpath)
library(DDSQLtools)
library(devtools)
load_all(path = "C:/Proyectos/DemoKin")
options(tibble.print_min=50)

# download data and save --------------------------------------------------

base_url <- 'https://population.un.org/dataportalapi/api/v1'
locations <- read.csv(paste0(base_url,'/locations/?format=csv'), sep='|', skip=1)
codes <- read.csv(paste0(base_url,'/indicators/?format=csv'), sep='|', skip=1) 
qx_code <- codes$Id[codes$ShortName == "qx1"]
dx_code <- codes$Id[codes$ShortName == "Deaths1"]
asfr_code <- codes$Id[codes$ShortName == "ASFR1"]
tfr_code <- codes$Id[codes$ShortName == "TFR5"]
e0_code <- codes$Id[codes$ShortName == "E0"]
mac_code <- codes$Id[codes$ShortName == "MAC5"]

# Guatemala, Argentina and LAC
my_locations <- c(32, 320, 904)

e0 <- map_df(my_locations, function(my_location){
  read.csv(paste0(base_url,
                  '/data/indicators/',e0_code,
                  '/locations/',my_location,
                  '/start/',1950,
                  '/end/',2020,
                  '/?format=csv'), sep='|', skip=1)}) %>% 
  dplyr::filter(Variant %in% "Median")
px <- map_df(my_locations, function(my_location){read.csv(paste0(base_url,
                      '/data/indicators/',qx_code,
                      '/locations/',my_location,
                      '/start/',1950,
                      '/end/',2020,
                      '/?format=csv'), sep='|', skip=1)}) %>% 
  dplyr::filter(Variant %in% "Median") %>% 
  dplyr::mutate(px = 1- Value) %>% 
  dplyr::select(Location, sex = Sex, year = TimeLabel, age = AgeStart, px)
dx <- map_df(my_locations, function(my_location){read.csv(paste0(base_url,
                     '/data/indicators/',dx_code,
                     '/locations/',my_location,
                     '/start/',1950,
                     '/end/',2020,
                     '/?format=csv'), sep='|', skip=1)}) %>% 
  dplyr::filter(Variant %in% "Median") %>% 
  dplyr::select(Location, sex = Sex, year = TimeLabel, age = AgeStart, Value)
tfr <- map_df(my_locations, function(my_location){read.csv(paste0(base_url,
                        '/data/indicators/',tfr_code,
                        '/locations/',my_location,
                        '/start/',1950,
                        '/end/',2020,
                        '/?format=csv'), sep='|', skip=1)}) %>% 
  dplyr::filter(Variant %in% "Median")
mac <- map_df(my_locations, function(my_location){read.csv(paste0(base_url,
                                                                  '/data/indicators/',mac_code,
                                                                  '/locations/',my_location,
                                                                  '/start/',1950,
                                                                  '/end/',2020,
                                                                  '/?format=csv'), sep='|', skip=1)}) %>% 
  dplyr::filter(Variant %in% "Median")
asfr <- map_df(my_locations, function(my_location){read.csv(paste0(base_url,
                        '/data/indicators/',asfr_code,
                        '/locations/',my_location,
                        '/start/',1950,
                        '/end/',2020,
                        '/?format=csv'), sep='|', skip=1)}) %>%  
  dplyr::filter(Variant %in% "Median") %>% 
  dplyr::select(Location, year = TimeLabel, age = AgeStart, fx = Value) %>% 
  dplyr::mutate(fx = fx/1000)

# basic demogr characteristics --------------------------------------------

plot_demographics <- bind_rows(e0 %>% filter(Sex == "Both sexes") %>% select(Location, year = TimeLabel, Value) %>% mutate(Indicator = "e0"),
          tfr %>% select(Location, year = TimeLabel, Value) %>% mutate(Indicator = "TFR"),
          mac %>% select(Location, year = TimeLabel, Value) %>% mutate(Indicator = "MAC"),
          dx  %>% filter(sex == "Both sexes", age >=10) %>% 
            summarise(Value = sum(Value * age)/sum(Value), .by = c(Location, year)) %>% 
            mutate(Indicator = "MAD")) %>% 
  ggplot(aes(year, Value, col = Location)) + 
  geom_line() +
  facet_wrap(~Indicator, scales = "free") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(plot = plot_demographics, filename = "T1 - Bereavement measures/plots/plot_demographics.pdf")

# build kin network --------------------------------------------------------

# build kin: focal is female
country_years <- expand.grid(country = c("Argentina", "Guatemala"), years = c(1950, 2015))
kin_net_country_years <- map_df(1:nrow(country_years) , function(i){
  # i = 4
  print(i)
  pf <- px %>% filter(sex == "Female", 
                      Location == country_years$country[i], 
                      year == as.character(country_years$year[i])) %>% pull(px)
  pm <- px %>% filter(sex == "Male", 
                      Location == country_years$country[i], 
                      year == as.character(country_years$year[i])) %>% pull(px)
  f <- asfr %>% filter(Location == country_years$country[i], 
                       year == as.character(country_years$year[i])) %>% pull(fx)
  f <- c(rep(0,15), f, rep(0,51))
  
  # Markov Chain and ex
  ages <- length(pf)
  age <- 0:(ages-1)
  Uf <- Um <- matrix(0, ages, ages)
  Uf[row(Uf)-1==col(Uf)] <- pf[-ages]
  Um[row(Um)-1==col(Um)] <- pm[-ages]
  ex_m <- data.frame(age = age, ex = solve(diag(1,ages)-Um) %>% colSums())
  ex_f <- data.frame(age = age, ex = solve(diag(1,ages)-Uf) %>% colSums())
  
  # Markov Chain joint-survival
  joint_ages <- expand.grid(sex_x = c("f"), sex_y = c("f", "m"), x=age, y=age, stringsAsFactors = FALSE)
  exy <- purrr::map_df(1:nrow(joint_ages), function(j){
    # j = 1500
    x <- joint_ages$x[j]
    y <- joint_ages$y[j]
    x_ages <- length((x+1):ages)
    y_ages <- length((y+1):ages)
    x_sex <- joint_ages$sex_x[j]
    y_sex <- joint_ages$sex_y[j]
    p_x <- p_y <- rep(0, ages)
    p_x[1:x_ages] <- if(x_sex == "m") pm[(x+1):ages] else pf[(x+1):ages]
    p_y[1:y_ages] <- if(y_sex == "m") pm[(y+1):ages] else pf[(y+1):ages]
    pxy <- p_x * p_y
    U <- matrix(0, ages, ages)
    U[row(U)-1==col(U)] <- pxy[-ages]
    U[ages, ages] <- pxy[ages]
    exy <- data.frame(x_0 = x, y_0 = y, 
                      t = age, 
                      x = x:(x+(ages-1)), y = y:(y+(ages-1)), 
                      x_sex = x_sex, y_sex = y_sex, 
                      exy = colSums(solve(diag(1,ages)-U)))
    exy$exy[dplyr::lag(exy$exy)==1] <- 0
    exy
    })

  # kin network
  px_kin <- bind_rows(
    data.frame(sex_kin = rep("m"), age_kin = age, px_kin = pm),
    data.frame(sex_kin = rep("f"), age_kin = age, px_kin = pf))
  ex_kin <- bind_rows(
    data.frame(sex_kin = rep("m"), age_kin = age, ex_kin = ex_m$ex),
    data.frame(sex_kin = rep("f"), age_kin = age, ex_kin = ex_f$ex))
  kin_out <- kin2sex(pf, pm, f, f, birth_female = .5, sex_focal = "f")$kin_full %>%
    mutate(kin = case_when(kin %in% c("os","ys") ~ "s",
                           kin %in% c("oa","ya") ~ "a",
                           kin %in% c("coa","cya") ~ "c",
                           kin %in% c("nos","nys") ~ "n",
                           T ~ kin)) %>%
    summarise(living = sum(living), dead = sum(dead), 
              .by = c(kin, age_kin, age_focal, sex_kin)) %>%
    left_join(data.frame(age_focal = age, px_focal = pf)) %>%
    left_join(data.frame(age_focal = age, ex_focal = ex_f$ex)) %>%
    left_join(px_kin) %>%
    left_join(ex_kin) %>% 
    left_join(exy %>% filter(t==0) %>% select(age_focal = x, age_kin = y, sex_kin = y_sex, exy)) %>%
    mutate(country = country_years$country[i], year = country_years$years[i])
  
  # return
  return(kin_out)
})
# save(kin_net_country_years, file = "data/kin_net_country_years.rda")
load("T1 - Bereavement measures/data/kin_net_country_years.rda")

# plot kin counts
plot_kin_net_country_years <- kin_net_country_years %>% 
  summarise(l = sum(living), d = sum(dead), .by = c(country, year, kin, age_focal)) %>% 
  filter(kin %in% c("d", "s", "m", "gm")) %>%
  mutate(kin = case_when(kin == "d" ~ "Children", 
                         kin == "s" ~ "Siblings", 
                         kin == "m" ~ "Parents", 
                         kin == "gm" ~ "Granparents", T ~ "Other")) %>% 
  ggplot(aes(age_focal, d, col=kin)) +
  geom_textline(aes(label = kin)) +
  scale_y_continuous(name = "") +
  scale_x_continuous(name = "Age of Focal",
                     breaks = seq(0,100,10), labels =  seq(0,100,10), 
                     limits = c(0,99), expand = c(0,0)) +
  theme_bw() +
  guides(color = "none") +
  facet_grid(rows = vars(country), cols = vars(year), scales = "free")
ggsave(plot = plot_kin_net_country_years, filename = "plots/kin_net_country_years.pdf")
  
# accumulated measures
summary_kin_net_country_years <- kin_net_country_years %>% 
  summarise(l = sum(living), d = sum(dead), .by = c(country, year, kin, age_focal)) %>% 
  filter(kin %in% c("d", "s", "m", "gm")) %>% 
  group_by(year, country, kin) %>% arrange(age_focal) %>% 
  mutate(D = cumsum(d)) %>% 
  ungroup %>% filter(age_focal %in% c(35, 65))

# functions ----------------------------------------------------------------

# main and other indicators
bereavement_indicators <- function(kin_net, x_focal, age_unexp, support_group){
  
  # load("T1 - Bereavement measures/data/kin_net_country_years.rda")
  # kin_net <- kin_net_country_years %>%
  #   filter(country == "Argentina", year == 1950)
  # x_focal <- 0
  # age_unexp <- 30
  # support_group <- c("s","m")
  
  # kin ever met
  ever_met <- kin_net %>%
    summarise(D = sum(dead[age_focal<x_focal]),
              l = sum(living[age_focal==x_focal]),
              .by = kin) %>% 
    mutate(L = D+l)
  
  # total death experience by age
  death_experienced <- kin_net  %>% 
    summarise(D = sum(dead), .by = c(age_focal)) %>% 
    mutate(D_cum = cumsum(D),
           D_porc = D/sum(D),
           D_porc_acum = cumsum(D_porc)) %>% 
    as.data.frame()
  
  # death in the future by kin
  death_to_experience <- kin_net  %>% 
    filter(age_focal >= x_focal) %>% 
    group_by(kin, age_focal) %>% 
    summarise(d = sum(dead)) %>% 
    group_by(kin) %>% 
    mutate(D = cumsum(d)) %>% 
    as.data.frame()
  
  # define age and get indicators
  ex_focal <- kin_net %>% filter(age_focal==x_focal) %>% pull(ex_focal) %>% unique
  lx_focal <- kin_net %>% distinct(age_focal, px_focal)
  lx_focal$cum_px_focal <- cumprod(lx_focal$px_focal)
  lx_focal$lx_focal <-c(1, lx_focal$cum_px_focal[-101])
  
  # types of kin
  kin_types <- unique(kin_net$kin)
  
  # T time expected in bereavement
  T2 <- kin_net %>% 
    filter(age_focal >= x_focal, age_kin > (age_focal - x_focal)) %>% 
    arrange(kin, age_focal) %>% 
    summarise(d = sum(dead), l = sum(living), 
              .by = c(kin, age_focal, px_focal, ex_focal)) %>%
    group_by(kin) %>% 
    mutate(px_y = cumprod(px_focal)/px_focal[age_focal==x_focal],
           ex_y = px_y*ex_focal,
           T2_x = d * ex_y) %>%
    summarise(T2 = sum(T2_x, na.rm = T))
  
  if(x_focal > 0){
    T1 <- kin_net %>%
      filter(age_kin<x_focal) %>% 
      summarise(T1 = sum(dead[age_focal < x_focal]),
                .by = kin) %>% 
      mutate(T1 = T1 * ex_focal)  
  }else{
    T1 <- T2 %>% mutate(T1 = 0) %>% select(-T2)
  }
  
  T3 <- kin_net %>% 
    filter(age_focal >= x_focal, age_kin <= (age_focal - x_focal)) %>% 
    arrange(kin, age_focal) %>% 
    summarise(d = sum(dead), l = sum(living), 
              .by = c(kin, age_focal, px_focal, ex_focal)) %>%
    group_by(kin) %>% 
    mutate(px_y = cumprod(px_focal)/px_focal[age_focal==x_focal],
           ex_y = px_y*ex_focal,
           T3_x = d * ex_y) %>%
    summarise(T3 = sum(T3_x, na.rm = T))
  
  T. <- inner_join(T1, T2) %>% inner_join(T3) %>% 
    mutate(T. = T1 + T2 + T3)
  
  # S time loss share
  S2 <- kin_net %>%
    filter(age_focal >= x_focal, age_kin > (age_focal - x_focal)) %>%
    left_join(lx_focal %>% select(age_focal, lx_focal)) %>% 
    arrange(kin, age_focal, x_focal) %>%
    select(kin, sex_kin, age_focal, age_kin, dead, lx_focal, exy) %>% 
    mutate(px_y = lx_focal/min(lx_focal[age_focal==x_focal]),
           exy = px_y * exy,
           S_x = dead * exy) %>%
    group_by(kin) %>%
    summarise(S2 = sum(S_x, na.rm = T))
  
  S3 <- kin_net %>%
    filter(age_focal >= x_focal, age_kin <= (age_focal - x_focal)) %>%
    left_join(lx_focal %>% select(age_focal, lx_focal)) %>% 
    arrange(kin, age_focal, x_focal) %>%
    select(kin, age_focal, age_kin, dead, lx_focal, exy) %>% 
    mutate(px_y = lx_focal/min(lx_focal[age_focal==x_focal]),
           exy = px_y * exy,
           S_x = dead * exy) %>%
    group_by(kin) %>%
    summarise(S3 = sum(S_x, na.rm = T))
  
  S1 <- S3 %>% mutate(S1 = 0) %>% select(-S3)
  
  S <- inner_join(S1, S2) %>% inner_join(S3) %>% 
    mutate(S = S1 + S2 + S3)
  
  # M portion time in bereavement
  M <- kin_net %>% 
    filter(age_focal < x_focal) %>% 
    summarise(M = 1-sum(dead*age_focal)/sum(dead)/x_focal, 
              .by = c(kin))
  
  # U Non-expected proportion of deaths: age definition
  U <- kin_net  %>% 
    filter(age_focal >= x_focal) %>% 
    summarise(U = sum(dead[age_kin<age_unexp])/sum(dead), .by = kin)
  
  # O overlapping/unnatural: ggm>gm, gm>m, m>s, m>d, a>d, d>gd, gd>ggd, a>c, c>d
  O <- kin_net  %>% 
    filter(age_focal>=x_focal) %>% 
    summarise(d = sum(dead), .by = c(kin, age_focal)) %>% 
    pivot_wider(names_from = kin, values_from = d) %>% 
    rowwise() %>% 
    mutate(
      gm_ggm_overlap = min(ggm, gm), gm_ggm_total = max(ggm, gm),
      m_gm_overlap = min(gm, m), m_gm_total = max(gm, m),
      a_gm_overlap = min(gm, a), a_gm_total = max(gm, a),
      s_m_overlap = min(s, m), s_m_total = max(s, m),
      d_m_overlap = min(d, m), d_s_total = max(d, m),
      gd_d_overlap = min(d, gd), gd_d_total = max(d, gd),
      ggd_gd_overlap = min(gd, ggd), ggd_gd_total = max(gd, ggd),
      c_a_overlap = min(a, c), c_a_total = max(c, a),
      n_s_overlap = min(n, s), n_s_total = max(n, s)) %>% 
    ungroup %>% 
    summarise(
      ggm = 0,
      gm = sum(gm_ggm_overlap)/sum(gm),
      m = sum(m_gm_overlap)/sum(m),
      a = sum(a_gm_overlap)/sum(a),
      s = sum(s_m_overlap)/sum(s),
      d = sum(d_m_overlap)/sum(d),
      gd = sum(gd_d_overlap)/sum(gd),
      ggd = sum(ggd_gd_overlap)/sum(ggd),
      c = sum(c_a_overlap)/sum(c),
      n = sum(n_s_overlap)/sum(n)
    ) %>% 
    pivot_longer(ggm:n, names_to = "kin", values_to = "O")
  
  ## Others
  # Death proportion over kin ever met
  P <- kin_net  %>% 
    summarise(D = sum(dead[age_focal<x_focal]),
              l = sum(living[age_focal==x_focal]),
              .by = kin) %>% 
    mutate(P = D/(D+l)) %>% 
    select(kin, P)
  
  # H Age with expected highest lost
  H <- kin_net %>% 
    filter(age_focal>=x_focal) %>% 
    summarise(D = sum(dead), .by = c(kin, age_focal)) %>% 
    group_by(kin) %>% arrange(-D) %>% slice(1) %>% 
    select(kin, H = age_focal)
  
  # L Loneliest age for lost
  L <- kin_net  %>% 
    filter(age_focal>=x_focal) %>% 
    summarise(L = sum(living[kin %in% support_group]), .by = c(age_focal)) %>% 
    arrange(L) %>% slice(1) %>% select(L = age_focal)
  L <- expand.grid(kin = kin_types, L = L$L)
    
  ## indicators
  bereav_indicators <- T. %>% 
    left_join(S) %>% 
    left_join(M) %>% 
    left_join(U) %>% 
    left_join(O) %>% 
    left_join(P) %>% 
    left_join(H) %>% 
    left_join(L)

  # out
  return(list(bereav_indicators = bereav_indicators,
              others = list(ever_met = ever_met,
                            death_experienced = death_experienced,
                            death_to_experience = death_to_experience))
         )
}

# extra information to characterize kin
extra_information <- function(kin_net){
  
  # load("T1 - Bereavement measures/data/kin_net_country_years.rda")
  # kin_net <- kin_net_country_years %>% filter(country == "Argentina", year == 2015)
  # x_focal <- 0
  
  # kin ever met
  ever_met <- inner_join(
    kin_net %>%
      summarise(D = sum(dead), .by = c(age_focal, kin)) %>% 
      arrange(kin, age_focal) %>% 
      mutate(D_cum = cumsum(D), 
             porc_D = D/sum(D),
             porc_D_cum = D_cum/sum(D),
             .by = kin) %>% 
      select(age_focal, kin, D, porc_D, D_cum, porc_D_cum), 
    kin_net %>%
      summarise(L = sum(living), .by = c(kin, age_focal))
    ) %>% 
    group_by(kin) %>% 
    mutate(E = L + lag(D_cum, default = 0))
  
  # out
  return(
    list(ever_met = ever_met)
  )
}

# compute indicators ------------------------------------------------------

load("T1 - Bereavement measures/data/kin_net_country_years.rda")

country_year <- kin_net_country_years %>% distinct(country, year)
country_year_extra_information <- map_df(1:nrow(country_year),
                                      # i = 4
                                      function(i){
                                        this_kin_net <- kin_net_country_years %>% 
                                          filter(country == country_year$country[i],
                                                 year == country_year$year[i])
                                        ex <- this_kin_net %>% 
                                          distinct(age_focal, ex_focal)
                                        extra_information_i <- extra_information(
                                          kin_net = this_kin_net)
                                        extra_information_i$ever_met %>% ungroup %>% 
                                          left_join(ex) %>% 
                                          mutate(country = country_year$country[i],
                                                 year = country_year$year[i], 
                                                 .before = 1)
})
save(country_year_extra_information, file = "T1 - Bereavement measures/data/country_year_extra_information.rda")

country_year_age <- kin_net_country_years %>%
  distinct(country, year) %>%
  cross_join(data.frame(age = c(0, 30, 60)))
country_year_age_indicators <- map_df(1:nrow(country_year_age),
                                         # i = 1
                                         function(i){
                                           this_kin_net <- kin_net_country_years %>% 
                                             filter(country == country_year_age$country[i],
                                                    year == country_year_age$year[i])
                                           bereavement_indicators_i <- bereavement_indicators(
                                             kin_net = this_kin_net, 
                                             x_focal = country_year_age$age[i], 
                                             age_unexp = 30, # to define 
                                             support_group = c("s","m") # to define dinamically
                                           )
                                           bereavement_indicators_i$bereav_indicators %>% 
                                             mutate(country = country_year_age$country[i],
                                                    year = country_year_age$year[i], 
                                                    age = country_year_age$age[i], 
                                                    .before = 1)
                                         })
save(country_year_age_indicators, 
     file = "T1 - Bereavement measures/data/country_year_age_indicators.rda")

# analysis paper----------------------------------------------------------------

# Argentina 2015
arg2015_extrainfo <- country_year_extra_information %>% ungroup %>% 
  filter(country == "Argentina", year == 2015)
arg2015_extrainfo %>% 
  summarise(D_cum = sum(D_cum), .by = age_focal) %>% 
  mutate(porc_D_cum = cumsum(D_cum/sum(D_cum))) %>% as.data.frame()
arg2015_bereavement <- country_year_age_indicators %>% ungroup %>% 
  filter(country == "Argentina", year == 2015) %>% 
  filter(kin %in% c("gm", "m", "d", "s")) %>%
  left_join(arg2015_extrainfo %>% distinct(age = age_focal, ex_focal))
table_arg2015_bereavement <- arg2015_bereavement %>% 
  select(age, ex_focal, kin, `T` = T., S, M, U, O) %>%  #, P, H, L) %>% 
  arrange(age, kin) %>% 
  kbl(align = "c", booktabs = T, format = "latex", digits = 2)

# argentina and guatemala 1950 and 2015
plot_arg_guatem_1950_2015 <- country_year_age_indicators %>% 
  ungroup %>% 
  filter(kin %in% c("gm", "m", "d", "s")) %>% 
  select(country, year, age, kin, `T` = T., S, M, U, O) %>% 
  pivot_longer(`T`:O, names_to = "Indicator", values_to = "Value") %>% 
  ggplot(aes(age, Value, col = country, linetype = factor(year))) +
  geom_line() + geom_point() +
  theme_bw() +
  scale_x_continuous(name ="", 
                     labels = c(0,30,60), 
                     breaks = c(0,30,60))+
  labs(color = "Year", linetype = "Year") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  facet_grid(rows = vars(Indicator), col = vars(kin), 
             scales = "free", switch = "y", as.table = T)
ggsave(plot = plot_arg_guatem_1950_2015, 
       filename = "T1 - Bereavement measures/plots/plot_arg_guatem_1950_2015.pdf") 

# plots -------------------------------------------------------------------

# plot time bereavement example
kin_lines <- data.frame(kin = c("grandmother", "mother", rep("daugter",4)),
                        x_focal     = c(0, 0,  20, 25,  35, 40),
                        x_kin       = c(40, 20, 0,  0,   0,  0),
                        x_focal_end = c(10, 60, 25, 100, 45, 100),
                        x_kin_end   = c(50, 80, 5,  75,  10, 60),
                        death = c(T, T, T, F, T, F))

kin_lines_plot <- ggplot(kin_lines) +
  geom_segment(aes(x = x_focal, y = x_kin, xend = x_focal_end, yend = x_kin_end, col=kin)) +
  coord_equal() +
  geom_point(data = . %>% filter(death),
             aes(x = x_focal_end, y = x_kin_end), shape = 3) +
  geom_vline(xintercept = 30, linetype=2, col =1) +
  coord_equal() +
  scale_x_continuous(expand = c(0,0), name = "Age Focal",
                     breaks = seq(0, 100, 10), labels = seq(0, 100, 10)) +
  scale_y_continuous(expand = c(0,0), name = "Age Kin", limits = c(0,100),
                     breaks = seq(0, 100, 10), labels = seq(0, 100, 10)) +
  theme_bw()
ggsave(plot = kin_lines_plot, filename = "plots/kin_lines_plot.pdf")


# plot freq all
kin_net  %>% 
  summarise(`d(x)` = sum(dead), .by = c(kin, age_focal)) %>% 
  arrange(kin, age_focal) %>% 
  mutate(`D(x)` = cumsum(`d(x)`), .by = c(kin)) %>%
  pivot_longer(`d(x)`:`D(x)`, names_to = "type", values_to = "D") %>% 
  ggplot(aes(age_focal, D, col=kin)) +
  geom_textline(aes(label = kin)) +
  scale_y_continuous(name = "") +
  scale_x_continuous(name = "Age of Focal",
                     breaks = seq(0,100,10), labels =  seq(0,100,10), 
                     limits = c(0,99), expand = c(0,0)) +
  theme_bw() +
  guides(color = "none") +
  facet_wrap(~type, scales = "free")

# plot freq youngest
kin_net  %>% 
  filter(kin %in% c("n", "d", "gd", "ggd")) %>% 
  summarise(`d(x)` = sum(dead), .by = c(kin, age_focal)) %>% 
  arrange(kin, age_focal) %>% 
  mutate(`D(x)` = cumsum(`d(x)`), .by = c(kin)) %>%
  pivot_longer(`d(x)`:`D(x)`, names_to = "type", values_to = "D") %>% 
  ggplot(aes(age_focal, D, col=kin)) +
  geom_textline(aes(label = kin)) +
  scale_y_continuous(name = "") +
  scale_x_continuous(name = "Age of Focal",
                     breaks = seq(0,100,10), labels =  seq(0,100,10), 
                     limits = c(0,99), expand = c(0,0)) +
  theme_bw() +
  guides(color = "none") +
  facet_wrap(~type, scales = "free")

# plot T3 example
load("data/kin_net_years.rda")
m_gm <- kin_net_years[[4]]  %>%
  filter(kin %in% c("gd","ggd")) %>%
  summarise(`d(x)` = sum(dead), .by = c(kin, age_focal))
m_gm %>% 
  ggplot(aes(age_focal, `d(x)`)) +
  geom_textline(aes(label = kin, col=kin), size = 5) +
  geom_area(data = m_gm %>% summarise(`d(x)` = min(`d(x)`), .by = "age_focal"), 
            aes(age_focal, `d(x)`), fill  ="lightgrey") +
  # geom_line(data = m_gm %>% summarise(`d(x)` = sum(`d(x)`), .by = "age_focal"), 
  #          aes(age_focal, `d(x)`), linetype = 2) +
  scale_y_continuous(name = "") +
  scale_x_continuous(name = "Age of Focal",
                     breaks = seq(0,100,10), labels =  seq(0,100,10),
                     limits = c(0,99), expand = c(0,0)) +
  theme_bw() +
  guides(color = "none")

load("data/kin_net_years.rda")
m_gm <- kin_net_years[[4]]  %>%
  filter(kin %in% c("gm","m")) %>%
  summarise(`d(x)` = sum(dead), .by = c(kin, age_focal))
P2 <- m_gm %>% 
  ggplot(aes(age_focal, `d(x)`)) +
  geom_textline(aes(label = kin, col=kin), size = 5) +
  geom_area(data = m_gm %>% summarise(`d(x)` = min(`d(x)`), .by = "age_focal"), 
            aes(age_focal, `d(x)`), fill  ="lightgrey") +
  scale_y_continuous(name = "") +
  scale_x_continuous(name = "Age of Focal",
                     breaks = seq(0,100,10), labels =  seq(0,100,10),
                     limits = c(0,99), expand = c(0,0)) +
  theme_bw() +
  guides(color = "none")
ggsave(plot = P2, filename = "plots/P2.pdf")

plot_kin_arg_2010 <- kin_net_years[[4]]  %>% 
  summarise(`d_k(x)` = sum(dead), .by = c(kin, age_focal)) %>% 
  arrange(kin, age_focal) %>% 
  mutate(`D_k(x)` = cumsum(`d_k(x)`), .by = c(kin)) %>%
  pivot_longer(`d_k(x)`:`D_k(x)`, names_to = "type", values_to = "D") %>% 
  ggplot(aes(age_focal, D, col=kin)) +
  geom_textline(aes(label = kin)) +
  scale_y_continuous(name = "") +
  scale_x_continuous(name = "Age of Focal",
                     breaks = seq(0,100,10), labels =  seq(0,100,10), 
                     limits = c(0,99), expand = c(0,0)) +
  theme_bw() +
  guides(color = "none") +
  facet_wrap(~type, scales = "free")
ggsave(plot = plot_kin_arg_2010, filename = "plots/kin_arg_2010.pdf")

load("data/arg_bereav_years.rda")
kable(arg_bereav_years %>% filter(year == 2010), digits = 2, caption = "Bereavement indicators in a female-dominant population under 2015 demographic conditions")

# end ---------------------------------------------------------------------
