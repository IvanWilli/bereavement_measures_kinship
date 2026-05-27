# compute kin bereavement measures for Argentina

library(tidyverse)
library(DemoTools)
library(geomtextpath)
library(devtools)
library(kableExtra)
library(DemoKin)
options(tibble.print_min=50)

# download data and save --------------------------------------------------

# Guatemala, Argentina and LAC
my_locations <- c(32, 320, 904)

lt_female_1950a2023 <- data.table::fread("D:/data+/WPP2024_Life_Table_Complete_Medium_Female_1950-2023.csv") %>% 
  filter(LocID %in% my_locations)
lt_male_1950a2023 <-   data.table::fread("D:/data+/WPP2024_Life_Table_Complete_Medium_Male_1950-2023.csv") %>% 
  filter(LocID %in% my_locations)
lt_female_2024a2100 <- data.table::fread("D:/data+/WPP2024_Life_Table_Complete_Medium_Female_2024-2100.csv") %>% 
  filter(LocID %in% my_locations)
lt_male_2024a2100 <-   data.table::fread("D:/data+/WPP2024_Life_Table_Complete_Medium_Male_2024-2100.csv") %>% 
  filter(LocID %in% my_locations)
asfr <- data.table::fread("D:/data+/WPP2024_Fertility_by_Age1.csv")

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
load("data/kin_net_country_years.rda")

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

# compute indicators ------------------------------------------------------

load("data/kin_net_country_years.rda")

# extra info
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
save(country_year_extra_information, 
     file = "data/country_year_extra_information.rda")

# indicators
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
                                           bereavement_indicators_i %>% 
                                             mutate(country = country_year_age$country[i],
                                                    year = country_year_age$year[i], 
                                                    age = country_year_age$age[i], 
                                                    .before = 1)
                                         })
save(country_year_age_indicators, 
     file = "data/country_year_age_indicators.rda")

# analysis paper----------------------------------------------------------------

load("data/kin_net_country_years.rda")
load("data/country_year_age_indicators.rda")
load("data/country_year_extra_information.rda")

# Argentina 2015
arg2015_extrainfo <- country_year_extra_information %>% ungroup %>% 
  filter(country == "Argentina", year == 2015)
arg2015_extrainfo %>% 
  summarise(D_cum = sum(D_cum), .by = age_focal) %>% 
  mutate(porc_D_cum = cumsum(D_cum/sum(D_cum))) %>% as.data.frame()
arg2015_bereavement <- country_year_age_indicators %>% ungroup %>% 
  filter(country == "Argentina", year == 2015) %>% 
  filter(kin %in% c("gm", "m", "d", "s")) %>%
  left_join(demokin_codes %>% select(kin = DemoKin, Labels_2sex)) %>% 
  left_join(arg2015_extrainfo %>% select(age = age_focal, kin, Living = L, D_cum, ex_focal))
table_arg2015_bereavement <- arg2015_bereavement  %>% 
  select(Age = age, `e(x)` = ex_focal, Kin = Labels_2sex, `L(x)` = Living, `D(x)` = D_cum, 
         `T(x)` = T., `S(x)` = S, `M(x)` = M_retr, `U(x)` = U_retr, `O(x)` = O_retr) %>%  #, P, H, L) %>% 
  arrange(Age, Kin) %>% 
  kbl(align = "c", booktabs = T, format = "latex", digits = 2)
table_arg2015_bereavement

# argentina and guatemala 1950 and 2015
arg_guatem_1950_2015_bereavement <- country_year_age_indicators %>% ungroup %>% 
  filter(kin %in% c("gm", "m", "d", "s")) %>% 
  left_join(demokin_codes %>% rename(kin = DemoKin)) %>% 
  left_join(country_year_extra_information %>% 
              select(country, year, age = age_focal, kin, Living = L, D_cum, ex_focal)) %>% 
  select(Country = country, Year = year,
         Age = age, `e(x)` = ex_focal, Kin = Labels_2sex, `L(x)` = Living, `D(x)` = D_cum, 
         `T(x)` = T., `S(x)` = S, `M(x)` = M_retr, `U(x)` = U_retr, `O(x)` = O_retr)

plot_arg_guatem_1950_2015 <- arg_guatem_1950_2015_bereavement %>%
  select(-`L(x)`, -`D(x)`) %>% 
  pivot_longer(`T(x)`:`O(x)`, names_to = "Indicator", values_to = "Value") %>% 
  ggplot(aes(Age, Value, col = Country, linetype = factor(Year))) +
  geom_line() + geom_point() +
  theme_bw() +
  scale_x_continuous(name ="Age of Focal", 
                     labels = c(0,30,60), 
                     breaks = c(0,30,60))+
  labs(color = "Country", linetype = "Year") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  facet_grid(rows = vars(Indicator), col = vars(Kin), 
             scales = "free", switch = "y", as.table = T)
plot_arg_guatem_1950_2015
ggsave(plot = plot_arg_guatem_1950_2015, 
       filename = "plots/plot_arg_guatem_1950_2015.pdf") 

table_arg_guatem_1950_2015 <- arg_guatem_1950_2015_bereavement %>% 
  arrange(Country, Age, Kin) %>% 
  kbl(align = "c", booktabs = T, format = "latex", digits = 2)

# overlapping plot m and ggm
m_gm <- kin_net_country_years  %>%
  filter(kin %in% c("m","gm")) %>%
  summarise(`d(x)` = sum(dead), .by = c(country, year, kin, age_focal))
m_gm_o <- m_gm  %>%
  summarise(overl = min(`d(x)`),
            total = sum(`d(x)`), .by = c(country, year, age_focal))
m_gm_mad <- m_gm  %>%
  summarise(mad = sum(`d(x)`*age_focal)/sum(`d(x)`), .by = c(country, year, kin))
plot_arg_guatem_1950_2015_overl <- ggplot() +
  geom_textline(data = m_gm_o,
            aes(age_focal, total), linetype = 2, label = "total") +
  geom_area(data = m_gm_o,
            aes(age_focal, overl), fill  ="lightgrey") +
  geom_textline(data = m_gm, 
                aes(age_focal, `d(x)`, label = kin, col=kin), size = 5) +
  geom_vline(data = m_gm_mad, aes(xintercept = mad, col = kin), linetype = 2) + 
  scale_y_continuous(name = "", expand = c(0, 0)) +
  scale_x_continuous(name = "Age of Focal",
                     breaks = seq(0,100,10), labels =  seq(0,100,10),
                     limits = c(0,99), expand = c(0,0)) +
  theme_bw() +
  guides(color = "none") +
  facet_grid(rows = vars(year), cols = vars(country), switch = "y")
ggsave(plot = plot_arg_guatem_1950_2015_overl, 
       filename = "plots/plot_arg_guatem_1950_2015_overl.pdf") 

# example plot
D_kins_arg_2015 <- m_gm %>% 
  filter(country == "Argentina", year == 2015) %>% 
  summarise(sum(`d(x)`), .by = kin) %>% pull(2)
D_overlap_arg_2015 <- m_gm_o %>% 
  filter(country == "Argentina", year == 2015) %>% 
  summarise(overl = sum(overl), total = sum(total)) %>% 
  as.numeric()

plot_arg_2015_overl_example <- ggplot() +
  geom_textline(data = m_gm_o %>% filter(country == "Argentina", year == 2015),
                aes(age_focal, total), linetype = 2, label = "total") +
  geom_area(data = m_gm_o%>% filter(country == "Argentina", year == 2015),
            aes(age_focal, overl), fill  ="lightgrey") +
  geom_textline(data = m_gm %>% filter(country == "Argentina", year == 2015), 
                aes(age_focal, `d(x)`, label = kin, col=kin), size = 5) +
  # geom_vline(data = m_gm_mad %>% filter(country == "Argentina", year == 2015), aes(xintercept = mad, col = kin), linetype = 2) + 
  scale_y_continuous(name = "", expand = c(0, 0)) +
  scale_x_continuous(name = "Age of Focal",
                     breaks = seq(0,100,10), labels =  seq(0,100,10),
                     limits = c(0,99), expand = c(0,0)) +
  annotate("text", x = 40, y = .01, label = round(D_overlap_arg_2015[1],1), col = 1) +
  annotate("text", x = 55, y = .08, label = round(D_overlap_arg_2015[2],1), col = 1) +
  annotate("text", x = 60, y = .03, label = round(D_kins_arg_2015[1],1), col = "#00BFC4") +
  annotate("text", x = 25, y = .07, label = round(D_kins_arg_2015[2],1), col = "#F8766D") +
  theme_bw() +
  guides(color = "none")
plot_arg_2015_overl_example
ggsave(plot = plot_arg_2015_overl_example, 
       filename = "plots/plot_arg_2015_overl_example.pdf") 

# plots -------------------------------------------------------------------

# plot time bereavement example
kin_lines <- data.frame(kin = c("grandmother", "mother", rep("children",4), "focal"),
                        x_focal     = c(0, 0,  20, 25,  35, 40,   0),
                        x_kin       = c(40, 20, 0,  0,   0,  0,   0),
                        x_focal_end = c(10, 60, 25, 100, 45, 100, 70),
                        x_kin_end   = c(50, 80, 5,  75,  10, 60,  70),
                        death = c(T, T, T, F, T, F, T))

kin_lines_plot <- ggplot(kin_lines %>% filter(kin != "focal")) +
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
  scale_color_manual(values = c(2, 4, 3)) +
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
kable(arg_bereav_years %>% filter(year == 2010), digits = 2, 
      caption = "Bereavement indicators in a female-dominant population under 2015 demographic conditions")

# end ---------------------------------------------------------------------

base_url <- 'https://population.un.org/dataportalapi/api/v1'
locations <- read.csv(paste0(base_url,'/locations/?format=csv'), sep='|', skip=1)
codes <- read.csv(paste0(base_url,'/indicators/?format=csv'), sep='|', skip=1) 
qx_code <- codes$Id[codes$ShortName == "qx1"]
dx_code <- codes$Id[codes$ShortName == "Deaths1"]
asfr_code <- codes$Id[codes$ShortName == "ASFR1"]
tfr_code <- codes$Id[codes$ShortName == "TFR5"]
e0_code <- codes$Id[codes$ShortName == "E0"]
mac_code <- codes$Id[codes$ShortName == "MAC5"]



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