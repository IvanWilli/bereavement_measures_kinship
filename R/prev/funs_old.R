# functions ----------------------------------------------------------------

# main and other indicators
bereavement_indicators <- function(kin_net, x_focal, age_unexp, support_group){
  
  # load("data/kin_net_country_years.rda")
  # kin_net <- kin_net_country_years %>%
  #   filter(country == "Argentina", year == 1950)
  # x_focal <- 30
  # age_unexp <- 30
  # support_group <- c("s","m")
  
  # define age and get indicators
  ex_focal <- kin_net %>% filter(age_focal==x_focal) %>% pull(ex_focal) %>% unique
  lx_focal <- kin_net %>% distinct(age_focal, px_focal)
  lx_focal$cum_px_focal <- cumprod(lx_focal$px_focal)
  lx_focal$lx_focal <-c(1, lx_focal$cum_px_focal[-nrow(lx_focal)])
  
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
  omega <- max(kin_net$age_focal)
  M_prosp <- kin_net %>% 
    filter(age_focal >= x_focal) %>% 
    summarise(M_prosp = 1 - sum(dead*(age_focal-x_focal))/sum(dead)/(omega-x_focal), 
              .by = c(kin))
  
  if(x_focal > 0){
    M_retr <- kin_net %>% 
      filter(age_focal < x_focal) %>% 
      summarise(M_retr = 1-sum(dead*age_focal)/sum(dead)/x_focal, 
                .by = c(kin))
  }else{
    M_retr <- M_prosp %>% mutate(M_retr = 0) %>% select(-M_prosp)
  }
  
  # U Non-expected proportion of deaths: age definition
  omega <- max(kin_net$age_focal)
  U_prosp <- kin_net  %>% 
    filter(age_focal > x_focal) %>% 
    summarise(U_prosp = sum(dead[age_kin<age_unexp], na.rm = T)/sum(dead, na.rm = T), .by = kin) 
  
  if(x_focal > 0){
    U_retr <- kin_net  %>% 
      filter(age_focal < x_focal) %>% 
      summarise(U_retr = sum(dead[age_kin<age_unexp], na.rm = T)/sum(dead, na.rm = T), .by = kin)
  }else{
    U_retr <- U_prosp %>% mutate(U_retr = 0) %>% select(-U_prosp)
  }
  
  # O overlapping/unnatural: ggm>gm, gm>m, m>s, m>d, a>d, d>gd, gd>ggd, a>c, c>d
  O_prosp <- kin_net  %>% 
    filter(age_focal>=x_focal) %>% 
    summarise(d = sum(dead), .by = c(kin, age_focal)) %>% 
    pivot_wider(names_from = kin, values_from = d) %>% 
    rowwise() %>% 
    mutate(
      gm_ggm_overlap = min(ggm, gm), gm_ggm_total = sum(ggm, gm),
      m_gm_overlap = min(gm, m),     m_gm_total = sum(gm, m),
      a_gm_overlap = min(gm, a),     a_gm_total = sum(gm, a),
      s_m_overlap = min(s, m),       s_m_total = sum(s, m),
      d_m_overlap = min(d, m),       d_m_total = sum(d, m),
      gd_d_overlap = min(d, gd),     gd_d_total = sum(d, gd),
      ggd_gd_overlap = min(gd, ggd), ggd_gd_total = sum(gd, ggd),
      c_a_overlap = min(a, c),       c_a_total = sum(c, a),
      n_s_overlap = min(n, s),       n_s_total = sum(n, s)) %>% 
    ungroup %>% 
    summarise(
      ggm = 0,
      gm = sum(gm_ggm_overlap)/sum(gm_ggm_total),
      m = sum(m_gm_overlap)/sum(m_gm_total),
      a = sum(a_gm_overlap)/sum(a_gm_total),
      s = sum(s_m_overlap)/sum(s_m_total),
      d = sum(d_m_overlap)/sum(d_m_total),
      gd = sum(gd_d_overlap)/sum(gd_d_total),
      ggd = sum(ggd_gd_overlap)/sum(ggd_gd_total),
      c = sum(c_a_overlap)/sum(c_a_total),
      n = sum(n_s_overlap)/sum(n_s_total)
    ) %>% 
    pivot_longer(ggm:n, names_to = "kin", values_to = "O_prosp")
  
  if(x_focal > 0){
    O_retr <- kin_net  %>% 
      filter(age_focal < x_focal) %>% 
      summarise(d = sum(dead), .by = c(kin, age_focal)) %>% 
      pivot_wider(names_from = kin, values_from = d) %>% 
      rowwise() %>% 
      mutate(
        gm_ggm_overlap = min(ggm, gm), gm_ggm_total = sum(ggm, gm),
        m_gm_overlap = min(gm, m),     m_gm_total = sum(gm, m),
        a_gm_overlap = min(gm, a),     a_gm_total = sum(gm, a),
        s_m_overlap = min(s, m),       s_m_total = sum(s, m),
        d_m_overlap = min(d, m),       d_m_total = sum(d, m),
        gd_d_overlap = min(d, gd),     gd_d_total = sum(d, gd),
        ggd_gd_overlap = min(gd, ggd), ggd_gd_total = sum(gd, ggd),
        c_a_overlap = min(a, c),       c_a_total = sum(c, a),
        n_s_overlap = min(n, s),       n_s_total = sum(n, s)) %>% 
      ungroup %>% 
      summarise(
        ggm = 0,
        gm = sum(gm_ggm_overlap)/sum(gm_ggm_total),
        m = sum(m_gm_overlap)/sum(m_gm_total),
        a = sum(a_gm_overlap)/sum(a_gm_total),
        s = sum(s_m_overlap)/sum(s_m_total),
        d = sum(d_m_overlap)/sum(d_m_total),
        gd = sum(gd_d_overlap)/sum(gd_d_total),
        ggd = sum(ggd_gd_overlap)/sum(ggd_gd_total),
        c = sum(c_a_overlap)/sum(c_a_total),
        n = sum(n_s_overlap)/sum(n_s_total)
      ) %>% 
      pivot_longer(ggm:n, names_to = "kin", values_to = "O_retr")
    
  }else{
    O_retr <- O_prosp %>% mutate(O_retr = 0) %>% select(-O_prosp)
  }
  
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
    arrange(L) %>% slice(1)
  L <- expand.grid(kin = kin_types, L = L$age_focal, L_living = L$L)
  
  ## indicators
  bereav_indicators <- T. %>% 
    left_join(S) %>% 
    left_join(M_retr) %>% 
    left_join(M_prosp) %>% 
    left_join(U_retr) %>% 
    left_join(U_prosp) %>%
    left_join(O_retr) %>% 
    left_join(O_prosp) %>% 
    left_join(P) %>% 
    left_join(H) %>% 
    left_join(L)
  
  # out
  return(bereav_indicators)
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
