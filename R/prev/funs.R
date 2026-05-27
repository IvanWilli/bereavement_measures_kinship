# functions ----------------------------------------------------------------

# main and other indicators
bereavement_indicators <- function(kin_net){
  
  # load("data/kin_net_country_years.rda")
  # kin_net <- kin_net_country_years %>%
  #   filter(country == "Argentina", year == 2025)

  x_focal <- 0
  e0_focal <- kin_net %>% filter(age_focal == x_focal) %>% distinct(ex_focal) %>% pull()
  
  # cumulative % loss and mean age of loss
  x_star_focal <- kin_net %>% filter(age_focal==x_focal) %>% pull(ex_focal) %>% round() %>% unique()
  
  D_e0 <- kin_net %>% 
    filter(age_focal <= e0_focal) %>%
    summarise(l = sum(living), d = sum(dead), .by = c(kin, age_focal)) %>% 
    # porcentaje de muerte hasta los 50 años (completó nuevos kins, no harbá más, excepto ggd)
    summarise(L = l[age_focal == e0_focal], 
              D = sum(d[age_focal < e0_focal]),
              Total = L + D,
              porc_d = D/Total,
              mean_age = sum((d[age_focal < e0_focal]) * (age_focal[age_focal < e0_focal]))/sum(d[age_focal < e0_focal]),
              .by = c(kin)) %>% 
    arrange(kin) %>% 
    select(kin, D, L, Total, porc_d, mean_age) %>% 
    mutate(e0 = e0_focal, .before = 1)
  
  # S time loss share
  lx_focal <- kin_net %>% distinct(age_focal, px_focal)
  lx_focal$cum_px_focal <- cumprod(lx_focal$px_focal)
  lx_focal$lx_focal <-c(1, lx_focal$cum_px_focal[-nrow(lx_focal)]) 
  S <- kin_net %>%
    filter(age_focal >= x_focal) %>%
    left_join(lx_focal %>% select(age_focal, lx_focal)) %>% 
    arrange(kin, age_focal, x_focal) %>%
    select(kin, sex_kin, age_focal, age_kin, dead, lx_focal, exy) %>% 
    mutate(px_y = lx_focal/min(lx_focal[age_focal==x_focal]),
           S_x = dead * px_y * exy) %>%
    group_by(kin) %>%
    summarise(S = sum(S_x, na.rm = T))
  
  # O overlapping/unnatural: ggm>gm, gm>m, m>s, m>d, a>d, d>gd, gd>ggd, a>c, c>d
  O <- kin_net  %>% 
    filter(age_focal>=x_focal, age_focal<e0_focal) %>% 
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
    pivot_longer(ggm:n, names_to = "kin", values_to = "O")
  
  ## indicators
  bereav_indicators <- D_e0 %>% 
    left_join(S) %>% 
    left_join(O)
  
  # out
  return(bereav_indicators)
}








# # dying alone
# surv_px_t <- function(lx, t){
#   OAG <- length(lx)
#   px_t <- lx[(1+t):OAG]/lx[1:(OAG-t)]
#   px_t <- c(px_t, rep(0, t))
# }
# 
# kin_net <- kin_net_country_years %>%
#   filter(country == "Argentina", year == 2015)
# this_kin_net <- kin_net %>%
#   filter(age_focal == 50, kin %in% c("d", "m", "s")) %>% 
#   select(kin, age_focal, px_focal, sex_kin, age_kin, px_kin, living)
# k0_t <- map_df(1:50, function(t){
#   # t = 10
#   data.frame(x = 0:100, 
#              px = this_kin_net$px_kin, 
#              living = this_kin_net$living,
#              sex_kin = this_kin_net$sex_kin,
#              kin = this_kin_net$kin) %>% 
#     mutate(lx = c(1, cumprod(px))[1:101],
#            px_t = surv_px_t(lx, t),
#            qx_t = 1 - px_t, .by = c(kin, sex_kin)) %>% 
#     summarise(mean_qx_t = sum(qx_t * living)/sum(living),
#               total_living = sum(living),
#               prob_k0 = round(mean_qx_t^total_living, 2),
#               .by = c(kin, sex_kin)) %>% 
#     summarise(mean_qx_t = sum(mean_qx_t * total_living)/sum(total_living),
#               total_living = sum(total_living),
#               prob_k0 = round(mean_qx_t ^ total_living, 2),
#               .by = c(kin)) %>% 
#     mutate(t = t, age_focal = 50 + t - 1, .before = 1)
# }) 
# # k0_t %>% 
# #   mutate(age_focal = 50 + t - 1) %>% 
# #   ggplot(aes(age_focal, prob_k0, col = kin)) +
# #   geom_line()
# k0 <- k0_t %>% 
#   left_join(this_kin_net %>% 
#               filter(sex_kin == "f") %>% 
#               distinct(age_kin, px_kin) %>% 
#               filter(age_kin >= 50) %>%
#               mutate(qx = 1 - px_kin, 
#                      lx = c(1, cumprod(px_kin))[1:51],
#                      prob_muere = lx * qx) %>% 
#               select(age_focal = age_kin,
#                      prob_muere)) 
# K0 <- k0 %>% 
#   summarise(k0 = sum(prob_muere * prob_k0), .by = c(kin))