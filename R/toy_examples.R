# toy examples for presentation

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(geomtextpath)

# Shared time ------------------------------------------------
kin_lines <- data.frame(
  kin = c("Grandparent", "Grandparent", "Parent", "Parent", rep("Children", 3), "Sibling"),
  x_focal     = c(-15, -15,     -15, -15,            20, 25, 30,            5),
  x_kin       = c(40,  35,       5,   10,            0,  0,  0,               0),
  x_focal_end = c(-5,  10,       30,   50,           21, 100, 90,           70),
  x_kin_end   = c(50,  60,       50,   75,           1,  75, 60,              65),
  death       = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
)
plot_kin_lines <- ggplot(kin_lines) +
  geom_segment(aes(x = x_focal, y = x_kin, xend = x_focal_end, yend = x_kin_end,
                   col = kin, linetype = kin)) +
  coord_equal() +
  geom_point(data = . %>% filter(death),
             aes(x = x_focal_end, y = x_kin_end), shape = 3) +
  geom_vline(xintercept = 0, linetype=2, col =1) +
  coord_equal() +
  geom_abline(intercept = 0, slope = 1, linetype = 2, col = "grey") +
  scale_x_continuous(expand = c(0,0), name = "Years since Focal born",
                     breaks = seq(-50, 100, 10), labels = seq(-50, 100, 10)) +
  scale_y_continuous(expand = c(0,0), name = "Age Kin", limits = c(0,100),
                     breaks = seq(0, 100, 10), labels = seq(0, 100, 10)) +
  scale_color_manual(values = c(2, 1, 3, 5, 6), name = "") +
  scale_linetype_manual(values = c(Grandparent = 1, Parent = 1, Children = 1, Focal = 2, Sibling = 1),
                        guide = "none") +
  
  theme_bw()
ggsave(filename = "plots/plot_kin_lines.pdf",
       plot = plot_kin_lines, 
       width = 8, height = 6)

# Overlapping ------------------------------------------------

# gammas
x <- seq(0, 90, by = 0.1)
dgamma_mean_sd <- function(x, mean, sd) {
  shape <- (mean / sd)^2
  scale <- sd^2 / mean
  dgamma(x, shape = shape, scale = scale)
}

# function to make data for each panel
make_panel_data <- function(mu_m, mu_gm, sd_m, sd_gm, panel_name) {
  
  wide <- tibble(
    x = x,
    M  = dgamma_mean_sd(x, mean = mu_m,  sd = sd_m),
    GM = dgamma_mean_sd(x, mean = mu_gm, sd = sd_gm),
    overlap = pmin(M, GM),
    total = M + GM,
    panel = panel_name
  )
  
  long <- wide %>%
    select(x, M, GM, panel) %>%
    pivot_longer(cols = c(M, GM), names_to = "kin", values_to = "density")
  
  list(wide = wide, long = long)
}

# Panel 1: original distributions
p1_dat <- make_panel_data(
  mu_m = 25, mu_gm = 10,
  sd_m = 11.5, sd_gm = 10,
  panel_name = " Pre-transition"
)

# Panel 2: shifted + more concentrated
p2_dat <- make_panel_data(
  mu_m = 60, mu_gm = 25,
  sd_m = 7, sd_gm = 6,
  panel_name = "Post-transition"
)

wide_all <- bind_rows(p1_dat$wide, p2_dat$wide)
long_all <- bind_rows(p1_dat$long, p2_dat$long)

sum(wide_all$overlap[wide_all$panel == " Pre-transition"])*2/sum( wide_all$total[wide_all$panel == " Pre-transition"])
sum(wide_all$overlap[wide_all$panel == "Post-transition"])*2/sum( wide_all$total[wide_all$panel == "Post-transition"])

Oplot_toy <- ggplot() +
  # filled densities
  # overlap
  geom_area(
    data = wide_all,
    aes(x = x, y = 2*overlap),
    fill = "grey50",
    alpha = 0.7
  ) +
  # individual density lines
  geom_textline(
    data = long_all,
    aes(x = x, y = density, color = kin, label = kin),
    linewidth = 1.1
  ) +
  # total line = sum of both
  geom_textline(
    data = wide_all,
    aes(x = x, y = total),
    color = "black",
    linetype = 2, label = "Total",
    linewidth = 1
  ) +
  facet_wrap(~panel, nrow = 1) +
  scale_color_manual(values = c(M = "steelblue4", GM = "firebrick3")) +
  scale_fill_manual(values = c(M = "steelblue2", GM = "tomato")) +
  labs(
    x = "Focal Age",
    y = "Deaths",
    color = NULL,
    fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )
Oplot_toy
ggsave(filename = "plots/Oplot_toy.pdf",
       plot = Oplot_toy, 
       width = 8, height = 6)

# Burden ------------------------------------------------

library(ggplot2)
library(patchwork)
library(ggplot2)
library(grid)

rects <- data.frame(
  xmin = c(15, 20, 25),
  xmax = c(60, 60, 60),
  ymin = c(0, 1, 2),
  ymax = c(1, 2, 3)
)

shade <- data.frame(
  xmin = 20,
  xmax = 60,
  ymin = 0,
  ymax = 1
)

p_soon_few <- ggplot() +
  geom_rect(
    data = shade,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey80",
    color = NA
  ) +
  geom_rect(
    data = rects,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  coord_cartesian(xlim = c(0, 62), ylim = c(0, 3.4), clip = "off") +
  theme_void() +
  annotate(
    "segment", x = 0, xend = 62, y = 0, yend = 0,
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.7
  ) +
  annotate(
    "segment", x = 0, xend = 0, y = 0, yend = 3.4,
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.7
  ) +
  annotate("text", x = c(15, 20, 45, 60), y = -0.12,
           label = c("15", "20", "45", "60"), size = 5) +
  annotate("text", x = -1.2, y = c(0.5, 1.5, 2.5),
           label = c("1", "2", "3"), size = 5) +
  annotate("text", x = -5, y = 1.7, label = "Children", angle = 90, size = 5) +
  annotate("text", x = 10, y = 2.5, label = "40/120=.33", size = 5) +
  labs(title = "Soon-Few") +
  theme(plot.title = element_text(hjust = 0.5, vjust =-5, size = 20, margin = margin(b = 0)),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 24))

shade <- data.frame(
  xmin = c(20, 25),
  xmax = 60,
  ymin = c(0, 1),
  ymax = c(1, 2)
)

p_soon_many <- ggplot() +
  geom_rect(
    data = shade,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey80",
    color = NA
  ) +
  geom_rect(
    data = rects,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  coord_cartesian(xlim = c(0, 62), ylim = c(0, 3.4), clip = "off") +
  theme_void() +
  annotate(
    "segment", x = 0, xend = 62, y = 0, yend = 0,
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.7
  ) +
  annotate(
    "segment", x = 0, xend = 0, y = 0, yend = 3.4,
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.7
  ) +
  annotate("text", x = c(15, 20, 25, 60), y = -0.12,
           label = c("15", "20", "25", "60"), size = 5) +
  annotate("text", x = 31, y = -0.35, label = "Age of Focal", size = 5) +
  annotate("text", x = -1.2, y = c(0.5, 1.5, 2.5),
           label = c("1", "2", "3"), size = 5) +
  annotate("text", x = -5, y = 1.7, label = "Children", angle = 90, size = 5) +
  annotate("text", x = 10, y = 2.5, label = "(40+35)/120=.63", size = 5) +
  labs(title = "Soon-Many") +
  theme(plot.title = element_text(hjust = 0.5, vjust =-5, size = 20, margin = margin(b = 0)),
        plot.margin = margin(t = 5.5, r = 5.5, b = 24, l = 24))

shade <- data.frame(
  xmin = 45,
  xmax = 60,
  ymin = 0,
  ymax = 1
)

p_later_few <- ggplot() +
  geom_rect(
    data = shade,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey80",
    color = NA
  ) +
  geom_rect(
    data = rects,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  coord_cartesian(xlim = c(0, 62), ylim = c(0, 3.4), clip = "off") +
  theme_void() +
  annotate(
    "segment", x = 0, xend = 62, y = 0, yend = 0,
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.7
  ) +
  annotate(
    "segment", x = 0, xend = 0, y = 0, yend = 3.4,
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.7
  ) +
  annotate("text", x = c(15, 20, 45, 60), y = -0.12,
           label = c("15", "20", "45", "60"), size = 5) +
  annotate("text", x = -1.2, y = c(0.5, 1.5, 2.5),
           label = c("1", "2", "3"), size = 5) +
  annotate("text", x = 10, y = 2.5, label = "15/120=.13", size = 5) +
  labs(title = "Late-Few") +
  theme(plot.title = element_text(hjust = 0.5, vjust =-5, size = 20, margin = margin(b = 0)))

shade <- data.frame(
  xmin = c(45, 50),
  xmax = 60,
  ymin = c(0, 1),
  ymax = c(1, 2)
)

p_later_many <- ggplot() +
  geom_rect(
    data = shade,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey80",
    color = NA
  ) +
  geom_rect(
    data = rects,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  coord_cartesian(xlim = c(0, 62), ylim = c(0, 3.4), clip = "off") +
  theme_void() +
  annotate(
    "segment", x = 0, xend = 62, y = 0, yend = 0,
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.7
  ) +
  annotate(
    "segment", x = 0, xend = 0, y = 0, yend = 3.4,
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.7
  ) +
  annotate("text", x = c(15, 20, 45, 50, 60), y = -0.12,
           label = c("15", "20", "45", "50", "60"), size = 5) +
  annotate("text", x = 31, y = -0.35, label = "Age of Focal", size = 5) +
  annotate("text", x = -1.2, y = c(0.5, 1.5, 2.5),
           label = c("1", "2", "3"), size = 5) +
  annotate("text", x = 10, y = 2.5, label = "(15+10)/120=.21", size = 5) +
  labs(title = "Late-Many") +
  theme(plot.title = element_text(hjust = 0.5, vjust =-5, size = 20, margin = margin(b = 0)),
        plot.margin = margin(t = 5.5, r = 5.5, b = 24, l = 5.5))

p1 <- p_soon_few
p2 <- p_later_few
p3 <- p_soon_many
p4 <- p_later_many

final_plot <- (p1 + p2) / (p3 + p4) + plot_layout(nrow = 2)

ggsave(
  filename = "plots/child_presence_quadrants.pdf",
  plot = final_plot, scale = 1.5,
  width = 8,
  height = 6
)
