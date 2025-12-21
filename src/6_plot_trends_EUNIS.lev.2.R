################################################################################
# Author: Gabriele Midolo
# Email: midolo@fzp.czu.cz
# Date: 17.12.2025
################################################################################

# Description: Get summary figure for EIV trends across EUNIS level 2 habitats

################################################################################

# 1. prepare data####

# load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

# set indicator names
ind.names <- data.frame(
  ind.name = c('EIV_L', 'EIV_T', 'EIV_M', 'EIV_N', 'EIV_R'),
  eiv_name_long = c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction')
)

# define paths
pth2preds <- './preds/pdp/'
pth2fig   <- './fig/'

pdp.data <- paste0(pth2preds, 'ESy2_alltrees.csv') %>%
  read_csv(show_col_types = FALSE) %>%
  left_join(ind.names, by = 'ind.name')

EIV_long_order <- c('Light', 'Temperature', 'Moisture', 'Nitrogen', 'Reaction')
ESy2_order <- c(
  'T1','T2','T3', # forests
  'R1','R2','R3','R4','R5','R6', # grasslands
  'S1','S2','S3','S4','S5','S6','S7','S9', # scrub
  'Q1','Q2','Q4','Q5','Q6' # wetland
)

# keep ESy2 ordering as present in the data
pdp.data$ESy2 <- factor(pdp.data$ESy2, levels = unique(pdp.data$ESy2))

# standardization to the baseline year
pdp.data_standardized <- pdp.data %>%
  group_by(ind.name, ESy2) %>%
  mutate(std_offset = mean[year == min(year, na.rm = TRUE)]) %>%
  mutate(
    std_mean          = mean - std_offset,
    std_lower_ci_95   = lower_ci_95 - std_offset,
    std_upper_ci_95   = upper_ci_95 - std_offset
  ) %>%
  ungroup() %>%
  mutate(eiv_name_long = factor(eiv_name_long, levels = EIV_long_order)) %>%
  select(eiv_name_long, ESy2, year, std_mean, std_lower_ci_95, std_upper_ci_95)

# get no. of plots (count)
count_data <- pdp.data %>%
  distinct(ind.name, ESy2, no.plots) %>%
  mutate(label = prettyNum(no.plots, big.mark = ',', scientific = FALSE)) %>%
  left_join(ind.names, by = "ind.name") %>%
  mutate(
    eiv_name_long = factor(eiv_name_long, levels = EIV_long_order),
    ESy2 = factor(ESy2, levels = ESy2_order)
  )

# calculate slope
slope_data <- pdp.data %>%
  mutate(decade = year / 10) %>%
  group_by(ind.name, ESy2) %>%
  do(tidy(lm(mean ~ decade, data = .))) %>%
  filter(term == 'decade') %>%
  rename(slope = estimate) %>%
  select(ind.name, ESy2, slope) %>%
  ungroup()

# categorize slopes
break_levels <- c(
  '\u2264 −0.04',
  '−0.04 to −0.02',
  '−0.02 to −0.01',
  '−0.01 to 0.01',
  '0.01 to 0.02',
  '0.02 to 0.04',
  '\u2265 0.04'
)
slope_data_categorized <- slope_data %>%
  left_join(ind.names, by = 'ind.name') %>%
  mutate(
    ESy2 = factor(ESy2, levels = levels(pdp.data$ESy2)),
    eiv_name_long = factor(eiv_name_long, levels = EIV_long_order),
    slope_category = case_when(
      slope <= -0.04                  ~ '\u2264 −0.04',
      slope >  -0.04 & slope <= -0.02 ~ '−0.04 to −0.02',
      slope >  -0.02 & slope <= -0.01 ~ '−0.02 to −0.01',
      slope >  -0.01 & slope <   0.01 ~ '−0.01 to 0.01',
      slope >=  0.01 & slope <   0.02 ~ '0.01 to 0.02',
      slope >=  0.02 & slope <   0.04 ~ '0.02 to 0.04',
      slope >=  0.04                  ~ '\u2265 0.04'
    ),
    slope_category = factor(slope_category, levels = break_levels)
  )

# Colors 
color_key <- hcl.colors(7, palette = 'RdYlBu', rev = TRUE)
names(color_key) <- break_levels


# 2. Export slope trends ####
library(flextable)
library(officer)

# prepare tab
tab_dat <- slope_data_categorized %>%
  left_join(
    count_data %>% rename(`no. plots` = label),
    by = c('ind.name', 'ESy2', 'eiv_name_long')
  ) %>%
  select(ESy2, eiv_name_long, `no. plots`, slope, slope_category) %>%
  mutate(
    eiv_name_long = factor(eiv_name_long, levels = EIV_long_order),
    slope = round(slope, 3),
    ESy2 = factor(ESy2, levels = ESy2_order)
  ) %>%
  arrange(ESy2, eiv_name_long) %>%
  left_join(
    read_delim('./data/EUNIS_ESy2_habitat.names.txt', col_names = F, show_col_types = F) %>%
      setNames(c('ESy2', 'ESy2_name'))
  ) %>%
  select(ESy2_name, everything()) %>%
  group_by(ESy2_name) %>%
  mutate(is_last_in_group = row_number() == n()) %>%
  ungroup() %>%
  left_join(pdp.data %>% 
              group_by(ESy2) %>% 
              filter(year == min(year)) %>% 
              select(ESy2, year) %>% 
              mutate(year = as.character(year)) %>% 
              rename(`Start Year` = year) %>%
              unique()) %>%
  left_join(pdp.data %>% 
              group_by(ESy2) %>% 
              filter(year == max(year)) %>% 
              select(ESy2, year) %>% 
              mutate(year = as.character(year)) %>% 
              rename(`End Year` = year) %>% 
              unique())

# prepare flextable (Added Start Year and End Year to col_keys)
ft <- tab_dat %>%
  flextable(col_keys = c("ESy2_name", "ESy2", "Start Year", "End Year", "eiv_name_long", "no. plots", "slope")) %>%
  set_header_labels(
    ESy2 = 'EUNIS code',
    ESy2_name = 'EUNIS name',
    eiv_name_long = 'EIV',
    `no. plots` = 'No. plots',
    slope = 'CM[EIV]\n/decade',
    `Start Year` = 'Start Year',
    `End Year` = 'End Year'
  ) %>%
  colformat_num(j = 'slope', digits = 4) %>%
  align(align = 'center', part = 'all') %>%
  set_table_properties(layout = "autofit", width = 1)

# apply colors based on slope category
for (lvl in names(color_key)) {
  ft <- ft %>%
    bg(i = ~ slope_category == lvl, j = 'slope', bg = color_key[[lvl]], part = 'body')
}

# vertical merge and horizontal lines
std_border <- fp_border(color = "black", width = 1)

ft <- ft %>% 
  merge_v(j = c('ESy2_name', 'ESy2')) %>%
  valign(j = c('ESy2_name', 'ESy2'), valign = 'center') %>%
  hline(i = ~ is_last_in_group == TRUE, part = 'body', border = std_border) %>%
  color(j = 'slope', color = 'black', part = 'body') %>%
  fix_border_issues()

# export
read_docx() %>%
  body_add_flextable(ft) %>%
  print(target = paste0(pth2fig, 'EUNIS.ESy2_slope_trends.docx'))


# 3. Plot CM[EIV] trends over time (pratial dependence) per EUNIS habitat####
y_breaks <- seq(-0.3, 0.3, 0.15)
y_limits_centrd <- range(y_breaks)

pdp.data_standardized$ESy2 <- factor(pdp.data_standardized$ESy2, levels = ESy2_order)

p <- ggplot(pdp.data_standardized, aes(year, std_mean)) +
  geom_rect(
    data = slope_data_categorized,
    aes(
      xmin = -Inf, xmax = Inf,
      ymin = y_limits_centrd[1],
      ymax = y_limits_centrd[2],
      fill = slope_category
    ),
    inherit.aes = FALSE,
    alpha = 0.8
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 3, color = 'grey70') +
  geom_ribbon(
    aes(ymin = std_lower_ci_95, ymax = std_upper_ci_95),
    alpha = 0.6, fill = 'grey70'
  ) +
  geom_line(linewidth = 0.6) +
  geom_text(
    data = count_data,
    aes(x = 1990, y = -0.22, label = label),
    hjust = .5, vjust = 1, size = 2.4
  ) +
  facet_grid(eiv_name_long ~ ESy2) +
  theme_bw() +
  scale_y_continuous(
    limits = y_limits_centrd,
    breaks = c(-0.20 ,0.00, 0.20),
    name = expression(CM[EIV]~'shift from baseline year')
  ) +
  scale_x_continuous(
    breaks = c(1960, 1990, 2020),
    name = NULL
  ) +
  scale_fill_manual(
    values = color_key,
    breaks = names(color_key),
    name = expression(Delta*CM[EIV]~'per decade slope')
  ) +
  theme(
    strip.text.y.left = element_text(angle = 270, vjust = 0.5, hjust = 0.5, face = 'bold'),
    strip.text.x = element_text(face = 'bold'),
    strip.background.y = element_blank(),
    strip.background.x = element_rect(fill = 'white'),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.1, 'lines'),
    # panel.border = element_blank()
    legend.background = element_rect(color = "#575757", linewidth = 0.4)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# print
print(p)

# export
ggsave(paste0(pth2fig, 'ESy2trends.pdf'), p, width = 11, height = 4.76, device = cairo_pdf)