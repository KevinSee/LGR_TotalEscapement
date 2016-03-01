# Author: Kevin See
# Purpose: Examine results of total adult escapement model over Lower Granite dam 
# Created: 2/19/2015
# Last Modified: 2/25/2016

library(R2jags)
library(ggplot2)
library(gridExtra)
library(dclone)
library(boot)
library(MCMCglmm)
library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggmcmc)

theme_set(theme_bw())

#-----------------------------------------------------------------
# Plotting functions
#-----------------------------------------------------------------
HPDI_boxplot = function(x, middle.def='Mode') {
  summ = data.frame(
    ymin = max(0, HPDinterval(as.mcmc(x))[,1]),
    lower = HPDinterval(as.mcmc(x), prob=0.5)[,1],
    upper = HPDinterval(as.mcmc(x), prob=0.5)[,2],
    ymax = HPDinterval(as.mcmc(x))[,2])
  if(middle.def=='Mode') summ$middle = posterior.mode(as.mcmc(x))
  if(middle.def=='Median') summ$middle = quantile(x, 0.5)
  if(middle.def=='Mean') summ$middle = mean(x)
  return(summ)
}

HPDI_outlier = function(x) {
  outs = subset(x, x < HPDinterval(as.mcmc(x))[,1] | x > HPDinterval(as.mcmc(x))[,2])
  if(diff(range(x))==0) outs = NA
  return(outs)
}

#-----------------------------------------------------------------
# load prepped data
#-----------------------------------------------------------------
# load the data for all years
load(file = 'DataPrepped.rda')
# load the data list prepped for JAGS
load(file = 'JAGS_DataList.rda')


#-----------------------------------------------------------------
# model diagnostics
#-----------------------------------------------------------------
names(jags_data_list)
i = 5
load(paste('ModelFits/LGD_', gsub('\\.', '_', names(jags_data_list)[i]), '.rda', sep=''))

my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('X.sig'))
my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('X.tot'))
my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('trap.rate.true'))
my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('reasc.avg'))
my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('reasc.true'))
my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('win.prop.avg'))
my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('win.prop.true'))
my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('win.prop.sigma'))
my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('wnc.avg'),
             par_labels = data.frame(Parameter = paste0('wnc.avg[', 1:2, ']'),
                                     Label = c('Perc Wild Morph', 'Perc Wild PBT')))
my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('wnc.true'),
             par_labels = data.frame(Parameter = paste0('wnc.true[', 1:2, ']'),
                                     Label = c('Perc Wild Morph', 'Perc Wild PBT')))

my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('acf'),
             par_labels = data.frame(Parameter = paste0('acf[', 1:4, ']'),
                                     Label = c('Day Rate', 'Re-Ascen Rate', 'Wild Rate (Morph)', 'Wild Rate (PBT)')))
my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('sigma'))
my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('trap.bump'))
my_ggs = ggs(as.mcmc(adult.pass.mod), family = c('^r$'))

ggs_density(my_ggs) +
  facet_wrap(~ Parameter, scales = 'free')
ggs_traceplot(my_ggs) +
  facet_wrap(~ Parameter, scales = 'free')
ggs_Rhat(my_ggs)
ggs_geweke(my_ggs)
ggs_caterpillar(my_ggs)
ggs_autocorrelation(my_ggs)
ggs_crosscorrelation(my_ggs)
ggs_running(my_ggs)
ggs_compare_partial(my_ggs)
# ggs_ppmean(my_ggs)

#----------------------------------------------------------------------
# combine results from all years
#----------------------------------------------------------------------
names(jags_data_list)
tot_post_list = vector('list', length(jags_data_list))
names(tot_post_list) = names(jags_data_list)
other_post_list = week_post_list = tot_post_list

for(i in 1:length(jags_data_list)) {
  load(paste('ModelFits/LGD_', gsub('\\.', '_', names(jags_data_list)[i]), '.rda', sep=''))
  attach.jags(adult.pass.mod)
  tot_post_list[[i]] = ldply(list('All.Fish' = X.tot.all, 
                                  'All.Wild.Fish.Morph' = X.tot.all.wild.morph, 
                                  'All.Wild.Fish.PBT' = X.tot.all.wild.pbt, 
                                  'Unique.Wild.Fish.Morph' = X.tot.new.wild.morph, 
                                  'Unique.Wild.Fish.PBT' = X.tot.new.wild.pbt, 
                                  'Daytime.Fish' = X.tot.day, 
                                  'Reascent.Fish' = X.tot.reasc, 
                                  'Night.Wild.Fish.Morph' = X.tot.night.wild.morph, 
                                  'Night.Wild.Fish.PBT' = X.tot.night.wild.pbt, 
                                  'Wild.Reascents.Morph' = X.tot.reasc.wild.morph, 
                                  'Wild.Reascents.PBT'= X.tot.reasc.wild.pbt), 
                             .id='Variable') %>% tbl_df() %>%
    rename(value = `1`)
  
  other_post_list[[i]] = ldply(list('Tot.Rand.Walk.SD' = X.sigma,
                                    'Avg.Day' = win.prop.avg,
                                    'SD.Day' = win.prop.sigma,
                                    'Day.Rate.AutoCorr' = acf[,1] %>% as.matrix(),
                                    'Avg.ReAsc' = reasc.avg,
                                    'SD.ReAsc' = reasc.sigma,
                                    'ReAsc.AutoCorr' = acf[,2] %>% as.matrix(),
                                    'Avg.Prop.Wild.Morph' = wnc.avg[,1] %>% as.matrix(),
                                    'Prop.Wild.Morph.AutoCorr' = acf[,3] %>% as.matrix(),
                                    'Avg.Prop.Wild.PBT' = wnc.avg[,2] %>% as.matrix(),
                                    'Prop.Wild.PBT.AutoCorr' = acf[,4] %>% as.matrix(),
                                    # 'Trap.Bump' = trap.bump,
                                    'Over.Dispersion' = r),
                               .id = 'Variable') %>% tbl_df() %>%
    rename(value = `1`)
  
  week_post_list[[i]] = ldply(list('New.Wild.Fish.PBT' = X.new.wild.pbt, 
                                   'All.Wild.Fish.PBT' = X.all.wild.pbt, 
                                   'New.Wild.Fish.Morph' = X.new.wild.morph, 
                                   'All.Wild.Fish.Morph' = X.all.wild.morph, 
                                   'All.Fish' = X.all, 
                                   'Daytime.Fish' = X.day,
                                   'Night.Wild.Fish' = X.night.wild.pbt,
                                   'Wild.Reascents' = X.reasc.wild.pbt, 
                                   'Trap.Rate' = trap.rate.true,
                                   'Day.Time.Rate' = true.prop, 
                                   'Night.Time.Rate' = 1 - true.prop, 
                                   'Re-Ascension.Rate' = reasc.true,
                                   'Wild.Morph.Rate' = wnc.true[,,1],
                                   'Wild.PBT.Rate' = wnc.true[,,2]), 
                              .id='Variable') %>% tbl_df() %>%
    gather(week_num, value, -Variable)
  detach.jags()
}

tot_post = ldply(tot_post_list, .id = 'Model') %>% tbl_df() %>%
  mutate(Model = as.character(Model),
         Species = sapply(strsplit(Model, '\\.'), function(x) x[1]),
         Year = sapply(strsplit(Model, '\\.'), function(x) x[2]),
         Year = as.integer(Year)) %>%
  select(Species, Year, Variable, value)

tot_summ = tot_post %>%
  group_by(Species, Year, Variable) %>%
  summarise(mean = mean(value),
            median = median(value),
            se = sd(value),
            cv = se / mean,
            low_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,1],
            upp_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,2]) %>%
  ungroup()

# compare total window counts, total trap estimates, model total
lgr_weekly %>%
  rename(Year = SpawnYear) %>%
  mutate(trap_est = ifelse(trap_est > 1e13, NA, trap_est)) %>%
  group_by(Species, Year) %>%
  summarise(tot_win = sum(win_cnt),
            tot_trap = floor(sum(trap_est, na.rm=T))) %>%
  left_join(tot_summ %>%
              filter(Variable == 'All.Fish') %>%
              select(Species, Year, median, low_ci, upp_ci) %>%
              mutate(median = round(median)))

other_post = ldply(other_post_list, .id = 'Model') %>% tbl_df() %>%
  mutate(Model = as.character(Model),
         Species = sapply(strsplit(Model, '\\.'), function(x) x[1]),
         Year = sapply(strsplit(Model, '\\.'), function(x) x[2]),
         Year = as.integer(Year)) %>%
  select(Species, Year, Variable, value)
other_summ =  other_post %>%
  group_by(Species, Year, Variable) %>%
  summarise(mean = mean(value),
            median = median(value),
            se = sd(value),
            cv = se / mean,
            low_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,1],
            upp_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,2]) %>%
  ungroup()

week_summ = ldply(week_post_list, .id = 'Model') %>% tbl_df() %>%
  mutate(Model = as.character(Model),
         Species = sapply(strsplit(Model, '\\.'), function(x) x[1]),
         Year = sapply(strsplit(Model, '\\.'), function(x) x[2]),
         Year = as.integer(Year)) %>%
  select(Species, Year, Variable, week_num, value) %>%
  group_by(Species, Year, Variable, week_num) %>%
  summarise(mean = mean(value),
            median = median(value),
            se = sd(value),
            cv = se / mean,
            low_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,1],
            upp_ci = HPDinterval(as.mcmc(value), prob = 0.95)[,2]) %>%
  ungroup() %>%
  mutate(week_num = as.integer(as.character(week_num))) %>%
  select(Species, Year, week_num, everything())

# pull out historical rate of daytime passage
hist_day_rate = ldply(jags_data_list, .id = 'Model', .fun = function(x) {
  hist_day_rate_month = with(x, day.fish / tot.fish) %>%
    rowMeans()
  hist_day_rate = hist_day_rate_month[x$month.vec]
  return(data.frame(week_num = 1:x$TotLadderWeeks, hist_day_rate))
}) %>% tbl_df() %>%
  rename(Year = SpawnYear)

# combine several model output and data input for plotting purposes
plot_df = left_join(week_summ,
                    lgr_weekly %>%
                      rename(Year = SpawnYear) %>%
                      group_by(Species, Year) %>%
                      mutate(week_num = 1:length(Year),
                             trap_est = ifelse(trap_est < 1e13, trap_est, NA)) %>%
                      ungroup()) %>%
  left_join(hist_day_rate)

#--------------------------------------------------
# make some plots
#--------------------------------------------------
# compare daytime estimates with window counts
plot_df %>%
  filter(Variable == 'Daytime.Fish') %>%
  mutate(inCI = ifelse(low_ci <= win_cnt & upp_ci >= win_cnt, T, F),
         bias = median - win_cnt,
         rel_err = ifelse(win_cnt > 0, abs(bias) / win_cnt, NA)) %>%
  select(Species, Year, week_num, win_cnt, low_ci, upp_ci, inCI, bias, rel_err) %>%
  group_by(Species, Year) %>%
  summarize(perc_inCI = round(sum(inCI) / n(), 2),
            avg_bias = round(mean(bias)),
            wt_avg_bias = round(weighted.mean(bias, w = win_cnt)),
            avg_abs_bias = round(mean(abs(bias))),
            wt_avg_abs_bias = round(weighted.mean(abs(bias), w = win_cnt)),
            avg_rel_err = round(mean(rel_err, na.rm = T), 3),
            wt_avg_rel_err = round(weighted.mean(rel_err, w = win_cnt, na.rm = T), 3))

# pick a species
spp = c('Chinook', 'Steelhead')[1]

# comparison of daytime estimates with window counts
plot_df %>%
  filter(Variable == 'Daytime.Fish',
         Species == spp) %>%
  ggplot(aes(x = win_cnt,
             y = median)) +
  geom_errorbar(aes(ymin = low_ci,
                    ymax = upp_ci)) +
  geom_point(size = 4) +
  geom_abline(intercept = 0,
              slope = 1,
              lty = 2,
              color = 'darkgray') +
  theme_bw() +
  geom_smooth(method = lm,
              formula = y ~ -1 + x,
              fullrange = T) +
  facet_wrap(~ Year, scales = 'free') +
  labs(x = 'Window Count',
       y = 'Model Estimate',
       title = paste(spp, 'Crossing During Day'))

# comparison of weekly estimates with trap estimates
plot_df %>%
  filter(Variable == 'All.Fish',
         Species == spp) %>%
  ggplot(aes(x = trap_est,
             y = median)) +
  geom_errorbar(aes(ymin = low_ci,
                    ymax = upp_ci)) +
  geom_errorbarh(aes(xmin = trap_est + qnorm(0.025) * trap_est_se,
                     xmax = trap_est + qnorm(0.975) * trap_est_se)) +
  geom_point(size = 4,
             aes(shape = trap_valid)) +
  scale_shape_manual(values = c('TRUE' = 19, 'FALSE' = 1),
                     name = 'Trap Valid') +
  geom_abline(intercept = 0,
              slope = 1,
              lty = 2,
              color = 'red') +
  
  geom_smooth(method = lm,
              formula = y ~ -1 + x,
              fullrange = T) +
  facet_wrap(~ Year, scales = 'free') +
  labs(x = 'Trap Estimate',
       y = 'Model Estimate',
       title = paste(spp, 'Crossing Each Week'))


# total escapement by week, with window, trap and model estimates
tot_p = plot_df %>%
  filter(Species == spp) %>%
  filter(Variable == 'All.Fish') %>%
  ggplot(aes(x = Start_Date, y = median)) +
  geom_ribbon(aes(ymin = low_ci, ymax = upp_ci), fill = 'lightgray') +
  geom_line(aes(y = win_cnt, color = 'Window')) +
  geom_point(aes(y = win_cnt, color = 'Window')) +
  geom_line(aes(y = trap_est, color = 'Trap')) +
  geom_point(aes(y = trap_est, 
                 color = 'Trap',
                 shape = trap_valid)) +
  scale_shape_manual(values = c('TRUE' = 19, 'FALSE' = 1),
                     name = 'Trap Valid') +
  geom_line(aes(color = 'Model')) +
  geom_point(aes(color = 'Model')) +
  scale_color_manual(values = c('Model' = 'black', 'Window' = 'blue', 'Trap' = 'red')) +
  theme_bw() +
  facet_wrap(~ Year, scales = 'free') +
  # scale_y_log10() +
  labs(x = 'Date', y = 'Estimate', color = 'Source', title = paste('All Fish -', spp))
tot_p

ggsave(paste0('Figures/TotEscape_TimeSeries_', spp, '.pdf'),
       tot_p,
       width = 9,
       height = 6,
       units = 'in')

plot_df %>%
  filter(Variable == 'All.Fish') %>%
  mutate(est_wind = floor(median - win_cnt),
         est_trap = floor(median - trap_est)) %>%
  group_by(Species, Year) %>%
  summarise(max_win_diff = max(est_wind, na.rm = T),
            mean_win_diff = round(mean(est_wind, na.rm = T)),
            min_win_diff = min(est_wind, na.rm = T),
            max_trap_diff = max(est_trap, na.rm = T),
            mean_trap_diff = round(mean(est_trap, na.rm = T)),
            min_trap_diff = min(est_trap, na.rm = T))
  
# look at rates of wild (morph & PBT) fish
plot_df %>%
  filter(Species == spp,
         Variable %in% c('Wild.Morph.Rate', 'Wild.PBT.Rate')) %>%
  mutate(wild_rate_morph = Wild.morph / trap_fish,
         wild_rate_pbt = Wild.PBT / trap_fish) %>%
  ggplot(aes(x = week_num, y = mean, color = Variable)) +
  geom_ribbon(aes(ymin = low_ci, ymax = upp_ci, group = Variable, fill = Variable), color = NA, alpha = 0.4) +
  geom_hline(data = filter(other_summ, 
                           Variable %in% c('Avg.Prop.Wild.Morph', 'Avg.Prop.Wild.PBT'),
                           Species == spp) %>%
               select(Year, median),
             aes(yintercept =  median), lty = 2, col = 'darkgreen') +
  geom_line(aes(group = Variable)) +
  geom_point(aes(group = Variable, size = win_cnt), position = position_dodge(width = 0.90)) +
  scale_color_manual(values = c('Wild.Morph.Rate' = 'blue4', 'Wild.PBT.Rate' = 'red'),
                     labels = c('Wild.Morph.Rate' = 'Morph', 'Wild.PBT.Rate' = 'PBT')) +
  scale_fill_manual(values = c('Wild.Morph.Rate' = 'lightskyblue', 'Wild.PBT.Rate' = 'mistyrose2'),
                    labels = c('Wild.Morph.Rate' = 'Morph', 'Wild.PBT.Rate' = 'PBT')) +
  geom_line(aes(y = wild_rate_morph), color = 'blue') +
  geom_point(aes(y = wild_rate_morph, size = win_cnt), color = 'blue') +
  geom_line(aes(y = wild_rate_pbt), color = 'hotpink') +
  geom_point(aes(y = wild_rate_pbt, size = win_cnt), color = 'hotpink') +
  theme_bw() +
  facet_wrap(~ Year) +
  labs(x = 'Week', y = 'Rate', title = paste('Rate of Wild Fish -', spp), size = 'Window Count')

# look at daytime passage rate
plot_df %>%
  filter(Species == spp,
         Variable == 'Day.Time.Rate') %>%
  mutate(day_rate = day_tags_W / tot_tags_W) %>%
  ggplot(aes(x = week_num, y = median)) +
  geom_ribbon(aes(ymin = low_ci, ymax = upp_ci), fill = 'lightgray') +
  geom_hline(data = filter(other_summ, 
                           Variable %in% c('Avg.Day'),
                           Species == spp) %>%
               select(Year, median),
             aes(yintercept =  median), lty = 2, col = 'darkgreen') +
  geom_line(aes(color = 'Model')) +
  geom_point(aes(color = 'Model')) +
  geom_line(aes(y = day_rate,
                color = 'Tags')) +
  geom_point(aes(y = day_rate, 
                 color = 'Tags',
                 size = tot_tags_W)) +
  geom_line(aes(y = hist_day_rate,
                color = 'Historical')) +
  geom_point(aes(y = hist_day_rate,
                 color = 'Historical')) +
  scale_color_manual(values = c('Historical' = 'orange',
                                'Tags' = 'blue',
                                'Model' = 'black')) +
  theme_bw() +
  facet_wrap(~ Year) +
  # coord_cartesian(ylim = c(0.7, 1.02)) +
  labs(x = 'Week', y = 'Rate', title = paste(spp, 'Daytime Passage Rate'), size = 'Wild Tags', color = 'Source')


# Look at posteriors of trap bump parameter
# other_post %>%
#   filter(Variable == 'Trap.Bump') %>%
#   mutate(Year = as.factor(Year)) %>%
#   ggplot(aes(x = value, color = Year, fill = Year)) +
#   geom_density(alpha = 0.2, lwd=1.5) +
#   theme_bw() +
#   scale_color_brewer(palette = 'Set1') +
#   scale_fill_brewer(palette = 'Set1') +
#   facet_wrap(~ Species, scales = 'free') +
#   labs(x = 'Posterior', title = 'Trap Bump')

# Look at posteriors of over dispersion parameter
# overdispersed if r is small, approximately Poisson if r is very large
other_post %>%
  filter(Variable == 'Over.Dispersion',
         value < 100) %>%
  mutate(Year = as.factor(Year)) %>%
  ggplot(aes(x = value, color = Year, fill = Year)) +
  geom_density(alpha = 0.2, lwd=1.5) +
  theme_bw() +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  facet_wrap(~ Species, scales = 'free') +
  labs(x = 'Posterior', title = 'Overdispersion')

# Look at posteriors of random walk variance
other_post %>%
  filter(Variable == 'Tot.Rand.Walk.SD') %>%
  mutate(Year = as.factor(Year)) %>%
  ggplot(aes(x = value, color = Year, fill = Year)) +
  geom_density(alpha = 0.2, lwd=1.5) +
  theme_bw() +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  facet_wrap(~ Species, scales = 'free') +
  labs(x = 'Posterior', title = 'Random Walk Std Dev.')

# Look at trap rate
plot_df %>%
  filter(Species == spp,
         Variable == 'Trap.Rate') %>%
  ggplot(aes(x = week_num, y = median)) +
  geom_ribbon(aes(ymin = Rate_MR - qnorm(0.975) * Rate_MR_se,
                  ymax = Rate_MR + qnorm(0.975) * Rate_MR_se), 
              fill = 'pink', alpha = 0.3) +
#   geom_ribbon(aes(ymin = trap_rate - qnorm(0.975) * trap_est_se,
#                   ymax = trap_rate + qnorm(0.975) * trap_est_se), 
#               fill = 'lightblue', alpha = 0.3) +
  geom_ribbon(aes(ymin = low_ci, ymax = upp_ci), fill = 'lightgray', alpha = 0.7) +
  geom_line(aes(y = Rate, color = 'Goal')) +
  # geom_point(aes(y = Rate, color = 'Goal')) +
  geom_line(aes(y = RateCalc, color = 'Time - Calculated')) +
  geom_point(aes(y = RateCalc, color = 'Time - Calculated', size = win_cnt)) +
  geom_line(aes(y = trap_rate, color = 'Mark Recapture')) +
  geom_point(aes(y = trap_rate, color = 'Mark Recapture', size = win_cnt)) +
#   geom_line(aes(y = trap_rate, color = 'Model Input')) +
#   geom_point(aes(y = trap_rate, color = 'Model Input', size = win_cnt)) +
  geom_line(aes(color = 'Estimate')) +
  geom_point(aes(color = 'Estimate')) +
  geom_point(aes(shape = trap_valid),
             color = 'red',
             size = 5) +
  scale_shape_manual(values = c('TRUE' = NA,
                                'FALSE' = 1)) +
  theme_bw() +
  scale_color_manual(values = c('Estimate' = 'black',
                                'Model Input' = 'blue',
                                'Time - Calculated' = 'orange',
                                'Mark Recapture' = 'red',
                                'Goal' = 'darkgreen'),
                     name = 'Rate') +
  facet_wrap(~ Year, scales = 'fixed') +
  labs(x = 'Week', 
       y = 'Trap Rate', 
       title = spp,
       size = 'Window Count', 
       shape = 'Valid Trap Rate')

# compare night-time ascension and re-acension totals
# for wild fish
plot_x_max = filter(tot_post,
                    Variable %in% c('Night.Wild.Fish.PBT', 'Wild.Reascents.PBT')) %>%
  # Variable %in% c('Night.Wild.Fish.Morph', 'Wild.Reascents.Morph')) %>%
  group_by(Species) %>%
  summarise(max_x = quantile(value, 0.999))
night_reasc_p = filter(tot_post,
                       Variable %in% c('Night.Wild.Fish.PBT', 'Wild.Reascents.PBT')) %>%
  # Variable %in% c('Night.Wild.Fish.Morph', 'Wild.Reascents.Morph')) %>%
  left_join(plot_x_max) %>%
  filter(value < max_x) %>%
  mutate(Variable = factor(Variable, levels = c('Night.Wild.Fish.PBT', 'Wild.Reascents.PBT'))) %>%
  ggplot(aes(x = value)) +
  geom_density(aes(color = Variable, fill = Variable), alpha=0.2, lwd=1.5) +
  facet_grid(Year ~ Species, scales = 'free') +
  scale_color_brewer(palette = 'Set1', name = 'Fish', labels = c('Night Passage', 'Re-ascents')) +
  scale_fill_brewer(palette = 'Set1', name = 'Fish', labels = c('Night Passage', 'Re-ascents')) +
  theme_bw() +
  theme(axis.text.y = element_blank()) +
  labs(x = 'Total Fish', 
       y = 'Density',
       title = 'Night Passage vs. Re-ascension\nWild Fish')
night_reasc_p

# by total fish
night_reasc_df = tot_post %>%
  filter(Variable %in% c('All.Fish', 'Daytime.Fish', 'Reascent.Fish')) %>%
  group_by(Variable) %>%
  mutate(iter = 1:n()) %>%
  ungroup() %>%
  spread(Variable, value) %>%
  mutate(Night.Fish = All.Fish - Daytime.Fish) %>%
  select(-(All.Fish:Daytime.Fish)) %>%
  gather(Variable, value, -(Species:iter)) %>%
  mutate(Variable = factor(Variable, levels = c('Night.Fish', 'Reascent.Fish')))

plot_x_max2 = filter(night_reasc_df,
                    Variable %in% c('Night.Fish', 'Reascent.Fish')) %>%
  group_by(Species) %>%
  summarise(max_x = quantile(value, 0.999))
night_reasc_p2 = filter(night_reasc_df,
                        Variable %in% c('Night.Fish', 'Reascent.Fish')) %>%
  left_join(plot_x_max2) %>%
  filter(value < max_x) %>%
  ggplot(aes(x = value)) +
  geom_density(aes(color = Variable, fill = Variable), alpha=0.2, lwd=1.5) +
  facet_grid(Year ~ Species, scale = 'free') +
  scale_color_brewer(palette = 'Set1', name = 'Fish', labels = c('Night Passage', 'Re-ascents')) +
  scale_fill_brewer(palette = 'Set1', name = 'Fish', labels = c('Night Passage', 'Re-ascents')) +
  theme_bw() +
  theme(axis.text.y = element_blank()) +
  labs(x = 'Total Fish', 
       y = 'Density',
       title = 'Night Passage vs. Re-ascension\nAll Fish')
night_reasc_p2

# for wild fish, by week
plot_df %>%
  filter(Species == spp,
         Variable == 'Wild.Reascents') %>%
  select(Species:week_num, Start_Date, ReAsc = median, ReAsc_lowCI = low_ci, ReAsc_uppCI = upp_ci) %>%
  full_join(plot_df %>%
              filter(Species == spp,
                     Variable == 'Night.Wild.Fish') %>%
              select(Species:week_num, Start_Date, Night = median, Night_lowCI = low_ci, Night_uppCI = upp_ci)) %>%
  full_join(lgr_weekly %>%
              filter(Species == spp) %>%
              select(Species, Year = SpawnYear, week_num, Start_Date, win_cnt, trap_est, trap_valid)) %>%
  ggplot(aes(x = week_num)) +
  geom_ribbon(aes(ymin = ReAsc_lowCI,
                  ymax = ReAsc_uppCI,
                  fill = 'reasc'), 
              alpha = 0.3) +
  geom_ribbon(aes(ymin = Night_lowCI,
                  ymax = Night_uppCI,
                  fill = 'night'), 
              alpha = 0.3) +
  geom_line(aes(y = ReAsc,
                color = 'reasc')) +
  geom_line(aes(y = Night,
                color = 'night')) +
  scale_fill_manual(values = c('reasc' = 'pink',
                               'night' = 'lightblue'),
                    labels = c('reasc' = 'Re-Ascension',
                               'night' = 'Night Passage'),
                    name = 'Variable') +
  scale_color_manual(values = c('reasc' = 'red',
                               'night' = 'blue'),
                     labels = c('reasc' = 'Re-Ascension',
                                'night' = 'Night Passage'),
                     name = 'Variable') +
  facet_wrap(~ Year, scales = 'free') +
  labs(x = 'Week',
       y = 'Fish',
       title = paste('Night-passage & Re-ascension - Wild', spp))

# difference between model estimate and window count / trap estimate, by week
plot_df %>%
  filter(Species == spp) %>%
  filter(Variable == 'All.Fish') %>%
  ggplot(aes(x = Start_Date)) +
  geom_ribbon(aes(ymin = low_ci - win_cnt,
                  ymax = upp_ci - win_cnt,
                  fill = 'Window'),
              alpha = 0.5) +
  geom_ribbon(aes(ymin = low_ci - trap_est,
                  ymax = upp_ci - trap_est,
                  fill = 'Trap'),
              alpha = 0.5) +
  geom_line(aes(y = median - win_cnt,
                color = 'Window')) +
  geom_line(aes(y = median - trap_est,
                color = 'Trap')) +
  geom_hline(yintercept = 0,
             lty = 2) +
  scale_fill_manual(name = 'Source',
                    values = c('Window' = 'lightblue',
                               'Trap' = 'pink')) +
  scale_color_manual(name = 'Source',
                     values = c('Window' = 'blue',
                                'Trap' = 'red')) +
  facet_wrap(~ Year, scales = 'free') +
  labs(x = 'Date',
       y = 'Model - Data',
       title = paste('Difference Between Model and Window Counts / Trap Estimates\nFor', spp))

#----------------------------------------------------------------------
# add IDFG estimates for comparison
#----------------------------------------------------------------------
idfg_lgr = na.omit(read.csv('Data/IDFG/sy09-15ThetasForKevin_20151117.csv')[-c(1:2),-1]) %>% tbl_df() %>%
  mutate(boot_iter = 1:nrow(.))
idfg_lgr_pt_est = read.csv('Data/IDFG/sy09-15ThetasForKevin_20151117.csv')[c(1:2),-1] %>% tbl_df()


idfg_pt_est = gather(idfg_lgr_pt_est, Name, estimate) %>% 
  bind_cols(data.frame(est_type = rep(c('PointEstimate', 'Median'), nrow(.)/2))) %>%
  mutate(Species = ifelse(grepl('sthd', Name), 'Sthd', 'Chnk'),
         Year = gsub('[[:alpha:]]', '', Name),
         Year = as.integer(Year),
         Type = ifelse(grepl('H$', Name), 'H', 
                       ifelse(grepl('HNC$', Name), 'HNC',
                              ifelse(grepl('W$', Name), 'W', NA))),
         Group = ifelse(Type == 'H', 'H', 'W')) %>%
  select(-Name) %>%
  arrange(est_type, Species, Year, Type)

idfg_summ = idfg_pt_est %>% 
  select(-Group) %>%
  mutate(Species = revalue(Species, c('Chnk' = 'Chinook',
                                      'Sthd' = 'Steelhead'))) %>%
  spread(Type, estimate) %>%
  mutate(Hatch_Morph = H,
         Hatch_PBT = H + HNC,
         Wild_Morph = W + HNC,
         Wild_PBT = W) %>%
  filter(est_type == 'PointEstimate') %>%
  select(-est_type, -(H:W))

# compare total window counts, total trap estimates, model total, and IDFG estimates
tot_summ %>%
  filter(Variable %in% c('Unique.Wild.Fish.Morph', 'Unique.Wild.Fish.PBT')) %>%
  select(Species, Year, Variable, estimate = median) %>%
  mutate(estimate = round(estimate),
         Variable = revalue(Variable, c('Unique.Wild.Fish.Morph' = 'ISEMP_Wild_Morph',
                                        'Unique.Wild.Fish.PBT' = 'ISEMP_Wild_PBT'))) %>%
  spread(Variable, estimate) %>%
  left_join(idfg_summ %>%
              select(-matches('Hatch')) %>%
              rename(IDFG_Wild_Morph = Wild_Morph,
                     IDFG_Wild_PBT = Wild_PBT)) %>%
  left_join(lgr_weekly %>%
              rename(Year = SpawnYear) %>%
              mutate(trap_rate = ifelse(trap_rate < 1e-12, NA, trap_rate)) %>%
              group_by(Species, Year) %>%
              summarise(win_Wild_Morph = round(sum(win_cnt * (Wild.morph / trap_fish), na.rm=T)),
                        win_Wild_PBT = round(sum(win_cnt * (Wild.PBT / trap_fish), na.rm=T)),
                        trap_Wild_Morph = round(sum(Wild.morph / trap_rate, na.rm=T)),
                        trap_Wild_PBT = round(sum(Wild.PBT / trap_rate, na.rm=T)))) %>%
  select(Species, Year, matches('Morph'), matches('PBT'))



# pull together IDFG bootstrap samples, reformat to match with ISEMP's posteriors
idfg_boot = gather(idfg_lgr, Name, value, -boot_iter) %>% 
  mutate(Species = ifelse(grepl('sthd', Name), 'Sthd', 'Chnk'),
         Species = revalue(Species, c('Sthd' = 'Steelhead', 'Chnk' = 'Chinook')),
         Year = gsub('[[:alpha:]]', '', Name),
         Year = as.integer(Year),
         Type = ifelse(grepl('H$', Name), 'H', 
                       ifelse(grepl('HNC$', Name), 'HNC',
                              ifelse(grepl('W$', Name), 'W', NA))),
         Group = ifelse(Type == 'H', 'H', 'W'),
         Model = paste(Species, Year, sep='.')) %>%
  select(-Name, -Model)

# Sample posteriors (and re-sample IDFG bootstrap samples)
# how many samples from posteriors?
n_samp = 2000

tot_comp = tot_post %>%
  group_by(Species, Year, Variable) %>%
  mutate(iter = 1:length(value)) %>%
  ungroup() %>%
  filter(Variable %in% c('Unique.Wild.Fish.Morph', 'Unique.Wild.Fish.PBT')) %>%
  mutate(Variable = revalue(Variable, c('Unique.Wild.Fish.Morph' = 'ISEMP_Morph', 'Unique.Wild.Fish.PBT' = 'ISEMP_PBT'))) %>%
  rename(Type = Variable) %>%
  spread(Type, value) %>%
  group_by(Species, Year) %>%
  sample_n(n_samp, replace = T) %>%
  mutate(iter = 1:n_samp) %>%
  ungroup() %>%
  full_join(idfg_boot %>%
              select(-Group) %>%
              spread(Type, value) %>%
              mutate(IDFG_Morph = W + HNC) %>%
              rename(IDFG_PBT = W) %>%
              select(-(H:HNC), -boot_iter) %>%
              group_by(Species, Year) %>%
              sample_n(n_samp, replace = T) %>%
              mutate(iter = 1:n_samp) %>%
              ungroup()) %>%
  select(-iter) %>%
  gather(Group, estimate, -(Species:Year)) %>%
  mutate(Type = ifelse(grepl('Morph$', Group), 'Morph', 'PBT'),
         Org = ifelse(grepl('^ISEMP', Group), 'ISEMP', 'IDFG'),
         yr_grp = paste(Year, Org, sep = '_')) %>%
  select(-Group)

pos_dodge = 0.9
idfg_isemp_p = tot_comp %>%
  filter(Type == 'PBT') %>%
  ggplot(aes(x = Year, y = estimate, fill = yr_grp)) +
  stat_summary(fun.y=HPDI_outlier, geom='point', size=1, color='gray', position=position_dodge(width=pos_dodge)) +
  stat_summary(fun.data=HPDI_boxplot, geom='errorbar', width=0.4, position=position_dodge(width=pos_dodge)) +
  stat_summary(fun.data=HPDI_boxplot, geom='boxplot', middle.def='Mode', position=position_dodge(width=pos_dodge)) +
  stat_summary(fun.y=median, geom='point', pch=1, position=position_dodge(width=pos_dodge)) +
  stat_summary(fun.y=mean, geom='point', pch=4, position=position_dodge(width=pos_dodge)) +
  scale_fill_brewer(palette='Paired', name='Estimate', labels=paste(rep(2010:2015, each=2), c('IDFG', 'ISEMP')), guide=guide_legend(nrow=2)) +
  geom_point(data = lgr_weekly %>%
               mutate(trap_rate = ifelse(trap_rate < 1e-12, NA, trap_rate)) %>%
               group_by(Species, Year = SpawnYear) %>%
               summarise(win_Wild_Morph = round(sum(win_cnt * (Wild.morph / trap_fish), na.rm=T)),
                         win_Wild_PBT = round(sum(win_cnt * (Wild.PBT / trap_fish), na.rm=T)),
                         trap_Wild_Morph = round(sum(Wild.morph / trap_rate, na.rm=T)),
                         trap_Wild_PBT = round(sum(Wild.PBT / trap_rate, na.rm=T))) %>%
               select(Species:Year, everything()) %>%
               gather(Source, estimate, -(Species:Year)) %>%
               mutate(Type = ifelse(grepl('PBT', Source), 'PBT', 'Morph')) %>%
               mutate(Source = gsub('_Wild_PBT', '', gsub('_Wild_Morph', '', Source)),
                      Source = revalue(Source, c('win' = 'Window', 'trap' = 'Trap'))) %>%
               filter(Type == 'PBT') %>%
               mutate(yr_grp = NA),
             aes(x = Year, 
                 y = estimate,
                 shape = Source),
             size = 5,
             color = 'darkred') +
  scale_shape_manual(values = c('Window' = 4,
                                'Trap' = 8)) +
  facet_wrap(~ Species, scales = 'free') +
  # facet_grid(Species ~ Type, scales = 'free') +
  theme_bw() +
  theme(legend.position='bottom', axis.text.x = element_text(size=7)) +
  labs(title='Total Escapement over Lower Granite Dam', y='Wild Fish', x='Run Year')
idfg_isemp_p

ggsave('Figures/IDFG_ISEMP_comparison.pdf',
       idfg_isemp_p,
       width = 9,
       height = 6,
       units = 'in')


# how do window counts compare for ISEMP and IDFG?
# we expect steelhead to differ somewhat by random chance (ISEMP uses unclipped, IDFG uses all)
lgr_weekly %>%
  group_by(Species, SpawnYear) %>%
  summarise(ISEMP_tot_win = sum(win_cnt)) %>%
  left_join(idfg_pt_est %>%
              mutate(Species = revalue(Species, c('Chnk' = 'Chinook',
                                                  'Sthd' = 'Steelhead'))) %>%
              filter(est_type == 'PointEstimate',
                     !(Species == 'Steelhead' & Type == 'H')) %>%
              group_by(Species, SpawnYear = Year) %>%
              summarise(IDFG_tot_win = sum(estimate))) %>%
  mutate(diff = ISEMP_tot_win - IDFG_tot_win)

tot_post %>%
  filter(Variable == 'All.Fish') %>%
  mutate(Year = as.factor(as.character(Year))) %>%
  ggplot(aes(x = Year, y = value, fill = Year)) +
  stat_summary(fun.y=HPDI_outlier, geom='point', size=1, color='gray', position=position_dodge(width=pos_dodge)) +
  stat_summary(fun.data=HPDI_boxplot, geom='errorbar', width=0.4, position=position_dodge(width=pos_dodge)) +
  stat_summary(fun.data=HPDI_boxplot, geom='boxplot', middle.def='Mode', position=position_dodge(width=pos_dodge)) +
  stat_summary(fun.y=median, geom='point', pch=1, position=position_dodge(width=pos_dodge)) +
  stat_summary(fun.y=mean, geom='point', pch=4, position=position_dodge(width=pos_dodge)) +
  scale_fill_brewer(palette = 'Set1') +
  geom_point(data = lgr_weekly %>%
               group_by(Species, Year = SpawnYear) %>%
               summarise(ISEMP_tot_win = sum(win_cnt),
                         ISEMP_tot_trap = sum(trap_est, na.rm = T)) %>%
               mutate(Year = as.factor(as.character(Year))) %>%
               gather(Source, value, -(Species:Year)),
             aes(shape = Source),
             size = 6) +
  scale_shape_manual(values = c('ISEMP_tot_win' = 19,
                                'ISEMP_tot_trap' = 17),
                     labels = c('ISEMP_tot_win' = 'Window',
                                'ISEMP_tot_trap' = 'Trap')) +
  facet_wrap(~ Species, scales = 'free')

# example plot
ex_plot_df = tot_post %>%
  group_by(Species, Year, Variable) %>%
  mutate(iter = 1:length(value)) %>%
  ungroup() %>%
  filter(Variable %in% c('All.Fish', 'Daytime.Fish', 'All.Wild.Fish.PBT', 'Unique.Wild.Fish.PBT'),
         value < 2.1e5)

ggplot(ex_plot_df, aes(x = value, group = Variable, color = Variable)) +
  geom_density(aes(fill=Variable), lwd=1.5, alpha=0.2) +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  theme_bw() +
  facet_grid(Year ~ Species, scales = 'free') +
  geom_vline(data = lgr_weekly %>%
               mutate(trap_est = ifelse(trap_est > 1e12, NA, trap_est)) %>%
               group_by(Species, Year = SpawnYear) %>%
               summarise(tot_win_cnt = sum(win_cnt, na.rm=T),
                         tot_trap_est = sum(trap_est, na.rm=T)) %>%
               gather(Method, estimate, -(Species:Year)) %>%
               mutate(Method = revalue(Method, c('tot_win_cnt' = 'Window Count',
                                                 'tot_trap_est' = 'Trap Estimate')),
                      Method = factor(Method, levels = c('Window Count', 'Trap Estimate'))),
             aes(xintercept = estimate,
                 linetype = Method),
             show_guide = T,
             color = 'gray') +
  scale_linetype_manual(values = c('Window Count' = 2,
                                   'Trap Estimate' = 3)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = 'Fish',
       y = 'Density',
       title = 'Posteriors')

# another example, showing estimates of PBT and morphological totals, as well as the IDFG point estimates
tot_plot_df = tot_post %>%
  group_by(Species, Year, Variable) %>%
  mutate(iter = 1:length(value)) %>%
  ungroup() %>%
  filter(Variable %in% c('All.Fish', 'Daytime.Fish', 'All.Wild.Fish.PBT', 'Unique.Wild.Fish.PBT', 'All.Wild.Fish.Morph', 'Unique.Wild.Fish.Morph'),
         value < 2.1e5) %>%
  mutate(Variable = factor(Variable, levels = c('All.Fish', 'Daytime.Fish', 'All.Wild.Fish.PBT', 'Unique.Wild.Fish.PBT', 'All.Wild.Fish.Morph', 'Unique.Wild.Fish.Morph')))

ggplot(tot_plot_df, aes(x = value, group = Variable, color = Variable)) +
  geom_density(aes(fill=Variable), lwd=1.5, alpha=0.2) +
  scale_color_brewer(palette = 'Paired') +
  scale_fill_brewer(palette = 'Paired') +
  theme_bw() +
  facet_grid(Year ~ Species, scales = 'free') +
  geom_vline(data = lgr_weekly %>%
               mutate(trap_est = ifelse(trap_est > 1e12, NA, trap_est)) %>%
               group_by(Species, Year = SpawnYear) %>%
               summarise(tot_win_cnt = sum(win_cnt, na.rm=T),
                         tot_trap_est = sum(trap_est, na.rm=T)) %>%
               gather(Method, estimate, -(Species:Year)) %>%
               mutate(Method = revalue(Method, c('tot_win_cnt' = 'Window Count',
                                                 'tot_trap_est' = 'Trap Estimate')),
                      Method = factor(Method, levels = c('Window Count', 'Trap Estimate'))),
             aes(xintercept = estimate,
                 linetype = Method),
             show_guide = T,
             color = 'gray') +
  geom_vline(data = idfg_summ %>%
               gather(Method, estimate, -(Species:Year)) %>%
               filter(grepl('Wild', Method)) %>%
               mutate(Method = paste0('IDFG_', Method)),
             aes(xintercept = estimate,
                 linetype = Method),
             show_guide = T) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = 'Fish',
       y = 'Density',
       title = 'Posteriors')

