# Author: Kevin See
# Purpose: pull together window counts and trap data from Lower Granite dam, queried by DART (http://www.cbr.washington.edu/dart)
# Created: 2/23/2016
# Last Modified: 6/17/2016
# Notes: 
# Window counts: http://www.cbr.washington.edu/dart/query/adult_daily
# Trap sample rates: http://www.cbr.washington.edu/dart/query/pitadult_valid - Sample Time/Rates
# Details about night passage and re-ascension fish: http://www.cbr.washington.edu/dart/query/pit_adult_window - Detection Details by TagID
# All fish in LGR trap: http://www.cbr.washington.edu/dart/query/pitadult_valid - Valid List - All in Trap Detail
# OR use the data directly from LGR trap database - provided by IDFG (downloaded 3/8/16)
# Rick Orme prepared mark-recapture data to estimate trap rate

#--------------------------------------------------------
library(plyr)
library(dplyr)
library(lubridate)
library(tidyr)
library(magrittr)
library(FSA)
library(Rcapture)
library(boot)
library(msm)
library(ggplot2)
library(readxl)

#--------------------------------------------------------
# pull together all window counts
wind_cnts = read.csv('Data/WindowCounts/adultdaily_1444247382_835.csv') %>% 
  bind_rows(read.csv('Data/WindowCounts/adultdaily_1444247422_27.csv')) %>% 
  bind_rows(read.csv('Data/WindowCounts/adultdaily_1444247479_476.csv')) %>% 
  bind_rows(read.csv('Data/WindowCounts/adultdaily_1444247495_810.csv')) %>% tbl_df()  %>%
  filter(Project == 'Lower Granite') %>%
  mutate(window_open = T,
         Date = ymd(Date)) %>%
  # drop Fall Chinook
  mutate(Chin = ifelse(month(Date) > 8 | (month(Date) == 8 & mday(Date) > 17), NA, Chin),
         JChin = ifelse(month(Date) > 8 | (month(Date) == 8 & mday(Date) > 17), NA, JChin)) %>%
  select(-Project, -Chinook.Run, -(Sock:TempC))

#--------------------------------------------------------
# for dividing by weeks
for(yr in unique(year(wind_cnts$Date))) {
  if(yr == 2009) strata_dates = data.frame(start_date = ymd(paste0(yr, '0701')) + weeks(0:51)) %>%
      mutate(end_date = start_date + days(7) - eseconds(1))
  if(yr > 2009) strata_dates = strata_dates %>%
      bind_rows(data.frame(start_date = ymd(paste0(yr, '0701')) + weeks(0:51)) %>%
                  mutate(end_date = start_date + days(7) - eseconds(1)))
}

strata_df = strata_dates %>%
  filter(!(month(end_date) == 6 & mday(end_date) > 23)) %>%
  bind_rows(strata_dates %>%
              filter(month(end_date) == 6,
                     mday(end_date) > 23) %>%
              mutate(end_date = ymd(paste0(year(start_date), '0701')) - eseconds(1))) %>%
  arrange(start_date)

week_strata = interval(strata_df$start_date, strata_df$end_date)

#--------------------------------------------------------
# pull together all trap rates from DART
# all Chinook trap rates are identical to steelhead, so only read in the steelhead ones
trap_rates_org = read.csv('Data/TrapSampleRates/pit_adult_valid_2009_3.csv') %>%
  bind_rows(read.csv('Data/TrapSampleRates/pit_adult_valid_2010_3.csv')) %>%
  bind_rows(read.csv('Data/TrapSampleRates/pit_adult_valid_2011_3.csv')) %>%
  bind_rows(read.csv('Data/TrapSampleRates/pit_adult_valid_2012_3.csv')) %>%
  bind_rows(read.csv('Data/TrapSampleRates/pit_adult_valid_2013_3.csv')) %>%
  bind_rows(read.csv('Data/TrapSampleRates/pit_adult_valid_2014_3.csv')) %>%
  bind_rows(read.csv('Data/TrapSampleRates/pit_adult_valid_2015_3.csv')) %>%
  rename(n.Samples = X.Samples,
         n.SbyC = X.SbyC,
         n.Close = X.Close) %>%
  mutate(Date = mdy(Date))

# deal with fact that some days have two rows in data (different rates during different times of day)
trap_rates = trap_rates_org %>%
  group_by(Date, Year, DOY) %>%
  summarise_each(funs(sum(., na.rm = T)), SampTime:SecondsInDay) %>%
  left_join(trap_rates_org %>%
              group_by(Date, Year, DOY) %>%
              summarise_each(funs(weighted.mean(., na.rm = T, w = TotalTime)), Rate:ActualRateInclusiveTime)) %>%
  mutate(trap_open = T,
         RateCalc = TotalTimeInclusive / SecondsInDay)


#--------------------------------------------------------
# night passage & reascension data
# read in data for individual tags
lgr_details = read.csv('Data/TagDetails/pitadultwindow_upper_tagid_GRA_2010_1_no_1_2010_365_2010.csv') %>%
  bind_rows(read.csv('Data/TagDetails/pitadultwindow_upper_tagid_GRA_2011_1_no_1_2011_365_2011.csv')) %>%
  bind_rows(read.csv('Data/TagDetails/pitadultwindow_upper_tagid_GRA_2012_1_no_1_2012_366_2012.csv')) %>%
  bind_rows(read.csv('Data/TagDetails/pitadultwindow_upper_tagid_GRA_2013_1_no_1_2013_365_2013.csv')) %>%
  bind_rows(read.csv('Data/TagDetails/pitadultwindow_upper_tagid_GRA_2014_1_no_1_2014_365_2014.csv')) %>%
  bind_rows(read.csv('Data/TagDetails/pitadultwindow_upper_tagid_GRA_2015_1_no_1_2015_281_2015.csv')) %>%
  bind_rows(read.csv('Data/TagDetails/pitadultwindow_upper_tagid_GRA_2009_3_no_1_2009_365_2009.csv')) %>% 
  bind_rows(read.csv('Data/TagDetails/pitadultwindow_upper_tagid_GRA_2010_3_no_1_2010_365_2010.csv')) %>%
  bind_rows(read.csv('Data/TagDetails/pitadultwindow_upper_tagid_GRA_2011_3_no_1_2011_365_2011.csv')) %>% 
  bind_rows(read.csv('Data/TagDetails/pitadultwindow_upper_tagid_GRA_2012_3_no_1_2012_366_2012.csv')) %>% 
  bind_rows(read.csv('Data/TagDetails/pitadultwindow_upper_tagid_GRA_2013_3_no_1_2013_365_2013.csv')) %>% 
  bind_rows(read.csv('Data/TagDetails/pitadultwindow_upper_tagid_GRA_2014_3_no_1_2014_365_2014.csv')) %>% 
  bind_rows(read.csv('Data/TagDetails/pitadultwindow_upper_tagid_GRA_2015_3_no_1_2015_281_2015.csv')) %>% 
  tbl_df() %>%
  filter(Ladder == 'GRA') %>%
  mutate(Date = ymd(Date),
         Detection.DateTime = ymd_hms(Detection.DateTime),
         Species = revalue(Species, c('1' = 'Chinook', '3' = 'Steelhead')),
         Period = factor(Period, levels = c('D', 'N')),
         ReAscent = ifelse(TagID.Ascent.Count > 1, T, F),
         SpawnYear = ifelse(Species == 'Chinook', year(Date),
                            ifelse(Date >= ymd(paste0(year(Date), '0701')), year(Date) + 1, year(Date)))) %>%
  filter(!(Species == 'Chinook' & (Date < ymd(paste0(year(Date), '0301')) | Date > ymd(paste0(year(Date), '0817')))))

# summarise rates by day
lgr_night_reasc_daily = select(lgr_details,
                              Species, SpawnYear, Date, TagID, Rear.Type, Period, ReAscent) %>%
  distinct() %>%
  arrange(Species, SpawnYear, TagID, Date) %>%
  group_by(Species, SpawnYear, Date) %>%
  summarise(tot_tags = n_distinct(TagID),
            day_tags = length(unique(TagID[Period=='D'])),
            reascent_tags = length(unique(TagID[ReAscent])),
            tot_tags_W = length(unique(TagID[Rear.Type == 'W'])),
            day_tags_W = length(unique(TagID[Period=='D' & Rear.Type == 'W'])),
            reascent_tags_W = length(unique(TagID[ReAscent & Rear.Type == 'W']))) %>%
  ungroup() %>%
  filter(SpawnYear >= 2010, SpawnYear <= max(SpawnYear[Species == 'Chinook'])) %>%
  mutate(week_num_org = NA)

# assign week numbers to each day
for(i in 1:length(week_strata)) {
  lgr_night_reasc_daily$week_num_org[with(lgr_night_reasc_daily, which(Date %within% week_strata[i]))] = i
}

# summarise weekly
lgr_night_reasc_weekly = lgr_night_reasc_daily %>%
  group_by(Species, SpawnYear, week_num_org) %>%
  summarise_each(funs(sum), tot_tags:reascent_tags_W) %>%
  ungroup()

# plot showing estimates of day-time rate for all tags vs. wild tags
ggplot(lgr_night_reasc_weekly, aes(x = day_tags / tot_tags, y = day_tags_W / tot_tags_W)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  facet_grid(Species ~ SpawnYear, scales = 'free') +
  labs(x = 'All Tags', y = 'Wild Tags', title = 'Day-time Rate')

# plot showing estimates of re-ascension rate for all tags vs. wild tags
ggplot(lgr_night_reasc_weekly, aes(x = reascent_tags / tot_tags, y = reascent_tags_W / tot_tags_W)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  facet_grid(Species ~ SpawnYear, scales = 'free') +
  labs(x = 'All Tags', y = 'Wild Tags', title = 'Re-ascension Rate')


#--------------------------------------------------------
# pull trap data directly from IDFG trap database
pbt.db.date = '20160308'

lgr_trap = read.csv(paste0('Data/IDFG/tblLGDMasterCombineExportJodyW_', pbt.db.date, '.csv')) %>% tbl_df() %>%
  rename(Tag.ID = LGDNumPIT) %>%
  mutate(Date = floor_date(ymd_hms(CollectionDate), unit =),
         LGDSpecies = ifelse(Tag.ID == '3D9.1C2D0A60A1', 1, LGDSpecies),
         LGDSpecies = ifelse(Tag.ID == '3D9.1C2D0BC097', 3, LGDSpecies),
         Species = ifelse(LGDSpecies == 1, 'Chinook', ifelse(LGDSpecies == 3, 'Steelhead', NA))) %>%
  # drop data from other species, and other runs of Chinook
  filter(Species %in% c('Chinook', 'Steelhead'),
         !(Species == 'Chinook' & Date > ymd(paste0(year(Date), '0817'))))

# summarise by spawnyear, species and day
lgr_trap_summ = lgr_trap %>%
  # filter out juveniles, keep only adults
  filter(LGDLifeStage == 'RF',
         # filter out sort by code fish
         PTAgisSxCGRAObs != 'Yes') %>%
  # use all Chinook in the trap, but only unclipped steelhead (to match with windown counts)
  filter(Species == 'Chinook' | (Species == 'Steelhead' & LGDMarkAD == 'AI')) %>%
  group_by(Species, SpawnYear, Date) %>%
  summarize(trap_fish = length(MasterID),
            Wild.morph = sum(LGDRear=='W'),
            Hatch.morph = sum(LGDRear=='H'),
            NA.morph = sum(is.na(LGDRear)),
            Wild.PBT = sum(grepl('W$', SRR)),
            HNC.PBT = sum(grepl('H$', SRR) & LGDMarkAD=='AI'),
            Hatch.PBT = sum(grepl('H$', SRR) & LGDMarkAD=='AD'),
            NA.PBT = sum(is.na(SRR)),
            n_invalid = sum(LGDValid != 1),
            tot_B_run_sthd = sum(LGDFLmm >= 780 & grepl('W$', SRR))) %>%
  ungroup() %>%
  mutate(tot_B_run_sthd = ifelse(Species == 'Steelhead', tot_B_run_sthd, NA)) %>%
  mutate(SpawnYear = gsub('^SY', '', SpawnYear),
         SpawnYear = as.integer(SpawnYear),
         SpawnYear = ifelse(!is.na(SpawnYear), SpawnYear, ifelse(Species == 'Chinook', year(Date), ifelse(Species == 'Steelhead' & Date >= ymd(paste0(year(Date), '0701')), year(Date) + 1, year(Date))))) %>%
  arrange(Species, Date)

#----------------------------------------------------------
# put window counts, trap counts and trap rates together in daily and weekly formats
#----------------------------------------------------------
# include Chinook jacks in window counts
lgr_daily = wind_cnts %>%
  mutate(Chin = ifelse(is.na(Chin), 0, Chin),
         JChin = ifelse(is.na(JChin), 0, JChin),
         AllChin = Chin + JChin) %>%
  gather(Species, win_cnt, -Date, -window_open) %>%
  filter(Species %in% c('AllChin', 'WStlhd')) %>%
  mutate(Species = revalue(Species, c('AllChin' = 'Chinook', 'WStlhd' = 'Steelhead')),
         win_cnt = ifelse(is.na(win_cnt), 0, win_cnt)) %>%
  full_join(lgr_trap_summ %>%
              select(-SpawnYear)) %>%
  mutate(SpawnYear = ifelse(Species == 'Chinook', year(Date),
                            ifelse(Date >= ymd(paste0(year(Date), '0701')), year(Date) + 1, year(Date)))) %>%
  filter(!(Species == 'Chinook' & (Date < ymd(paste0(year(Date), '0301')) | Date > ymd(paste0(year(Date), '0817')))),
         SpawnYear >= 2010,
         SpawnYear <= 2015) %>%
  # pull in trap rates
  left_join(trap_rates %>%
              select(Date, trap_open, matches('Rate'), one_of(c('TotalTimeInclusive', 'SecondsInDay')))) %>%
  # estimate total escapement for each day by trap fish / trap rate
  mutate(trap_est = ifelse(RateCalc > 0, 
                           ifelse(Species == 'Steelhead', Wild.morph / RateCalc, 
                                  ifelse(Species == 'Chinook', (Wild.morph + Hatch.morph) / RateCalc, NA)), NA),
         week_num_org = NA,
  # create filters for whether window and/or trap is open
         window_open = ifelse(is.na(window_open), F, window_open),
         trap_open = ifelse(is.na(trap_open) | RateCalc == 0, F, trap_open)) %>%
  select(Species, Date, Year, SpawnYear, matches('week_num'), everything())
# assign week number
for(i in 1:length(week_strata)) {
  lgr_daily$week_num_org[with(lgr_daily, which(Date %within% week_strata[i]))] = i
}
lgr_daily %<>%
  left_join(lgr_night_reasc_daily) %>%
  group_by(Species, SpawnYear) %>%
  summarise(min_week = min(week_num_org)) %>%
  left_join(lgr_daily) %>%
  group_by(Species, SpawnYear) %>%
  mutate(week_num = week_num_org - min_week + 1) %>%
  select(-min_week)

lgr_daily %>%
  group_by(Species, SpawnYear) %>%
  summarise(min_week = min(week_num),
            min_org_week = min(week_num_org),
            max_week = max(week_num),
            min_date = min(Date),
            max_date = max(Date))


# transform to weekly summaries
lgr_weekly = lgr_daily %>%
  group_by(Species, SpawnYear, week_num_org, week_num) %>%
  summarise_each(funs(sum(., na.rm=T)), win_cnt, trap_fish:tot_B_run_sthd) %>%
  left_join(lgr_daily %>%
              group_by(Species, SpawnYear, week_num_org, week_num) %>%
              summarise_each(funs(mean(., na.rm=T)), Rate:ActualRateInclusiveTime)) %>%
  left_join(lgr_daily %>%
              group_by(Species, SpawnYear, week_num_org, week_num) %>%
              summarise(Start_Date = min(Date),
                        RateCalc = sum(TotalTimeInclusive, na.rm=T) / (172800 * length(SecondsInDay)),
                        window_open = ifelse(sum(window_open) > 0, T, F),
                        trap_open = ifelse(sum(trap_open) > 0, T, F))) %>%
  ungroup() %>%
  left_join(lgr_night_reasc_weekly) %>%
  select(Species:week_num, Start_Date, everything()) %>%
    mutate(trap_est = ifelse(trap_open & RateCalc > 0, 
                             ifelse(Species == 'Steelhead', Wild.PBT / RateCalc, 
                                    ifelse(Species == 'Chinook', (Wild.PBT + Hatch.PBT) / RateCalc, NA)), NA))

# do daily and weekly data have same number of fish at window and trap in total?
identical(group_by(lgr_weekly, Species, SpawnYear) %>%
            summarise(tot_cnts = sum(win_cnt),
                      tot_trap = sum(trap_fish)),
          group_by(lgr_daily, Species, SpawnYear) %>%
            summarise(tot_cnts = sum(win_cnt, na.rm = T),
                      tot_trap = sum(trap_fish, na.rm=T)))


# compare trap estimates and window counts
lgr_weekly %>%
  ggplot(aes(x = win_cnt, y = trap_est, color = as.factor(SpawnYear), fill = as.factor(SpawnYear))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  geom_smooth(method = lm, formula = y ~ -1 + x, fullrange = T) +
  # facet_grid(SpawnYear ~ Species, scales = 'free') +
  # facet_grid(Species ~ SpawnYear, scales = 'free') +
  facet_wrap(~ Species + SpawnYear, scales = 'free', nrow = 2) +
  scale_color_brewer(palette = 'Dark2') +
  scale_fill_brewer(palette = 'Dark2') +
  labs(x = 'Window Count', 
       y = 'Trap Estimate', 
       color = 'Spawn Year', 
       fill = 'Spawn Year',
       title = 'All Chinook and Unclipped Steelhead') +
  theme_bw() +
  theme(legend.position = 'bottom')

#--------------------------------------------------------
# for mark-recapture estimate of trap rate, read in data prepped by Rick Orme
trap_mr_org = read_excel('/Users/kevin/Dropbox/ISEMP/Analysis_Projects/LowerGraniteDam_Escapement/TotalEscapement/Data/LGR PIT obs for trap rate MCR 3-17-16 .xlsx', 
                         1, 
                         skip = 5) %>%
  mutate(Species = revalue(Species, c('steelhead' = 'Steelhead')))

# pull out Chinook and steelhead, try to group by week_num_org to match rest of LGR data
trap_mr = trap_mr_org %>%
  rename(Date = `Trap date`,
         M = `TRAP adj dn time`,
         C = `LAD adj dn time`,
         SpawnYear = `Spawn yr`,
         Year = `Obs Year YYYY`,
         Tag_Code = `Tag Code`,
         tag_yr = `unique tag-year`,
         rear_type = `rear type`,
         SRR = `SRR Code`,
         StartDate = `sample week min date`,
         week_num_org2 = `both spec. samp Week`) %>%
  mutate(R = ifelse(M == 1 & C == 1, 1, 0)) %>%
  select(SpawnYear, Year, Tag_Code, rear_type, SRR, Date, M, C, R) %>%
  filter(!is.na(Date),
         !is.na(M)) %>%
  arrange(Date) %>%
  mutate(week_num_org = NA)

# assign week number
for(i in 1:length(week_strata)) {
  trap_mr$week_num_org[with(trap_mr, which(Date %within% week_strata[i]))] = i
}

# summarise by week
mr_n_fish = trap_mr %>%
  group_by(Year, week_num_org) %>%
  summarise_each(funs(sum), M:R)

# put together list of capture histories
caphist_list = dlply(trap_mr, .(week_num_org), capHistConvert, cols2use = c('M', 'C', 'R'), in.type = 'individual', out.type = 'frequency')

# estimate trap rate based on mark-recapture
trap_rate_mr = ldply(caphist_list, function(x) {
  mod = try(closedp.t(x, dfreq=T), silent = T)
  if(class(mod)[1] == 'try-error') return(data.frame(N = NA, N_se = NA, p = NA, p_se = NA))
  mod_N = mod$results['Mt', c('abundance', 'stderr')]
  p_mean = inv.logit(coef(mod$glm[['Mt']])[2:3])
  p_se = deltamethod(list(~ 1 / (1 + exp(-x1)), ~ 1 / (1 + exp(-x2))), mean = coef(mod$glm[['Mt']])[2:3], cov = vcov(mod$glm[['Mt']])[2:3,2:3])
  return(data.frame(N = mod_N[1], N_se = mod_N[2], p = p_mean[1], p_se = p_se[1]))
}, .id = 'week_num_org') %>% tbl_df()



left_join(mr_n_fish, trap_rate_mr) %>%
  # filter(week_num_org %in% 253:258) %>%
  select(-(N:N_se)) %>%
  mutate(p_cv = p_se / p) %>%
  filter(p_cv > 0.5)

#-------------------------------------------------
# include the mark recapture estimate of trap rate
lgr_weekly %<>%
  left_join(trap_rate_mr %>%
              mutate(p_cv = p_se / p) %>%
              filter(p_cv < 0.7) %>%
              select(week_num_org, 
                     Rate_MR = p, 
                     Rate_MR_se = p_se)) %>%
  mutate(trap_rate = ifelse(!is.na(Rate_MR) & Rate_MR > 0, Rate_MR, RateCalc),
         # trap_rate = RateCalc,
         trap_est = ifelse(!is.na(trap_rate),
                           ifelse(Species == 'Steelhead', Wild.morph / trap_rate, 
                                  ifelse(Species == 'Chinook', (Wild.morph + Hatch.morph) / trap_rate, NA)), NA),
         trap_rate_se = ifelse(is.na(Rate_MR_se), 0.001, Rate_MR_se),
         trap_est_se = sqrt(trap_rate_se^2 * (-trap_fish * trap_rate^(-2))^2),
         # set up parameters describing trap rate as a beta distribution
         trap_alpha = ((1 - trap_rate) / trap_rate_se^2 - 1 / trap_rate) * trap_rate^2,
         trap_beta = trap_alpha * (1 / trap_rate - 1),
         trap_alpha = ifelse(trap_rate_se == 0, NA, trap_alpha),
         trap_beta = ifelse(trap_rate_se == 0, NA, trap_beta),
         trap_alpha = ifelse(trap_open, trap_alpha, 1e-12),
         trap_beta = ifelse(trap_open, trap_beta, 1)) %>%
  # check to see if some trap estimates seem valid
  mutate(Prob_Less = pbinom(trap_fish, win_cnt, trap_rate, lower.tail=T),
         Prob_More = pbinom(trap_fish, win_cnt, trap_rate, lower.tail=F),
         lower_trap_lim = qbinom(0.005, win_cnt, trap_rate, lower.tail=T) / trap_rate,
         upper_trap_lim = qbinom(0.995, win_cnt, trap_rate, lower.tail=T) / trap_rate,
         trap_valid = ifelse(trap_open & lower_trap_lim <= win_cnt & upper_trap_lim >= win_cnt & trap_fish < win_cnt, T, F),
         trap_valid = ifelse(trap_fish < 10, F, trap_valid),
         trap_valid = ifelse(abs(trap_est - win_cnt) / win_cnt > 10, F, trap_valid))

# look at trap rate
spp = c('Chinook', 'Steelhead')[1]
lgr_weekly %>%
  filter(Species == spp) %>%
  select(SpawnYear, trap_open, trap_valid, Start_Date, week_num, Rate, ActualRateInclusiveTime, RateCalc, Rate_MR, Rate_MR_se) %>%
  gather(Method, Estimate, -(SpawnYear:week_num), -Rate_MR_se) %>%
  mutate(Method = revalue(Method, c('Rate' = 'Set Rate', 'ActualRateInclusiveTime' = 'Naive Calc. Rate', 'RateCalc' = 'Obs. Trap Rate', 'Rate_MR' = 'M-R Rate'))) %>%
  rename(Rate_se = Rate_MR_se,
         Trap_Open = trap_open) %>%
  mutate(Rate_se = ifelse(Method == 'M-R Rate', Rate_se, NA)) %>%
  arrange(SpawnYear, week_num, Method) %>%
  ggplot(aes(x = Start_Date, y = Estimate, color = Method)) +
  geom_line() +
  # geom_point(aes(shape = Trap_Open)) +
  geom_point(aes(shape = trap_valid)) +
  scale_shape_manual(values = c('TRUE' = 19, 'FALSE' = 1)) +
  geom_errorbar(aes(ymin = Estimate + qnorm(0.025) * Rate_se,
                    ymax = Estimate + qnorm(0.975) * Rate_se)) +
  facet_wrap(~ SpawnYear, scales = 'free_x') +
  theme_bw()


# compare window counts and trap rates now
spp = c('Chinook', 'Steelhead')[2]
lgr_weekly %>%
  filter(Species == spp) %>%
  # filter(trap_est < 1e13) %>%
  ggplot(aes(x = win_cnt, y = trap_est, color = as.factor(SpawnYear), fill = as.factor(SpawnYear))) +
  geom_point(aes(shape = trap_valid)) +
  scale_shape_manual(values = c('TRUE' = 19, 'FALSE' = 1)) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  geom_smooth(data = lgr_weekly %>%
                filter(Species == spp,
                       trap_valid),
              method = lm, formula = y ~ -1 + x, fullrange = T) +
  facet_wrap(~ SpawnYear, scales = 'free') +
  labs(x = 'Window', y = 'Trap', color = 'Spawn Year', fill = 'Spawn Year', shape = 'Trap Valid',
       title = spp) +
  theme_bw() +
  theme(title = element_text(size = 20))

filter(lgr_weekly,
       Species == spp) %>%
  select(SpawnYear, window_open, trap_open, trap_valid, week_num, win_cnt, trap_est) %>%
  gather(Method, Estimate, -(SpawnYear:week_num)) %>%
  mutate(Method = revalue(Method, c('win_cnt' = 'Window', 'trap_est' = 'Trap')),
         Open = ifelse(Method == 'Window', window_open, trap_open)) %>%
  select(SpawnYear, week_num, Method, Open, trap_valid, Estimate) %>%
  arrange(SpawnYear, week_num, Method) %>%
  ggplot(aes(x = week_num, y = Estimate, color = Method)) +
  geom_point(aes(shape = trap_valid)) +
  scale_shape_manual(values = c('TRUE' = 19, 'FALSE' = 1)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ SpawnYear, scales = 'free') +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'Week', title = spp, color = 'Source', shape = 'Trap Valid')

filter(lgr_weekly,
       Species == spp) %>%
  select(SpawnYear, trap_open, win_cnt, week_num, Rate, RateCalc, Rate_MR, Rate_MR_se) %>%
  arrange(SpawnYear, week_num) %>%
  filter(!is.na(Rate_MR)) %>%
  ggplot(aes(x = RateCalc, y = Rate_MR)) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_point(aes(shape = trap_open, size = log(win_cnt))) +
  scale_shape_manual(values = c('TRUE' = 19, 'FALSE' = 1)) +
  geom_errorbar(aes(ymin = Rate_MR + qnorm(0.025) * Rate_MR_se,
                    ymax = Rate_MR + qnorm(0.975) * Rate_MR_se)) +
  facet_wrap(~ SpawnYear, scales = 'free') +
  labs(x = 'Obs Rate',
       y = 'M-R Rate',
       shape = 'Trap Open',
       size = 'Log Window Count')


#-------------------------------------------------
# read in historical data about night passage
night_data = read_excel('Data/Night Count Percents at LGD.xlsx', skip=17)

night_passage = night_data[17:21, 7:9] %>%
  mutate_each(funs(as.integer)) %>%
  mutate(Month = month(4:8, label = T)) %>%
  gather(Year, day_fish, -Month) %>%
  full_join(night_data[17:21, 11:18] %>%
              mutate_each(funs(as.integer)) %>%
              mutate(Month = month(4:8, label = T)) %>%
              gather(Year, night_fish, -Month)) %>%
  arrange(Year, Month) %>%
  mutate(day_fish = ifelse(is.na(day_fish) & Year >= 2002, 0, day_fish),
         night_fish = ifelse(is.na(night_fish), 0, night_fish),
         tot_fish = day_fish + night_fish,
         Species = 'Chinook') %>%
  bind_rows(night_data[1:10, 2:9] %>%
              mutate_each(funs(as.integer)) %>%
              mutate(Month = month(3:12, label = T)) %>%
              gather(Year, day_fish, -Month) %>%
              left_join(night_data[1:10, 11:18] %>%
                          mutate_each(funs(as.integer)) %>%
                          mutate(Month = month(3:12, label = T)) %>%
                          gather(Year, night_fish, -Month)) %>%
              mutate(day_fish = ifelse(is.na(day_fish), 0, day_fish),
                     night_fish = ifelse(is.na(night_fish), 0, night_fish),
                     tot_fish = day_fish + night_fish,
                     Species = 'Steelhead')) %>%
  select(Species, Month, Year, everything()) %>%
  arrange(Species, Month, Year)

# night_passage %>%
#   filter(Species == 'Chinook') %>%
#   select(Month:day_fish) %>%
#   spread(Year, day_fish)

#------------------------------------------
# save compiled data
save(lgr_weekly, night_passage, week_strata, file = 'DataPrepped.rda')

#------------------------------------------
# compile everything into a list to pass to JAGS
# average number of fish per week from historical data set
hist_month_means = expand.grid(list(Month = month(1:12, label = T),
                                    Species = c('Chinook', 'Steelhead'))) %>%
  full_join(night_passage %>%
              group_by(Species, Month) %>%
              summarise(mean_tot = mean(tot_fish, na.rm = T) / (30/7)) %>%
              ungroup())



jags_data_list = dlply(lgr_weekly, .(Species, SpawnYear), function(x) {
  month_vec = month(x$Start_Date)
  spp = unique(x$Species)
  
  # weighting of historical day time passage, vs observed tags
  hist_tags_all = hist_month_means %>%
    filter(Species == spp) %>%
    select(mean_tot) %>% as.matrix() %>% as.vector()
  hist_tags_all = hist_tags_all[month_vec]
  
  tot.fish = expand.grid(list(Month = month(1:12, label = T),
                              Year = filter(night_passage,
                                            Species == spp,
                                            !is.na(tot_fish)) %>%
                                select(Year) %>%
                                distinct() %>% as.matrix() %>% as.vector())) %>%
    left_join(night_passage %>%
                filter(Species == spp)) %>%
    select(Month, Year, tot_fish) %>%
    spread(Year, tot_fish, fill = 0) %>%
    select(-Month) %>%
    as.matrix()
  
  day.fish = expand.grid(list(Month = month(1:12, label = T),
                              Year = filter(night_passage,
                                            Species == spp,
                                            !is.na(tot_fish)) %>%
                                select(Year) %>%
                                distinct() %>% as.matrix() %>% as.vector())) %>%
    left_join(night_passage %>%
                filter(Species == spp)) %>%
    select(Month, Year, day_fish) %>%
    spread(Year, day_fish, fill = 0) %>%
    select(-Month) %>%
    as.matrix()
  
  obs_tags = mutate(x, 
                    obs_tags = ifelse(is.na(tot_tags_W), 0, tot_tags_W)) %>%
    select(obs_tags) %>% as.matrix() %>% as.vector()
  hist_tags = hist_tags_all * median(obs_tags / hist_tags_all, na.rm=T)
  # calculate weighting of historic to observed tags
  rho = hist_tags / (hist_tags + obs_tags)
  rho[is.na(rho)] = 0

  
  data_list = list('TotLadderWeeks' = nrow(x),
                   'Y.window' = x %>% mutate(win_cnt = ifelse(!window_open, NA, win_cnt)) %>%
                     select(win_cnt) %>% as.matrix() %>% as.vector(),
                   'Y.trap' = x %>% mutate(y_trap = ifelse(Species == 'Chinook', trap_fish, Wild.morph),
                                           y_trap = ifelse(!trap_open | !trap_valid, NA, y_trap)) %>%
                     select(y_trap) %>% as.matrix() %>% as.vector(),
                   'trap.fish' = x %>% select(trap_fish) %>% as.matrix() %>% as.vector(),
                   'wild.morph' = x %>% select(Wild.morph) %>% as.matrix() %>% as.vector(),
                   'wild.pbt' = x %>% select(Wild.PBT) %>% as.matrix() %>% as.vector(),
                   'trap.rate' = x %>% mutate(trap_rate = ifelse(!trap_open | !trap_valid | is.na(trap_fish), 0, trap_rate)) %>%
                     select(trap_rate) %>% as.matrix() %>% as.vector(),
                   'trap.alpha' = x %>% select(trap_alpha) %>% as.matrix() %>% as.vector(),
                   'trap.beta' = x %>% select(trap_beta) %>% as.matrix() %>% as.vector(),
                   'Tot.tags' = x %>% mutate(Tot_tags = ifelse(is.na(tot_tags_W), 0, tot_tags_W)) %>%
                     select(Tot_tags) %>% as.matrix() %>% as.vector(),
                   'ReAsc.tags' = x %>% mutate(ReAsc_tags = ifelse(reascent_tags_W > tot_tags_W, tot_tags_W, reascent_tags_W),
                                               ReAsc_tags = ifelse(!trap_open, NA, ReAsc_tags)) %>%
                     select(ReAsc_tags) %>% as.matrix() %>% as.vector(),
                   'DC.tags' = x %>% mutate(Day_tags = ifelse(day_tags_W > tot_tags_W, tot_tags_W, day_tags_W),
                                            Day_tags = ifelse(!trap_open, NA, Day_tags)) %>%
                     select(Day_tags) %>% as.matrix() %>% as.vector(),
                   'B.trap' = x %>% select(tot_B_run_sthd) %>% as.matrix() %>% as.vector(),
                   'month.vec' = month_vec,
                   'n.hist.yrs' = ncol(tot.fish),
                   'tot.fish' = tot.fish,
                   'day.fish' = day.fish,
                   'rho' = rho
  )
  return(data_list)
})

# save data that's ready for JAGS
save(jags_data_list, hist_month_means, file = 'JAGS_DataList.rda')
