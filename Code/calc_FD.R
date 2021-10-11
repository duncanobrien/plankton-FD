### Calculate Planktonic FD ###

#Preamble
require(readxl) #read.xlsx function
require(tidyverse) #dplyr, ggplot etc.
require(janitor) #row_to_names function
require()
source("Code/tidy_FD_fn.R")


##########################################################################################
## Read in Plankton Abundance Data ##
##########################################################################################
source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Kinneret/Data/kinneret_plankton_data.R")
source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Kinneret/Data/kinneret_environmental_data.R")
#phyo, zoo + combined raw concentrations at week, month and year groupings

source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Zurich/Data/zurich_plankton_data.R")
source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Zurich/Data/LZ_environmental_data.R")

source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Madison/Data/madison_plankton_data.R")
source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Madison/Data/madison_environmental_data.R")

source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Windermere/Data/windermere_plankton_data.R")
source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Windermere/Data/windermere_environmental_data.R")

##########################################################################################
## Read in Plankton Trait Data ##
##########################################################################################
phyto.kin.traits.dat <- read_xlsx("Data/fuzzy_phytoplankton_traits.xlsx",sheet = 3) %>%
  slice(-c(1)) %>% #drop redundant first col
  janitor::row_to_names(row_number = 1) %>%
  select(-c(Species,Notes)) %>%
  slice(-3) %>% #drop 2-Aphanizomenon oval as no data
  mutate(across(c(lgth_1:col_T,n_fix:sil_T),as.numeric)) %>%
  mutate(across(c(mob,troph),as.factor))

zoo.kin.traits.dat <- read_xlsx("Data/fuzzy_zooplankton_traits.xlsx",sheet = 1) %>%
  slice(-c(1)) %>%
  janitor::row_to_names(row_number = 1) %>%
  select(-c(Notes)) %>%
  drop_na()%>% #drop no data species
  mutate(across(c(lgth_1:omniherb),as.numeric)) %>% #ensure numeric and not character
  mutate(across(c(rv_mech,f_mode),as.factor))

phyto.LZ.traits.dat <- read_xlsx("Data/fuzzy_phytoplankton_traits.xlsx",sheet = 6) %>%
  slice(-c(1)) %>%
  janitor::row_to_names(row_number = 1) %>%
  select(-c(Notes)) %>%
  #slice(-c(4,19,22,23,24,45,52,53,60,63,70,71)) %>% #drop as no data
  drop_na()%>% #drop no data species
  mutate(across(c(lgth_1:col_T,n_fix:sil_T),as.numeric)) %>% #ensure numeric and not character
  mutate(across(c(mob,troph),as.factor))

zoo.LZ.traits.dat <- read_xlsx("Data/fuzzy_zooplankton_traits.xlsx",sheet = 4) %>%
  slice(-c(1)) %>%
  janitor::row_to_names(row_number = 1) %>%
  dplyr::select(-c(Notes)) %>%
  drop_na()%>% #drop no data species
  mutate(across(c(lgth_1:omniherb),as.numeric)) %>% #ensure numeric and not character
  mutate(across(c(rv_mech,f_mode),as.factor))

phyto.mad.traits.dat <- read_xlsx("Data/fuzzy_phytoplankton_traits.xlsx",sheet = 5) %>%
  slice(-c(1)) %>%
  janitor::row_to_names(row_number = 1) %>%
  select(-c(Notes)) %>%
  #slice(-c(1,7,8,17,19,23,32)) %>% #drop as no data
  drop_na()%>% #drop no data species
  mutate(across(c(lgth_1:col_T,n_fix:sil_T),as.numeric)) %>% #ensure numeric and not character
  mutate(across(c(mob,troph),as.factor))

zoo.mad.traits.dat <- read_xlsx("Data/fuzzy_zooplankton_traits.xlsx",sheet = 3) %>%
  slice(-c(1)) %>%
  janitor::row_to_names(row_number = 1) %>%
  select(-c(Notes)) %>%
  slice(-c(8,9,11)) %>% #drop as no data
  drop_na()%>% #drop no data species
  mutate(across(c(lgth_1:omniherb),as.numeric)) %>% #ensure numeric and not character
  mutate(across(c(rv_mech,f_mode),as.factor))

phyto.wind.traits.dat <- read_xlsx("Data/fuzzy_phytoplankton_traits.xlsx",sheet = 4) %>%
  slice(-c(1)) %>%
  janitor::row_to_names(row_number = 1) %>%
  select(-c(Notes)) %>%
  mutate(across(c(lgth_1:col_T,n_fix:sil_T),as.numeric)) %>% #ensure numeric and not character
  mutate(across(c(mob,troph),as.factor))

##########################################################################################
## Estimate Phytoplankton FD ##
##########################################################################################
FD_metrics <- c("MNND","Rao","FRic","FDis","FDiv","FEve","SRic","SDiv","RARDis","RARScr")

# Kinneret #
phyto.kin.fuzFDs.mth <-tidyFD(plank.kin.combo.mth[,4:48], phyto.kin.traits.dat, trophic.lvl = "phyto",
                              traittype = "fuzzy", method = FD_metrics, correction="cailliez")
phyto.kin.fuzFDs.mth <- cbind(date = as.numeric(plank.kin.combo.mth$Date),phyto.kin.fuzFDs.mth)
write.csv(phyto.kin.fuzFDs.mth,"Data/raw_FD/FD_kin_phyto_mth_raw.csv",row.names = FALSE)
phyto.kin.fuzFDs.mth <- read.csv("Data/raw_FD/FD_kin_phyto_mth_raw.csv")

phyto.kin.fuzFDs.yr <-tidyFD(plank.kin.combo.yr[,2:46], phyto.kin.traits.dat, trophic.lvl = "phyto",
                             traittype = "fuzzy", method = FD_metrics, correction="cailliez")
phyto.kin.fuzFDs.yr <- cbind(date = as.numeric(plank.kin.combo.yr$Date),phyto.kin.fuzFDs.yr)
write.csv(phyto.kin.fuzFDs.yr,"Data/raw_FD/FD_kin_phyto_yr_raw.csv",row.names = FALSE)
phyto.kin.fuzFDs.yr <- read.csv("Data/raw_FD/FD_kin_phyto_yr_raw.csv")

# Lower Zurich #
phyto.LZ.fuzFDs.mth <-tidyFD(plank.LZ.combo.mth[,4:178], phyto.LZ.traits.dat, trophic.lvl = "phyto",
                             traittype = "fuzzy", method = FD_metrics, correction="cailliez", ndim =9)
phyto.LZ.fuzFDs.mth <- cbind(date = as.numeric(plank.LZ.combo.mth$date),phyto.LZ.fuzFDs.mth)
write.csv(phyto.LZ.fuzFDs.mth,"Data/raw_FD/FD_LZ_phyto_mth_raw.csv",row.names = FALSE)
phyto.LZ.fuzFDs.mth <- read.csv("Data/raw_FD/FD_LZ_phyto_mth_raw.csv")

phyto.LZ.fuzFDs.yr <-tidyFD(plank.LZ.combo.yr[,2:176], phyto.LZ.traits.dat, trophic.lvl = "phyto",
                            traittype = "fuzzy", method = FD_metrics, correction="cailliez")
phyto.LZ.fuzFDs.yr <- cbind(date = as.numeric(plank.LZ.combo.yr$date),phyto.LZ.fuzFDs.yr)
write.csv(phyto.LZ.fuzFDs.yr,"Data/raw_FD/FD_LZ_phyto_yr_raw.csv",row.names = FALSE)
phyto.LZ.fuzFDs.yr <- read.csv("Data/raw_FD/FD_LZ_phyto_yr_raw.csv")

# Mendota #
phyto.mad.fuzFDs.mth <-tidyFD(plank.mad.combo.mth[,4:200], phyto.mad.traits.dat, trophic.lvl = "phyto",
                              traittype = "fuzzy", method = FD_metrics, correction="cailliez")
phyto.mad.fuzFDs.mth <- cbind(date = as.numeric(plank.mad.combo.mth$date),phyto.mad.fuzFDs.mth)
write.csv(phyto.mad.fuzFDs.mth,"Data/raw_FD/FD_mad_phyto_mth_raw.csv",row.names = FALSE)
phyto.mad.fuzFDs.mth <- read.csv("Data/raw_FD/FD_mad_phyto_mth_raw.csv")

phyto.mad.fuzFDs.yr <-tidyFD(plank.mad.combo.yr[,2:198], phyto.mad.traits.dat, trophic.lvl = "phyto",
                             traittype = "fuzzy", method = FD_metrics, correction="cailliez")
phyto.mad.fuzFDs.yr <- cbind(date = as.numeric(plank.mad.combo.yr$date),phyto.mad.fuzFDs.yr)
write.csv(phyto.mad.fuzFDs.yr,"Data/raw_FD/FD_mad_phyto_yr_raw.csv",row.names = FALSE)
phyto.mad.fuzFDs.yr <- read.csv("Data/raw_FD/FD_mad_phyto_yr_raw.csv")

# Windermere #
phyto.wind.fuzFDs.mth <-tidyFD(phyto_env.windmthdata[,4:20], phyto.wind.traits.dat, trophic.lvl = "phyto",
                               traittype = "fuzzy", method = FD_metrics, correction="cailliez")
phyto.wind.fuzFDs.mth <- cbind(date = as.numeric(phyto_env.windmthdata$Date),phyto.wind.fuzFDs.mth)
write.csv(phyto.wind.fuzFDs.mth,"Data/raw_FD/FD_wind_phyto_mth_raw.csv",row.names = FALSE)
phyto.wind.fuzFDs.mth <- read.csv("Data/raw_FD/FD_wind_phyto_mth_raw.csv")

phyto.wind.fuzFDs.yr <-tidyFD(phyto_env.windyrdata[,2:18], phyto.wind.traits.dat, trophic.lvl = "phyto",
                              traittype = "fuzzy", method = FD_metrics, correction="cailliez")
phyto.wind.fuzFDs.yr <- cbind(date = as.numeric(phyto_env.windyrdata$Date),phyto.wind.fuzFDs.yr)
write.csv(phyto.wind.fuzFDs.yr,"Data/raw_FD/FD_wind_phyto_yr_raw.csv",row.names = FALSE)
phyto.wind.fuzFDs.yr <- read.csv("Data/raw_FD/FD_wind_phyto_yr_raw.csv")

##########################################################################################
## Estimate Zooplankton FD ##
##########################################################################################
FD_metrics <- c("MNND","Rao","FRic","FDis","FDiv","FEve","SRic","SDiv","RARDis","RARScr")

# Kinneret #
zoo.kin.fuzFDs.mth <-tidyFD(plank.kin.combo.mth[,49:79], zoo.kin.traits.dat, trophic.lvl = "zoo",
                            traittype = "fuzzy", method = FD_metrics, correction="cailliez")
zoo.kin.fuzFDs.mth <- cbind(date = as.numeric(plank.kin.combo.mth$Date),zoo.kin.fuzFDs.mth)
write.csv(zoo.kin.fuzFDs.mth,file = "Data/raw_FD/FD_kin_zoo_mth_raw.csv",row.names = F)
zoo.kin.fuzFDs.mth <- read.csv("Data/raw_FD/FD_kin_zoo_mth_raw.csv")

zoo.kin.fuzFDs.yr <-tidyFD(plank.kin.combo.yr[,47:77], zoo.kin.traits.dat, trophic.lvl = "zoo",
                           traittype = "fuzzy", method = FD_metrics, correction="cailliez")
zoo.kin.fuzFDs.yr <- cbind(date = as.numeric(plank.kin.combo.yr$Date),zoo.kin.fuzFDs.yr)
write.csv(zoo.kin.fuzFDs.yr,file = "Data/raw_FD/FD_kin_zoo_yr_raw.csv",row.names = F)
zoo.kin.fuzFDs.yr <- read.csv("Data/raw_FD/FD_kin_zoo_yr_raw.csv")

# Lower Zurich #
zoo.LZ.fuzFDs.mth <-tidyFD(plank.LZ.combo.mth[,179:205], zoo.LZ.traits.dat, trophic.lvl = "zoo",
                           traittype = "fuzzy", method = FD_metrics, correction="cailliez")
zoo.LZ.fuzFDs.mth <- cbind(date = as.numeric(plank.LZ.combo.mth$date),zoo.LZ.fuzFDs.mth)
write.csv(zoo.LZ.fuzFDs.mth,file = "Data/raw_FD/FD_LZ_zoo_mth_raw.csv",row.names = F)
zoo.LZ.fuzFDs.mth <- read.csv("Data/raw_FD/FD_LZ_zoo_mth_raw.csv")

zoo.LZ.fuzFDs.yr <-tidyFD(plank.LZ.combo.yr[,177:203], zoo.LZ.traits.dat, trophic.lvl = "zoo",
                          traittype = "fuzzy", method = FD_metrics, correction="cailliez")
zoo.LZ.fuzFDs.yr <- cbind(date = as.numeric(plank.LZ.combo.yr$date),zoo.LZ.fuzFDs.yr)
write.csv(zoo.LZ.fuzFDs.yr,file = "Data/raw_FD/FD_LZ_zoo_yr_raw.csv",row.names = F)
zoo.LZ.fuzFDs.yr <- read.csv("Data/raw_FD/FD_LZ_zoo_yr_raw.csv")

# Mendota #
zoo.mad.fuzFDs.mth <-tidyFD(plank.mad.combo.mth[,201:224], zoo.mad.traits.dat, trophic.lvl = "zoo",
                            traittype = "fuzzy", method = FD_metrics, correction="cailliez")
zoo.mad.fuzFDs.mth <- cbind(date = as.numeric(plank.mad.combo.mth$date),zoo.mad.fuzFDs.mth)
write.csv(zoo.mad.fuzFDs.mth,file = "Data/raw_FD/FD_mad_zoo_mth_raw.csv",row.names = F)
zoo.mad.fuzFDs.mth <- read.csv("Data/raw_FD/FD_mad_zoo_mth_raw.csv")

zoo.mad.fuzFDs.yr <-tidyFD(plank.mad.combo.yr[,199:222], zoo.mad.traits.dat, trophic.lvl = "zoo",
                           traittype = "fuzzy", method = FD_metrics, correction="cailliez")
zoo.mad.fuzFDs.yr <- cbind(date = as.numeric(plank.mad.combo.yr$date),zoo.mad.fuzFDs.yr)
write.csv(zoo.mad.fuzFDs.yr,file = "Data/raw_FD/FD_mad_zoo_yr_raw.csv",row.names = F)
zoo.mad.fuzFDs.yr <- read.csv("Data/raw_FD/FD_mad_zoo_yr_raw.csv")
