setwd("D:\\PosDoc\\Coronavirus\\USomicron\\USomicron")
library(dplyr)
library(xlsx)
library(readxl)
library(tidyverse)     ## data wrangling + ggplot2
library(colorspace)    ## adjust colors
library(rcartocolor)   ## Carto palettes
library(ggforce)       ## sina plots
library(ggdist)        ## halfeye plots
library(ggridges)      ## ridgeline plots
library(ggbeeswarm)    ## beeswarm plots
library(gghalves)      ## off-set jitter
library(systemfonts)   ## custom fonts
library(latex2exp)
library(zoo)
library(data.table)
# Parameters -------------------------------------------------------------------
set.seed(1432)
#Total cost of vaccine clinic setup 
cost_setup = 3022840

#Expenditure on advertisement and awareness campaigns
cost_advertisement = 242986305.11 #242,986,305.11
#Total cost of vaccine storage and transport 
cost_storage_and_transport = 7205179.89
#Cost of vaccine administration (all other costs) 
cost_administration = 1752152923.18 #1,752,152,923.18
#Other vaccination expenses (mobile and homebound vaccinations) 
#Expenses that will benefit all vaccinations. ie. CIR, DOITT developed applications
expenses_benefits = 31000000
# Total cost of vaccines
vaccines_cost = 282663374.8 #282,663,374.8


#Indirect costs
pcpi_nyc = 69288 ## Per-capita personal income NYC
# https://tradingeconomics.com/united-states/employment-rate#:~:text=Employment%20Rate%20in%20the%20United%20States%20averaged%2059.22%20percent%20from,percent%20in%20April%20of%202020.
perc_vac_emp = 0.665 ## proportion of vaccinated people who is employed (18-64 yo)


wdl_vac = 0.5 #work days lost due to visit for vaccination
pm_adverse_1 = 0.517 #proportion of adverse reaction First dose Moderna
pm_adverse_2 = 0.748 #proportion of adverse reaction First dose Moderna
pp_adverse_1 = 0.48 #proportion of adverse reaction First dose Moderna
pp_adverse_2 = 0.642 #proportion of adverse reaction First dose Moderna
pjj_adverse = 0.76 #proportion of adverse reaction First dose Moderna
wdl_adverse_1 = 1.66 #(sd = 1.48)working days lost due to adverse reactions first dose
wdl_adverse_2 = 1.39 #(sd = 0.82)working days lost due to adverse reactions

#Direct costs
cost_outpatient_appointment = 893.0 #outpatient appointment (symptomatic cases)
n_outpatient_visits = 0.5 #(total number of mild cases / 2) per mild case - ASSUMED
cost_transp_outpatients = 44.49 #for each visit
cost_hosp_nICU = 25188.0 #cost of hospital non-ICU admission
cost_hosp_ICU = 70098.0 #cost of ICU admission
n_ED_visits = 1 #for each severe non-hospitalized case- ASSUMED
cost_ED_care = 2200 #cost ED care
n_EMS_calls = 2.5 #per hospitalized case
cost_transp_EMS = 1005.0
r = 0.03 #discount rate
cost_lifelost = 441325 ##240676 #455484 #per year of life lost #average of statistical life in US is
# between US$ 9-10 mi with life expectancy of 79 years - REVISE
max_cost_life = 10300000
# Cost of Illness

symp_isolation = 10 #days out of work
# hospitalization - take it from JAMA paper and add another 4 days
duration_hosp_niCU = c(6,6,6,6,6,3) #data for each strain for non-ICU
duration_hosp_ICU = c(15,15,15,15,15,7)  #data for each strain for ICU
days_beforeafter = 3.5+4


basedate = as.Date("2020-09-01")
basedate_vac = as.Date("2020-12-12")
enddate = as.Date("2022-03-31")
population = 332968798


idx_1 = 1
idx_2 = 2
#  Reading file function -----------------------------------------------------------------

# Let's create a function to read the incidence file

  read_file_incidence <- function(index,type,strain = c(1,2,3,4,5,6),st2 = "usa",beta = "106",ag="all"){
    
    data.cases1 = read.table(paste0("../data/results_prob_0_",beta,"_herd_immu_10_",index,"_",st2,"/simlevel_",type,"_inc_",ag,".dat"),',',h = T) 
    data.cases1 = data.cases1[,-1]
    
    data.cases2 = read.table(paste0("../data/results_prob_0_",beta,"_herd_immu_10_",index,"_",st2,"/simlevel_",type,"2_inc_",ag,".dat"),',',h = T) 
    data.cases2 = data.cases2[,-1]
    
    data.cases3 = read.table(paste0("../data/results_prob_0_",beta,"_herd_immu_10_",index,"_",st2,"/simlevel_",type,"3_inc_",ag,".dat"),',',h = T) 
    data.cases3 = data.cases3[,-1]
    
    data.cases4 = read.table(paste0("../data/results_prob_0_",beta,"_herd_immu_10_",index,"_",st2,"/simlevel_",type,"4_inc_",ag,".dat"),',',h = T) 
    data.cases4 = data.cases4[,-1]
    
    data.cases5 = read.table(paste0("../data/results_prob_0_",beta,"_herd_immu_10_",index,"_",st2,"/simlevel_",type,"5_inc_",ag,".dat"),',',h = T) 
    data.cases5 = data.cases5[,-1]
    
    data.cases6 = read.table(paste0("../data/results_prob_0_",beta,"_herd_immu_10_",index,"_",st2,"/simlevel_",type,"6_inc_",ag,".dat"),',',h = T) 
    data.cases6 = data.cases6[,-1]
    
    l = list(data.cases1,data.cases2,data.cases3,data.cases4,data.cases5,data.cases6)
    
    return(l[strain])
  }

# And a function to bootstrap 
  
  fc <- function(d, i){
    return(mean(d[i],na.rm=T))
  }

  # discount formula for YLL
  formula = function(x) min((cost_lifelost/r) - (1/((1+r)^x))*(cost_lifelost/r),max_cost_life)
  #formula = function(x) (cost_lifelost/r) - (1/((1+r)^x))*(cost_lifelost/r)
  
# Vaccination costs (direct) ----------------------------------------------

# Vaccination costs (indirect) -------------------------------------------------------
# calculate the loss of workdays to go to get vaccine from 15-64 y.o.
## Let's read the data that was provided
data_vac = read.table("../data/results_prob_0_106_herd_immu_10_1_usa/vaccine_working.dat", h=F)
head(data_vac)
tail(data_vac)

n_days_vac = rowSums(data_vac)
n_days_vac = n_days_vac*wdl_vac*perc_vac_emp

# Now, for adverse reaction, we need to calculate the proportion
# of each vaccine that was administered

#Let's read the data

data = read.table("../data/results_prob_0_106_herd_immu_10_1_usa/vaccine_working.dat", h=F)

# we need to add Janssen to second dose, to be able to use it for fully vaccinated
# in the next section


# number of working days lost due to vaccination
n_days_ad_jensen = (data[,3]*wdl_adverse_1)*pjj_adverse*perc_vac_emp
n_days_ad_moderna = 
  (data[,5]+data[,8])*wdl_adverse_2*pm_adverse_2*perc_vac_emp+
  data[,2]*wdl_adverse_1*pm_adverse_1*perc_vac_emp
n_days_ad_pfizer = 
  (data[,4]+data[,7])*wdl_adverse_2*pp_adverse_2*perc_vac_emp+
  data[,1]*wdl_adverse_1*pp_adverse_1*perc_vac_emp

# total
n_days_work_lost = n_days_vac+n_days_ad_pfizer+n_days_ad_jensen+n_days_ad_moderna

indirect_vaccination_cost = n_days_work_lost*pcpi_nyc/365*population/100000

# COVID-19 costs ----------------------------------------------------------

# We need to read the real data, re-scale the hospitalization, and rework the
# other outcomes

data.cases = read.csv("../data/data_us_april.csv")

LHS_tab = read.csv("../data/LHSmatrix.csv")

# Illness and Hospitalization (direct) ---------------------------------------------------------

# we want to see the hospitalization scaling factor

#total hospitalization per 100,000 population from the beginning of vaccination
total_hosp = data.cases %>% mutate(rhos = frollmean(hos,7)) %>% filter(date >= basedate_vac,date <= enddate) %>% pull(rhos) %>% sum()/population*100000
total_hosp

#Let's see this number in the simulation

hos_sim = read_file_incidence(idx_1,"hos")
icu_sim = read_file_incidence(idx_1,"icu")

hos = Reduce('+', hos_sim) # adding the strains
icu = Reduce('+', icu_sim)# adding the strains
nn = nrow(hos)
  
v_date = basedate+seq(0,nn-1) #creating a vector with the dates of simulation


sum.sim.hos = sum(frollmean(rowMeans(hos),7)[v_date >= basedate_vac & v_date <= enddate])
sum.sim.icu = sum(frollmean(rowMeans(icu),7)[v_date >= basedate_vac & v_date <= enddate])

factor_hos = total_hosp/(sum.sim.hos+sum.sim.icu)
factor_hos

asymp = Reduce('+',read_file_incidence(idx_1,"asymp")) # adding the strains
inf = Reduce('+',read_file_incidence(idx_1,"inf")) # adding the strains
mild = Reduce('+',read_file_incidence(idx_1,"mild")) # adding the strains

# Let's set inf to be severe non-hospitalized
inf_nh = inf - hos - icu
# number of extra hospitalization after scaling
n_extra = (sum.sim.hos+sum.sim.icu)*factor_hos - (sum.sim.hos+sum.sim.icu)

#total number
sum.sim.asymp = sum(asymp[v_date >= basedate_vac & v_date <= enddate,])/ncol(asymp)
sum.sim.mild = sum(mild[v_date >= basedate_vac & v_date <= enddate,])/ncol(mild)
sum.sim.sev = sum(inf_nh[v_date >= basedate_vac & v_date <= enddate,])/ncol(inf_nh)

factor_non_hos = (sum(c(sum.sim.asymp,sum.sim.mild,sum.sim.sev))-n_extra)/sum(c(sum.sim.asymp,sum.sim.mild,sum.sim.sev))

# Now we want to bootstrap the mean of those matrices

sum.sim.mild = colSums(mild)
#sum.sim.mild = boot::boot(sum.sim,fc,500)$t[,1]

sum.sim.inf = colSums(inf_nh)
#sum.sim.inf = boot::boot(sum.sim,fc,500)$t[,1]

sum.sim.hos = colSums(hos)
#sum.sim.hos = boot::boot(sum.sim,fc,500)$t[,1]

sum.sim.icu = colSums(icu)
#sum.sim.icu = boot::boot(sum.sim,fc,500)$t[,1]

# we increased the hospitalizations by some amount, therefore,
# we decrease the other infections proportionately

sum.sim.mild = sum.sim.mild*factor_non_hos
sum.sim.inf = sum.sim.inf*factor_non_hos

sum.sim.hos = sum.sim.hos*factor_hos
sum.sim.icu = sum.sim.icu*factor_hos


#cost mild infection of hospital
set.seed(1432)
cost_symp = (boot::boot(sum.sim.mild,fc,1000)$t[,1])*
  n_outpatient_visits*LHS_tab$Outpatient
#cost severe non-hospitalized infection
set.seed(1432)
cost_inf = (boot::boot(sum.sim.inf,fc,1000)$t[,1])*LHS_tab$ED.care
#cost for hospitalizations
set.seed(1432)
cost_hos = (boot::boot(sum.sim.hos,fc,1000)$t[,1])*(LHS_tab$Hospitalization.without.ICU+n_EMS_calls*LHS_tab$EMS.transportation)
set.seed(1432)
cost_icu = (boot::boot(sum.sim.icu,fc,1000)$t[,1])*(LHS_tab$Hospitalization.with.ICU+n_EMS_calls*LHS_tab$EMS.transportation)

cost_hospital = (cost_symp+cost_inf+cost_hos+cost_icu)*population/100000
#cost_hospital

###
# For the SCENARIO without vaccination
###

hos_sim = read_file_incidence(idx_2,"hos")
icu_sim = read_file_incidence(idx_2,"icu")

hos = Reduce('+', hos_sim) # adding the strains
icu = Reduce('+', icu_sim)# adding the strains

nn = nrow(hos)
v_date = basedate+seq(0,nn-1) #creating a vector with the dates of simulation

sum.sim.hos2 = sum(rowMeans(hos)[v_date >= basedate_vac & v_date <= enddate])
sum.sim.icu2 = sum(rowMeans(icu)[v_date >= basedate_vac & v_date <= enddate])


asymp = Reduce('+',read_file_incidence(idx_2,"asymp")) # adding the strains
inf = Reduce('+',read_file_incidence(idx_2,"inf")) # adding the strains
mild = Reduce('+',read_file_incidence(idx_2,"mild")) # adding the strains

# Let's set inf to be severe non-hospitalized
inf_nh = inf - hos - icu
# number of extra hospitalization after scaling
n_extra = (sum.sim.hos2+sum.sim.icu2)*factor_hos - (sum.sim.hos2+sum.sim.icu2)

#total number
sum.sim.asymp2 = sum(asymp[v_date >= basedate_vac & v_date <= enddate,])/ncol(asymp)
sum.sim.mild2 = sum(mild[v_date >= basedate_vac & v_date <= enddate,])/ncol(mild)
sum.sim.sev2 = sum(inf_nh[v_date >= basedate_vac & v_date <= enddate,])/ncol(inf_nh)

factor_non_hos_2 = (sum(c(sum.sim.asymp2,sum.sim.mild2,sum.sim.sev2))-n_extra)/sum(c(sum.sim.asymp2,sum.sim.mild2,sum.sim.sev2))


# Bootstraping
sum.sim.mild2 = colSums(mild)
#sum.sim.mild2 = boot::boot(sum.sim,fc,500)$t[,1]

sum.sim.inf2 = colSums(inf_nh)
#sum.sim.inf2 = boot::boot(sum.sim,fc,500)$t[,1]

sum.sim.hos2 = colSums(hos)
#sum.sim.hos2 = boot::boot(sum.sim,fc,500)$t[,1]

sum.sim.icu2 = colSums(icu)
#sum.sim.icu2 = boot::boot(sum.sim,fc,500)$t[,1]

# we increased the hospitalizations by some amount, therefore,
# we decrease the other infections proportionately

sum.sim.mild2 = sum.sim.mild2*factor_non_hos_2
sum.sim.inf2 = sum.sim.inf2*factor_non_hos_2

sum.sim.hos2 = sum.sim.hos2*factor_hos
sum.sim.icu2 = sum.sim.icu2*factor_hos

#cost mild infection of hospital
set.seed(1432)
cost_symp = (boot::boot(sum.sim.mild2,fc,1000)$t[,1])*
  n_outpatient_visits*LHS_tab$Outpatient
#cost severe non-hospitalized infection
set.seed(1432)
cost_inf = (boot::boot(sum.sim.inf2,fc,1000)$t[,1])*LHS_tab$ED.care
#cost for hospitalizations
set.seed(1432)
cost_hos = (boot::boot(sum.sim.hos2,fc,1000)$t[,1])*(LHS_tab$Hospitalization.without.ICU+n_EMS_calls*LHS_tab$EMS.transportation)
set.seed(1432)
cost_icu = (boot::boot(sum.sim.icu2,fc,1000)$t[,1])*(LHS_tab$Hospitalization.with.ICU+n_EMS_calls*LHS_tab$EMS.transportation)

cost_hospital2 = (cost_symp+cost_inf+cost_hos+cost_icu)*population/100000
#cost_hospital2



# Illness and Hospitalization (indirect) ---------------------------------------------------------
# 
# LOS: 
#   https://www.medrxiv.org/content/10.1101/2022.01.11.22269045v1
# 
# https://stacks.cdc.gov/view/cdc/113758
# https://stacks.cdc.gov/view/cdc/114452


# We need to look at specific age groups. Therefore, it is better if we create a function
# for computing the total number of days lost in hospital and illness

days_of_leave <- function(agegroup = "all",idx = 1){
  hos_sim = read_file_incidence(idx,"hos",c(1,2,3,4,5,6),"usa","106",agegroup)
  icu_sim = read_file_incidence(idx,"icu",c(1,2,3,4,5,6),"usa","106",agegroup)
  
  inf_sim = read_file_incidence(idx,"inf",c(1,2,3,4,5,6),"usa","106",agegroup)
  asymp_sim = read_file_incidence(idx,"asymp",c(1,2,3,4,5,6),"usa","106",agegroup)
  mild_sim = read_file_incidence(idx,"mild",c(1,2,3,4,5,6),"usa","106",agegroup)
  
  #for symptoms
  inf_nh = Reduce('+',inf_sim) - Reduce('+',hos_sim) -Reduce('+',icu_sim) 
  
  fh = c(factor_non_hos,factor_non_hos_2)
  #we are considering Asymp, because the code is calibrated for reported cases
  symp = colSums(inf_nh+Reduce('+',mild_sim)+Reduce('+',asymp_sim))*fh[idx]
  
  v_symp = symp#boot::boot(symp,fc,500)$t[,1]
  
  h1 = lapply(hos_sim, colSums)
  #h1b = lapply(h1, boot::boot,statistic=fc,R=500)
  #h11b = lapply(h1b, function(x) x[["t"]])
  #h11b = do.call(cbind,h11b)
  h11b = do.call(cbind,h1)
  
  hc1 = lapply(icu_sim, colSums)
  #hc1b = lapply(hc1, boot::boot,statistic=fc,R=500)
  #hc11b = lapply(hc1b, function(x) x[["t"]])
  #hc11b = do.call(cbind,hc11b)
  hc11b = do.call(cbind,hc1)
  
  h11b = factor_hos*h11b
  hc11b = factor_hos*hc11b
  
  dh1 = h11b %*% (duration_hosp_niCU+days_beforeafter)
  dc1 = hc11b %*% (duration_hosp_ICU+days_beforeafter)
  dsymp = v_symp*symp_isolation
  
  total_days = dh1+dc1+dsymp
  return(total_days)
}

ag = c("ag3","ag4","ag5","ag6")
total_days_1 = lapply(ag, days_of_leave, idx = 1)
total_days_2 = lapply(ag, days_of_leave, idx = 2)

total1_hos = Reduce("+",total_days_1)[,1]
total2_hos = Reduce("+",total_days_2)[,1]


#
cost_indirect_ill1 = total1_hos*perc_vac_emp*pcpi_nyc/365*population/100000
cost_indirect_ill2 = total2_hos*perc_vac_emp*pcpi_nyc/365*population/100000

# Results -----------------------------------------------------------------

pal_plots = c("#885687","#9ebcda","#756bb1","#bcbddc")


# Direct -----------------------------------------------------------------
set.seed(1432)
ind_vac = boot::boot(indirect_vaccination_cost,fc,1000)$t[,1]

#Initial Value Investment

bb_vac = cost_hospital#+cost_indirect_ill1
total_cost_vaccination = bb_vac#+indirect_vaccination_cost
total_cost_no_vac = cost_hospital2#+cost_indirect_ill2

###

# Final value of investment
FVI = total_cost_no_vac - total_cost_vaccination
bb = FVI

# Let's read the file containing the amount of people that died at age x in each sim
# plot --------------------------------------------------------------------


## theme for horizontal charts
theme_flip <-
  theme(
    axis.text.x = element_text(face = "plain", family = "Helvetica", size = 22),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(), 
    axis.title.x = element_text(face = "bold", family = "Helvetica", size = 22),
    panel.grid.major.x = element_blank(),#element_line(color = "grey90", size = .6),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect("white"),
    axis.ticks.x = element_line(color = "black", size = 1.0),
    axis.ticks.length.x = unit(0.1,"cm"),
    legend.position = "none", 
    legend.text = element_text(family = "Helvetica", size = 18),
    legend.title = element_text(face = "bold", size = 18),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

dd = data.frame(ROI=bb)
df_point = data.frame(x = c(as.factor(1),as.factor(1)),y = quantile(dd$ROI,c(0.025,0.975),name=F))

mean(dd$ROI)
df_point

ggplot(dd, aes(x = forcats::fct_rev(as.factor(1)), y = ROI,color = as.factor(1),fill = as.factor(1)))+
  geom_boxplot(
    width = .08, fill = "white",
    size = 1.1, outlier.shape = NA) +
  ggdist::stat_halfeye(
    adjust = .33,
    width = .3, 
    color = NA,
    position = position_nudge(x = .07)
  ) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .2, 
    alpha = .5, size = 1.5,position = position_nudge(x = .08)
  ) +
  stat_summary(fun = mean,geom="point", color= "#6a51a3", fill = "#6a51a3",shape = 21, size = 2.5)+
  geom_point(data=df_point, aes(x = x,y=y),color= "#6a51a3", fill = "#6a51a3",shape = 23, size = 2.5)+
  scale_y_continuous(expand=expansion(mult=c(0,0)))+
  scale_x_discrete(expand=expansion(mult=c(0,0)))+
  labs(y="Savings",y=NULL)+
  coord_flip() +
  scale_color_manual(values = pal_plots[3], guide = "none") +
  scale_fill_manual(values = pal_plots[3], guide = "none") +
  theme_flip


ggsave(
  "../data/savings_direct.pdf",
  device = "pdf",
  width = 7.5,
  height = 4,
  dpi = 300
)




# Indirect -----------------------------------------------------------------

#Initial Value Investment
set.seed(1432)
bb_vac = boot::boot(cost_indirect_ill1,fc,1000)$t[,1]
total_cost_vaccination = bb_vac+indirect_vaccination_cost
set.seed(1432)
total_cost_no_vac = boot::boot(cost_indirect_ill2,fc,1000)$t[,1]

###

# Final value of investment
FVI = total_cost_no_vac - total_cost_vaccination
bb = FVI



# Let's read the file containing the amount of people that died at age x in each sim
# plot --------------------------------------------------------------------


## theme for horizontal charts
theme_flip <-
  theme(
    axis.text.x = element_text(face = "plain", family = "Helvetica", size = 22),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(), 
    axis.title.x = element_text(face = "bold", family = "Helvetica", size = 22),
    panel.grid.major.x = element_blank(),#element_line(color = "grey90", size = .6),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect("white"),
    axis.ticks.x = element_line(color = "black", size = 1.0),
    axis.ticks.length.x = unit(0.1,"cm"),
    legend.position = "none", 
    legend.text = element_text(family = "Helvetica", size = 18),
    legend.title = element_text(face = "bold", size = 18),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

dd = data.frame(ROI=bb)
df_point = data.frame(x = c(as.factor(1),as.factor(1)),y = quantile(dd$ROI,c(0.025,0.975),name=F))

mean(dd$ROI)
df_point

ggplot(dd, aes(x = forcats::fct_rev(as.factor(1)), y = ROI,color = as.factor(1),fill = as.factor(1)))+
  geom_boxplot(
    width = .08, fill = "white",
    size = 1.1, outlier.shape = NA) +
  ggdist::stat_halfeye(
    adjust = .33,
    width = .3, 
    color = NA,
    position = position_nudge(x = .07)
  ) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .2, 
    alpha = .5, size = 1.5,position = position_nudge(x = .08)
  ) +
  stat_summary(fun = mean,geom="point", color= "#6a51a3", fill = "#6a51a3",shape = 21, size = 2.5)+
  geom_point(data=df_point, aes(x = x,y=y),color= "#6a51a3", fill = "#6a51a3",shape = 23, size = 2.5)+
  scale_y_continuous(expand=expansion(mult=c(0,0)))+
  scale_x_discrete(expand=expansion(mult=c(0,0)))+
  labs(y="Savings",y=NULL)+
  coord_flip() +
  scale_color_manual(values = pal_plots[1], guide = "none") +
  scale_fill_manual(values = pal_plots[1], guide = "none") +
  theme_flip


ggsave(
  "../data/savings_indirect.pdf",
  device = "pdf",
  width = 7.5,
  height = 4,
  dpi = 300
)
