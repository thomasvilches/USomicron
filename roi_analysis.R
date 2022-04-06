setwd("~/PosDoc/Coronavirus/USomicron/Code/")
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
pcpi_nyc = 74472 ## Per-capita personal income NYC
perc_vac_emp = 0.7176 ## proportion of vaccinated people that is employed (18-64 yo)
wdl_vac = 0.5 #work days lost due to visit for vaccination
pm_adverse_1 = 0.517 #proportion of adverse reaction First dose Moderna
pm_adverse_2 = 0.748 #proportion of adverse reaction First dose Moderna
pp_adverse_1 = 0.48 #proportion of adverse reaction First dose Moderna
pp_adverse_2 = 0.642 #proportion of adverse reaction First dose Moderna
pjj_adverse = 0.76 #proportion of adverse reaction First dose Moderna
wdl_adverse_1 = 1.66 #(sd = 1.48)working days lost due to adverse reactions first dose
wdl_adverse_2 = 1.39 #(sd = 0.82)working days lost due to adverse reactions

#Direct costs
cost_outpatient_appointment = 893.0/10 #outpatient appointment (symptomatic cases)
n_outpatient_visits = 0.5 #(total number of mild cases / 2) per mild case - ASSUMED
cost_transp_outpatients = 44.49 #for each visit
cost_hosp_nICU = 25188.0 #cost of hospital non-ICU admission
cost_hosp_ICU = 70098.0 #cost of ICU admission
n_ED_visits = 1 #for each severe non-hospitalized case- ASSUMED
cost_ED_care = 530.0 #cost ED care
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

# This is using the data about expenses provided by NYC
direct_vaccination_cost = cost_setup+cost_administration+cost_advertisement+
  cost_storage_and_transport+expenses_benefits+vaccines_cost

# Vaccination costs (indirect) -------------------------------------------------------
# calculate the loss of workdays to go to get vaccine from 15-64 y.o.
## Let's read the data that was provided
data_vac = read_excel("data/NYC_Daily_COVID-19_Vax_by_UHF_AgeGroup_2022-02-08_1500.xlsx")
head(data_vac)
tail(data_vac)

#data_vac %>% mutate(DATE = as.Date(DATE)) %>% ggplot()+geom_line(aes(x = DATE,y = N_FULLY_VACCINATED,color = AGE_GROUP))

#clean the number 999
data_vac = data_vac[data_vac$UHF != 999,]
### let's take the sum of facilities
data = data_vac %>% group_by(AGE_GROUP) %>% summarise(partially = sum(N_PARTIALLY_VACCINATED), fully = sum(N_FULLY_VACCINATED))
data
# The age groups are not exactly what we need, but let's use it anyway. 
# This is conservative (overestimates the vaccination cost)

# We want to calculate the number of working days that were spent due to vaccination

working_group = c("15to24","25to44","45to64")
n_vacs_first = sum(data$partially[data$AGE_GROUP %in% working_group])
n_vacs_second = sum(data$fully[data$AGE_GROUP %in% working_group])
# We need to do the same for booster dose

data_booster = read.csv("data/Demo_additional_dose_age_2022-02-08_1500.csv",sep = ";")
head(data_booster)
tail(data_booster)

data_booster$DATE = as.Date(data_booster$DATE)

data_booster_total = data_booster %>% filter(DATE >= as.Date("2020-12-14"), RESIDENCY == "NYC") %>%
  group_by(AGE_GROUP) %>% summarise(total_b = sum(N_ADDITIONAL_VACCINATED))

data_booster_total

working_group = c("18to24","25to34","35to44","45to54","55to64")
n_vacs_booster = sum(data_booster_total$total_b[data_booster_total$AGE_GROUP %in% working_group])
n_days_vac = (n_vacs_first+n_vacs_second+n_vacs_booster)*wdl_vac*perc_vac_emp

# Now, for adverse reaction, we need to calculate the proportion
# of each vaccine that was administered

#Let's read the data

data = read.csv("data/Dose_admin_bymonth_2022-02-08_1500.csv",sep = ";")

# we need to add Janssen to second dose, to be able to use it for fully vaccinated
# in the next section
df = data %>% group_by(VAC_CODE) %>% 
  summarise(total = sum(ALL_DOSES),
            first = sum(ADMIN_DOSE1),
            second  = sum(ADMIN_DOSE2)+sum(ADMIN_SINGLE), booster = sum(ADMIN_ADDITIONAL))

# The total number is in df
vaccines = c("Janssen","Moderna","Pfizer")
pp = df %>% filter(VAC_CODE %in% vaccines)

first = pp$first/sum(pp$first)
second = pp$second/sum(pp$second)
boost = pp$booster/sum(pp$booster)

# number of working days lost due to vaccination
n_days_ad_jensen = (n_vacs_second*second[1]*wdl_adverse_1+n_vacs_booster*boost[1]*wdl_adverse_2)*pjj_adverse*perc_vac_emp
n_days_ad_moderna = 
  (n_vacs_second*second[2]+n_vacs_booster*boost[2])*wdl_adverse_2*pm_adverse_2*perc_vac_emp+
  n_vacs_first*first[2]*wdl_adverse_1*pm_adverse_1*perc_vac_emp
n_days_ad_pfizer = 
  (n_vacs_second*second[3]+n_vacs_booster*boost[3])*wdl_adverse_2*pp_adverse_2*perc_vac_emp+
  n_vacs_first*first[3]*wdl_adverse_1*pp_adverse_1*perc_vac_emp

# total
n_days_work_lost = n_days_vac+n_days_ad_pfizer+n_days_ad_jensen+n_days_ad_moderna

indirect_vaccination_cost = n_days_work_lost*pcpi_nyc/365

# COVID-19 costs ----------------------------------------------------------

# We need to read the real data, re-scale the hospitalization, and rework the
# other outcomes

data.cases = read.csv("../data/data_us_april.csv")


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
cost_symp = (sum.sim.mild)*
  n_outpatient_visits*(cost_outpatient_appointment)
#cost severe non-hospitalized infection
cost_inf = (sum.sim.inf)*cost_ED_care
#cost for hospitalizations
cost_hos = (sum.sim.hos)*(cost_hosp_nICU+n_EMS_calls*cost_transp_EMS)
cost_icu = (sum.sim.icu)*(cost_hosp_ICU+n_EMS_calls*cost_transp_EMS)

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
cost_symp = (sum.sim.mild2)*
     n_outpatient_visits*(cost_outpatient_appointment)
#cost severe non-hospitalized infection
cost_inf = (sum.sim.inf2)*cost_ED_care
#cost for hospitalizations
cost_hos = (sum.sim.hos2)*(cost_hosp_nICU+n_EMS_calls*cost_transp_EMS)
cost_icu = (sum.sim.icu2)*(cost_hosp_ICU+n_EMS_calls*cost_transp_EMS)

cost_hospital2 = (cost_symp+cost_inf+cost_hos+cost_icu)*population/100000
#cost_hospital2


bb = boot::boot(cost_hospital,fc,R=500)
mean(bb$t[,1])
boot::boot.ci(bb,0.95)$bca

bb = boot::boot(cost_hospital2,fc,R=500)
mean(bb$t[,1])
boot::boot.ci(bb,0.95)$bca
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


# Years of Life Lost ------------------------------------------------------
#fist of all, let's find the scaling factor for deaths.

deaths = read_file_incidence(idx_1,"ded")
ded = Reduce('+', deaths)
nn = nrow(ded)
v_date = basedate+seq(0,nn-1) #creating a vector with the dates of simulation
sum.sim.ded = sum(ded[v_date >= basedate_vac & v_date <= enddate,])/ncol(ded)
total_deaths = data.cases %>% filter(date_of_interest >= basedate_vac,date_of_interest <= enddate) %>% pull(inc_deaths) %>% sum()/population*100000
total_deaths
factor_deaths = total_deaths/sum.sim.ded

# Let's read the file containing the amount of people that died at age x in each sim

age_of_death= read.table("data/results_prob_0_121_1_newyorkcity/year_of_death.dat",h=F)
dim(age_of_death)

life_exp = read.csv("data/life_exp.csv",sep = ";",h=F)$V2[1:nrow(age_of_death)]


vector_cost = unlist(lapply(life_exp,formula))
#plot(vector_cost)

vyll = as.vector(vector_cost %*% as.matrix(age_of_death))

vyll1 = vyll*factor_deaths*population/100000#boot::boot(vyll*factor_deaths*population/100000,fc,500)

age_of_death= read.table("data/results_prob_0_121_2_newyorkcity/year_of_death.dat",h=F)
dim(age_of_death)

vyll = as.vector(vector_cost %*% as.matrix(age_of_death))

vyll2 = vyll*factor_deaths*population/100000#boot::boot(vyll*factor_deaths*population/100000,fc,500)


cost_yll1 = vyll1#(vyll1$t[,1])
cost_yll2 = vyll2#(vyll2$t[,1])

#mean(vyll1$t)
#boot::boot.ci(vyll1)

#mean(vyll2$t)
#boot::boot.ci(vyll2)


# YLL ---------------------------------------------------------------------


age_of_death= read.table("data/results_prob_0_121_1_newyorkcity/year_of_death.dat",h=F)
dim(age_of_death)

#plot(vector_cost)

vyll = as.vector(life_exp %*% as.matrix(age_of_death))

vyll1 = boot::boot(vyll*factor_deaths*population/100000,fc,500)
mean(vyll1$t)
boot::boot.ci(vyll1)

n_people1 = colSums(age_of_death)
average_yll1 = vyll/n_people1

quantile(average_yll1,c(0.025,0.975,0.5),na.rm = T)

age_of_death= read.table("data/results_prob_0_121_2_newyorkcity/year_of_death.dat",h=F)
dim(age_of_death)
n_people2 = colSums(age_of_death)
#plot(vector_cost)

vyll = as.vector(life_exp %*% as.matrix(age_of_death))

vyll2 = boot::boot(vyll*factor_deaths*population/100000,fc,500)
mean(vyll2$t)
boot::boot.ci(vyll2)


df = data.frame(value=c(vyll1$t[,1],vyll2$t[,1]),scen = c(rep("vac",length(vyll1$t[,1])),rep("novac",length(vyll2$t[,1]))))

kruskal.test(value~scen, data = df)

np1 = boot::boot(n_people1*factor_deaths*population/100000,fc,500)$t[,1]
np2 = boot::boot(n_people2*factor_deaths*population/100000,fc,500)$t[,1]


cc = (cost_yll2-cost_yll1)/(np2-np1)
mean(cc)
quantile(cc,c(0.025,0.975,0.5),na.rm=T)
# Results -----------------------------------------------------------------

#Initial Value Investment
IVI = direct_vaccination_cost
bb_vac = boot::boot(cost_hospital+cost_indirect_ill1+cost_yll1,fc,R=500)$t[,1]
total_cost_vaccination = indirect_vaccination_cost+bb_vac
total_cost_no_vac = boot::boot(cost_hospital2+cost_indirect_ill2+cost_yll2,fc,R=500)$t[,1]

# Final value of investment
FVI = total_cost_no_vac - total_cost_vaccination
ROI = (FVI - IVI)/IVI

#Initial Value Investment
IVI_society = direct_vaccination_cost+indirect_vaccination_cost
total_cost_vaccination_society = bb_vac
total_cost_no_vac_society = total_cost_no_vac

# Final value of investment
FVI_society = total_cost_no_vac_society - total_cost_vaccination_society
ROI_society = (FVI_society - IVI_society)/IVI_society

df = data.frame(cost = c(total_cost_vaccination,total_cost_no_vac),
                scen= c(rep("Vaccination",length(total_cost_vaccination)),rep("No Vaccination",length(total_cost_no_vac))))

df %>% group_by(scen) %>% summarise(m = mean(cost),ci1 = quantile(cost,0.025,name=F),ci2 = quantile(cost,0.975,name=F)) %>%
  mutate(m=format(m,big.mark=","),ci1=format(ci1,big.mark=","),ci2=format(ci2,big.mark=","))



pal_plots = c("#885687","#9ebcda","#756bb1","#bcbddc")


# FVI-IVI -----------------------------------------------------------------


#Initial Value Investment
IVI = direct_vaccination_cost
bb_vac = cost_hospital+cost_indirect_ill1+cost_yll1
total_cost_vaccination = indirect_vaccination_cost+bb_vac
total_cost_no_vac = cost_hospital2+cost_indirect_ill2+cost_yll2

###
xx = total_cost_no_vac-bb_vac

# Final value of investment
FVI = total_cost_no_vac - total_cost_vaccination
xx = FVI - IVI
bb = boot::boot(xx,fc,R=500)
mean(bb$t[,1])
bbb=boot::boot.ci(bb,0.95)
bbb$bca
# Better plot -------------------------------------------------------------

factor_cost = total_cost_no_vac/total_cost_vaccination

total_cost_no_vac2 = total_cost_no_vac
total_cost_vaccination2 = total_cost_vaccination*4 #this is for plot

df_plot = data.frame(cost = c(total_cost_vaccination2,total_cost_no_vac2),
                scen= c(rep("Vaccination",length(total_cost_vaccination2)),rep("No Vaccination",length(total_cost_no_vac2))))
# 

theme_set(theme_void(base_family = "Helvetica"))

theme_update(
  axis.text.x = element_text(color = "black", face = "bold", size = 26, 
                             margin = margin(t = 6)),
  axis.text.y = element_text(color = "black", size = 22, hjust = 1, 
                             margin = margin(r = 6), family = "Helvetica"),
  axis.line.x = element_line(color = "black", size = 1),
  panel.grid.major.y = element_line(color = "grey90", size = .6),
  plot.background = element_rect(fill = "white", color = "white"),
  plot.margin = margin(rep(20, 4))
)


## custom colors
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]

## theme for horizontal charts
theme_flip <-
  theme(
    axis.text.x = element_text(face = "plain", family = "Helvetica", size = 22),
    axis.text.y = element_text(face = "bold", family = "Helvetica", size = 30),
    axis.title.x = element_text(face = "bold", family = "Helvetica", size = 22),
    panel.grid.major.x = element_blank(),#element_line(color = "grey90", size = .6),
    panel.grid.major.y = element_blank(),
    axis.ticks.x = element_line(color = "black", size = 1.0),
    axis.ticks.length.x = unit(0.1,"cm"),
    legend.position = "none", 
    legend.text = element_text(family = "Helvetica", size = 18),
    legend.title = element_text(face = "bold", size = 18),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

  label_t = c(TeX(r"($C_v$)"),TeX(r"($C_b$)"))
  
  ggplot(df_plot, aes(x = forcats::fct_rev(scen), y = cost/4/(1e9), 
                   color = scen, fill = scen)) +
    geom_boxplot(
      width = .1, fill = "white",
      size = 1.5, outlier.shape = NA) +
    ggdist::stat_halfeye(
      adjust = .33,
      width = .5, 
      color = NA,
      position = position_nudge(x = .09)
    ) +
    gghalves::geom_half_point(
      side = "l", 
      range_scale = .3, 
      alpha = .5, size = 1.5,position = position_nudge(x = .03)
    ) +
    scale_y_continuous(sec.axis = sec_axis(~.*4,name = "Cost (billion US$)"))+
    scale_x_discrete(breaks = c("Vaccination","No Vaccination"),labels = label_t, expand=expansion(mult=c(0,0)))+
    labs(y="Cost (billion US$)",y=NULL)+
    coord_flip() +
    scale_color_manual(values = my_pal[c(1,4)], guide = "none") +
    scale_fill_manual(values = my_pal[c(1,4)], guide = "none") +
    theme_flip
  
  
  
  
  ggsave(
    "../figures/cost.png",
    device = "png",
    width = 7,
    height = 6,
    dpi = 300,
  )
  
  ggsave(
    "../figures/cost.pdf",
    device = "pdf",
    width = 7,
    height = 6,
    dpi = 300,
  )
  
  
  
  
  
  ## general theme
  theme_set(theme_void(base_family = "Helvetica"))
  
  theme_update(
    axis.text.x = element_text(color = "black", face = "bold", size = 26, 
                               margin = margin(t = 6)),
    #axis.text.y = element_text(color = "black", size = 22, hjust = 1, 
    #margin = margin(r = 6), family = "Arial"),
    axis.line.x = element_line(color = "black", size = 1),
    panel.grid.major.y = element_line(color = "grey90", size = .6),
    plot.background = element_rect(fill = "white", color = "white"),
    plot.margin = margin(rep(20, 4))
  )
  ## theme for horizontal charts
  theme_flip <-
    theme(
      axis.text.x = element_text(face = "plain", family = "Helvetica", size = 22),
      #axis.text.y = element_text(face = "bold", family = "Arial", size = 26),
      axis.title.x = element_text(face = "bold", family = "Helvetica", size = 22),
      panel.grid.major.x = element_blank(),#element_line(color = "grey90", size = .6),
      panel.grid.major.y = element_blank(),
      axis.ticks.x = element_line(color = "black", size = 1.0),
      axis.ticks.length.x = unit(0.1,"cm"),
      legend.position = "none", 
      legend.text = element_text(family = "Helvetica", size = 18),
      legend.title = element_text(face = "bold", size = 18),
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    )

dd = as.data.frame(ROI)
df_point = data.frame(x = c(as.factor(1),as.factor(1)),y = quantile(ROI,c(0.025,0.975),name=F))

ggplot(dd, aes(x = forcats::fct_rev(as.factor(1)), y = ROI,color = as.factor(1),fill = as.factor(1))) +
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
  stat_summary(fun = mean,geom="point", color= "#88419d", fill = "#88419d",shape = 21, size = 2.5)+
  geom_point(data=df_point, aes(x = x,y=y),color= "#88419d", fill = "#88419d",shape = 23, size = 2.5)+
  scale_y_continuous(limits= c(54.2,63.8),breaks = seq(55,63,2),expand=expansion(mult=c(0,0)))+
  scale_x_discrete(expand=expansion(mult=c(0,0)))+
  labs(y="Return of investment",y=NULL)+
  coord_flip() +
  scale_color_manual(values = pal_plots[1], guide = "none") +
  scale_fill_manual(values = pal_plots[1], guide = "none") +
  theme_flip



ggsave(
  "../figures/ROI.pdf",
  device = "pdf",
  width = 7,
  height = 4,
  dpi = 300
)


dd %>% summarise(m = mean(ROI),ci1 = quantile(ROI,0.025,name=F),ci2 = quantile(ROI,0.975,name=F))


table = data.frame(Vaccination=c(mean(direct_vaccination_cost),mean(indirect_vaccination_cost), mean(cost_hospital), mean(cost_indirect_ill1),mean(cost_yll1)),
                   NoVaccination=c(NA,NA, mean(cost_hospital2), mean(cost_indirect_ill2),mean(cost_yll2)))
rownames(table) = c("Direct Vaccination Cost","Indirect Vaccination Cost","Direct cost of Illness","Illness indirect cost","YLL cost")
table %>% mutate(Vaccination = format(Vaccination, big.mark=","),
                 NoVaccination = format(NoVaccination, big.mark=","))


colSums(table,na.rm=T) %>% format( big.mark=",")


# Data for plot -----------------------------------------------------------

dd = data.frame(ROI=ROI_society)
df_point = data.frame(x = c(as.factor(1),as.factor(1)),y = quantile(ROI_society,c(0.025,0.975),name=F))


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
  scale_y_continuous(limits= c(26.5,32.5),breaks = seq(27,32,1),expand=expansion(mult=c(0,0)))+
  scale_x_discrete(expand=expansion(mult=c(0,0)))+
  labs(y="Return of investment",y=NULL)+
  coord_flip() +
  scale_color_manual(values = pal_plots[3], guide = "none") +
  scale_fill_manual(values = pal_plots[3], guide = "none") +
  theme_flip


ggsave(
  "../figures/ROI_soc.pdf",
  device = "pdf",
  width = 7,
  height = 4,
  dpi = 300
)


dd %>% summarise(m = mean(ROI),ci1 = quantile(ROI,0.025,name=F),ci2 = quantile(ROI,0.975,name=F))


table = data.frame(Vaccination=c(mean(direct_vaccination_cost),mean(indirect_vaccination_cost), mean(cost_hospital), mean(cost_indirect_ill1),mean(cost_yll1)),
                   NoVaccination=c(NA,NA, mean(cost_hospital2), mean(cost_indirect_ill2),mean(cost_yll2)))
rownames(table) = c("Direct Vaccination Cost","Indirect Vaccination Cost","Direct cost of Illness","Illness indirect cost","YLL cost")
table %>% mutate(Vaccination = format(Vaccination, big.mark=","),
                 NoVaccination = format(NoVaccination, big.mark=","))


colSums(table,na.rm=T) %>% format( big.mark=",")


# YLL figure --------------------------------------------------------------

mean(cost_yll2)/mean(cost_yll1)

df_plot = data.frame(cost = c(cost_yll1,cost_yll2),
                     scen= c(rep("Vaccination",length(cost_yll1)),rep("No Vaccination",length(cost_yll2))))
# 

df_plot

theme_set(theme_void(base_family = "Helvetica"))

theme_update(
  axis.text.x = element_text(color = "black", face = "bold", size = 26, 
                             margin = margin(t = 6)),
  axis.text.y = element_text(color = "black", size = 22, hjust = 1, 
                             margin = margin(r = 6), family = "Helvetica"),
  axis.line.x = element_line(color = "black", size = 1),
  panel.grid.major.y = element_line(color = "grey90", size = .6),
  plot.background = element_rect(fill = "white", color = "white"),
  plot.margin = margin(rep(20, 4))
)

#15946e
## custom colors
#my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]
my_pal <- c(rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7)],"#15946e")

## theme for horizontal charts
theme_flip <-
  theme(
    axis.text.x = element_text(face = "plain", family = "Helvetica", size = 22),
    axis.text.y = element_text(face = "bold", family = "Helvetica", size = 30),
    axis.title.x = element_text(face = "bold", family = "Helvetica", size = 22),
    panel.grid.major.x = element_blank(),#element_line(color = "grey90", size = .6),
    panel.grid.major.y = element_blank(),
    axis.ticks.x = element_line(color = "black", size = 1.0),
    axis.ticks.length.x = unit(0.1,"cm"),
    legend.position = "none", 
    legend.text = element_text(family = "Helvetica", size = 18),
    legend.title = element_text(face = "bold", size = 18),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

label_t = c(TeX(r"($C_v$)"),TeX(r"($C_b$)"))

df_point = data.frame(x = c("Vaccination","Vaccination"),y = c(32,35.8))

ggplot(df_plot %>% filter(scen == "Vaccination"), aes(x = forcats::fct_rev(scen), y = cost/(1e9), 
                    color = scen, fill = scen)) +
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
  stat_summary(fun = mean,geom="point", color= "#88419d", fill = "#88419d",shape = 21, size = 2.5)+
  geom_point(data=df_point, aes(x = x,y=y),color= "#88419d", fill = "#88419d",shape = 23, size = 2.5)+
  scale_y_continuous(limits = c(30.5,37.5),breaks = seq(31,37), expand=expansion(mult=c(0,0)))+
  scale_x_discrete(breaks = c("Vaccination"),labels = label_t, expand=expansion(mult=c(0,0)))+
  labs(y="Cost (billion US$)",y=NULL)+
  coord_flip() +
  scale_color_manual(values = pal_plots[1], guide = "none") +
  scale_fill_manual(values = pal_plots[1], guide = "none") +
  theme_flip




ggsave(
  "../figures/VSLvac.pdf",
  device = "pdf",
  width = 7,
  height = 4,
  dpi = 300,
)

df_point = data.frame(x = c("No Vaccination","No Vaccination"),y = c(140.7,157.7))


ggplot(df_plot %>% filter(scen == "No Vaccination"), aes(x = forcats::fct_rev(scen), y = cost/(1e9), 
                    color = scen, fill = scen)) +
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
  scale_y_continuous(limits = c(137.5,162.5),breaks = seq(140,160,4),expand=expansion(mult=c(0,0)))+
  scale_x_discrete(breaks = c("No Vaccination"),labels = label_t, expand=expansion(mult=c(0,0)))+
  labs(y="Cost (billion US$)",y=NULL)+
  coord_flip() +
  scale_color_manual(values = pal_plots[3], guide = "none") +
  scale_fill_manual(values = pal_plots[3], guide = "none") +
  theme_flip




ggsave(
  "../figures/VSLnovac.pdf",
  device = "pdf",
  width = 7,
  height = 4,
  dpi = 300,
)



# Let's read the file containing the amount of people that died at age x in each sim

# 
# age_groups1 = c(0,18,25,35,45,55,65,75)
# age_groups2 = c(17,24,34,44,54,64,74,99)
# names = paste(age_groups1,"-",age_groups2)

# # 
# age_groups1 = c(0,10,20,30,40,50,60,70,80,90)
# age_groups2 = c(9,19,29,39,49,59,69,79,89,99)
# names = paste(age_groups1,"-",age_groups2)
# 
# 
# age_of_death= read.table("data/results_prob_0_121_1_newyorkcity/year_of_death.dat",h=F)
# dim(age_of_death)
# 
# life_exp = read.csv("data/life_exp.csv",sep = ";",h=F)$V2[1:nrow(age_of_death)]
# 
# 
# ff <- function(x){
#   aux = life_exp[(age_groups1[x]+1):(age_groups2[x]+1)]*age_of_death[(age_groups1[x]+1):(age_groups2[x]+1),]
#   
#   return(colSums(aux))
# }
# 
# bb <- function(x){
#   a = boot::boot(as.vector(x),fc,500)$t[,1]
#   c = c(mean(a),quantile(a,0.025,name=F),quantile(a,0.975,name=F))
# }
# 
# listb = lapply(seq(1:length(age_groups1)),ff)
# 
# ll = lapply(listb,bb)
# m1 = matrix(unlist(ll), ncol = 3,byrow = T)
# 
# 
# 
# age_of_death= read.table("data/results_prob_0_121_2_newyorkcity/year_of_death.dat",h=F)
# dim(age_of_death)
# 
# ff <- function(x){
#   aux = life_exp[(age_groups1[x]+1):(age_groups2[x]+1)]*age_of_death[(age_groups1[x]+1):(age_groups2[x]+1),]
#   
#   return(colSums(aux))
# }
# 
# bb <- function(x){
#   a = boot::boot(as.vector(x),fc,500)$t[,1]
#   c = c(mean(a),quantile(a,0.025,name=F),quantile(a,0.975,name=F))
# }
# 
# listb = lapply(seq(1:length(age_groups1)),ff)
# 
# ll = lapply(listb,bb)
# m2 = matrix(unlist(ll), ncol = 3,byrow = T)
# 
# 
# df = data.frame(group = c(names,names),mm = rbind(m1,m2),scen = c(rep("Vaccination",length(names)),rep("No Vaccination",length(names))))
# head(df)
# 
# df = df %>% rename(mean=mm.1,ci1=mm.2,ci2=mm.3) %>% mutate(mean = mean*factor_deaths*population/100000,
#                                                            ci1 = ci1*factor_deaths*population/100000,ci2=ci2*factor_deaths*population/100000)
# head(df)
# 
# theme_set(theme_bw())
# theme_flip <-
#   theme(
#     axis.text.x = element_text(color = "black",face = "bold", family = "Helvetica", size = 16,angle = 45,hjust=01.0,vjust=1.0),
#     axis.text.y = element_text(color = "black",face = "bold", family = "Helvetica", size = 20),
#     axis.title.y = element_text(face = "bold", family = "Helvetica", size = 24,angle = 90),
#     panel.grid.major.x = element_blank(),#element_line(color = "grey90", size = .6),
#     panel.grid.major.y = element_line(color = "grey90", size = .6),
#     axis.ticks.x = element_line(color = "black", size = 1.0),
#     axis.ticks.length.x = unit(0.1,"cm"),
#     legend.position = c(0.28,0.88),
#     legend.box.background = element_rect(color = "black", size = .6),
#     legend.text = element_text(face="bold",family = "Helvetica", size = 14),
#     legend.title = element_text(face = "bold", size = 18),
#     panel.border = element_rect(colour = "black", fill=NA, size=1)
#   )
# 
# ggplot(data = df)+
#   geom_col(aes(x=group,y=mean, color=scen,fill = scen,fill = after_scale(colorspace::lighten(fill, .5))),position=position_dodge2(),size=1.0)+
#   geom_errorbar(aes(x=group,ymin=ci1,ymax=ci2, color=scen),width=0.5,size = 1.25,position=position_dodge(0.9))+
#   scale_color_manual(values = my_pal[c(1,4)],name = NULL) +
#   scale_fill_manual(values = my_pal[c(1,4)],name = NULL) +
#   scale_y_continuous(expand=expansion(mult=c(0,0.05)), name = "Years of life lost")+
#   scale_x_discrete(name = NULL)+
#   theme_flip
# 
# 
# 
# ggsave(
#   "../figures/YLL.png",
#   device = "png",
#   width = 5.5,
#   height = 5,
#   dpi = 300
# )
# 
# 
# ggsave(
#   "../figures/YLL.pdf",
#   device = "pdf",
#   width = 5.5,
#   height = 5,
#   dpi = 300
# )
# 
# 
# df
# write.csv(df,"../figures/YLL.csv",row.names = F)
# 
# 
# ### Average
# 
# 
# ffm <- function(x){
#   aux = life_exp[(age_groups1[x]+1):(age_groups2[x]+1)]*age_of_death[(age_groups1[x]+1):(age_groups2[x]+1),]
#   
#   mm = sum(aux)/sum(age_of_death[(age_groups1[x]+1):(age_groups2[x]+1),])
#   
#   return(mm)
# }
# 
# 
# listb = lapply(seq(1:length(age_groups1)),ffm)
# 
# mm = matrix(unlist(listb), ncol = 1,byrow = T)
# row.names(mm) = names


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

dd = data.frame(ROI=bb$t[,1])
df_point = data.frame(x = c(as.factor(1),as.factor(1)),y = quantile(dd$ROI,c(0.025,0.975),name=F))


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
  "../figures/ROI_soc.pdf",
  device = "pdf",
  width = 7,
  height = 4,
  dpi = 300
)

