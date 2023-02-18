setwd("D:/PosDoc/Coronavirus/USomicron/USomicron/")
library(dplyr)
library(stringr)
library(ggplot2)
library(lubridate)
library(data.table)
theme_set(theme_bw())
enddate=as.Date("2022-12-12")#as.Date("2022-01-31")
startvacdate = as.Date("2020-12-12")
basedate = as.Date("2020-09-01")


population = 332968798
# Deaths per state --------------------------------------------------------

temp <- tempfile()
download.file("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv",temp)
data.cases <- read.csv(temp,h=T,stringsAsFactors = F)
unlink(temp)
rm(temp)
head(data.cases)
ignore_states = c("Virgin Islands","Guam","Puerto Rico","Northern Mariana Islands","American Samoa")
data.cases %>% filter(!(state %in% ignore_states)) %>% select(state) %>% unique() %>% pull(state) %>% length()

data.cases = data.cases %>% filter(!(state %in% ignore_states)) %>% 
  mutate(date=as.Date(date))

data.cases.test = data.cases %>% arrange(date)  %>% group_by(state) %>% 
  mutate(inc_cases = diff(c(0,cases)),inc_deaths = diff(c(0,deaths))) %>%
  filter(date<=enddate, date>= startvacdate) %>%
  mutate(inc_deaths = ifelse(inc_deaths<0,0,inc_deaths),inc_cases = ifelse(inc_cases<0,0,inc_cases))


ggplot()+geom_line(data = data.cases.test,aes(x=date,y=inc_cases, color = state))
#ggplot()+geom_line(data = data.cases.test,aes(x=date,y=inc_deaths, color = state))

total.deaths = data.cases.test %>% pull(inc_deaths) %>% sum()
total.deaths

# df = data.cases.test %>% group_by(state) %>% 
#   summarise(deaths.state = sum(inc_deaths),ratio = sum(inc_deaths)/total.deaths )%>%
#   rename(State = state, Total.Deaths=deaths.state,Ratio = ratio)


df = data.cases.test %>% group_by(state) %>% 
  summarise(deaths.state = sum(inc_deaths),cases.state = sum(inc_cases),ratio = sum(inc_deaths)/sum(inc_cases))%>%
  rename(State = state, Total.Deaths=deaths.state,Ratio = ratio)

sum(df$Ratio)


# Vaccination primary -----------------------------------------------------


#########################################3

library(readxl)
library(dplyr)
library(ggplot2)

### https://data.cdc.gov/Vaccinations/COVID-19-Vaccination-Demographics-in-the-United-St/km4m-vcsb

data = read.csv("../COVID-19_Vaccination_Demographics_in_the_United_States_National.csv")
head(data)

fields_data = unique(data$Demographic_category)
fields_data


# [5:11,12:17, 18:24, 25:39, 40:49, 50:64, 65:74, 75:100]
fields_data = fields_data[grepl("Ages_",fields_data,fixed = T)]
fields_data
fields_data = fields_data[-c(7,2,8,11)]
fields_data


field_correct = c("Ages_5-11_yrs","Ages_12-17_yrs","Ages_18-24_yrs",
                  "Ages_25-39_yrs","Ages_40-49_yrs","Ages_50-64_yrs",
                  "Ages_65-74_yrs","Ages_75+_yrs")

data = data %>% filter(Demographic_category %in% field_correct) %>% mutate(Date = as.Date(Date,"%m/%d/%Y")) %>%
  filter(Date <= enddate)
head(data)


ggplot()+geom_line(data = data, aes(x = Date,y = Administered_Dose1,color = Demographic_category))

nr = length(unique(data$Date))
nr
m1 = matrix(0,nr,length(field_correct))
data_ = data[order(data$Date),]
for(i in 1:length(field_correct)){
  
  d = data_$Date[data_$Demographic_category == field_correct[i]]
  x = data_$Administered_Dose1[data_$Demographic_category == field_correct[i]]
  
  x = c(0,x)
  dff_ = diff(x)
  
  m1[,i] = dff_
}

colnames(m1) = field_correct

head(m1)
m = round(m1/population*100000)
head(m)

write.table(m, "../dose_1_us_dez.dat", col.names = F, row.names = F)


df = stack(as.data.frame(m1))
df$idx = rep(seq(1,nrow(m1)),ncol(m1))
df$dose = rep("First Dose",nrow(df))

head(df)

ggplot()+geom_line(data = df,aes(x = idx ,y=values,color =ind))


m1 = matrix(0,nr,length(field_correct))
data_ = data[order(data$Date),]
for(i in 1:length(field_correct)){
  
  d = data_$Date[data_$Demographic_category == field_correct[i]]
  x = data_$Series_Complete_Yes[data_$Demographic_category == field_correct[i]]
  
  x = c(0,x)
  dff_ = diff(x)
  
  m1[,i] = dff_
}

colnames(m1) = field_correct

head(m1)
m = round(m1/population*100000)
head(m)

write.table(m,"../dose_2_us_dez.dat",col.names = F,row.names = F)


df2 = stack(as.data.frame(m1))
df2$idx = rep(seq(1,nrow(m1)),ncol(m1))
df2$dose = rep("Second dose",nrow(df2))

head(df)

df_ = rbind(df,df2)
ggplot()+geom_line(data = df_,aes(x = idx ,y=values,color =ind))+facet_grid(.~dose)



sum(df2$values)/population
sum(df$values)/population




df_2 = data[data$Date == as.Date("2021-12-01"),]
df_2
sum(df_2$Administered_Dose1)/population

v = min(data$Date)+seq(0,nrow(m1)-1)
sum(m1[v == as.Date("2021-11-19"),])



# Cases, Deaths,  Vaccination and Hospitalization -------------------------
#################################################################
#setwd("~/PosDoc/Coronavirus/USomicron/")
### https://covid.cdc.gov/covid-data-tracker/#vaccination-trends

#need to remove lines from header


## There is a problem here because on December 01 they removed the previous information about first and second booster
# So, we will go back to one single vector for booster and people may take one or two of them


data_boost = read.csv("../trends_in_number_of_covid19_vaccinations_in_the_us_firstbooster.csv")

data_boost = data_boost %>% #filter(Date.Type == "Admin") %>% 
  select(Date = `ï..Date`,Daily.Count.People.Receiving.a.First.Booster.Dose, Daily.Count.of.People.Ages.50...Receiving.a.Second.Booster.Dose) %>%
  rename(Booster = Daily.Count.People.Receiving.a.First.Booster.Dose, SecBooster = Daily.Count.of.People.Ages.50...Receiving.a.Second.Booster.Dose) %>% 
  mutate(Date = as.Date(Date), Booster = as.numeric(Booster), SecBooster = as.numeric(SecBooster)) %>%
  filter(Date <= enddate)


ggplot()+geom_col(data=data_boost,aes(x=Date,y=Booster))+scale_y_continuous()

data_boost$Booster100 = round(data_boost$Booster*100000/population)
data_boost$SecBooster100 = round(data_boost$SecBooster*100000/population)


min(data_boost$Date)-as.Date("2020-09-01")

write.table(c(0,data_boost$Booster100),"../booster1_aug.dat",row.names = F,col.names = F)

write.table(c(0,data_boost$SecBooster100),"../booster2_aug.dat",row.names = F,col.names = F)



ggplot()+geom_col(data=data_boost,aes(x=Date,y=Booster))+scale_y_continuous()

data_boost$Booster100 = round(data_boost$Booster*100000/population)

min(data_boost$Date)-as.Date("2020-09-01")

write.table(c(0,data_boost$Booster100),"../booster1_aug.dat",row.names = F,col.names = F)

sum(data_boost$Booster100)
aux = data_boost[data_boost$Date >= as.Date("2021-12-01") & data_boost$Date <= enddate,]
aa = mean(aux$Booster)
aa*100000/population
ggplot()+geom_col(data=aux,aes(x=Date,y=Booster))+geom_hline(yintercept = aa)






##############################3
# 

### https://data.cdc.gov/Vaccinations/COVID-19-Vaccination-Trends-in-the-United-States-N/rh2h-3yt2

# We decided to add up the Booster and Bivalent and use it after bivalent came up as a unique vector
# Whoever is eligible, will take it


data_boost = read.csv("../COVID-19_Vaccination_Trends_in_the_United_States_National_and_Jurisdictional.csv")

not_want = c("MH", "PR", "GU", "MP", "US", "PW", "FM", "HI", "AS")

df = data_boost %>% filter(!(Location %in% not_want), date_type == "Admin") %>%
  mutate(Date = mdy(Date)) %>%
  group_by(Date) %>%
  summarise(booster = sum(Booster_Daily), sec_boost_50 = sum(Second_Booster_50Plus_Daily),
            biv_boost = sum(Bivalent_Booster_Daily), state = paste0(Location, collapse = "; "))

df %>% View

df %>% group_by(Date) %>% summarise(tt = sum(booster))

df$m7biv = frollmean(df$biv_boost, 7, align = "right")

ggplot(df)+geom_line(aes(Date, biv_boost, color = "Bivalent"))+
  geom_line(aes(Date, sec_boost_50, color = "Sec booster 50+"))+
  geom_line(aes(Date, booster, color = "Booster"))+
  geom_line(aes(Date, m7biv, color = "mean"), size = 1.5)+
  scale_x_date(limits = as.Date(c("2022-06-01", "2022-12-01")))+
  scale_y_continuous(limits = c(0, 1000000))+
  scale_color_manual(values = c("black", "blue", "red", "purple"))



df$soma = df$booster+df$biv_boost
df$m7biv = frollmean(df$soma, 7, align = "right")

ggplot(df)+#geom_line(aes(Date, biv_boost, color = "Bivalent"))+
  geom_line(aes(Date, sec_boost_50, color = "Sec booster 50+"))+
  geom_line(aes(Date, soma, color = "Sum"))+
  geom_line(aes(Date, m7biv, color = "mean"), size = 1.5)+
  scale_x_date(limits = as.Date(c("2022-06-01", "2022-12-01")))+
  scale_y_continuous(limits = c(0, 1000000))+
  scale_color_manual(values = c("black", "blue", "red", "purple"))+
  geom_vline(aes(xintercept = as.Date("2022-08-30")))




df = df %>% mutate(soma100 = round(soma*100000/population), sec_boost50_100 = round(sec_boost_50*100000/population))


write.table(c(0,df$soma100),"../booster1_dez.dat",row.names = F,col.names = F)

write.table(c(0,df$sec_boost50_100),"../booster2_dez.dat",row.names = F,col.names = F)




#############################
#### hospitalization
### https://healthdata.gov/Hospital/COVID-19-Reported-Patient-Impact-and-Hospital-Capa/g62h-syeh


#############################
#### hospitalization
### https://healthdata.gov/Hospital/COVID-19-Reported-Patient-Impact-and-Hospital-Capa/g62h-syeh
library(zoo)

data.hos = read.csv("../COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_State_Timeseries.csv")
head(data.hos)
data.hos = data.hos[,c("state","date","previous_day_admission_adult_covid_confirmed","previous_day_admission_pediatric_covid_confirmed")]
head(data.hos)
na.fill(data.hos,as.numeric(0))

data.hos
unique(data.hos$state)

data.hos2 = data.hos %>% mutate(date = as.Date(date,"%Y/%m/%d")) %>% filter(!(state %in% c("PR","VI","AS")) ) %>% 
  group_by(date) %>% 
  summarise(hos_adult = sum(previous_day_admission_adult_covid_confirmed,na.rm=TRUE),
            hos_pediatric = sum(previous_day_admission_pediatric_covid_confirmed,na.rm=TRUE),
            hos = sum(previous_day_admission_adult_covid_confirmed+previous_day_admission_pediatric_covid_confirmed,na.rm=TRUE))%>%
  filter(date>=as.Date("2020-09-01"),date<=enddate)

sum(data.hos2$hos)
max(data.hos2$date)

ggplot()+geom_col(data=data.hos2,aes(x=date,y=hos))

# for entire US -----------------------------------------------------------


temp <- tempfile()
download.file("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv",temp)
data.cases <- read.csv(temp,h=T,stringsAsFactors = F)
unlink(temp)
rm(temp)


head(data.cases)
ignore_states = c("Virgin Islands","Guam","Puerto Rico","Northern Mariana Islands","American Samoa")
data.cases %>% filter(!(state %in% ignore_states)) %>% select(state) %>% unique() %>% pull(state) %>% length()

data.cases = data.cases %>% filter(!(state %in% ignore_states)) %>% 
  mutate(date=as.Date(date))

data.cases.test = data.cases %>% arrange(date)  %>% group_by(state) %>% 
  mutate(inc_cases = diff(c(0,cases)),inc_deaths = diff(c(0,deaths))) %>%
  filter(date<=enddate, date>= as.Date("2020-09-01")) %>%
  mutate(inc_deaths = ifelse(inc_deaths<0,0,inc_deaths),inc_cases = ifelse(inc_cases<0,0,inc_cases)) %>%
  group_by(date) %>% summarise(inc_cases = sum(inc_cases),inc_deaths = sum(inc_deaths))


head(data.cases.test)
ggplot()+geom_line(data = data.cases.test,aes(x=date,y=inc_cases))

ggplot()+geom_line(data = data.cases.test,aes(x=date,y=inc_deaths))

data.cases = cbind(data.cases.test,data.hos2[-1])
head(data.cases)



write.csv(data.cases,"../data_us_feb2023.csv",row.names=F)





# for entire US 2 -----------------------------------------------------------


temp <- tempfile()
download.file("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv",temp)
data.cases <- read.csv(temp,h=T,stringsAsFactors = F)
unlink(temp)
rm(temp)


head(data.cases)
ignore_states = c("Virgin Islands","Guam","Puerto Rico","Northern Mariana Islands","American Samoa")
data.cases %>% filter(!(state %in% ignore_states)) %>% select(state) %>% unique() %>% pull(state) %>% length()

data.cases = data.cases %>% filter(!(state %in% ignore_states)) %>% 
  mutate(date=as.Date(date))

data.cases.test = data.cases %>% arrange(date)  %>% group_by(state) %>% 
  mutate(inc_cases = diff(c(0,cases)),inc_deaths = diff(c(0,deaths))) %>%
  filter(date<=as.Date("2022-03-31"), date>= as.Date("2020-10-01")) %>%
  mutate(inc_deaths = ifelse(inc_deaths<0,0,inc_deaths),inc_cases = ifelse(inc_cases<0,0,inc_cases)) %>%
  group_by(date) %>% summarise(inc_cases = sum(inc_cases),inc_deaths = sum(inc_deaths))


head(data.cases.test)
ggplot()+geom_col(data = data.cases.test,aes(x=date,y=inc_cases))

ggplot()+geom_col(data = data.cases.test,aes(x=date,y=inc_deaths))

data.cases = cbind(data.cases.test,data.hos2[-1])
head(data.cases)

write.csv(data.cases.test,"data_us_seyed.csv",row.names=F)




