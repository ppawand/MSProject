####################################################
# library
####################################################
library(readxl)
library(tidyverse)
library(lme4)
library(effects)
library(ggpubr)
library(car)
library(MuMIn) # to calculate R2
library(emmeans)

##########################################
# theme
####################################################

theme <- theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill = NA, linewidth = 1),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.border = element_rect(linewidth  = 1))

################################################
#Temperature

temp <- read_excel("~/Desktop/MSProject/MS_Paper/data/TM.final.xlsx", sheet = 2)
temp <- temp %>%
  pivot_longer(c(C, W, R, WR),
               names_to = "trt",
               values_to = "temperature",
               values_drop_na = TRUE) %>%
  filter(temperature > 10)

boxplot(temp$temperature)

# removing outliers
Q1 <- quantile(temp$temperature, .25)
Q3 <- quantile(temp$temperature, .75)
IQR <- IQR(temp$temperature)
temp <- subset(temp, temp$temperature > (Q1 - 1.5*IQR) & temp$temperature < (Q3 + 1.5*IQR))
boxplot(temp$temperature)

# adding date and month columns
temp <- temp %>%
  mutate(datee = format(date, format = "%F"),
         month = format(date, format = "%b"),
         day = format(date, format = "%d")) %>% 
  filter(month != "Nov")

####################################################
#  Soil Moisture

moist<- read_excel("~/Desktop/MSProject/MS_Paper/data/TM.final.xlsx", sheet = 3)
moist <- moist %>%
  pivot_longer(c(C, W, R, WR),
               names_to = "trt",
               values_to = "moisture",
               values_drop_na = TRUE)%>%
  filter(moisture > 0)

# removing outliers
Q1m <- quantile(moist$moisture, .25)
Q3m <- quantile(moist$moisture, .75)
IQRm <- IQR(moist$moisture)
moist <- subset(moist, moist$moisture > (Q1m - 1.5*IQRm) & moist$moisture < (Q3m + 1.5*IQRm))
boxplot(moist$moisture)

# adding date and month columns
moist <- moist %>% 
  mutate(datee = format(date, format = "%F"),
         month = format(date, format = "%b"),
         day = format(date, format = "%d")) %>% 
  filter(month != "Nov")

####################################################

# daily average temperature
temp.average <- temp %>% 
  group_by(wp, plot, trt, irrigation, datee, day, month) %>% 
  summarize(mean.temp = mean(temperature),
            max.temp = max(temperature),
            min.temp = min(temperature),
            dtr = max.temp - min.temp)

# daily average moisture 
moist.average <- moist %>% 
  group_by(wp, plot, trt, irrigation, datee, day, month) %>% 
  summarize(mean.moist = mean(moisture))

####################################################
# merging both files

daily.TM <- merge(temp.average, moist.average)
daily.TM <- daily.TM[-1, ] # removing first row - not a full day average
daily.TM$warming <-ifelse(daily.TM$trt %in% c("C", "R"), "Ambient", "OTC")
daily.TM$residue <-ifelse(daily.TM$trt %in% c("C", "W"), "noResidue", "Residue")
daily.TM <- daily.TM %>% 
  relocate(warming, residue, .before = irrigation)
daily.TM[, c(1:6, 9)] <- lapply(daily.TM[, c(1:6, 9)], factor)
daily.TM$trt <- factor(daily.TM$trt, levels = c("C", "W", "R", "WR"))
daily.TM$month <- factor(daily.TM$month, levels = c("Jun", "Jul","Aug", "Sep", "Oct", "Nov"))
daily.TM$datee <- strptime(daily.TM$datee, format  = "%Y-%m-%d")
daily.TM$datee <- as.POSIXct(daily.TM$datee)


# averaging temperatures by hours of the day


hourly.avg.temp <- temp %>% 
  mutate(hours = lubridate::hour(date)) %>% 
  group_by(wp, plot, irrigation, trt, month, hours) %>% 
  summarize(temp.avg = mean(temperature))

hourly.avg.temp$warming <-ifelse(hourly.avg.temp$trt %in% c("C", "R"), "Ambient", "OTC")
hourly.avg.temp$residue <-ifelse(hourly.avg.temp$trt %in% c("C", "W"), "noResidue", "Residue")


# separating night and day temperature. For this we use 7 AM - 9 PM (14 hours) 
# as daylight hours
day.temp <- hourly.avg.temp %>%
  filter(hours %in% c(7:21)) %>%
  filter(month != "Oct")


night.temp <- hourly.avg.temp %>%
  filter(hours %in% !c(7:21)) %>%
  filter(month != "Oct")
  
# fitting models for night and day temperature separately
day_temp_model <- lmer(temp.avg ~ irrigation*warming*residue + (1|wp) + (1|month) +
                         (1|plot),
                       data = day.temp)

Anova(day_temp_model)
summary(day_temp_model)
plot(day_temp_model)
qqnorm(residuals(temp_model))
hist(residuals(temp_model))
r.squaredGLMM(temp_model)


night_temp_model <- lmer(temp.avg ~ irrigation*warming*residue + (1|wp) + (1|month) +
                         (1|plot),
                       data = night.temp)
Anova(night_temp_model)
summary(night_temp_model)
plot(night_temp_model)
qqnorm(residuals(night_temp_model))
hist(residuals(night_temp_model))

# plot the hourly temperature

hourly.temp <- hourly.avg.temp %>% 
  group_by(irrigation, trt, hours) %>% 
  summarize(avg.temp = mean(temp.avg),
            sd = sd(temp.avg),
            n = n(),
            se = sd/sqrt(n))

(hourly.stemp.plot <- ggplot(hourly.temp, aes(x = hours, y =  avg.temp, color = trt)) +
  geom_line(linewidth = 1) +
  #geom_vline(xintercept = c(7, 21), lty = "dotted", color = "red") +
  geom_errorbar(aes(ymin = avg.temp - se,
                    ymax = avg.temp + se)) +
  facet_wrap(~irrigation) +
  scale_x_continuous(breaks = seq(0, 23, 2),
                     labels = c("12 AM", "2 AM", "4 AM", "6 AM", " 8 AM",
                                "10 AM", "12 PM", "2 PM", "4 PM", "6 PM",
                                "8 PM", "10 PM")) +
  labs(x = "Time of the day", y = expression("Soil Temperature ("*~degree*C*")")) +
  scale_color_manual(values = c("#131E3A", "#0018F9","#980019", "#4C0099"),
                     labels = c("C", "R", "OTC", "OTC + R")) +
  theme +
  theme(axis.text.x=element_text(angle=90,hjust=1)))
  


###########################################################
# calculating weekly mean

weekly.TM <- daily.TM %>%
  group_by(irrigation, trt, week = cut(datee, "week", start.on.monday = FALSE)) %>%
  summarise(temp = mean(mean.temp),
            moist = mean(mean.moist),
            max.tempr = mean(max.temp),
            min.tempr = mean(min.temp),
            mean.dtr = mean(dtr))

weekly.TM$week <- strptime(weekly.TM$week, format  = "%Y-%m-%d")
weekly.TM$week <- as.POSIXct(weekly.TM$week)

#####################################################

# temporal variation of temperature (graph)
# mean temperature
(fig1 <- ggplot(weekly.TM, aes(week, temp, color = trt))+
  geom_line(linewidth = 1) + 
   scale_color_manual(values = c("#131E3A", "#980019", "#0018F9", "#4C0099"),
                      labels = c("C", "OTC", "R", "OTC + R")) +
  xlab( "Months") +
  ylab(expression("Soil Temperature ("*~degree*C*")")) + 
  facet_wrap(~irrigation) + 
   scale_y_continuous(limits = c(15, 35)) +
   theme)


# daily temperature range
(dt <- ggplot(weekly.TM)+
    geom_line(aes(x = week, y = mean.dtr, color = trt), linewidth = 1) + 
    scale_color_manual(values = c("#131E3A", "#980019", "#0018F9", "#4C0099"),
                       labels = c("C", "OTC", "R", "OTC + R")) +
    xlab( "Months") +
    ylab(expression("Daily temperature range ("*~degree*C*")")) + 
    facet_wrap(~ irrigation) + theme)
  
#######################################################
# precipitation
precipitation <- read_excel("~/Desktop/MSProject/MS_Paper/data/weather.data.xlsx")
weekly.rain <- precipitation %>%
  group_by(week = cut(date.time, "week", start.on.monday = FALSE)) %>%
  summarise(rain = sum(precipitation))

weekly.rain$week <- strptime(weekly.rain$week, format  = "%Y-%m-%d")
weekly.rain$week <- as.POSIXct(weekly.rain$week)


####################################################
# temporal variation of moisture (graph)
(fig2 <- ggplot(weekly.TM, aes(week, color = trt)) + 
  geom_line(aes(y = moist), linewidth = 1) +  
   geom_bar(data = weekly.rain, aes(x = week, y = rain/100), width = 0.01, 
                                    stat = "identity", 
                                    fill = "deepskyblue3" , 
                                    color = "deepskyblue3") +
   scale_y_continuous(sec.axis = sec_axis (~ .*100, name = "Rainfall (mm)")) +
   scale_color_manual(values = c("#131E3A", "#980019", "#0018F9", "#4C0099"),
                      labels = c("C", "OTC", "R", "OTC + R")) +
  xlab( "Months") +
  ylab(expression("VWC ("*m^3/m^3*")")) + 
  facet_wrap(~irrigation) + theme)

################################################################

# summarizing by month
monthly.avg.TM <- daily.TM %>% 
  group_by(wp, plot, warming, residue, irrigation, month) %>% 
  summarize(moisture = mean(mean.moist),
            temp = mean(mean.temp),
            mean.dtr = mean(dtr))

# removing Oct temperature to include only peak season data for analysis
monthly.avg.TM <- monthly.avg.TM %>%
  filter(month != "Oct")

# temperature model

temp_model <- lmer(temp ~ warming * residue * irrigation +
                           (1|wp) + (1|month) + (1|plot),
                         data = monthly.avg.TM)


Anova(temp_model)
summary(temp_model)
plot(temp_model)
qqnorm(residuals(temp_model))
hist(residuals(temp_model))
r.squaredGLMM(temp_model)


# pairwise comparison
# interaction was not significant. we ran pairwise comparison for the
# significant main effects.
temp.res.emm <- emmeans(temp_model, ~ residue)
contrast(temp.res.emm, "consec", simple = "each", combine = TRUE)

temp.irr.emm <- emmeans(temp_model, ~ irrigation)
contrast(temp.irr.emm, "consec", simple = "each", combine = TRUE)


####################################################
# temp by trt (graph)

temp.plot.avg <- 
  as.data.frame(emmeans(temp_model,  ~ warming * residue * irrigation))

(fig3 <- ggplot(monthly.avg.TM, aes(x = warming, color = residue))  +
    geom_boxplot(aes(y = temp, color  = residue), alpha = 0.5) +
    geom_point(aes(y = temp, color  = residue), alpha = 0.5,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3)) +
    geom_point(data = temp.plot.avg, 
             aes(x = warming, y = emmean, group = residue),
             size = 5, shape = 15,
             position = position_dodge(width = 0.75)) +
    geom_errorbar(data = temp.plot.avg,
                aes(x = warming, ymin = lower.CL, 
                    ymax = upper.CL, group = residue), width = 0.2,
                stat = "identity", linewidth = 1,
                position = position_dodge(width = 0.75)) +
    scale_color_manual(labels = c("No Residue", "Residue"),
                     values = c("#000000", "#009E73")) +
    facet_grid( ~ irrigation) +
    labs(x="Treatments", y = expression("Soil Temperature ("*~degree*C*")")) +
    theme)
##################################################
# daily temperature range model

dtr_model <- lmer(mean.dtr ~ warming * residue * irrigation + (1|wp) + 
                    (1|month) + (1|plot),
                  data = monthly.avg.TM)

Anova(dtr_model)
summary(dtr_model)
plot(dtr_model)
qqnorm(residuals(dtr_model))
hist(residuals(dtr_model))
r.squaredGLMM(dtr_model)

# 

# pairwise comparison
# interaction was not significant. we ran pairwise comparison for the
# significant main effects.
dtr.res.emm <- emmeans(dtr_model, ~ residue)
contrast(dtr.res.emm, "consec", simple = "each", combine = TRUE)

dtr.otc.emm <- emmeans(dtr_model, ~ warming)
contrast(dtr.otc.emm, "consec", simple = "each", combine = TRUE)


dtr.plot.avg <- 
  as.data.frame(emmeans(dtr_model,  ~ warming * residue * irrigation))

(dtr.graph <- ggplot(monthly.avg.TM, aes(x = warming, color = residue))  +
    geom_boxplot(aes(y = mean.dtr, color  = residue), alpha = 0.5) +
    geom_point(aes(y = mean.dtr, color  = residue), alpha = 0.5,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3)) +
    geom_point(data = dtr.plot.avg, 
               aes(x = warming, y = emmean, group = residue),
               size = 5, shape = 15,
               position = position_dodge(width = 0.75)) +
    geom_errorbar(data = dtr.plot.avg,
                  aes(x = warming, ymin = lower.CL, 
                      ymax = upper.CL, group = residue), width = 0.2,
                  stat = "identity", linewidth = 1,
                  position = position_dodge(width = 0.75)) +
    scale_color_manual(labels = c("No Residue", "Residue"),
                       values = c("#000000", "#009E73")) +
    labs(x="Treatments", y = expression("Daily temperature range ("*~degree*C*")")) +
    facet_grid( ~ irrigation) +
    theme +
    theme(axis.title.x = element_blank()))




##################################################
# moisture model

moist_model <- lmer(moisture ~ warming * residue * irrigation + (1|wp) + 
                      (1|month) + (1|plot),
                    data = monthly.avg.TM)

Anova(moist_model)
summary(moist_model)
plot(moist_model)
qqnorm(residuals(moist_model))
hist(residuals(moist_model))
r.squaredGLMM(moist_model)

#pairwise comparison
moist.emm <- emmeans(moist_model, ~ warming:residue:irrigation)
contrast(moist.emm, "consec", simple = "each", combine = TRUE)

####################################################
# moisture by trt

moist.plot.avg <- 
  as.data.frame(emmeans(moist_model,  ~ warming * residue * irrigation))

(moist.graph <- ggplot(monthly.avg.TM, aes(x = warming, color = residue))  +
    geom_boxplot(aes(y = moisture, color  = residue), alpha = 0.5) +
    geom_point(aes(y = moisture, color  = residue), alpha = 0.5,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3)) +
    geom_point(data = moist.plot.avg, 
               aes(x = warming, y = emmean, group = residue),
               size = 5, shape = 15,
               position = position_dodge(width = 0.75)) +
    geom_errorbar(data = moist.plot.avg,
                  aes(x = warming, ymin = lower.CL, 
                      ymax = upper.CL, group = residue), width = 0.2,
                  stat = "identity", linewidth = 1,
                  position = position_dodge(width = 0.75)) +
    scale_color_manual(labels = c("No Residue", "Residue"),
                       values = c("#000000", "#009E73")) +
    labs(x="Treatments", y = expression("VWC ("*m^3/m^3*")")) +
    facet_grid( ~ irrigation) +
    theme)



#################################################

# air temperature
air.temp <- read_excel("~/Desktop/MSProject/MS_Paper/data/hobo.final.xlsx")

# removing outliers
boxplot(air.temp$temperature)
# no outliers
# adding date and month column
air.temp <- air.temp %>%
  mutate(datee = format(date, format = "%F"),
         month = format(date, format = "%b"),
         day = format(date, format = "%d"))

air.temp$warming <-ifelse(air.temp$trt %in% c("C", "R"), "Ambient", "OTC")

air.temp$residue <-ifelse(air.temp$trt %in% c("C", "W"), "noResidue", "Residue")

# daily average temperature
air.temp.average <- air.temp %>% 
  group_by(wp, trt, irrigation, datee, day, month) %>% 
  summarize(mean.temp = mean(temperature),
            mean.light = mean(light))

air.temp.average$warming <-ifelse(air.temp.average$trt %in% c("C", "R"), "Ambient", "OTC")

air.temp.average$residue <-ifelse(air.temp.average$trt %in% c("C", "W"), "noResidue", "Residue")

air.temp.average <- air.temp.average %>%  
  relocate(warming, residue, .before = irrigation)

air.temp.average[, c(1:5, 8)] <- lapply(air.temp.average[, c(1:5, 8)], factor)
air.temp.average$trt <- factor(air.temp.average$trt, levels = c("C", "W", "R", "WR"))
air.temp.average$month <- factor(air.temp.average$month, levels = c("Jul","Aug", "Sep", "Oct"))
air.temp.average$datee <- strptime(air.temp.average$datee, format  = "%Y-%m-%d")
air.temp.average$datee <- as.POSIXct(air.temp.average$datee)

#########################################################

# averaging temperatures by hours of the day


hourly.avg.airtemp <- air.temp %>% 
  mutate(hours = lubridate::hour(date)) %>% 
  group_by(warming, hours) %>% 
  summarize(temp.avg = mean(temperature),
            sd = sd(temperature),
            n = n(),
            se = sd/sqrt(n))


 (hourly.atemp.plot <- ggplot(hourly.avg.airtemp, aes(x = hours, y =  temp.avg, color = warming)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = temp.avg - se,
                    ymax = temp.avg + se)) +
  scale_x_continuous(breaks = seq(0, 23, 2),
                     labels = c("12 AM", "2 AM", "4 AM", "6 AM", " 8 AM",
                                "10 AM", "12 PM", "2 PM", "4 PM", "6 PM",
                                "8 PM", "10 PM")) +
  labs(x = "Time of the day", y = expression("Air Temperature ("*~degree*C*")")) +
  scale_color_manual(values = c("#131E3A", "#980019"),
                     labels = c("C","OTC")) +
  theme +
  theme(axis.text.x=element_text(angle=90,hjust=1)))



# temporal variation of air temperature (graph)
# weekly mean

weekly.air.temperature <- air.temp.average %>%
  group_by(irrigation, trt, week = cut(datee, "week", start.on.monday = FALSE)) %>%
  summarise(air.temp = mean(mean.temp))
           
weekly.air.temperature$week <- strptime(weekly.air.temperature$week, format  = "%Y-%m-%d")
weekly.air.temperature$week <- as.POSIXct(weekly.air.temperature$week)

(fig6 <- ggplot(weekly.air.temperature, aes(week, air.temp, color = trt)) +
  geom_line(linewidth = 1) + 
  scale_color_manual(values = c("#131E3A", "#980019", "#0018F9", "#4C0099"),
                       labels = c("C", "OTC", "R", "OTC + R")) +
  labs(x = "Months", y = (expression("Air temperature ("*~degree*C*")"))) + 
  facet_wrap(~irrigation) + 
    theme +
    theme(strip.background = element_blank(),
          strip.text = element_blank()))



########################################################
# air temperature model

monthly.air.temp <- air.temp.average %>% 
  group_by(wp, warming, residue, irrigation, month) %>% 
  summarize(temp = mean(mean.temp),
            light = (mean(mean.light))) %>%
  filter(month != "Oct")


air_temp_model <- lmer(temp ~ warming * residue * irrigation +
                         (1|wp) + (1|month),
     data = monthly.air.temp)

Anova(air_temp_model)
summary(air_temp_model)
plot(air_temp_model)
qqnorm(residuals(air_temp_model))
hist(residuals(air_temp_model))
r.squaredGLMM(air_temp_model)

# pairwise comparison
air.temp.emm <-emmeans(air_temp_model,  ~ irrigation)
contrast(air.temp.emm,"consec", simple = "each", combine = TRUE)



####################################################
# air temp by trt (graph)

airtemp.plot.avg <- 
  as.data.frame(emmeans(air_temp_model,  ~ warming * residue * irrigation))

(fig7 <- ggplot(monthly.air.temp, aes(x = warming, color = residue))  +
    geom_boxplot(aes(y = temp, color  = residue), alpha = 0.5) +
    geom_point(aes(y = temp, color  = residue), alpha = 0.5,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3)) +
    geom_point(data = airtemp.plot.avg, 
               aes(x = warming, y = emmean, group = residue),
               size = 5, shape = 15,
               position = position_dodge(width = 0.75)) +
    geom_errorbar(data = airtemp.plot.avg,
                  aes(x = warming, ymin = lower.CL, 
                      ymax = upper.CL, group = residue), width = 0.1,
                  stat = "identity", linewidth = 0.8,
                  position = position_dodge(width = 0.75)) +
    scale_color_manual(labels = c("No Residue", "Residue"),
                       values = c("#000000", "#009E73")) +
    labs(x="Treatments", y = expression("Air temperature ("*~degree*C*")")) +
    facet_grid( ~ irrigation) +
    theme +
    theme(strip.background = element_blank(),
          strip.text = element_blank()))


####################################################


####################################################



# Plots

# figure 1

ggarrange(fig1 + rremove("xlab") + rremove("x.text"), 
          fig6 + rremove("xlab"),
          nrow = 2,
          common.legend = T,
          labels = "auto",
          vjust = 0.5)


# figure 2
dtr.graph

# figure 3
ggarrange(fig2 + rremove("xlab"), 
          moist.graph + rremove("xlab"),
          nrow = 2,
          common.legend = F,
          labels = "auto")

# supplemental figure S1


ggarrange(hourly.stemp.plot + rremove ("xlab"), 
          hourly.atemp.plot,
          nrow = 2,
          common.legend = F, 
          labels = "auto")

# Supplemental figure S2
dt

# supplemental figure S3
ggarrange(fig3 + rremove("xlab") + rremove("x.text"),
          fig7 + rremove("xlab"),
          nrow = 2,
          common.legend = T,
          labels = "auto")




