####################################################
# library
####################################################
library(readxl)
library(tidyverse)
library(lme4)
library(cowplot)
library(scales)
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
        strip.background = element_rect(fill = NA, size = 1),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.border = element_rect(size  = 1))

################################################
#Temperature

temp <- read_excel("~/Desktop/MSProject/data/TM.final.xlsx", sheet = 2)
temp <- temp %>%
  pivot_longer(c(C, W, R, WR),
               names_to = "trt",
               values_to = "temperature",
               values_drop_na = TRUE) %>%
  filter(temperature > 10)

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

moist<- read_excel("~/Desktop/MSProject/data/TM.final.xlsx", sheet = 3)
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
  group_by(block, plot, trt, irrigation, datee, day, month) %>% 
  summarize(mean.temp = mean(temperature),
            max.temp = max(temperature),
            min.temp = min(temperature),
            dtr = max.temp - min.temp)

# daily average moisture 
moist.average <- moist %>% 
  group_by(block, plot, trt, irrigation, datee, day, month) %>% 
  summarize(mean.moist = mean(moisture))

####################################################
# merging both files

daily.TM <- merge(temp.average, moist.average)
daily.TM <- daily.TM[-1, ] # removing first row - not a full day average
daily.TM$warming <-ifelse(daily.TM$trt %in% c("C", "R"), "Ambient", "OTC")
daily.TM$residue <-ifelse(daily.TM$trt %in% c("C", "W"), "noResidue", "Residue")

daily.TM[, c(1:4, 7, 13, 14)] <- lapply(daily.TM[, c(1:4, 7, 13, 14)], factor)
daily.TM$trt <- factor(daily.TM$trt, levels = c("C", "W", "R", "WR"))
daily.TM$month <- factor(daily.TM$month, levels = c("Jun", "Jul","Aug", "Sep", "Oct", "Nov"))
daily.TM$datee <- strptime(daily.TM$datee, format  = "%Y-%m-%d")
daily.TM$datee <- as.POSIXct(daily.TM$datee)
str(daily.TM)


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
(temp_graph <- ggplot(weekly.TM, aes(week, temp, color = trt))+
  geom_line(size = 1) + 
   scale_color_manual(values = c("#131E3A", "#980019", "#0018F9", "#660033"),
                      labels = c("C", "OTC", "R", "OTC + R")) +
  xlab( "Months") +
  ylab(expression("Soil Temperature ("*~degree*C*")")) + 
  facet_wrap(~irrigation) + 
   scale_y_continuous(limits = c(15, 35)) +
   theme)

# max and min temp
(max_min <- ggplot(weekly.TM)+
    geom_line(aes(x = week, y = max.tempr), color = "black") +
    geom_line(aes(x = week, y = min.tempr), color = "green") +
    xlab( "Months") +
    ylab(expression("Soil Temperature ("*~degree*C*")")) + 
    facet_grid(irrigation~trt) + theme)

# daily temperature range
(dt <- ggplot(weekly.TM)+
    geom_line(aes(x = week, y = mean.dtr, color = trt), size = 1) + 
    scale_color_manual(values = c("#131E3A", "#980019", "#0018F9", "#660033"),
                       labels = c("C", "OTC", "R", "OTC + R")) +
    ylab(expression("Daily temperature range ("*~degree*C*")")) + 
    facet_wrap(~ irrigation) + theme +
    theme(axis.title.x = element_blank()))

#######################################################
# precipitation
precipitation <- read_excel("~/Desktop/MSProject/data/weather.data.xlsx")
str(precipitation)
weekly.rain <- precipitation %>%
  group_by(week = cut(date.time, "week", start.on.monday = FALSE)) %>%
  summarise(rain = sum(precipitation))

weekly.rain$week <- strptime(weekly.rain$week, format  = "%Y-%m-%d")
weekly.rain$week <- as.POSIXct(weekly.rain$week)


####################################################
# temporal variation of moisture (graph)
(vwc_graph <- ggplot(weekly.TM, aes(week, color = trt)) + 
  geom_line(aes(y = moist), size = 1) +  
   geom_bar(data = weekly.rain, aes(x = week, y = rain/100), width = 0.01, 
                                    stat = "identity", 
                                    fill = "deepskyblue3" , 
                                    color = "deepskyblue3") +
   scale_y_continuous(sec.axis = sec_axis (~ .*100, name = "Rainfall (mm)")) +
   scale_color_manual(values = c("#131E3A", "#980019", "#0018F9", "#660033"),
                      labels = c("C", "OTC", "R", "OTC + R")) +
  xlab( "Months") +
  ylab(expression("Volumetric water content ("*m^3/m^3*")")) + 
  facet_wrap(~irrigation) + theme)

################################################################

# temperature model

monthly.avg.TM <- daily.TM %>% 
  group_by(block, plot, warming, residue, irrigation, month) %>% 
  summarize(moisture = mean(mean.moist),
            temp = mean(mean.temp),
            mean.dtr = mean(dtr))


# temperature model

temp_model <- lmer(temp ~ warming * residue * irrigation +
                           (1|block) + (1|month),
                         data = monthly.avg.TM)

Anova(temp_model, type = 3)
summary(temp_model)
plot(temp_model)
qqnorm(residuals(temp_model))
hist(residuals(temp_model))
r.squaredGLMM(temp_model)

# pairwise comparison
emmeans(temp_model, specs = pairwise ~ residue) # residue decreased temp by 0.6
emmeans(temp_model, specs = pairwise ~ irrigation) # irrigation decreased temp by 0.8

# temp vs moisture
temp2 <- lmer(temp ~ moisture + (1|block) + (1|month),
              data = monthly.avg.TM )
Anova(temp2, type =3)
r.squaredGLMM(temp2)

####################################################
# temp by trt (graph)

temp.avg <-
  as.data.frame(emmeans(temp_model, ~warming * residue * irrigation)) 

(temp.avg.graph <- ggplot(temp.avg, aes(warming, emmean))  +
  geom_point(aes(color = residue), size = 3, shape = 15, alpha = 0.8,
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL, color = residue), width = 0.05,alpha = 0.8,
                position = position_dodge(width = 0.2)) +
  scale_color_manual(labels = c("No Residue", "Residue"),
                     values = c("red3", "blue3")) +
  labs(x="Treatments", y = expression("Soil Temperature ("*~degree*C*")")) +
  facet_grid( ~ irrigation) +
    scale_y_continuous(limits = c(18, 32)) +
  theme)

##################################################
# moisture model


moist_model <- lmer(moisture ~ warming * residue * irrigation +
                      (1|block) + (1|month),
           data = monthly.avg.TM)

Anova(moist_model, type = 3)
summary(moist_model)
plot(moist_model)
qqnorm(residuals(moist_model))
hist(residuals(moist_model))
r.squaredGLMM(moist_model)

#pairwise comparison
emmeans(moist_model, specs = pairwise ~ warming:residue:irrigation)

####################################################
# moisture by trt

moist.avg <-
  as.data.frame(emmeans(moist_model, ~warming * residue * irrigation))

(moist_avg_graph <- ggplot(moist.avg, aes(warming, emmean)) +
  geom_point(aes(color = residue), size = 3,shape =15, alpha = 0.8,
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL, color = residue), width = 0.05,alpha= 0.8,
                position = position_dodge(width = 0.2)) +
  labs(x="Treatments", y = expression("Volumetric water content ("*m^3/m^3*")")) +
  scale_color_manual(labels = c("No Residue", "Residue"),
                     values = c("red3", "blue3")) +
    facet_grid( ~ irrigation) +
    scale_y_continuous(limits = c(0.06, 0.25)) +
    theme)


#################################################

# air temperature
air.temp <- read_excel("~/Desktop/MSProject/data/hobo.final.xlsx")

# removing outliers
boxplot(air.temp$temperature)
# no outliers
# adding date and month column
air.temp <- air.temp %>%
  mutate(datee = format(date, format = "%F"),
         month = format(date, format = "%b"),
         day = format(date, format = "%d"))

# daily average temperature
air.temp.average <- air.temp %>% 
  group_by(block, trt, irrigation, datee, day, month) %>% 
  summarize(mean.temp = mean(temperature),
            mean.light = mean(light))

air.temp.average$warming <-ifelse(air.temp.average$trt %in% c("C", "R"), "Ambient", "OTC")

air.temp.average$residue <-ifelse(air.temp.average$trt %in% c("C", "W"), "noResidue", "Residue")

air.temp.average[, c(1:3, 6,9,10)] <- lapply(air.temp.average[, c(1:3, 6, 9, 10)], factor)
air.temp.average$trt <- factor(air.temp.average$trt, levels = c("C", "W", "R", "WR"))
air.temp.average$month <- factor(air.temp.average$month, levels = c("Jul","Aug", "Sep", "Oct"))
air.temp.average$datee <- strptime(air.temp.average$datee, format  = "%Y-%m-%d")
air.temp.average$datee <- as.POSIXct(air.temp.average$datee)

#########################################################
# temporal variation of air temperature (graph)
#facet by trt
# weekly mean

weekly.air.temperature <- air.temp.average %>%
  group_by(irrigation, trt, week = cut(datee, "week", start.on.monday = FALSE)) %>%
  summarise(air.temp = mean(mean.temp))
           
weekly.air.temperature$week <- strptime(weekly.air.temperature$week, format  = "%Y-%m-%d")
weekly.air.temperature$week <- as.POSIXct(weekly.air.temperature$week)

(air.temp.trend <- ggplot(weekly.air.temperature, aes(week, air.temp, color = trt)) +
  geom_line(size = 1) + 
    scale_color_manual(values = c("#131E3A", "#980019", "#0018F9", "#660033"),
                       labels = c("C", "OTC", "R", "OTC + R")) +
  labs(x = "Months", y = (expression("Air temperature ("*~degree*C*")"))) + 
  facet_wrap(~irrigation) + 
    
    theme)

########################################################
# air temperature model

monthly.air.temp <- air.temp.average %>% 
  group_by(block, warming, residue, irrigation, month) %>% 
  summarize(temp = mean(mean.temp),
            light = (mean(mean.light)))


air_temp_model <- lmer(temp ~ warming * residue * irrigation 
                       + (1|block) + (1|month),
     data = monthly.air.temp)

Anova(air_temp_model, type = 3)
summary(air_temp_model)
plot(air_temp_model)
qqnorm(residuals(air_temp_model))
hist(residuals(air_temp_model))
r.squaredGLMM(air_temp_model)

# pairwise comparison
emmeans(air_temp_model, specs = pairwise ~ warming) # OTC increased temp by 2.2

# air temperature model continuous predictors

airtemp.scaled <- as.data.frame(scale(monthly.air.temp[, 6:7], center = F))
airtemp.scaled$block <-monthly.air.temp$block
airtemp.scaled$month <-monthly.air.temp$month

air_temp_model1 <- lmer(temp ~ light + (1|block) + (1|month),
                        data = airtemp.scaled)

Anova(air_temp_model1, type = 3)
summary(air_temp_model)
plot(air_temp_model1)
qqnorm(residuals(air_temp_model1))
hist(residuals(air_temp_model1))
r.squaredGLMM(air_temp_model1)



####################################################
# air temp by trt (graph)
airtemp.avg <- 
  as.data.frame(emmeans(air_temp_model, ~warming * residue * irrigation))

(airtemp.avg.graph <- ggplot(airtemp.avg, aes(warming, emmean)) +
  geom_point(aes(color = residue), size = 3, shape = 15, alpha = 0.8,
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL, color = residue), width = 0.05, alpha = 0.8,
                position = position_dodge(width = 0.2)) +
  labs(x="Treatments", y = expression("Air temperature ("*~degree*C*")")) +
  scale_color_manual(labels = c("No Residue", "Residue"),
                     values = c("red3", "blue3")) +
  facet_wrap(~irrigation) + 
    scale_y_continuous(limits = c(10,45)) +
  theme)
  


####################################################
# light between OTC and control


lightOTC <- monthly.air.temp %>% group_by(block, warming) %>%
  summarise(mean.light = mean(light))

light_model <- lmer(mean.light ~ warming + (1|block),
                    data = lightOTC)


Anova(light_model, type=3)

summary(light_model)
# pairwise comparision
emmeans(light_model, specs = pairwise ~ warming)

# OTC decreased light intensity by 5.8%


# Plots

cowplot :: plot_grid(temp_graph, temp.avg.graph, ncol = 2, labels = "auto")
cowplot :: plot_grid(air.temp.trend, airtemp.avg.graph, ncol = 2, labels = "auto")

cowplot :: plot_grid(vwc_graph, moist_avg_graph, ncol = 2, labels = "auto")


plot_grid(fig1, fig3, ncol = 2, labels = "auto")
plot_grid(fig2, fig4, ncol = 2, labels = "auto")
plot_grid(fig6, fig7, ncol = 2, labels = "auto")
fig5
fig8
plot_grid(fig5, fig8, ncol = 2, labels = "auto")

