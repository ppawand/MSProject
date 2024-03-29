####################################################
# library
####################################################
library(readxl)
library(tidyverse)
library(lme4)
library(effects)
library(cowplot)
library(scales)
library(car)
library(MuMIn)

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

temp <- read_excel("~/Desktop/MS data/EM50/TM.final.xlsx", sheet = 2)
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

moist<- read_excel("~/Desktop/MS data/EM50/TM.final.xlsx", sheet = 3)
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
daily.TM$warming <-ifelse(daily.TM$trt %in% c("C", "R"), "Ambient", "Warmed")
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
(fig1 <- ggplot(weekly.TM, aes(week, temp, color = irrigation))+
  geom_line() + scale_color_manual(labels = c("Dryland", "Irrigated"),
                                   values = c("red3", "blue3")) +
  xlab( "Months") +
  ylab(expression("Soil Temperature ("*~degree*C*")")) + 
  facet_wrap(~trt) + theme)

# max and min temp
(max_min <- ggplot(weekly.TM)+
    geom_line(aes(x = week, y = max.tempr), color = "black") +
    geom_line(aes(x = week, y = min.tempr), color = "green") +
    xlab( "Months") +
    ylab(expression("Soil Temperature ("*~degree*C*")")) + 
    facet_grid(irrigation~trt) + theme)

# daily temperature range
(dt <- ggplot(weekly.TM)+
    geom_line(aes(x = week, y = mean.dtr, color = irrigation)) + 
    scale_color_manual(labels = c("Dryland", "Irrigated"),
                                     values = c("red3", "blue3")) +
    xlab( "Months") +
    ylab(expression("Daily temperature range ("*~degree*C*")")) + 
    facet_wrap(~trt) + theme)

####################################################
# temporal variation of moisture (graph)
(fig2 <- ggplot(weekly.TM, aes(week, moist, color = irrigation)) + 
  geom_line() +  scale_color_manual(labels = c("Dryland", "Irrigated"),
                                    values = c("red3", "blue3")) +
  xlab( "Months") +
  ylab(expression("Volumetric water content ("*m^3/m^3*")")) + 
  facet_wrap(~trt) + theme)

################################################################

# temperature model
monthly.avg.TM <- daily.TM %>% 
  group_by(block, plot, warming, residue, irrigation, month) %>% 
  summarize(moisture = mean(mean.moist),
            temp = mean(mean.temp),
            mean.dtr = mean(dtr))


# best temperature model
bestFitModel(x = monthly.avg.TM, dependent_var = monthly.avg.TM$temp)

temp_model <- lmer(temp ~ warming + residue + irrigation +
                           (1|block) + (1|month), REML = F,
                         data = monthly.avg.TM)

anova(temp_model)
Anova(temp_model)
summary(temp_model)
plot(temp_model)
qqnorm(residuals(temp_model))
hist(residuals(temp_model))
r.squaredGLMM(temp_model)

####################################################
# temp by trt (graph)

temp_plot <- daily.TM %>% 
  group_by(block, plot, warming, residue, irrigation) %>%
  summarise(temp = mean(mean.temp))

temp.avg <- temp_plot %>% 
  group_by(warming, residue, irrigation) %>%
  summarize(mean.T = mean(temp),
            n = n(),
            se = sd(temp)/sqrt(n))

(fig3 <- ggplot(temp.avg, aes(warming, mean.T))  +
  geom_point(aes(color = residue), size = 3, 
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean.T - se,
                    ymax = mean.T + se, color = residue), width = 0.05,
                position = position_dodge(width = 0.2)) +
  scale_color_manual(labels = c("No Residue", "Residue"),
                     values = c("red3", "blue3")) +
  labs(x="Treatments", y = expression("Soil Temperature ("*~degree*C*")")) +
  facet_grid( ~ irrigation) +
  theme)

##################################################
# moisture model

bestFitModel(x = monthly.avg.TM, dependent_var = monthly.avg.TM$moisture)

moist_model <- lmer(moisture ~ warming * residue * irrigation + (1|block) + (1|month), REML = F,
           data = monthly.avg.TM)

anova(moist_model)
Anova(moist_model)
summary(moist_model)
plot(moist_model)
qqnorm(residuals(moist_model))
hist(residuals(moist_model))
r.squaredGLMM(moist_model)

####################################################
# moisture by trt

moist_plot <- daily.TM %>% 
  group_by(block, plot, residue, warming, irrigation) %>%
  summarise(moist = mean(mean.moist))

moist.avg <- moist_plot %>% 
  group_by(residue, warming, irrigation) %>%
  summarize(mean.moist = mean(moist),
            n = n(),
            se = sd(moist)/sqrt(n))

(fig4 <- ggplot(moist.avg, aes(warming, mean.moist)) +
  geom_point(aes(color = residue), size = 3, 
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean.moist - se,
                    ymax = mean.moist + se, color = residue), width = 0.05,
                position = position_dodge(width = 0.2)) +
  labs(x="Treatments", y = expression("Volumetric water content ("*m^3/m^3*")")) +
  scale_color_manual(labels = c("No Residue", "Residue"),
                     values = c("red3", "blue3")) +
    facet_grid( ~ irrigation) +
    theme)

####################################################

# soil temp vs moisture
moiturevstemp <- lmer(moisture ~ temp + (1|block) + (1|month), data = monthly.avg.TM, REML = F)
anova(moiturevstemp) # not significant
Anova(moiturevstemp)

#################################################

# air temperature
air.temp <- read_excel("~/Desktop/MS data/hobo/hobo.final.xlsx")

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

air.temp.average$warming <-ifelse(air.temp.average$trt %in% c("C", "R"), "Ambient", "Warmed")

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

(fig6 <- ggplot(weekly.air.temperature, aes(week, air.temp, color = irrigation)) +
  geom_line() + scale_color_manual(labels = c("Dryland", "Irrigated"),
                                   values = c(c("red3", "blue3"))) +
  labs(x = "Months", y = (expression("Air temperature ("*~degree*C*")"))) + 
  facet_wrap(~trt) + theme +
  scale_y_continuous(limits = c(10, 40)))

########################################################
# air temperature model

monthly.air.temp <- air.temp.average %>% 
  group_by(block, warming, residue, irrigation, month) %>% 
  summarize(temp = mean(mean.temp),
            light = (mean(mean.light/1000)))

bestFitModel(x = monthly.air.temp, dependent_var = monthly.air.temp$temp)

air_temp_model <- lmer(temp ~ warming + residue + irrigation +
       warming:irrigation + light + (1|block) + (1|month), REML = F,
     data = monthly.air.temp)

anova(air_temp_model)
Anova(air_temp_model)
summary(air_temp_model)
plot(air_temp_model)
qqnorm(residuals(air_temp_model))
hist(residuals(air_temp_model))
r.squaredGLMM(air_temp_model)

####################################################
# air temp by trt (graph)

air_temp_plot <- air.temp.average %>% 
  group_by(block, warming, residue, irrigation) %>%
  summarise(temp = mean(mean.temp))

airtemp.avg <- air_temp_plot %>% 
  group_by(warming, residue, irrigation) %>%
  summarize(mean.T = mean(temp),
            n = n(),
            se = sd(temp)/sqrt(n))

(fig7 <- ggplot(airtemp.avg, aes(warming, mean.T)) +
  geom_point(aes(color = residue), size = 3, 
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean.T - se,
                    ymax = mean.T + se, color = residue), width = 0.05,
                position = position_dodge(width = 0.2)) +
  labs(x="Treatments", y = expression("Air temperature ("*~degree*C*")")) +
  scale_color_manual(labels = c("No Residue", "Residue"),
                     values = c("red3", "blue3")) +
  facet_wrap(~irrigation) + 
  theme)
  
####################################################
# air temp vs light

avgair <- air.temp.average %>% group_by(block, trt, irrigation) %>%
  summarise(mean.T = mean(mean.temp),
            light = mean(mean.light/1000))

e <- effect(term = "light", mod = air_temp_model) %>% as.data.frame()

(fig8 <- ggplot(avgair, aes(light, mean.T)) +
  geom_point(aes(shape = trt, fill = irrigation)) + 
  geom_line(data = e, aes(light, fit), color = "black", size = 1, linetype = 1) +
  theme +  scale_fill_manual(labels = c("Dryland", "Irrigated"),
                             values = c("red3", "blue3")) +
  labs( x= "light (Lux)", y = expression("Air temperature ("*~degree*C*")")) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_x_continuous(labels = unit_format(unit = "e+03")) +
  scale_y_continuous(limits = c(20, 40)))

####################################################
# light between OTC and control


lightOTC <- monthly.air.temp %>% group_by(block, warming) %>%
  summarise(mean.light = mean(light))


light_model <- lmer(mean.light ~ warming + (1|block),
                    data = lightOTC, REML = F)

anova(light_model)
Anova(light_model)

summary(light_model)
# OTC decreased light intensity by 5.8%


# Plots
plot_grid(fig1, fig3, ncol = 2, labels = "auto")
plot_grid(fig2, fig4, ncol = 2, labels = "auto")
plot_grid(fig6, fig7, ncol = 2, labels = "auto")
fig5
fig8
plot_grid(fig5, fig8, ncol = 2, labels = "auto")

