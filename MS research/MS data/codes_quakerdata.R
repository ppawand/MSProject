
##############################
# library
##############################
library(readxl)
library(tidyverse)
library(dplyr)
library(lme4)
library(effects)
library(ggplot2)
library(scales)
library(cowplot)
library(MuMIn)
library(rsq) # calculate R2 for mixed models with inverse Gaussian distribution
##############################
# theme
##############################

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




####################################################
# Summarizing available nutrients
####################################################
# average temperature and moisture of Oct (month of soil sampling)
vol.MT <- daily.TM %>% 
  filter(month == "Oct") %>% 
  group_by(block, plot, subplot, trt) %>% 
  summarize(mean.T = mean(mean.temp),
            mean.M = mean(mean.moist))

nutrients <- read_excel("~/Desktop/MS data/nutrients.xlsx")
nutrients[, 2:6] <- lapply(nutrients[, 2:6], factor)
nutrients$trt <- factor(nutrients$trt, levels = c("C", "W", "R", "WR"))
str(nutrients)

average.nutrients <- nutrients %>%
  group_by(trt) %>%
  summarize(mnitN = mean(nitN),
            sdnitN = sd(nitN),
            mNH4 = mean(NH4),
            sdNH4 = sd(NH4),
            mP = mean(P),
            sdP = sd(P),
            mK = mean(K),
            sdK = sd(K),
            mMg = mean(Mg),
            sdMg = sd(Mg),
            mCa = mean (Ca),
            sdCa = sd(Ca),
            mS =mean(S),
            sdS = sd(S),
            mB = mean(B),
            sdB = sd(B),
            mZn = mean(Zn),
            sdZn = sd(Zn),
            mFe = mean(Fe),
            sdFe = sd(Fe),
            mCU = mean(Cu),
            sdCu = sd(Cu),
            mMn = mean(Mn),
            sdMn = sd(Mn),
            mNa = mean(Na),
            sdNa = sd(Na),
            mepH = mean(pH),
            sdpH = sd(pH),
            meanCEC= mean(CEC),
            sdCEC = sd(CEC))



####################################
# SOM
##################################
nutrients <- read_excel("~/Desktop/MS data/nutrients.xlsx")
nutrients$warming <-ifelse(nutrients$trt %in% c("C", "R"), "Control", "OTC")
nutrients$residue <-ifelse(nutrients$trt %in% c("C", "W"), "noResidue", "Residue")
nutrients[, c(2:5, 25, 26)] <- lapply(nutrients[, c(2:5, 25, 26)], factor)
nutrients$trt <- factor(nutrients$trt, levels = c("C", "W", "R", "WR"))

nutrients$log.om <- log10(nutrients$OM)

###############################

## soil organic carbon model log transformed

bestFitModel2(x = nutrients, dependent_var = nutrients$log.om)

log_model_c <- lmer(log.om ~ warming + residue + irrigation + (1|block) , 
                REML = F, data = nutrients)

Anova(log_model_c)
anova(log_model_c)
plot(log_model_c)
qqnorm(residuals(log_model_c))
hist(residuals(log_model_c))
r.squaredGLMM(log_model_c)

nutrients %>% group_by(irrigation) %>%
  summarise(om <- mean(OM))

# OM model with continuous predictors (log transformed)

log_model_c1 <- lmer(log.om ~ vol.water + temp + (1|block) , 
                REML = F, data = nutrients)
anova(log_model_c1)
Anova(log_model_c1)
plot(log_model_c1)
qqnorm(residuals(log_model_c1))
hist(residuals(log_model_c1))
r.squaredGLMM(log_model_c1)

####################################################
# glmer with no transformation
bestGeneralizedModel(x = nutrients, dependent_var = nutrients$OM)

gmodel_c <- glmer(OM ~ warming + residue + irrigation + (1|block) , 
                family = Gamma (link = "log") , data = nutrients)

anova(gmodel_c)
Anova(gmodel_c)
plot(gmodel_c)
r.squaredGLMM(gmodel_c)

# continuous predictors only
gmodel_c1 <- glmer(OM ~ vol.water + temp + (1|block) , 
                 family = Gamma(link = "log"), data = nutrients)
anova(gmodel_c1)
Anova(gmodel_c1)
plot(gmodel_c1)
r.squaredGLMM(gmodel_c1)

####################################################
# OM vs temp
om_plot <- nutrients %>% group_by(block, plot, residue, warming, irrigation) %>%
  summarise(mean.om = mean(OM),
            temp = mean(temp))

pred_c <- effect(mod = gmodel_c1, term = "temp") %>% as.data.frame()

ggplot(data = om_plot) +
  geom_point(aes(temp, mean.om, color = irrigation), size = 2) +
  geom_line(data = pred_c, aes(x = temp, y = fit), size = 1) +
  geom_ribbon(data = pred_c, 
              aes(ymin = lower, ymax = upper, x = temp), alpha = 0.1) +
  scale_color_manual(labels = c("Dryland", "Irrigated"),
                    values = c("brown4", "navyblue")) +
  scale_y_continuous(limits = c(0, 1.5)) +
  labs( y = "Organic matter (%)", 
        x = expression("Soil temperature ("*~degree*C*")")) +
  theme


# OM by trt plot


OM.average <- om_plot %>%
  group_by(residue,warming, irrigation) %>%
  summarize(mean.soc = mean(mean.om),
            n = n(),
            se = sd(mean.om)/sqrt(n))

(plot1 <- ggplot(OM.average, aes(warming, mean.soc)) + 
    geom_point(aes(color = residue), size = 3, 
               position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = mean.soc - se,
                      ymax = mean.soc + se, color = residue), width = 0.05,
                  position = position_dodge(width = 0.2)) +
  labs(x="Treatments", y = "Soil Organic Matter (%)") +
  scale_color_manual(labels = c("No Residue", "Residue"),
                     values = c("red3", "blue3")) +
    scale_y_continuous(limits = c(0.2, 0.8)) +
  facet_wrap(~irrigation) +
  theme)

#########################################################
# MBC
mbc <- read_excel("~/Desktop/MS data/MBC/mbc.final.xlsx")
mbc$warming <-ifelse(mbc$trt %in% c("C", "R"), "Control", "OTC")
mbc$residue <-ifelse(mbc$trt %in% c("C", "W"), "noResidue", "Residue")
mbc$log.mbc <- log10(mbc$mbc)
mbc$P <- mbc$P * 1.121 # converiting into kg/ha
mbc[, c(2:5, 15, 16)] <- lapply(mbc[, c(2:5, 15, 16)], factor)
mbc$trt <- factor(mbc$trt, levels = c("C", "W", "R", "WR"))


####################################################
# MBC model(log transformed)

bestFitModel2(x = mbc, dependent_var = mbc$log.mbc)

model_mbc <- lmer(log.mbc ~ warming + residue + irrigation + warming:residue +
             (1|block), REML = F,
            data = mbc)

anova(model_mbc)
Anova(model_mbc)
plot(model_mbc)
qqnorm(residuals(model_mbc))
hist(residuals(model_mbc))
r.squaredGLMM(model_mbc)

mbc %>% group_by(irrigation) %>%
  summarise(om <- mean(mbc))

# predictor only model (log transformed)

model_mbc1 <- lmer(log.mbc ~ OM + vol.water + temp + pH + nitN + P + NH4 + (1|block),
                  REML = F, data = mbc)
anova(model_mbc1)
Anova(model_mbc1)
plot(model_mbc1)
qqnorm(residuals(model_mbc1))
hist(residuals(model_mbc1))
r.squaredGLMM(model_mbc1)

###############################################

# glmer with no transformation

bestGeneralizedModel(x = mbc, dependent_var = mbc$mbc)


gmodel_mbc <- glmer(mbc ~ warming + residue + irrigation + warming:residue + 
                      (1|block),
                    family = Gamma(link = "log"),
                    data = mbc)
anova(gmodel_mbc)
Anova(gmodel_mbc)
plot(gmodel_mbc)
r.squaredGLMM(gmodel_mbc)

# continuous predictors only

mbc.scaled <- as.data.frame(scale(mbc[, 6:14], center = F))

mbc.scaled$block <-mbc$block

gmodel_mbc1 <- glmer(mbc ~ OM + pH + temp + nitN + P + NH4 + (1|block),
                   family = Gamma(link = "log"), data = mbc.scaled)
anova(gmodel_mbc1)
Anova(gmodel_mbc1)
summary(gmodel_mbc1)
plot(gmodel_mbc1)
r.squaredGLMM(gmodel_mbc1)

###############################################
# mbc by trt

mbc.plot <- mbc %>% group_by(block, plot, residue, warming, irrigation) %>%
  summarise(mean.mbc = mean(mbc),
            nitrate = mean(nitN),
            phosphorus = mean(P))

mbc.average <- mbc.plot %>% 
  group_by(warming, residue, irrigation) %>% 
  summarize(mbc = mean(mean.mbc),
            n = n(),
            se = sd(mean.mbc)/sqrt(n))

(plot3 <- ggplot(mbc.average, aes(warming, mbc)) +
  geom_point(aes(color = residue), size = 3, 
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mbc - se,
                    ymax = mbc + se, color = residue), width = 0.05,
                position = position_dodge(width = 0.2)) +
  labs(x="Treatments", y = ("Microbial Biomass Carbon (mg/kg)" )) +
  scale_color_manual(labels = c("No Residue", "Residue"),
                     values = c("red3", "blue3")) +
    facet_wrap(~ irrigation) +
    theme)

# mbc and nitrate
ggplot(mbc.plot, aes(nitrate, mean.mbc, color = irrigation)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", alpha = 0.1) +
  scale_color_manual(labels = c("Dryland", "Irrigated"),
                     values = c("brown4", "navyblue")) +
  labs(x = "Nitrate (ppm)", y = "Microbial biomass(mg/kg)") +
  scale_x_continuous( limits = c(0, 11.5)) +
  scale_y_continuous(limits =  c(0, 200)) +
  scale_shape_manual(values = c(21,22, 23, 24)) +
  theme

# mbc vs P
ggplot(mbc.plot, aes(phosphorus, mean.mbc, color = irrigation)) + 
  geom_point(size = 2) +
  geom_smooth(method = "lm", alpha = 0.1) +
  scale_color_manual(labels = c("Dryland", "Irrigated"),
                     values = c("brown4", "navyblue")) +
  labs(x = "Phospohrus (kg/ha)", y = "Microbial biomass(mg/kg)") +
  scale_x_continuous(limits = c(20, 60)) +
  scale_y_continuous(limits =  c(0, 200)) +
  scale_shape_manual(values = c(21,22, 23, 24)) +
  theme

###############################################

# soil respiration
respiration <- read_excel("~/Desktop/MS data/Respiration data/respiration.xlsx")
respiration$warming <-ifelse(respiration$trt %in% c("C", "R"), "Control", "OTC")
respiration$residue <-ifelse(respiration$trt %in% c("C", "W"), "noResidue", "Residue")
respiration[, c(2:6, 14,15)] <- lapply(respiration[, c(2:6, 14, 15)], factor)
respiration$trt <- factor(respiration$trt, labels = c("C", "W", "R", "WR"))

###############################################

# respiration model log transformed
bestFitModel(respiration, dependent_var = respiration$log.flux)

resp_model <- lmer(log.flux ~ warming + residue + irrigation +
             residue:irrigation + (1|block) + (1|month), REML = F,
           data = respiration)

anova(resp_model)
Anova(resp_model)
summary(resp_model)
plot(resp_model)
qqnorm(residuals(resp_model))
hist(residuals(resp_model))
r.squaredGLMM(resp_model)

respiration %>% group_by(trt) %>%
  summarise(mean = mean(flux))

# predictor only model log transformed
resp_model1 <- lmer(log.flux ~ temp + moist + OM + pH + (1|block) + (1|month),
                    REML  = F,
                    data = respiration)

anova(resp_model1)
Anova(resp_model1)
summary(resp_model1)
plot(resp_model1)
qqnorm(residuals(resp_model1))
hist(residuals(resp_model1))
r.squaredGLMM(resp_model1)

###############################################

# glmer with no transformation

bestGeneralizedModel2(respiration, dependent_var = respiration$flux)

gresp_model <- glmer(flux ~ warming + residue + irrigation +
                     residue:irrigation + (1|block) + (1|month),
                     family = Gamma(link = "log"),
                     control = glmerControl(optimizer = "bobyqa"),
                   data = respiration)

anova(gresp_model)
Anova(gresp_model)
plot(gresp_model)
r.squaredGLMM(gresp_model)

# continuous predictors only
gresp_model1 <- glmer(flux ~ temp + moist + OM + (1|block) + (1|month), 
                      family = Gamma(link = "log"),
                      control = glmerControl(optimizer = "bobyqa"),
                      data = respiration)

anova(gresp_model1)
Anova(gresp_model1)
plot(gresp_model1)
r.squaredGLMM(gresp_model1)

###############################################
# respiration by trt
resp_plot <- respiration %>% 
  group_by(block, plot, residue, warming, irrigation) %>%
  summarise(mean.resp = mean(flux))

resp.average <- resp_plot %>%
  group_by(warming, residue, irrigation) %>% 
  summarize(mean = mean(mean.resp),
            n = n(),
            se = sd(mean.resp)/sqrt(n))

(plot4 <- ggplot(resp.average, aes(warming, mean)) +
  geom_point(aes(color = residue), size = 3,
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se, color = residue), width = 0.05,
                position = position_dodge(width = 0.2)) +
  labs(x="Treatments", y = expression("Soil Flux ("*mu~"mol" ~CO[2]~ m^-2~s^-1*")")) +
  scale_color_manual(labels = c("No residue", "Residue"),
                     values = c("red3", "blue3")) +
    facet_wrap(~irrigation) +
  theme)




####################################################
# yield
yield_data <-  read_excel("~/Desktop/MS data/yield/yield.xlsx")
yield_data$warming <-ifelse(yield_data$trt %in% c("C", "R"), "Control", "OTC")
yield_data$residue <-ifelse(yield_data$trt %in% c("C", "W"), "noResidue", "Residue")
yield_data[, c(2:5,14,15)] <- lapply(yield_data[, c(2:5,14,15)], factor)
yield_data$trt <- factor(yield_data$trt, labels = c("C", "W", "R", "WR"))

####################################################
# above ground biomass model
bestFitModel2(yield_data, dependent_var = yield_data$above.ground.biomass)

aboveground_model <- lmer(above.ground.biomass ~ warming + residue + irrigation +
                            (1|block), REML = F,
                          data = yield_data)

anova(aboveground_model)
Anova(aboveground_model)

summary(aboveground_model)
plot(aboveground_model)
qqnorm(residuals(aboveground_model))
hist(residuals(aboveground_model))
r.squaredGLMM(aboveground_model)

yield_data %>% group_by(irrigation) %>%
  summarise(mean = mean(above.ground.biomass))

# model with continuous predictors 

aboveground_model1 <- lmer(above.ground.biomass ~ nitN + P + NH4 +
                            (1|block), REML = F,
                          data = yield_data)

anova(aboveground_model1)
Anova(aboveground_model1)
summary(aboveground_model1)
r.squaredGLMM(aboveground_model1)

####################################################

# aboveground biomass by trt
agb <- yield_data %>%
  group_by(warming,residue, irrigation) %>% 
  summarize(mean = mean(above.ground.biomass),
            n = n(),
            se = sd(above.ground.biomass)/sqrt(n))

(plot5 <- ggplot(agb, aes(warming, mean)) +
  geom_point(aes(color = residue), size = 3, 
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se, color = residue), width = 0.05,
                position = position_dodge(width = 0.2)) +
  labs(x="Treatments", y = "Above ground biomass (kg/ha)") +
  scale_color_manual(labels = c("No Residue", "Residue"),
                     values = c("red3", "blue3")) +
    scale_y_continuous(limits = c(1000, 5000)) +
  facet_wrap(~irrigation) +
  theme)
  
####################################################

# seed.cotton yield model
bestFitModel2(yield_data, dependent_var = yield_data$seed.cotton)

yield_model <- lmer(seed.cotton ~ warming + residue + irrigation +
             warming:irrigation + (1|block), REML = F,
           data = yield_data)

anova(yield_model)
Anova(yield_model)

summary(yield_model)
plot(yield_model)
qqnorm(residuals(yield_model))
hist(residuals(yield_model))
r.squaredGLMM(yield_model)

# model with continuous predictors only
yield_model1 <- lmer(seed.cotton ~ nitN + P + NH4 + (1|block),
                   REML = F,
                   data = yield_data)
anova(yield_model1)
Anova(yield_model1)

plot(yield_model1)
qqnorm(residuals(yield_model1))
hist(residuals(yield_model1))
r.squaredGLMM(yield_model1)

yield_data %>% group_by(irrigation) %>%
  summarize(mean = mean(seed.cotton))

####################################################

# seed cotton yield  by trt
seed.cotton <- yield_data %>%
  group_by(warming, residue, irrigation) %>% 
  summarize(mean = mean(seed.cotton),
            n = n(),
            se = sd(seed.cotton)/sqrt(n))

(plot6 <- ggplot(seed.cotton, aes(warming, mean)) +
  geom_point(aes(color = residue), size = 3,
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se, color = residue), width = 0.05,
                position = position_dodge(width = 0.2)) +
  labs(x="Treatments", y = "Seed Cotton Yield (kg/ha)") +
  scale_color_manual(labels = c("No Residue", "Residue"),
                     values = c("red3", "blue3")) +
    scale_y_continuous(limits = c(1000, 3000)) +
    facet_wrap(~irrigation) +
      theme)


####################################################
# # below ground biomass model

bestFitModel2(yield_data, dependent_var = yield_data$rootbiomass)

bgb_model <- lmer(rootbiomass ~ warming + residue + irrigation +
                    (1|block), REML = F,
                  data = yield_data)

anova(bgb_model)
Anova(bgb_model)

summary(bgb_model)
plot(bgb_model)
qqnorm(residuals(bgb_model))
hist(residuals(bgb_model))
r.squaredGLMM(bgb_model)

# model with continuous predictors only
bgb_model1 <- lmer(rootbiomass ~ nitN + P + NH4 + (1|block),
                   REML = F,
                   data = yield_data)
anova(bgb_model1)
Anova(bgb_model1)
summary(bgb_model1)
plot(bgb_model1)
qqnorm(residuals(bgb_model1))
hist(residuals(bgb_model1))
r.squaredGLMM(bgb_model1)


yield_data %>% group_by(trt) %>%
  summarise(mean = mean(rootbiomass))

yield_data %>% group_by(irrigation) %>%
  summarise(mean = mean(rootbiomass))

####################################################

# below ground biomass by trt

bgb <- yield_data %>% group_by(warming, residue, irrigation) %>%
  summarize (mean = mean(rootbiomass),
             n = n(),
             se = sd(rootbiomass)/ sqrt(n))

(plot7 <- ggplot(bgb, aes(warming, mean))  +
  geom_point(aes(color = residue), size = 3, 
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se, color = residue), width = 0.05,
                position = position_dodge(width = 0.2)) +
  scale_color_manual(labels = c("No Residue", "Residue"),
                     values = c("red3", "blue3")) +
  labs(x="Treatments", y = "Root biomass (g)") +
  facet_grid( ~ irrigation) +
  theme)



####################################################

# bolls per plant model

bestFitModel2(yield_data, dependent_var = yield_data$boll.per.plant)

boll_model <- lmer(boll.per.plant ~ warming + residue + irrigation +
                    (1|block), REML = F,
                  data = yield_data)

anova(boll_model)
Anova(boll_model)

summary(boll_model)
plot(boll_model)
qqnorm(residuals(boll_model))
hist(residuals(boll_model))
r.squaredGLMM(boll_model)

####################################################
# bolls per plants by trt
bolls <- yield_data %>%
  group_by(residue, warming, irrigation) %>% 
  summarize(mean = mean(boll.per.plant),
            n = n(),
            se = sd(boll.per.plant)/n)

(plot8 <- ggplot(bolls, aes(warming, mean)) +
  geom_point(aes(color = residue), size = 3,
             position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean - se,
                    ymax = mean + se, color = residue), width = 0.05,
                position = position_dodge(width = 0.2)) +
  labs(x="Treatments", y = "Number of bolls per plant") +
  scale_color_manual(labels = c("No residue", "Residue"),
                     values = c("red3", "blue3")) +
    facet_wrap(~ irrigation) +
  theme)

































































