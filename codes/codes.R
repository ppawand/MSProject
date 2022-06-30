
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
library(rsq) # to calculate R2 for mixed model with inverse Gaussian distribution

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


setwd("~/Library/CloudStorage/OneDrive-TexasTechUniversity/data analysis")

####################################################
# Summarizing available nutrients
####################################################

nutrients <- read_excel("nutrients.xlsx")
nutrients$warming <-ifelse(nutrients$trt %in% c("C", "R"), "Ambient", "Warmed")
nutrients$residue <-ifelse(nutrients$trt %in% c("C", "W"), "noResidue", "Residue")

# Reordering column names, so warming and residue are earlier
nutrients <- nutrients %>% 
  relocate(warming, residue, .before = irrigation)
nutrients[, 2:7] <- lapply(nutrients[, 2:7], factor)

# averaging nutrients by treatments

average.nutrients <- nutrients %>%
  group_by(warming, residue, irrigation) %>% 
  summarize(mnitN = mean(nitN),
            sdnitN = sd(nitN),
            mNH4 = mean(NH4),
            sdNH4 = sd(NH4),
            mOM = mean(OM), 
            sdOM = sd(OM),
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

nutrients$log.om <- log10(nutrients$OM)

## soil organic carbon model log transformed

bestFitModel2(x = nutrients, dependent_var = nutrients$log.om)

# refitting best model
log_model_c <- lmer(log.om ~ warming + residue + irrigation + (1|block),
                    data = nutrients)

Anova(log_model_c)
anova(log_model_c)
plot(log_model_c)
qqnorm(residuals(log_model_c))
hist(residuals(log_model_c))
r.squaredGLMM(log_model_c)


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

# refitting with different distribution
gmodel_c1 <- glmer(OM ~ warming + residue + irrigation + (1|block) , 
                family = Gamma (link = "inverse") , data = nutrients)
gmodel_c2 <- glmer(OM ~ warming + residue + irrigation + (1|block) , 
                  family = poisson (link = "log") , data = nutrients)
gmodel_c3 <- lmer(OM ~ warming + residue + irrigation + (1|block), 
                  data = nutrients)
gmodel_c4 <- glmer(OM ~ warming + residue + irrigation + (1|block) , 
                   family = inverse.gaussian (link ="log") , data = nutrients) # NvG: Note there are different link functions you can use, not just log!

AIC(gmodel_c1, gmodel_c2, gmodel_c3, gmodel_c4)

# Winner is the last one! 

Anova(gmodel_c4)
plot(gmodel_c4)
qqnorm(residuals(gmodel_c4))

# continuous predictors only 

gmodel_carbon1 <- glmer(OM ~ vol.water + temp + (1|block) , 
                 family = Gamma(link = "log"), data = nutrients)
# let's try the inverse gaussian again
gmodel_carbon2 <- glmer(OM ~ vol.water + temp + (1|block) , 
                   family = inverse.gaussian (link ="log"), data = nutrients)


AIC(gmodel_carbon1, gmodel_carbon2)

# Again - the one with the inverse.gaussian is better (gmodel_c2)

anova(gmodel_carbon2)
Anova(gmodel_carbon2)
plot(gmodel_carbon2)


####################################################
# OM vs temp
om_plot <- nutrients %>% group_by(block, plot, residue, warming, irrigation) %>%
  summarise(mean.om = mean(OM),
            temp = mean(temp))

pred_c <- effect(mod = gmodel_c1, term = "temp") %>% as.data.frame()

(om_temp_graph <- ggplot(data = om_plot) +
  geom_point(aes(temp, mean.om, color = irrigation), size = 2) +
  geom_line(data = pred_c, aes(x = temp, y = fit), size = 1) +
  geom_ribbon(data = pred_c, 
              aes(ymin = lower, ymax = upper, x = temp), alpha = 0.1) +
  scale_color_manual(labels = c("Dryland", "Irrigated"),
                     values = c("brown4", "navyblue")) +
  scale_y_continuous(limits = c(0, 1.5)) +
  labs( y = "Organic matter (%)", 
        x = expression("Soil temperature ("*~degree*C*")")) +
  theme)

# OM by trt plot

OM.average <- om_plot %>%
  group_by(residue,warming, irrigation) %>%
  summarize(mean.soc = mean(mean.om),
            n = n(),
            se = sd(mean.om)/sqrt(n))

(om_graph <- ggplot(OM.average, aes(warming, mean.soc)) + 
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

mbc <- read_excel("mbc.final.xlsx")
mbc$warming <-ifelse(mbc$trt %in% c("C", "R"), "Ambient", "Warmed")
mbc$residue <-ifelse(mbc$trt %in% c("C", "W"), "noResidue", "Residue")
mbc$log.mbc <- log10(mbc$mbc)
mbc$P <- mbc$P * 1.121 # converiting into kg/ha
mbc[, c(2:5, 15, 16)] <- lapply(mbc[, c(2:5, 15, 16)], factor)
mbc$trt <- factor(mbc$trt, levels = c("C", "W", "R", "WR"))


####################################################
# MBC model(log transformed)
#best mbc model
bestFitModel2(x = mbc, dependent_var = mbc$log.mbc)
# refitting best model
model_mbc <- lmer(log.mbc ~ warming + residue + irrigation + warming:residue +
             (1|block), REML = F,
            data = mbc)
# singular fit (overfitting)
Anova(model_mbc)
plot(model_mbc)
qqnorm(residuals(model_mbc))
hist(residuals(model_mbc)) # still not normally distributed
r.squaredGLMM(model_mbc)


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

# refitting with different distribution
gmodel_mbc1 <- glmer(mbc ~ warming + residue + irrigation + warming:residue + 
                      (1|block),
                    family = Gamma(link = "log"),
                    data = mbc)
gmodel_mbc2 <- glmer(mbc ~ warming + residue + irrigation + warming:residue + 
                       (1|block),
                     family = Gamma(link = "inverse"),
                     data = mbc)

gmodel_mbc3 <- glmer(mbc ~ warming + residue + irrigation + warming:residue + 
                       (1|block),
                     family = inverse.gaussian(link = "log"),
                     data = mbc)

AIC(gmodel_mbc1, gmodel_mbc2, gmodel_mbc3)

# first model is winner.

anova(gmodel_mbc1)
Anova(gmodel_mbc1)
plot(gmodel_mbc1)
r.squaredGLMM(gmodel_mbc1)

# continuous predictors only

mbc.scaled <- as.data.frame(scale(mbc[, 6:14], center = F))

mbc.scaled$block <-mbc$block

gmodel_mbiomass1 <- glmer(mbc ~ OM + pH + temp + nitN + P + NH4 + (1|block),
                   family = Gamma(link = "log"), data = mbc.scaled)

gmodel_mbiomass2 <- glmer(mbc ~ OM + pH + temp + nitN + P + NH4 + (1|block),
                          family = inverse.gaussian(link = "log"), 
                          data = mbc.scaled)
AIC(gmodel_mbiomass1, gmodel_mbiomass2)

# first is winnner

Anova(gmodel_mbiomass1)
plot(gmodel_mbiomass1)
r.squaredGLMM(gmodel_mbiomass1)

###############################################
# mbc by trt

mbc.plot <- mbc %>% group_by(block, plot, residue, warming, irrigation) %>%
  summarise(mean.mbc = mean(mbc))

mbc.average <- mbc.plot %>% 
  group_by(warming, residue, irrigation) %>% 
  summarize(mbc = mean(mean.mbc),
            n = n(),
            se = sd(mean.mbc)/sqrt(n))

(mbc_graph <- ggplot(mbc.average, aes(warming, mbc)) +
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
(mbc_N_graph <- ggplot(mbc, aes(nitN, mbc, color = irrigation)) + 
   geom_point(aes(shape = trt, fill = irrigation)) +
  geom_smooth(method = "lm", se = F) +
  scale_color_manual(labels = c("Dryland", "Irrigated"),
                     values = c("red3", "blue3")) +
  scale_fill_manual(labels = c("Dryland", "Irrigated"),
                     values = c("red3", "blue3")) +
  theme + 
  labs(x = "Nitrate (ppm)", y = "Microbial biomass(mg/kg)") +
  scale_x_continuous( limits = c(0, 12)) +
  scale_y_continuous(limits =  c(0, 200)) +
  scale_shape_manual(values = c(21,22, 23, 24)))

# mbc vs P
(mbc_p_graph <- ggplot(mbc, aes(P, mbc, color = irrigation)) + 
  geom_point(aes(shape = trt, fill = irrigation)) +
  geom_smooth(method = "lm", se = F) +
  scale_color_manual(labels = c("Dryland", "Irrigated"),
                     values = c("red3", "blue3")) +
  scale_fill_manual(labels = c("Dryland", "Irrigated"),
                    values = c("red3", "blue3")) +
  theme + 
  labs(x = "Phospohrus (kg/ha)", y = "Microbial biomass(mg/kg)") +
  scale_x_continuous(limits = c(20, 60)) +
  scale_y_continuous(limits =  c(0, 200)) +
  scale_shape_manual(values = c(21,22, 23, 24)))

###############################################

# soil respiration
respiration <- read_excel("respiration.xlsx")
respiration$warming <-ifelse(respiration$trt %in% c("C", "R"), "Ambient", "Warmed")
respiration$residue <-ifelse(respiration$trt %in% c("C", "W"), "noResidue", "Residue")
respiration[, c(2:6, 14,15)] <- lapply(respiration[, c(2:6, 14, 15)], factor)
respiration$trt <- factor(respiration$trt, labels = c("C", "W", "R", "WR"))

###############################################

# respiration model log transformed
bestFitModel(respiration, dependent_var = respiration$log.flux)

# refitting best model
resp_model <- lmer(log.flux ~ warming + residue + irrigation +
             residue:irrigation + (1|block) + (1|month), REML = F,
           data = respiration)
# singular fit (overfitting)

Anova(resp_model)
summary(resp_model)
plot(resp_model)
qqnorm(residuals(resp_model))
hist(residuals(resp_model))
r.squaredGLMM(resp_model)


# predictor only model log transformed
resp_model1 <- lmer(log.flux ~ temp + moist + OM + pH + (1|block) + (1|month),
                    REML  = F,
                    data = respiration)
# singular fit (overfitting)


Anova(resp_model1)
summary(resp_model1)
plot(resp_model1)
qqnorm(residuals(resp_model1))
hist(residuals(resp_model1))
r.squaredGLMM(resp_model1)

###############################################
# lets try without transformation

# glmer with no transformation
# best model
bestGeneralizedModel2(respiration, dependent_var = respiration$flux)

# refitting with different distribution
gresp_model1 <- glmer(flux ~ warming + residue + irrigation +
                     residue:irrigation + (1|block) + (1|month),
                     family = Gamma(link = "log"),
                     control = glmerControl(optimizer = "bobyqa"),
                   data = respiration)
gresp_model2 <- glmer(flux ~ warming + residue + irrigation +
                       residue:irrigation + (1|block) + (1|month),
                     family = Gamma(link = "inverse"),
                     control = glmerControl(optimizer = "bobyqa"),
                     data = respiration)
gresp_model3 <- glmer(flux ~ warming + residue + irrigation +
                       residue:irrigation + (1|block) + (1|month),
                     family = inverse.gaussian(link = "log"),
                     control = glmerControl(optimizer = "bobyqa"),
                     data = respiration)

AIC(gresp_model1, gresp_model2, gresp_model3)

# model 2 is the winner

Anova(gresp_model2)
plot(gresp_model2)


# continuous predictors only
m_gresp1 <- glmer(flux ~ temp + moist  + (1|block) + (1|month), 
                      family = Gamma(link = "log"),
                      control = glmerControl(optimizer = "bobyqa"),
                      data = respiration)
m_gresp2 <- glmer(flux ~ temp + moist  + (1|block) + (1|month), 
                      family = Gamma(link = "inverse"),
                      control = glmerControl(optimizer = "bobyqa"),
                      data = respiration)
m_gresp3 <- glmer(flux ~ temp + moist  + (1|block) + (1|month), 
                      family = inverse.gaussian(link = "log"),
                      data = respiration)

AIC(m_gresp1, m_gresp2, m_gresp3)
# tird model is the best

Anova(m_gresp3)
plot(m_gresp3)


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

(flux_graph <- ggplot(resp.average, aes(warming, mean)) +
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
yield_data <-  read_excel("yield.xlsx")
yield_data$warming <-ifelse(yield_data$trt %in% c("C", "R"), "Ambient", "Warmed")
yield_data$residue <-ifelse(yield_data$trt %in% c("C", "W"), "noResidue", "Residue")
yield_data[, c(2:5,14,15)] <- lapply(yield_data[, c(2:5,14,15)], factor)
yield_data$trt <- factor(yield_data$trt, labels = c("C", "W", "R", "WR"))

####################################################
# above ground biomass model
bestFitModel2(yield_data, dependent_var = yield_data$above.ground.biomass)

#refitting best model
aboveground_model <- lmer(above.ground.biomass ~ warming + residue + irrigation +
                            (1|block), REML = F,
                          data = yield_data)


Anova(aboveground_model)
summary(aboveground_model)
plot(aboveground_model)
qqnorm(residuals(aboveground_model))
hist(residuals(aboveground_model))
r.squaredGLMM(aboveground_model)


# model with continuous predictors 

aboveground_model1 <- lmer(above.ground.biomass ~ nitN + P + NH4 +
                            (1|block), REML = F,
                          data = yield_data)


Anova(aboveground_model1)
summary(aboveground_model1)
r.squaredGLMM(aboveground_model1)


# aboveground biomass by trt
agb <- yield_data %>%
  group_by(warming,residue, irrigation) %>% 
  summarize(mean = mean(above.ground.biomass),
            n = n(),
            se = sd(above.ground.biomass)/sqrt(n))

(agb_graph <- ggplot(agb, aes(warming, mean)) +
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

# refitting best model
yield_model <- lmer(seed.cotton ~ warming + residue + irrigation +
             warming:irrigation + (1|block), REML = F,
           data = yield_data)


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

Anova(yield_model1)
plot(yield_model1)
qqnorm(residuals(yield_model1))
hist(residuals(yield_model1))
r.squaredGLMM(yield_model1)


####################################################

# seed cotton yield  by trt
seed.cotton <- yield_data %>%
  group_by(warming, residue, irrigation) %>% 
  summarize(mean = mean(seed.cotton),
            n = n(),
            se = sd(seed.cotton)/sqrt(n))

(seedcotton_graph <- ggplot(seed.cotton, aes(warming, mean)) +
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

Anova(bgb_model1)
plot(bgb_model1)
qqnorm(residuals(bgb_model1))
hist(residuals(bgb_model1))
r.squaredGLMM(bgb_model1)


####################################################

# below ground biomass by trt

bgb <- yield_data %>% group_by(warming, residue, irrigation) %>%
  summarize (mean = mean(rootbiomass),
             n = n(),
             se = sd(rootbiomass)/ sqrt(n))

(bgb_graph <- ggplot(bgb, aes(warming, mean))  +
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



# bolls per plants by trt
bolls <- yield_data %>%
  group_by(residue, warming, irrigation) %>% 
  summarize(mean = mean(boll.per.plant),
            n = n(),
            se = sd(boll.per.plant)/n)

(bolls_graph <- ggplot(bolls, aes(warming, mean)) +
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

































































