
##############################
# library
##############################
library(readxl)
library(tidyr)
library(dplyr)
library(lme4)
library(car)
library(effects)
library(ggplot2)
library(scales)
library(multcomp)
library(emmeans)

##############################
# theme
##############################

theme <- theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill = NA, linewidth = 1),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.border = element_rect(linewidth = 1))

####################################################
# Summarizing available nutrients
####################################################

nutrients <- read_excel("~/Desktop/MSProject/MS_Paper/data/nutrients.xlsx")
nutrients$warming <-ifelse(nutrients$trt %in% c("C", "R"), "Ambient", "OTC")
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
# summarizing OM by plots. ( averaging two pseudo replicates from within the same 
# to to get one value per plot)

om.plot <- nutrients %>% 
  group_by(wp, plot, warming, residue, irrigation) %>% 
  summarise(OM = mean(OM),
            vol.water = mean(vol.water),
            temp = mean(temp))

# OM model
om_model <- lmer(OM ~ warming * residue * irrigation + (1|wp),
                 data = om.plot)

Anova(om_model)
plot(om_model)
qqnorm(residuals(om_model))
hist(residuals(om_model))
r.squaredGLMM(om_model)

# pairwise comparison
# interaction was not significant. we ran pairwise comparison for the
# significant main effects.
om.emm <- emmeans(om_model, ~ irrigation)
contrast(om.emm, "consec", simple = "each", combine = TRUE)


# OM model with other numeric predictors 

om_model1<- lmer(OM ~ vol.water + temp + (1|wp) , data = om.plot)
Anova(om_model1)
plot(om_model1)
qqnorm(residuals(om_model1))
hist(residuals(om_model1))
r.squaredGLMM(om_model1)

# soil temperature has a significant effects on organic matter content.

####################################################
# OM vs temp


pred_c <- effect(mod = om_model1, term = "temp") %>% as.data.frame()

(om_temp_graph <- ggplot(data = om.plot) +
  geom_jitter(aes(temp, OM, fill = irrigation, color = irrigation,
              shape = residue), size = 5, alpha = 0.7) +
  geom_line(data = pred_c, aes(x = temp, y = fit), size = 1) +
  geom_ribbon(data = pred_c, 
              aes(ymin = lower, ymax = upper, x = temp), alpha = 0.2) +
  scale_fill_manual(labels = c("Dryland", "Irrigated"),
                     values = c("brown4", "navyblue")) +
  scale_color_manual(labels = c("Dryland", "Irrigated"),
                      values = c("brown4", "navyblue")) + 
  scale_y_continuous(limits = c(0, 1.2)) +
  scale_shape_manual(labels = c("No residue", "Residue"),
                       values = c(21,22)) +
  labs( y = "Soil organic matter (%)", 
        x = expression("Soil temperature ("*~degree*C*")")) +
  theme)



# OM by trt plot

OM.average <- as.data.frame(emmeans(om_model, ~warming * residue * irrigation))

(om_graph <- ggplot(om.plot, aes(x = warming, color = residue))+
    geom_boxplot(aes(y = OM, fill = residue), alpha = 0.01, color = "gray") +
    geom_point(aes(y = OM, color  = residue), alpha = 0.2, size = 2,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3)) +
    geom_point(data = OM.average,
               aes(x = warming, y = emmean, group = residue),
               size = 5, shape = 15,
               position = position_dodge(width = 0.75)) +
    geom_errorbar(data = OM.average,
                  aes(x = warming, ymin = emmean - SE,
                      ymax = emmean + SE, group = residue), width = 0.2,
                  position = position_dodge(width = 0.75)) +
  labs(x="Treatments", y = "Soil organic matter (%)") +
  facet_wrap(~irrigation) +
  scale_color_manual(labels = c("No residue", "Residue"),
                     values = c("#000000", "#009E73")) +
  scale_fill_manual(labels = c("No residue", "Residue"),
                       values = c("#000000", "#009E73")) +
  theme)


# om by irrigation

OM.irrigation <- as.data.frame(emmeans(om_model, ~irrigation))

(om_irr_plot <- ggplot(om.plot, aes(x = irrigation))+
    geom_boxplot(aes(y = OM), color = "gray") +
    geom_jitter(aes(y = OM), alpha = 0.2, size = 2, width = 0.1) +
    geom_point(data = OM.irrigation,
               aes(x = irrigation, y = emmean),
               size = 5, shape = 15) +
    geom_errorbar(data = OM.irrigation,
                  aes(x = irrigation, ymin = emmean - SE,
                      ymax = emmean + SE), linewidth = 1, width = 0.1) +
    labs(x="Treatments", y = "Soil organic matter (%)") +
    annotate(geom="text", x= 1, y= 0.9,
             label="Irrigation: P = 0.005", size = 4.5) +
    theme)



#########################################################
# MBC

mbc <- read_excel("~/Desktop/MSProject/MS_Paper/data/mbc.final.xlsx")
mbc$warming <-ifelse(mbc$trt %in% c("C", "R"), "Ambient", "OTC")
mbc$residue <-ifelse(mbc$trt %in% c("C", "W"), "noResidue", "Residue")

# Reordering column names, so warming and residue are earlier
mbc <- mbc %>% 
  relocate(warming, residue, .before = irrigation)

mbc[, 2:7] <- lapply(mbc[, 2:7], factor)
mbc$P <- mbc$P * 1.121 # converting into kg/ha

####################################################
# summarizing MBC by plots.(averaging two pseudo replicates from within the same 
# to to get one value per plot)

mbc.plot <- mbc %>% 
  group_by(wp, plot, trt, warming, residue, irrigation) %>% 
  summarise(mbc = mean(mbc),
            vol.water = mean(vol.water),
            temp = mean(temp),
            pH = mean(pH),
            om = mean(OM),
            P = mean(P),
            nitrate = mean(nitN))
            

# OM model
mbc_model<- lmer(mbc~warming * residue * irrigation + (1|wp),
                 data = mbc.plot,
                 control = lmerControl(check.conv.singular = 
                                                 .makeCC(action = "ignore", 
                                                         tol = 1e-4)))

Anova(mbc_model)
plot(mbc_model)
qqnorm(residuals(mbc_model))
hist(residuals(mbc_model))
# residuals are not normally distributed.
# refitting the generalized liner mixed effects model with different distribution
# log linked gamma distribution
mbc_model1 <- glmer(mbc ~ warming * residue * irrigation + (1|wp),
                    family = Gamma(link = "log"),
                    data = mbc.plot)
# inverse gaussian distribution
mbc_model2 <- glmer(mbc ~ warming * residue * irrigation + (1|wp),
                    family = inverse.gaussian(link = "log"),
                    data = mbc.plot)

#comparing two distribution
AIC( mbc_model1, mbc_model2)

# Gamma distribution looks better
Anova(mbc_model1)
plot(mbc_model1)
r.squaredGLMM(mbc_model1)

# pairwise comparison
# we ran pairwise comparison for the significant interaction and
# significant main effects.
mbc.emm <- emmeans(mbc_model1, ~ irrigation)
contrast(mbc.emm, "consec", simple = "each", combine = TRUE)

mbc.emm.int <- emmeans(mbc_model1, ~ residue:warming)
contrast(mbc.emm.int, "consec", simple = "each", combine = TRUE)


# continuous predictors only

mbc_numeric <- lmer(mbc ~ om + pH + temp + vol.water + nitrate + P + (1|wp),
                    data = mbc.plot,
                    control=lmerControl(check.conv.singular =
                                          .makeCC(action = "ignore",
                                                  tol = 1e-4)))
Anova(mbc_numeric)
plot(mbc_numeric)
qqnorm(residuals(mbc_numeric))
hist(residuals(mbc_numeric))
r.squaredGLMM(mbc_numeric)

###############################################
# mbc by trt

mbc.average <- 
  as.data.frame(emmeans(mbc_model1 , ~warming * residue * irrigation, type = "response"))

(mbc_graph <- ggplot(mbc.plot, aes(x = warming, color = residue)) +
    geom_boxplot(aes(y = mbc, fill = residue), alpha = 0.01, color = "gray") +
    geom_point(aes(y = mbc, color  = residue), alpha = 0.2, size = 2,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3)) +
    geom_point(data = mbc.average,
               aes(x = warming, y = response, group = residue),
               size = 5, shape = 15,
               position = position_dodge(width = 0.75)) +
    geom_errorbar(data = mbc.average, 
                aes(x = warming, ymin = response- SE,
                    ymax = response + SE, group = residue), width = 0.2,
                position = position_dodge(width = 0.75)) +
    facet_wrap(~ irrigation) + 
    labs(x="Treatments", y = ("Microbial biomass carbon (mg/kg)" )) +
    scale_color_manual(labels = c("No residue", "Residue"),
                     values = c("#000000", "#009E73")) +
    scale_fill_manual(labels = c("No residue", "Residue"),
                       values = c("#000000", "#009E73")) +
    theme +
    theme(strip.background = element_blank(),
          strip.text = element_blank()))

# mbc and om

pred_mbc_om <- effect(mod = mbc_numeric, term = "om") %>% as.data.frame()

(mbc_om_graph <- ggplot(data = mbc.plot) +
    geom_jitter(aes(om, mbc, color = irrigation, fill = irrigation, 
                    shape = residue), size = 5, alpha = 0.7) +
    geom_line(data = pred_mbc_om, aes(x = om, y = fit), size = 1) +
    geom_ribbon(data = pred_mbc_om, 
                aes(ymin = lower, ymax = upper, x = om), alpha = 0.2) +
    scale_color_manual(labels = c("Dryland", "Irrigated"),
                      values = c("brown4", "navyblue")) +
    scale_fill_manual(labels = c("Dryland", "Irrigated"),
                       values = c("brown4", "navyblue")) +
    labs(x = "Soil organic matter (%)", y = "Microbial biomass carbon (mg/kg)") +
    scale_shape_manual(labels = c("No residue", "Residue"),
                       values = c(21,22)) +
  theme)


# mbc and nitrate

pred_mbc_nit <- effect(mod = mbc_numeric, term = "nitrate") %>% as.data.frame()

(mbc_nitrate_graph <- ggplot(data = mbc.plot) +
    geom_jitter(aes(nitrate, mbc, color = irrigation, fill = irrigation, shape = residue), size = 5, alpha = 0.7) +
    geom_line(data = pred_mbc_nit, aes(x = nitrate, y = fit), size = 1) +
    geom_ribbon(data = pred_mbc_nit, 
                aes(ymin = lower, ymax = upper, x = nitrate), alpha = 0.2) +
    scale_color_manual(labels = c("Dryland", "Irrigated"),
                       values = c("brown4", "navyblue")) +
    scale_fill_manual(labels = c("Dryland", "Irrigated"),
                      values = c("brown4", "navyblue")) +
    labs(x = "Nitrate N (ppm)", y = "Microbial biomass carbon (mg/kg)") +
    scale_shape_manual(labels = c("No residue", "Residue"),
                       values = c(21,22)) +
    scale_x_continuous(breaks = seq(2, 12, by = 2), limits = c(2, 12)) +
    theme)


# mbc by otc and residue

mbc.otc.residue <- 
  as.data.frame(emmeans(mbc_model1 , ~ warming * residue, type = "response"))

(mbc_otc_res_plot <- ggplot(mbc.plot, aes(x = warming, color = residue)) +
    geom_boxplot(aes(y = mbc, fill = residue), alpha = 0.01, color = "gray") +
    geom_point(aes(y = mbc, color  = residue), alpha = 0.2, size = 2,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
    geom_point(data = mbc.otc.residue,
               aes(x = warming, y = response, group = residue),
               size = 5, shape = 15,
               position = position_dodge(width = 0.75)) +
    geom_errorbar(data = mbc.otc.residue, 
                  aes(x = warming, ymin = response - SE,
                      ymax = response + SE, group = residue), linewidth = 1,
                  width = 0.1,
                  position = position_dodge(width = 0.75)) +
    labs(x="Treatments", y = ("Microbial biomass carbon (mg/kg)" )) +
    scale_color_manual(labels = c("No residue", "Residue"),
                       values = c("#000000", "#009E73")) +
    scale_fill_manual(labels = c("No residue", "Residue"),
                       values = c("#000000", "#009E73")) +
    annotate(geom="text", x= 0.45, y= 185,
             label="OTC: ns \nResidue: P = 0.001 \nOTC x Residue: P = 0.006",
             size = 4.5, hjust = 0) +
    theme)

# mbc by irrigation

mbc.irrigation <- 
  as.data.frame(emmeans(mbc_model1 , ~irrigation, type = "response"))

(mbc_irri_plot <- ggplot(mbc.plot, aes(x = irrigation)) +
    geom_boxplot(aes(y = mbc), color = "gray") +
    geom_jitter(aes(y = mbc), alpha = 0.2, size = 2, width = 0.1) +
    geom_point(data = mbc.irrigation,
               aes(x = irrigation, y = response),
               size = 5, shape = 15) +
    geom_errorbar(data = mbc.irrigation, 
                  aes(x = irrigation, ymin = response - SE,
                      ymax = response + SE), linewidth = 1, width = 0.1) +
    labs(x="Treatments", y = ("Microbial biomass carbon (mg/kg)" )) +
    annotate(geom="text", x= 1, y= 180,
             label="Irrigation: P = 0.03", size = 4.5) +
    theme)


###############################################

# soil respiration
respiration <- read_excel("~/Desktop/MSProject/MS_Paper/data/respiration.xlsx")
respiration$warming <-ifelse(respiration$trt %in% c("C", "R"), "Ambient", "OTC")
respiration$residue <-ifelse(respiration$trt %in% c("C", "W"), "noResidue", "Residue")

# Reordering column names, so warming and residue are earlier
respiration <- respiration %>% 
  relocate(warming, residue, .before = irrigation)

respiration[, c(2:8)] <- lapply(respiration[, c(2:8)], factor)
respiration$trt <- factor(respiration$trt, labels = c("C", "W", "R", "WR"))

###############################################
# respiration model
resp_model <- lmer(flux ~ warming * residue * irrigation
                   + (1|wp) + (1|month),
                   data  = respiration,
                   control=lmerControl(check.conv.singular =
                                         .makeCC(action = "ignore",
                                                 tol = 1e-4)))

Anova(resp_model)
plot(resp_model)
qqnorm(residuals(resp_model))
hist(residuals(resp_model))

# residuals are not normally distributed.
# refitting the generalized liner mixed effects model with different distribution

gresp_model1 <- glmer(flux ~ warming * residue * irrigation +
                        (1|wp) + (1|month),
                     family = Gamma(link = "log"),
                   data = respiration,
                   control=glmerControl(check.conv.singular =
                                         .makeCC(action = "ignore",
                                                 tol = 1e-4)))
gresp_model2 <- glmer(flux ~ warming * residue * irrigation +
                        (1|wp) + (1|month),
                     family = Gamma(link = "inverse"),
                     data = respiration,
                     control=glmerControl(check.conv.singular =
                                           .makeCC(action = "ignore",
                                                   tol = 1e-4)))
gresp_model3 <- glmer(flux ~ warming * residue * irrigation +
                        (1|wp) + (1|month),
                     family = inverse.gaussian(link = "log"),
                     data = respiration,
                     control=glmerControl(check.conv.singular =
                                           .makeCC(action = "ignore",
                                                   tol = 1e-4)))

AIC(gresp_model1, gresp_model2, gresp_model3)

# model 2 (inverse link gamma distribution looks better
Anova(gresp_model2)
plot(gresp_model2)
r.squaredGLMM(gresp_model2)
# pairwise comparison
# we ran pairwise comparison for the significant interaction and
# significant main effects.
resp.emm <- emmeans(gresp_model2, ~ residue)
contrast(resp.emm, "consec", simple = "each", combine = TRUE)

resp.emm.int <- emmeans(gresp_model2, ~ warming:irrigation)
contrast(resp.emm.int, "consec", simple = "each", combine = TRUE)
 
# continuous predictor only model transformed
resp_model1 <- lmer(flux ~ temp + moist + OM + pH + (1|wp) + (1|month),
                    data = respiration)



Anova(resp_model1)
plot(resp_model1)
qqnorm(residuals(resp_model1))
hist(residuals(resp_model1))


# non normal distribution of residuals

# generalized linear models for continuous predictors only
m_gresp1 <- glmer(flux ~ temp + moist  + OM + pH + (1|wp) + (1|month), 
                      family = Gamma(link = "log"),
                      control = glmerControl(optimizer = "bobyqa"),
                      data = respiration)
m_gresp2 <- glmer(flux ~ temp + moist  +OM + pH + (1|wp) + (1|month), 
                      family = Gamma(link = "inverse"),
                      control = glmerControl(optimizer = "bobyqa"),
                      data = respiration)

AIC(m_gresp1, m_gresp2)

# second model is better

Anova(m_gresp2)
plot(m_gresp2)
r.squaredGLMM(m_gresp2)


###############################################
# respiration by trt
resp.average <- 
  as.data.frame(emmeans(gresp_model2, ~warming * residue * irrigation, type = "response")) 

(flux_graph <- ggplot(respiration, aes(x = warming, color = residue)) +
    geom_boxplot(aes(y = flux, fill = residue), alpha = 0.01, color = "gray") +
    geom_point(aes(y = flux, color  = residue), alpha = 0.2, size = 2,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3)) +
    geom_point(data = resp.average,
               aes(x = warming, y = response, group = residue),
               size = 5, shape = 15,
               position = position_dodge(width = 0.75)) +
    geom_errorbar(data = resp.average, 
                  aes(x = warming, ymin = response - SE,
                      ymax = response + SE, group = residue), width = 0.2,
                  position = position_dodge(width = 0.75)) +
    facet_wrap(~ irrigation) + 
    labs(x="Treatments", 
           y = expression("Soil respiration (" * mu ~ "mol" ~ CO[2] ~ m^{paste("\u2212", 2)} ~ s^{paste("\u2212", 1)} * ")"))+
    scale_color_manual(labels = c("No residue", "Residue"),
                       values = c("#000000", "#009E73")) +
    scale_fill_manual(labels = c("No residue", "Residue"),
                       values = c("#000000", "#009E73")) +
    scale_y_continuous(limits = c(0, 15)) +
    theme +
    theme(strip.background = element_blank(),
          strip.text = element_blank()))

# respiration by residue

resp.residue <- 
  as.data.frame(emmeans(gresp_model2, ~residue, type = "response" ))

(flux_res_plot <- ggplot(respiration, aes(x = residue)) +
    geom_boxplot(aes(y = flux), color = "gray") +
    geom_jitter(aes(y = flux), alpha = 0.2, size = 2, width = 0.1) +
    geom_point(data = resp.residue,
               aes(x = residue, y = response),
               size = 5, shape = 15) +
    geom_errorbar(data = resp.residue, 
                  aes(x = residue, ymin = response - SE,
                      ymax = response + SE), linewidth = 1, width = 0.1) +
    labs(x="Treatments", 
         y = expression("Soil respiration (" * mu ~ "mol" ~ CO[2] ~ m^{paste("\u2212", 2)} ~ s^{paste("\u2212", 1)} * ")"))+
    scale_y_continuous(limits = c(0, 15)) +
    annotate(geom="text", x= 1, y= 14,
             label= "Residue: P < 0.001",
             size = 4.5) +
    scale_x_discrete(labels = c("No residue", "Residue")) +
    theme)

# respiration by otc

resp.otc <- 
  as.data.frame(emmeans(gresp_model2, ~warming, type = "response" ))

(flux_otc_plot <- ggplot(respiration, aes(x = warming)) +
    geom_boxplot(aes(y = flux), color = "gray") +
    geom_jitter(aes(y = flux), alpha = 0.2, size = 2, width = 0.1) +
    geom_point(data = resp.otc,
               aes(x = warming, y = response),
               size = 5, shape = 15) +
    geom_errorbar(data = resp.otc, 
                  aes(x = warming, ymin = response - SE,
                      ymax = response + SE), linewidth = 1, width = 0.1) +
    labs(x="Treatments", 
         y = expression("Soil respiration (" * mu ~ "mol" ~ CO[2] ~ m^{paste("\u2212", 2)} ~ s^{paste("\u2212", 1)} * ")"))+
    scale_y_continuous(limits = c(0, 15)) +
    annotate(geom="text", x= 1, y= 15,
             label= "OTC: P = 0.03",
             size = 4.5) +
    theme)


# respiration by otc and irrigation
resp.otc.irrigation <- 
  as.data.frame(emmeans(gresp_model2, ~warming * irrigation, type = "response" )) 

(flux_otc_irri_plot <- ggplot(respiration, aes(x = warming, color = irrigation)) +
    geom_boxplot(aes(y = flux, fill = irrigation), alpha = 0.01, color = "gray") +
    geom_point(aes(y = flux, color  = irrigation), alpha = 0.2, size = 2,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
    geom_point(data = resp.otc.irrigation,
               aes(x = warming, y = response, group = irrigation),
               size = 5, shape = 15,
               position = position_dodge(width = 0.75)) +
    geom_errorbar(data = resp.otc.irrigation, 
                  aes(x = warming, ymin = response - SE,
                      ymax = response + SE, group = irrigation),linewidth = 1, width = 0.1,
                  position = position_dodge(width = 0.75)) +
    labs(x="Treatments", 
         y = expression("Soil respiration (" * mu ~ "mol" ~ CO[2] ~ m^{paste("\u2212", 2)} ~ s^{paste("\u2212", 1)} * ")"))+
    scale_color_manual(labels = c("Dryland", "Irrigated"),
                       values = c("brown4", "navyblue")) +
    scale_y_continuous(limits = c(0, 10)) +
    annotate(geom="text", x= 0.45, y= 9.4,
             label="OTC: P = 0.03 \nIrrigation: ns \nOTC x Irrigation: P = 0.02",
             size = 4.5, hjust = 0) +
    theme)


fig_4 <- ggarrange(om_irr_plot + rremove("xlab"),
          mbc_irri_plot + rremove("xlab"),
          flux_res_plot + rremove("xlab"),
          flux_otc_plot + rremove("xlab"),
          labels = "auto")

fig_5 <- ggarrange(mbc_otc_res_plot + rremove("xlab"),
          flux_otc_irri_plot + rremove("xlab"),
          common.legend = F,
          labels = "auto")


####################################################
# yield
yield_data <-  read_excel("~/Desktop/MSProject/MS_Paper/data/yield.xlsx")
yield_data$warming <-ifelse(yield_data$trt %in% c("C", "R"), "Ambient", "OTC")
yield_data$residue <-ifelse(yield_data$trt %in% c("C", "W"), "noResidue", "Residue")

# Reordering column names, so warming and residue are earlier
yield_data <- yield_data %>% 
  relocate(warming, residue, .before = irrigation)

yield_data[, c(2:7)] <- lapply(yield_data[, c(2:7)], factor)
yield_data$trt <- factor(yield_data$trt, labels = c("C", "W", "R", "WR"))

####################################################
# above ground biomass model

aboveground_model <- lmer(above.ground.biomass ~ warming * residue * irrigation +
                            (1|wp),
                          data = yield_data)


Anova(aboveground_model)
plot(aboveground_model)
qqnorm(residuals(aboveground_model))
r.squaredGLMM(aboveground_model)
# pairwise comparison for significant main effects & interaction effects

agb.emm <- emmeans(aboveground_model, ~ irrigation)
contrast(agb.emm, "consec", simple = "each", combine = TRUE)

# model with continuous predictors 

aboveground_model1 <- lmer(above.ground.biomass ~ nitN + P + OM +
                            (1|wp),
                          data = yield_data)


Anova(aboveground_model1)
plot(aboveground_model1)
qqnorm(residuals(aboveground_model1))
hist(residuals(aboveground_model1))



# above ground biomass by trt
agb <- as.data.frame(emmeans(aboveground_model, ~warming * residue * irrigation))

(agb_graph <- ggplot(yield_data, aes(x = warming, color = residue)) +
    geom_point(aes(y = above.ground.biomass), alpha = 0.2, size = 2,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
    geom_point(data = agb,
               aes(x = warming, y = emmean, group = residue),
               size = 5, shape = 15,
               position = position_dodge(width = 0.75)) +
    geom_errorbar(data = agb,
                  aes(x = warming, ymin = emmean- SE,
                      ymax = emmean + SE, group = residue), width = 0.2,
                  position = position_dodge(width =0.75)) +
    facet_wrap(~ irrigation) + 
    labs(x="Treatments", y = ("Aboveground biomass (Kg/ha)"))+
    scale_color_manual(labels = c("No residue", "Residue"),
                       values = c("#000000", "#009E73")) +
    theme)

 
# aboveground_biomass and OM

pred_agb_om <- effect(mod =aboveground_model1 , term = "OM") %>% as.data.frame()

(agb_om_graph <- ggplot(data = yield_data) +
    geom_jitter(aes(OM, above.ground.biomass, color = irrigation, fill = irrigation, 
                    shape = residue), size = 5, alpha = 0.7) +
    geom_line(data = pred_agb_om, aes(x = OM, y = fit), size = 1) +
    geom_ribbon(data = pred_agb_om, 
                aes(ymin = lower, ymax = upper, x = OM), alpha = 0.2) +
    scale_color_manual(labels = c("Dryland", "Irrigated"),
                       values = c("brown4", "navyblue")) +
    scale_fill_manual(labels = c("Dryland", "Irrigated"),
                      values = c("brown4", "navyblue")) +
    labs(x = "Soil organic Matter (%)", y = "Above Ground Biomass (Kg/ha)") +
    scale_shape_manual(labels = c("No residue", "Residue"),
                       values = c(21,22)) +
    theme)


# aboveground biomass by irrigation
agb_irr <- as.data.frame(emmeans(aboveground_model, ~irrigation))

(agb_irr_graph <- ggplot(yield_data, aes(x = irrigation)) +
    geom_boxplot(aes(y = above.ground.biomass), color = "gray") +
    geom_jitter(aes(y = above.ground.biomass), alpha = 0.2, size = 2, width = 0.1) +
    geom_point(data = agb_irr,
               aes(x = irrigation, y = emmean),
               size = 5, shape = 15) +
    geom_errorbar(data = agb_irr,
                  aes(x = irrigation, ymin = emmean - SE,
                      ymax = emmean + SE),linewidth = 1, width = 0.1) +
    labs(x="Treatments", y = ("Aboveground biomass (Kg/ha)"))+
    annotate(geom="text", x= 1, y=5200, label="Irrigation: P < 0.001", size = 4.5) +
    theme)


####################################################

# seed.cotton yield model

yield_model <- lmer(seed.cotton ~ warming * residue * irrigation + (1|wp),
           data = yield_data)


Anova(yield_model)
plot(yield_model)
qqnorm(residuals(yield_model))
hist(residuals(yield_model))
r.squaredGLMM(yield_model)

# pairwise comparison for significant main effects & interaction effects

yield.emm.irri <- emmeans(yield_model, ~ irrigation)
contrast(yield.emm.irri, "consec", simple = "each", combine = TRUE)

yield.emm.residue <- emmeans(yield_model, ~ residue)
contrast(yield.emm.residue, "consec", simple = "each", combine = TRUE)

# model with continuous predictors only
yield_model1 <- lmer(seed.cotton ~ nitN + P + OM + (1|wp),
                   data = yield_data)

Anova(yield_model1)
plot(yield_model1)
qqnorm(residuals(yield_model1))
hist(residuals(yield_model1))


####################################################

# seed cotton yield  by trt
seed.cotton <-
  as.data.frame(emmeans(yield_model, ~warming * residue * irrigation))

(yield_graph <- ggplot(yield_data, aes(x = warming, color = residue)) +
    geom_point(aes(y = seed.cotton, color  = residue), alpha = 0.2, size = 2,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
    geom_point(data = seed.cotton,
               aes(x = warming, y = emmean, group = residue),
               size = 5, shape = 15, 
               position = position_dodge(width = 0.75)) +
    geom_errorbar(data = seed.cotton,
                  aes(x = warming, ymin = emmean - SE,
                      ymax = emmean + SE, group = residue), width = 0.2,
                  position = position_dodge(width = 0.75)) +
    facet_wrap(~ irrigation) + 
    labs(x="Treatments", y = ("Seed cotton yield (kg/ha)")) +
    scale_color_manual(labels = c("No residue", "Residue"),
                       values = c("#000000", "#009E73")) +
    theme +
    theme(strip.background = element_blank(),
          strip.text = element_blank()))


# yield by irrigation

yield_irr <- as.data.frame(emmeans(yield_model, ~irrigation))

(yield_irr_graph <- ggplot(yield_data, aes(x = irrigation)) +
    geom_boxplot(aes(y = seed.cotton), color = "gray") +
    geom_jitter(aes(y = seed.cotton), alpha = 0.2, size = 2, width = 0.1) +
    geom_point(data = yield_irr,
               aes(x = irrigation, y = emmean),
               size = 5, shape = 15) +
    geom_errorbar(data = yield_irr,
                  aes(x = irrigation, ymin = emmean - SE,
                      ymax = emmean + SE), linewidth = 1, width = 0.1) +
    labs(x="Treatments", y = ("Seed cotton yield (kg/ha)"))+
    annotate(geom="text", x= 1, y=3500, label="Irrigation: P = 0.008", size = 4.5) +
    theme)

# yield by residue

yield_residue <- as.data.frame(emmeans(yield_model, ~residue))

(yield_residue_graph <- ggplot(yield_data, aes(x = residue)) +
    geom_boxplot(aes(y = seed.cotton),color = "gray") +
    geom_jitter(aes(y = seed.cotton), alpha = 0.2, size = 2, width = 0.1) +
    geom_point(data = yield_residue,
               aes(x = residue, y = emmean),
               size = 5, shape = 15) +
    geom_errorbar(data = yield_residue,
                  aes(x = residue, ymin = emmean - SE,
                      ymax = emmean + SE), linewidth = 1, width = 0.1) +
    labs(x="Treatments", y = ("Seed cotton yield (kg/ha)"))+
    scale_x_discrete(labels = c("No residue", "Residue")) +
    annotate(geom="text", x= 1, y=3500, label="Residue: P = 0.03", size = 4.5) +
    theme)



####################################################
# # below ground biomass model


bgb_model <- lmer(rootbiomass ~ warming * residue * irrigation +
                    (1|wp), data = yield_data)


Anova(bgb_model)
plot(bgb_model)
qqnorm(residuals(bgb_model))
hist(residuals(bgb_model))
r.squaredGLMM(bgb_model)


# pairwise comparison for significant predictors

bgb.emm <- emmeans(bgb_model, ~ irrigation)
contrast(bgb.emm, "consec", simple = "each", combine = TRUE)

# model with continuous predictors only
bgb_model1 <- lmer(log(rootbiomass) ~ nitN + P + OM + (1|wp),
                   data = yield_data)

Anova(bgb_model1)
plot(bgb_model1)
qqnorm(residuals(bgb_model1))
hist(residuals(bgb_model1))



####################################################

# below ground biomass by trt

bgb <- as.data.frame(emmeans(bgb_model, ~warming * residue * irrigation))

(bgb_graph <- ggplot(yield_data, aes(x = warming, color = residue)) +
    geom_point(aes(y = rootbiomass, color  = residue), alpha = 0.2, size = 2,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
    geom_point(data = bgb,
               aes(x = warming, y = emmean, group = residue),
               size = 5, shape = 15,
               position = position_dodge(width = 0.75)) +
    geom_errorbar(data = bgb,
                  aes(x = warming, ymin = emmean - SE,
                      ymax = emmean + SE, group = residue), width = 0.2,
                  position = position_dodge(width = 0.75)) +
    facet_wrap(~ irrigation) + 
    labs(x="Treatments", y = ("Belowground Biomass (g)"))+
    scale_color_manual(labels = c("No residue", "Residue"),
                       values = c("#000000", "#009E73")) +
    theme +
    theme(strip.background = element_blank(),
          strip.text = element_blank()))

# belowground biomass by irrigation

bgb_irr <- as.data.frame(emmeans(bgb_model, ~irrigation))



(bgb_irr_graph <- ggplot(yield_data, aes(x = irrigation)) +
    geom_boxplot(aes(y = rootbiomass), color = "gray") +
    geom_jitter(aes(y = rootbiomass), alpha = 0.2, size = 2, width = 0.1) +
    geom_point(data = bgb_irr,
               aes(x = irrigation, y = emmean),
               size = 5, shape = 15) +
    geom_errorbar(data = bgb_irr,
                  aes(x = irrigation, ymin = emmean - SE,
                      ymax = emmean + SE), linewidth = 1, width = 0.1) +
    labs(x="Treatments", y = ("Belowground biomass (Kg/ha)"))+
    annotate(geom="text", x= 1, y = 4.5, label="Irrigation: P = 0.002", size = 4.5) +
    theme)

# belowground biomass by residue

bgb_residue <- as.data.frame(emmeans(bgb_model, ~residue))

(bgb_residue_graph <- ggplot(yield_data, aes(x = residue)) +
    geom_boxplot(aes(y = rootbiomass), color = "gray") +
    geom_jitter(aes(y = rootbiomass), alpha = 0.2, size = 2, width = 0.1) +
    geom_point(data = bgb_residue,
               aes(x = residue, y = emmean),
               size = 5, shape = 15) +
    geom_errorbar(data = bgb_residue,
                  aes(x = residue, ymin = emmean - SE,
                      ymax = emmean + SE), linewidth = 1, width = 0.1) +
    labs(x="Treatments", y = ("Belowground biomass (Kg/ha)"))+
    scale_x_discrete(labels = c("No residue", "Residue")) +
    annotate(geom="text", x= 1, y = 4.5, label="Residue: P = 0.05", size = 4.5) +
    theme)


fig_7 <- ggarrange(agb_irr_graph + rremove("xlab"),
          bgb_irr_graph + rremove("xlab"),
          yield_irr_graph + rremove("xlab"),
          yield_residue_graph + rremove("xlab"),
          align = "v",
          labels = "auto")

ggsave(filename = "Figure 7.pdf", plot = fig_7, width = 10, height = 8, units = "in",
       dpi = 350)


#######################################################
# plots

setwd("~/Desktop/MSProject/MS_Paper/graphs/")

ggsave(filename = "Figure 4.pdf", plot = fig_4, width = 10, height = 8, units = "in",
       dpi = 350)

ggsave(filename = "Figure 5.pdf", plot = fig_5, width = 10, height = 5, units = "in",
       dpi = 350)
# Supplemental Figure S3

ggarrange(om_graph + rremove("xlab") + rremove("x.text") ,
          mbc_graph + rremove("xlab") + rremove("x.text"),
          flux_graph + rremove("xlab"),
          nrow = 3,
          common.legend = T,
          align = "v",
          labels = "auto",
          vjust = 0.3)

# figure 6
fig_6 <- ggarrange(om_temp_graph,
          mbc_om_graph,
          mbc_nitrate_graph,
          common.legend = T,
          allign = "v",
          labels = "auto",
          vjust = 0.3)

ggsave(filename = "Figure 6.pdf", plot = fig_6, width = 10, height = 8, units = "in",
       dpi = 350)

# Supplemental Figure S4

ggarrange(agb_graph + rremove("xlab") + rremove("x.text") ,
          bgb_graph + rremove("xlab") + rremove("x.text"),
          yield_graph + rremove("xlab"),
          nrow = 3,
          common.legend = T,
          align = "v",
          labels = "auto",
          vjust = 0.3,
          hjust = -2.5)






























































