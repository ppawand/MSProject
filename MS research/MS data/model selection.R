########################
# Soil Organic matter
#######################

nutrients <- read_excel("~/Desktop/MS data/nutrients.xlsx")
nutrients[, 2:6] <- lapply(nutrients[, 2:6], factor)
nutrients$trt <- factor(nutrients$trt, levels = c("C", "W", "R", "WR"))
nutrients$log.om <- log10(nutrients$OM)

# SOM model selection

model_c <- lmer(log.om ~ trt + irrigation + (1|block/plot) , REML = F, 
                data = nutrients)
model_c1 <- lmer(log.om ~ trt + irrigation + water.content + (1|block/plot) , REML = F,
                data = nutrients)
anova(model_c, model_c1)

model_c2 <- lmer(log.om ~ trt + irrigation + water.content + temp + (1|block/plot) ,
                 REML = F,
                 data = nutrients)
anova(model_c1, model_c2)

model_c3 <- lmer(log.om ~ trt + irrigation + water.content + trt:irrigation + (1|block/plot) ,
                 REML = F,
                 data = nutrients)
anova(model_c1, model_c3)

model_c4 <- lmer(log.om ~ trt + irrigation + water.content + pH + (1|block/plot) ,
                 REML = F,
                 data = nutrients)

anova(model_c1, model_c4)


# best fitted model
anova(model_c1)
summary(model_c1)
plot(model_c1)
qqnorm(residuals(model_c1))
hist(residuals(model_c1))

#####################
# Microbial Biomass
#####################
mbc <- read_excel("~/Desktop/MS data/MBC/mbc.final.xlsx")
mbc[, 2:6] <- lapply(mbc[, 2:6], factor)
mbc$trt <- factor(mbc$trt, levels = c("C", "W", "R", "WR"))


# MBC model selection

model_mbc <-lmer(mbc ~ trt + irrigation + (1|block/plot), REML = F,
                 data = mbc)
model_mbc1 <-lmer(mbc ~ trt + irrigation + water.content + (1|block/plot), REML = F,
                 data = mbc)
anova(model_mbc, model_mbc1)

model_mbc2 <-lmer(mbc ~ trt + irrigation + temp + (1|block/plot), REML = F,
                 data = mbc)
anova(model_mbc, model_mbc2)

model_mbc3 <-lmer(mbc ~ trt + irrigation + trt:irrigation + (1|block/plot), REML = F,
                 data = mbc)
anova(model_mbc, model_mbc3)

model_mbc4 <-lmer(mbc ~ trt + irrigation + pH + (1|block/plot), REML = F,
                 data = mbc)
anova(model_mbc, model_mbc4)
model_mbc5 <-lmer(mbc ~ trt + irrigation + OM + (1|block/plot), REML = F,
                 data = mbc)
anova(model_mbc, model_mbc5)

model_mbc6 <-lmer(mbc ~ trt + irrigation + P+ (1|block/plot), REML = F,
                 data = mbc)
anova(model_mbc, model_mbc6)

# best fit model
summary(model_mbc)
anova(model_mbc)
plot(model_mbc)
qqnorm(residuals(model_mbc))
hist(residuals(model_mbc))

# respiration model selection
respiration <- read_excel("~/Desktop/MS data/Respiration data/respiration.xlsx")
respiration[, 2:7] <- lapply(respiration[, 2:7], factor)
respiration$trt <- factor(respiration$trt, labels = c("C", "W", "R", "WR"))

# model selection
resp_model <- lmer(log.flux ~ trt + irrigation +
                     (1|block/plot) + (1|month), REML = F, data  = respiration)
resp_model1 <-lmer(log.flux ~ trt + irrigation + month +
                     (1|block/plot), REML = F,  data  = respiration)
anova(resp_model, resp_model1)

resp_model2 <-lmer(log.flux ~ trt + irrigation + month +
                     (1|block/plot) + (1|month), REML = F, data  = respiration)
anova(resp_model1, resp_model2)

resp_model3 <-lmer(log.flux ~ trt + irrigation + month + temp +
                     (1|block/plot), REML = F, data  = respiration)
anova(resp_model1, resp_model3)

resp_model4 <-lmer(log.flux ~ trt + irrigation + month + moist +
                     (1|block/plot), REML = F, data  = respiration)
anova(resp_model1, resp_model4)

resp_model5 <-lmer(log.flux ~ trt + irrigation + month + pH +
                     (1|block/plot), REML = F, data  = respiration)
anova(resp_model1, resp_model5)

resp_model6 <-lmer(log.flux ~ trt + irrigation + month + OM +
                     (1|block/plot), REML = F, data  = respiration)
anova(resp_model1, resp_model6)

resp_model7 <-lmer(log.flux ~ trt + irrigation + month + trt:irrigation +
                     (1|block/plot), REML = F, data  = respiration)
anova(resp_model1, resp_model7)

resp_model8 <-lmer(log.flux ~ trt + irrigation + month + trt:month +
                     (1|block/plot), REML = F, data  = respiration)
anova(resp_model1, resp_model8)

resp_model9 <-lmer(log.flux ~ trt + irrigation + month + trt:month + irrigation:month +
                     (1|block/plot), REML = F, data  = respiration)
anova(resp_model8, resp_model9)

# best fit model
summary(resp_model9)
anova(resp_model9)
plot(resp_model9)
qqnorm(residuals(resp_model9))
hist(residuals(resp_model9))





