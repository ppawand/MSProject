
# bestFitModel function for model selection for time series data

bestFitModel <- function (x = monthly.avg.TM, dependent_var = x$temp) {
  
  m1 <- lmer(dependent_var ~ warming * residue * irrigation + (1|block) + (1|month), REML = F,
             data = x)
  
  m2 <- lmer(dependent_var ~ warming + residue + irrigation + warming:residue +
               warming:irrigation + residue:irrigation + (1|block) + (1|month), REML = F,
             data = x)
  
  m3 <- lmer(dependent_var ~ warming + residue + irrigation + warming:residue +
               warming:irrigation + (1|block) + (1|month), REML = F,
             data = x)
  
  m4 <- lmer(dependent_var ~ warming + residue + irrigation + warming:residue +
               residue:irrigation + (1|block) + (1|month), REML = F,
             data = x)
  
  
  m5 <- lmer(dependent_var ~ warming + residue + irrigation +
               warming:irrigation + residue:irrigation + (1|block) + (1|month), REML = F,
             data = x)
  
  
  m6 <- lmer(dependent_var ~ warming + residue + irrigation + warming:residue +
               (1|block) + (1|month), REML = F,
             data = x)
  
  m7 <- lmer(dependent_var ~ warming + residue + irrigation +
               warming:irrigation + (1|block) + (1|month), REML = F,
             data = x)
  
  m8 <- lmer(dependent_var ~ warming + residue + irrigation +
               residue:irrigation + (1|block) + (1|month), REML = F,
             data = x)
  
  m9 <- lmer(dependent_var ~ warming + residue + irrigation +
               (1|block) + (1|month), REML = F,
             data = x)
  aic <- as.data.frame(AIC(m1, m2, m3, m4, m5, m6, m7, m8, m9))
  
  best.model <- aic[which.min(aic$AIC), ]
  
  return(best.model)
}




# bestFitModel2 function for model selection (non time series data)

bestFitModel2 <- function (x = nutrients, dependent_var = x$log.om) {
  
  m1 <- lmer(dependent_var ~ warming * residue * irrigation + (1|block), REML = F,
             data = x)
  
  m2 <- lmer(dependent_var ~ warming + residue + irrigation + warming:residue +
               warming:irrigation + residue:irrigation + (1|block), REML = F,
             data = x)
  
  m3 <- lmer(dependent_var ~ warming + residue + irrigation + warming:residue +
               warming:irrigation + (1|block), REML = F,
             data = x)
  
  m4 <- lmer(dependent_var ~ warming + residue + irrigation + warming:residue +
               residue:irrigation + (1|block), REML = F,
             data = x)
  
  
  m5 <- lmer(dependent_var ~ warming + residue + irrigation +
               warming:irrigation + residue:irrigation + (1|block), REML = F,
             data = x)
  
  
  m6 <- lmer(dependent_var ~ warming + residue + irrigation + warming:residue +
               (1|block), REML = F,
             data = x)
  
  m7 <- lmer(dependent_var ~ warming + residue + irrigation +
               warming:irrigation + (1|block), REML = F,
             data = x)
  
  m8 <- lmer(dependent_var ~ warming + residue + irrigation +
               residue:irrigation + (1|block), REML = F,
             data = x)
  
  m9 <- lmer(dependent_var ~ warming + residue + irrigation +
               (1|block), REML = F,
             data = x)
  aic <- as.data.frame(AIC(m1, m2, m3, m4, m5, m6, m7, m8, m9))
  
  best.model <- aic[which.min(aic$AIC), ]
  
  return(best.model)
}


# bestGeneralizedModel function for model selection (non time series data)

bestGeneralizedModel <- function (x = nutrients, dependent_var = x$OM) {
  
  m1 <- glmer(dependent_var ~ warming * residue * irrigation + (1|block), 
             family = Gamma(link = "log"),
             data = x)
  
  m2 <- glmer(dependent_var ~ warming + residue + irrigation + warming:residue +
               warming:irrigation + residue:irrigation + (1|block),
              family = Gamma(link = "log"),
              data = x)
  
  m3 <- glmer(dependent_var ~ warming + residue + irrigation + warming:residue +
               warming:irrigation + (1|block), 
              family = Gamma(link = "log"),
              data = x)
  
  m4 <- glmer(dependent_var ~ warming + residue + irrigation + warming:residue +
               residue:irrigation + (1|block),
              family = Gamma(link = "log"),
              data = x)
  
  
  m5 <- glmer(dependent_var ~ warming + residue + irrigation +
               warming:irrigation + residue:irrigation + (1|block),
              family = Gamma(link = "log"),
              data = x)
  
  
  m6 <- glmer(dependent_var ~ warming + residue + irrigation + warming:residue +
               (1|block),
              family = Gamma(link = "log"),
              data = x)
  
  m7 <- glmer(dependent_var ~ warming + residue + irrigation +
               warming:irrigation + (1|block),
              family = Gamma(link = "log"),
              data = x)
  
  m8 <- glmer(dependent_var ~ warming + residue + irrigation +
               residue:irrigation + (1|block),
              family = Gamma(link = "log"),
              data = x)
  
  m9 <- glmer(dependent_var ~ warming + residue + irrigation +
               (1|block),
              family = Gamma(link = "log"),
              data = x)
  aic <- as.data.frame(AIC(m1, m2, m3, m4, m5, m6, m7, m8, m9))
  
  best.model <- aic[which.min(aic$AIC), ]
  
  return(best.model)
}


# bestGeneralizedModel function for model selection for time series data

bestGeneralizedModel2 <- function (x = nutrients, dependent_var = x$OM) {
  
  m1 <- glmer(dependent_var ~ warming * residue * irrigation + (1|block) + (1|month), 
              family = Gamma(link = "log"),
              control = glmerControl(optimizer = "bobyqa"),
              data = x)
  
  m2 <- glmer(dependent_var ~ warming + residue + irrigation + warming:residue +
                warming:irrigation + residue:irrigation + (1|block) + (1|month),
              family = Gamma(link = "log"),
              control = glmerControl(optimizer = "bobyqa"),
              data = x)
  
  m3 <- glmer(dependent_var ~ warming + residue + irrigation + warming:residue +
                warming:irrigation + (1|block) + (1|month), 
              family = Gamma(link = "log"),
              control = glmerControl(optimizer = "bobyqa"),
              data = x)
  
  m4 <- glmer(dependent_var ~ warming + residue + irrigation + warming:residue +
                residue:irrigation + (1|block) + (1|month),
              family = Gamma(link = "log"),
              control = glmerControl(optimizer = "bobyqa"),
              data = x)
  
  
  m5 <- glmer(dependent_var ~ warming + residue + irrigation +
                warming:irrigation + residue:irrigation + (1|block) + (1|month),
              family = Gamma(link = "log"),
              control = glmerControl(optimizer = "bobyqa"),
              data = x)
  
  
  m6 <- glmer(dependent_var ~ warming + residue + irrigation + warming:residue +
                (1|block)+ (1|month),
              family = Gamma(link = "log"),
              control = glmerControl(optimizer = "bobyqa"),
              data = x)
  
  m7 <- glmer(dependent_var ~ warming + residue + irrigation +
                warming:irrigation + (1|block) + (1|month),
              family = Gamma(link = "log"),
              control = glmerControl(optimizer = "bobyqa"),
              data = x)
  
  m8 <- glmer(dependent_var ~ warming + residue + irrigation +
                residue:irrigation + (1|block) + (1|month),
              family = Gamma(link = "log"),
              control = glmerControl(optimizer = "bobyqa"),
              data = x)
  
  m9 <- glmer(dependent_var ~ warming + residue + irrigation +
                (1|block) + (1|month),
              family = Gamma(link = "log"),
              control = glmerControl(optimizer = "bobyqa"),
              data = x)
  aic <- as.data.frame(AIC(m1, m2, m3, m4, m5, m6, m7, m8, m9))
  
  best.model <- aic[which.min(aic$AIC), ]
  
  return(best.model)
}
