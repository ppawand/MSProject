library(ggplot2)
library(readxl)
library(dplyr)

data <- read_xlsx("~/desktop/photosyn.xlsx")

data <- data %>% group_by(plot, trt, date) %>% 
  summarize(mean = mean(A), 
            n = n(),
            se = sd(A)/sqrt(n))

ggplot(data, aes(factor(date,levels = c("July", "August", "September")), mean, fill = trt)) + 
  geom_col(position = position_dodge(preserve = "single"),
           width = 0.9) +
  theme_bw() + 
  facet_wrap(~ plot) + scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank(),
        panel.background = element_rect(color = NA),
        plot.background = element_rect(color = NA)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.5,
                position = position_dodge(preserve = "single", width = 0.9)) +
  xlab("Months") +
  ylab(bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')'))
 
data2 <- data %>% group_by(plot, trt, date) %>% 
  summarize(mean = mean(E), 
            n = n(),
            se = sd(E)/sqrt(n))

ggplot(data2, aes(factor(date,levels = c("July", "August", "September")), mean, fill = trt)) + 
  geom_col(position = position_dodge(preserve = "single"),
           width = 0.9) +
  theme_bw() + 
  facet_wrap(~ plot) + scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank(),
        panel.background = element_rect(color = NA),
        plot.background = element_rect(color = NA)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.5,
                position = position_dodge(preserve = "single", width = 0.9)) +
  xlab("Months") +
  ylab(bquote('Transpiration (mol' ~ m^-2~s^-1*')'))

data3 <- data %>% group_by(plot, trt, date) %>% 
  summarize(mean = mean(gs), 
            n = n(),
            se = sd(E)/sqrt(n))

ggplot(data3, aes(factor(date,levels = c("July", "August", "September")), mean, fill = trt)) + 
  geom_col(position = position_dodge(preserve = "single"),
           width = 0.9) +
  theme_bw() + 
  facet_wrap(~ plot) + scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank(),
        panel.background = element_rect(color = NA),
        plot.background = element_rect(color = NA)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.5,
                position = position_dodge(preserve = "single", width = 0.9)) +
  xlab("Months") +
  ylab(bquote('Stomatal Conductance (mol' ~ m^-2~s^-1*')'))
        