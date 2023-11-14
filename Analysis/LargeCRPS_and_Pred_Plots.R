####
## LargeCRPS_and_Pred_Plots.R
## Script to plot the predictive accuracy metrics for large data
####

## Data and Directory
load("~/R/Research/MiniBatch/BestFiles/Analysis/LargeResults.RData")
setwd("~/R/Research/MiniBatch/Bestfiles/Analysis")

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(gridExtra)
library(ggplotify)

# Data Preparation --------------------------------------------------------

## CRPS Data Prep
CRPS_dat1 <- tibble(CRPS = as.numeric(CRPS_AVG),
                    parameter = rep(rownames(CRPS_AVG),ncol(CRPS_AVG)),
                    model = c(rep(colnames(CRPS_AVG),each=nrow(CRPS_AVG))),
                    type='Continuous')
CRPS_dat1$model <- as_factor(CRPS_dat1$model)

CRPS_dat2 <- tibble(CRPS = as.numeric(CRPS_AVG_discrete),
                    parameter = rep(rownames(CRPS_AVG_discrete),ncol(CRPS_AVG_discrete)),
                    model = c(rep(colnames(CRPS_AVG_discrete),each=nrow(CRPS_AVG_discrete))),
                    type='Discrete')
CRPS_dat2$model <- as_factor(CRPS_dat2$model)

CRPS_dat <- bind_rows(CRPS_dat1, CRPS_dat2) 
CRPS_dat <- CRPS_dat %>% 
  mutate(model=str_remove_all(CRPS_dat$model, '_Gibbs')) 
CRPS_dat$model <- factor(CRPS_dat$model, 
                         levels = c("MB64", "MB32", "MB16", "MB8", 
                                    "MB4", "MB2", "BG", "NN", "Full"))
CRPS_dat <- CRPS_dat %>% 
  mutate(model1=as.numeric(model)) %>% 
  mutate(type=as_factor(type))

## Prediction Data Prep
load("Large_CRPS_for_preds.RData")

## CRPS for Predictions
pred_CRPS_dat <- tibble(CRPS = c(PRED_CRPS_AVG, PRED_CRPS_AVG_discrete),
                        type = c(rep('Continuous', 8), rep('Discrete', 8)),
                        model = c(rep(c("NN", "BG", "MB2", "MB4", 
                                        "MB8", "MB16", "MB32", "MB64"), 2))) %>% na.omit()
pred_CRPS_dat$model <- factor(pred_CRPS_dat$model, 
                              levels = c("MB64", "MB32", "MB16", "MB8", 
                                         "MB4", "MB2", "BG", "NN", "Full"))
pred_CRPS_dat <- pred_CRPS_dat %>%
  mutate(model1=as.numeric(model)) %>% 
  mutate(type=as_factor(type))

## RMSE for Predictions
RMSE_dat <- tibble(RMSE = c(AVG_RMSE, AVG_RMSE_discrete), 
                   type = c(rep('Continuous', 8), rep('Discrete', 8)),
                   model = c(rep(c("NN", "BG", "MB2", "MB4", 
                                   "MB8", "MB16", "MB32", "MB64"), 2))) %>% na.omit()
RMSE_dat$model <- factor(RMSE_dat$model, 
                         levels = c("MB64", "MB32", "MB16", "MB8", 
                                    "MB4", "MB2", "BG", "NN"))
RMSE_dat <- RMSE_dat %>%
  mutate(model1=as.numeric(model)) %>% 
  mutate(type=as_factor(type))

## ES Data Prep
ES_dat <- tibble(ES=c(ES_AVG, ES_AVG_discrete),
                 type=c(rep('Continuous',8), rep('Discrete',8)),
                 model=c(names(ES_AVG), names(ES_AVG_discrete)))
ES_dat$model <- factor(ES_dat$model, 
                       levels = c("MB64", "MB32", "MB16", "MB8", 
                                  "MB4", "MB2", "BG", "NN", "Full"))
ES_dat <- ES_dat %>%
  mutate(model1=as.numeric(model)) %>% 
  mutate(type=as_factor(type))

## ES For Predictions
#PRED_ES_dat <- tibble(ES=c(PRED_ES_AVG, PRED_ES_AVG_discrete),
#                 type=c(rep('Continuous',8), rep('Discrete',8)),
#                 model=c(rep(c("NN", "BG", "MB2", "MB4", 
#                               "MB8", "MB16", "MB32", "MB64"), 2))) %>% na.omit()
#PRED_ES_dat$model <- factor(PRED_ES_dat$model, 
#                            levels = c("MB64", "MB32", "MB16", "MB8", 
#                                       "MB4", "MB2", "BG", "NN", "Full"))
#PRED_ES_dat <- PRED_ES_dat %>%  
#  mutate(model1=as.numeric(model)) %>% 
#  mutate(type=as_factor(type))

## General Prep
plot_names <- c("1.56", "3.13", "6.25", "12.5", 
                "25", "50", "Barker", "NN")

# CRPS Plots --------------------------------------------------------------

## Beta0
b0 <- CRPS_dat %>% 
  filter(parameter == 'beta0') %>% 
  ggplot(aes(model1, CRPS, col=type)) +
  geom_point(show.legend=F) + geom_line(lwd=1, show.legend=F) +
  geom_vline(xintercept=6.5) +
  ggtitle(expression('a) '~beta[0])) +
  labs(y='CRPS', x='Model and % of Data Used') +
  scale_x_continuous(breaks=(1:ncol(CRPS_AVG)),
                     labels=plot_names)

## Beta1
b1 <- CRPS_dat %>% 
  filter(parameter=='beta1') %>% 
  arrange(CRPS) %>% 
  ggplot(aes(model1, CRPS, col=type)) +
  geom_point(show.legend=F) + geom_line(lwd=1, show.legend=F) +
  geom_vline(xintercept=6.5) +
  ggtitle(expression('b) '~beta[1])) +
  labs(y='CRPS', x='Model and % of Data Used') +
  scale_x_continuous(breaks=(1:ncol(CRPS_AVG)),
                     labels=plot_names)

## Beta2
b2 <- CRPS_dat %>% 
  filter(parameter=='beta2') %>% 
  arrange(CRPS) %>% 
  ggplot(aes(model1, CRPS, col=type)) +
  geom_point(show.legend=F) + geom_line(lwd=1, show.legend=F) +
  geom_vline(xintercept=6.5) +
  ggtitle(expression('c) '~beta[2])) +
  labs(y='CRPS', x='Model and % of Data Used') +
  scale_x_continuous(breaks=(1:ncol(CRPS_AVG)),
                     labels=plot_names)

## Beta3
b3 <- CRPS_dat %>% 
  filter(parameter=='beta3') %>% 
  arrange(CRPS) %>% 
  ggplot(aes(model1, CRPS, col=type)) +
  geom_point(show.legend=F) + geom_line(lwd=1, show.legend=F) +
  geom_vline(xintercept=6.5) +
  ggtitle(expression('d) '~beta[2])) +
  labs(y='CRPS', x='Model and % of Data Used') +
  scale_x_continuous(breaks=(1:ncol(CRPS_AVG)),
                     labels=plot_names)

## Sigma2
s2 <- CRPS_dat %>% 
  filter(parameter=='sigma2') %>% 
  arrange(CRPS) %>% 
  ggplot(aes(model1, CRPS, col=type)) +
  geom_point(show.legend=F) + geom_line(lwd=1, show.legend=F) +
  geom_vline(xintercept=6.5) +
  ggtitle(expression('e) '~sigma^2)) +
  labs(y='CRPS', x='Model and % of Data Used') +
  scale_x_continuous(breaks=(1:ncol(CRPS_AVG)),
                     labels=plot_names)

## omega*sigma2*alpha
val <- CRPS_dat %>% 
  filter(parameter=='value') %>% 
  arrange(CRPS) %>% 
  ggplot(aes(model1, CRPS, col=type)) +
  geom_point(show.legend=F) + geom_line(lwd=1, show.legend=F) +
  geom_vline(xintercept=6.5) +
  ggtitle(expression('f) '~sigma^2~omega~'/'~phi)) +
  labs(y='CRPS', x='Model and % of Data Used') +
  scale_x_continuous(breaks=(1:ncol(CRPS_AVG)),
                     labels=plot_names)


# ES Plots ----------------------------------------------------------------

## Plot the Energy Score for each Model
ES.plot <- ES_dat %>% 
  ggplot(aes(model1, ES, col=type)) +
  geom_point(show.legend=F) + geom_line(lwd=1, show.legend=F) +
  geom_vline(xintercept=6.5) +
  ggtitle("f)  Multivariate Energy Score") +
  labs(y='Energy Score', x='Model and % of Data Used', col='Prior Type') +
  scale_x_continuous(breaks=(1:8),
                     labels=plot_names)

# RMSE Plot ---------------------------------------------------------------

## RMSE Plotted with sd(y)
RMSE.plot <- RMSE_dat %>% 
  ggplot(aes(model1, RMSE, col=type)) +
  geom_point(show.legend=F) + geom_line(lwd=1, show.legend=F) +
  geom_hline(yintercept=AVG_sd) + 
  geom_vline(xintercept=6.5) +
  ggtitle('b)  Predictive Error') +
  labs(y='RMSE', x='Model and % of Data Used') +   
  ylim(0, 14) +
  scale_x_continuous(breaks=(1:8),
                     labels= plot_names)

## RMSE Plotted on its own
RMSE_only.plot <- RMSE_dat %>% 
  ggplot(aes(model1, RMSE, col=type)) +
  geom_point(show.legend=F) + 
  geom_line(lwd=1, show.legend=F) +
  geom_vline(xintercept=6.5) +
  ggtitle('b)  Predictive Error') +
  labs(y='RMSE', x='Model and % of Data Used') +
  scale_x_continuous(breaks=(1:8),
                     labels= plot_names)


# CRPS for Prediction Plots -----------------------------------------------

## Plot of CRPS values on Predictions
pred_CRPS.plot <- pred_CRPS_dat %>% 
  ggplot(aes(model1, CRPS, col=type)) +
  geom_point(show.legend=F) + geom_line(lwd=1, show.legend=F) +
  geom_vline(xintercept=6.5) +
  ggtitle('a)  CRPS For Predictions') +
  labs(y='CRPS', x='Model and % of Data Used') +  
  ylim(0, 3.5) +
  scale_x_continuous(breaks=(1:8),
                     labels= plot_names)

## Plot of CRPS values on Predictions
pred_CRPS_only.plot <- pred_CRPS_dat %>% 
  ggplot(aes(model1, CRPS, col=type)) +
  geom_point(show.legend=F) + geom_line(lwd=1, show.legend=F) +
  geom_vline(xintercept=6.5) +
  ggtitle('a)  CRPS For Predictions') +
  labs(y='CRPS', x='Model and % of Data Used') +   
  scale_x_continuous(breaks=(1:8),
                     labels= plot_names)


# ES for Prediction Plots -------------------------------------------------

## Plot the Energy Score for the Predictions
#pred_ES.plot <- PRED_ES_dat %>% 
#  ggplot(aes(model1, ES, col=type)) +
#  geom_point() + geom_line(lwd=1, show.legend=F) +
#  geom_vline(xintercept=6.5) +
#  ggtitle("Multivariate Energy Score") +
#  labs(y='Energy Score', x='Model and % of Data Used', col='Prior Type') +
#  scale_x_continuous(breaks=(1:8),
#                     labels=plot_names)

# Legend ------------------------------------------------------------------

legend.plot <- RMSE_dat %>% 
  ggplot(aes(model1, RMSE, col=type)) +
  geom_point() + xlim(0,0) +
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank()) +
  xlab(element_blank()) + ylab(element_blank()) + 
  scale_color_manual(name='Prior Type',
                     breaks=c('Continuous', 'Discrete'),
                     values=c('Continuous'='#F8776D', 'Discrete'='#00BFC4'))

# Plot Grids --------------------------------------------------------------

## All CRPS Plots
grid.arrange(b0, b1, b2,
             b3, s2, val, nrow=2,
             top="CRPS for Parameters and Predictions",
             right=as.grob(legend.plot))

## All RMSE Plots
grid.arrange(RMSE.plot, 
             RMSE_only.plot,
             right=as.grob(legend.plot))

## All CRPS Prediction Plots
grid.arrange(pred_CRPS.plot, 
             pred_CRPS_only.plot,
             right=as.grob(legend.plot))

## All CRPS and RMSE Prediction Plots
grid.arrange(pred_CRPS.plot, 
             RMSE.plot,
             pred_CRPS_only.plot,
             RMSE_only.plot,
             right=as.grob(legend.plot))


# All Plots ---------------------------------------------------------------

## Parameter CRPS
b0
b1
b2
b3
s2
val

## Energy Score Plots
ES.plot

## Energy Score for Predictions
#pred_ES.plot

## RMSE Plots
RMSE.plot
RMSE_only.plot

## CRPS for Preds
pred_CRPS.plot
pred_CRPS_only.plot

## General Plots
legend.plot


# Creating PDFs -----------------------------------------------------------

## CRPS plots
pdf(file = "./Plots/Large_CRPS.pdf",
    width = 11, height = 6)

grid.arrange(b0, b1, b2, 
             b3, s2, val, nrow=2,
             top="CRPS for Parameters ",
             right=as.grob(legend.plot))

dev.off()

## Prediction plots
pdf(file = "./Plots/Large_RMSE.pdf",
    width = 9, height = 3.25)

grid.arrange(pred_CRPS_only.plot,
             RMSE_only.plot, nrow=1,
             top="CRPS and RMSE for Predictions",
             right=as.grob(legend.plot))

dev.off()
