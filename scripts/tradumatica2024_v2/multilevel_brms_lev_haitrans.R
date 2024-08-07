library(brms)
#library(cmdstanr)
library(tidyverse)
library(sjPlot)
#install.packages("tidybayes")
library(brmstools)
library(tidybayes)
library(sjstats)
library(gapminder)
#install.packages("sjstats")
#install.packages('assertthat')
#install.packages('remotes')
#install.packages('brmsfit')
#library(brmsfit)
#install.packages('ggpubr')
#install.packages('rstatix')
#install.packages('datarium')
#install.packages('report')
#2remotes::install_github('m-clark/mixedup')
#install.packages('mixedup')
library(assertthat)
library(report)
library(ggpubr)
library(rstatix)
library(mixedup)
library(broom.mixed)
library(ggplot2)
library(dplyr)
library(ggh4x)
library(ggdist)
library("bayestestR")



set.seed(1234) #reproducibility
bayes_seed <- 1234
#set_cmdstan_path(path = '/home/mrios/.cmdstan/cmdstan-2.32.2')
eyetracking <- read.csv('/home/mrios/workspace/imminent_R/Tradumatica_v2/lev_prod_results.csv', header = TRUE, sep = ",")
eyetracking
eyetracking %>% sample_n_by(condition, text, size = 1)
#delete NA
#eyetracking <- na.omit(eyetracking, cols=c('productivity_Matecat'))
eyetracking <- eyetracking[!(is.na(eyetracking$pemt_speed)), ]
eyetracking <- eyetracking[!(is.na(eyetracking$no_of_searches_in_external_resources)), ]

colnames(eyetracking)
#summary stats
stats <- eyetracking %>%
  group_by(condition, text) %>%
  get_summary_stats(dlev_corpus, type = "mean_sd")
stats
write.csv(stats,'/home/mrios/workspace/imminent_R/Tradumatica_v2/prod_haitransTradumatica2024/summary_stats_lev.csv')




#
###
#MFD_ST ~ 1 + condition + text 
###

####
#model participant
#MFD_ST ~ 1 + condition + text + (1 + condition | participant)
####
fit2 <- brm(formula = dlev_corpus ~ 1 + condition + text + no_of_searches_in_external_resources + translation_years + PEMT_years + (1 + condition | participant), 
            data = eyetracking,
            warmup = 1000, iter = 10000, chains = 4, cores = 6,
            control=list(adapt_delta=0.9), seed = bayes_seed
            )
fit2
#r <- 0.05 * sd(eyetracking$pemt_speed)
#r
#tab_model(fit2)
text_summ <-summary(fit2)
sink("/home/mrios/workspace/imminent_R/Tradumatica_v2/prod_haitransTradumatica2024/multilevel_brms2_lev_summary.txt")
text_summ
sink()
#print(fit2)
#sjPlot::tab_model(fit2)
fit2
p_summary <- posterior_summary(fit2)
p_summary
write.csv(p_summary, '/home/mrios/workspace/imminent_R/Tradumatica_v2/prod_haitransTradumatica2024/multilevel_brms2_psummary.csv')
#write.csv(text_summ, '/home/mrios/workspace/imminent_R/quality/multilevel_brms1_summary.txt')
randeff <- ranef(fit2)
randeff
write.csv(randeff, '/home/mrios/workspace/imminent_R/Tradumatica_v2/prod_haitransTradumatica2024/multilevel_brms2_randeff.csv')
#plot(fit2)
pp_check(fit2)

fixef(fit2)
coef(fit2) 
extract_random_effects(fit2)
print(extract_random_effects(fit2), n=42)

conditional_effects(fit2)



post_fit2 <- describe_posterior(fit2)
post_fit2
write.csv(post_fit2, '/home/mrios/workspace/imminent_R/Tradumatica_v2/prod_haitransTradumatica2024/multilevel_brms2_lev_describepost.csv')




#pemt vs lev
#####

#tidy(fit2)
#describe_posterior(fit2)
fit3 <- brm(formula = dlev_corpus ~ 1 + condition + pemt_speed + (1 + pemt_speed | condition), 
            data = eyetracking,
            warmup = 1000, iter = 10000, chains = 4, cores = 3,
            control=list(adapt_delta=0.99), seed = bayes_seed
)
fit3
#r <- 0.05 * sd(eyetracking$pemt_speed)
#r
#tab_model(fit2)
text_summ <-summary(fit3)
sink("/home/mrios/workspace/imminent_R/Tradumatica_v2/prod_haitransTradumatica2024/multilevel_brms3_lev_summary.txt")
text_summ
sink()
#print(fit2)
#sjPlot::tab_model(fit2)
fit3

pp_check(fit3)

fixef(fit3)
coef(fit3) 
extract_random_effects(fit3)
print(extract_random_effects(fit3), n=42)

ceff <- conditional_effects(fit3, effects="pemt_speed:condition", re_formula = NULL)

plot(ceff, plot = FALSE) [[1]] +
  labs(x = "PEMT speed", y = "Edit distance")



post_fit3 <- describe_posterior(fit3)
post_fit3
write.csv(post_fit3, '/home/mrios/workspace/imminent_R/Tradumatica_v2/prod_haitransTradumatica2024/multilevel_brms3_lev_describepost.csv')

#quality

fit4 <- brm(formula = dlev_corpus ~ 1 + condition + quality_score + (1 + quality_score | condition), 
            data = eyetracking,
            warmup = 1000, iter = 10000, chains = 4, cores = 3,
            control=list(adapt_delta=0.99), seed = bayes_seed
)
fit4
#r <- 0.05 * sd(eyetracking$pemt_speed)
#r
#tab_model(fit2)
text_summ <-summary(fit4)
sink("/home/mrios/workspace/imminent_R/Tradumatica_v2/prod_haitransTradumatica2024/multilevel_brms4_lev_summary.txt")
text_summ
sink()
#print(fit2)
#sjPlot::tab_model(fit2)
fit4

pp_check(fit4)

fixef(fit4)
coef(fit4) 
extract_random_effects(fit4)
print(extract_random_effects(fit4), n=42)

ceff <- conditional_effects(fit4, effects="quality_score:condition", re_formula = NULL)
ceff
plot(ceff, plot = FALSE) [[1]] +
  labs(x = "Quality", y = "Edit distance")

post_fit4 <- describe_posterior(fit4)
post_fit4
write.csv(post_fit4, '/home/mrios/workspace/imminent_R/Tradumatica_v2/prod_haitransTradumatica2024/multilevel_brms4_lev_describepost.csv')
