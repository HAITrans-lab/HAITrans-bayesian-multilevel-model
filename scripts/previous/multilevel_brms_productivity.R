library(brms)
library(cmdstanr)
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
eyetracking <- read.csv('/home/mrios/workspace/test_R/table_productivity.csv', header = TRUE, sep = ",")
eyetracking
eyetracking %>% sample_n_by(condition, text, size = 1)
#delete NA
#eyetracking <- na.omit(eyetracking, cols=c('productivity_Matecat'))
eyetracking <- eyetracking[!(is.na(eyetracking$productivity_Matecat)), ]

colnames(eyetracking)
#summary stats
stats <- eyetracking %>%
  group_by(condition, text) %>%
  get_summary_stats(productivity_Matecat, type = "mean_sd")
stats
write.csv(stats,'/home/mrios/workspace/test_R/prod/summary_stats_st.csv')

##
#condition model
##
fit0 <- brm(formula = productivity_Matecat ~ 1 + condition, #(1 + condition | participant)
            data = eyetracking,
            warmup = 1000, iter = 5000, chains = 4, cores = 6,
            seed = bayes_seed,
            backend = "cmdstanr"
)
fit0
text_summ <-summary(fit0)
sink("/home/mrios/workspace/test_R/prod/multilevel_brms0_summary.txt")
text_summ
#tidy(fit1)
sink()
#sjPlot::tab_model(fit1)
fit0
p_summary <- posterior_summary(fit0)
write.csv(p_summary, '/home/mrios/workspace/test_R/prod/multilevel_brms0_psummary.csv')

plot(fit0)
pp_check(fit0)

fixef(fit0)
#coef(fit0)  #coef(fit1)$condition TODO
#coef(fit1)$text[conditions]
conditional_effects(fit0)



loo0 <- loo(fit0) 
loo0

tidy(fit0)
tidy(fit0)
describe_posterior(fit0)
bayestestR::hdi(fit0)
###
#MFD_ST ~ 1 + condition + text 
###

fit0b <- brm(formula = productivity_Matecat ~ 1 + condition + text, #(1 + condition | participant)
            data = eyetracking,
            warmup = 1000, iter = 5000, chains = 4, cores = 6,
            seed = bayes_seed,
            backend = "cmdstanr"
)
fit0b
text_summ <-summary(fit0b)
sink("/home/mrios/workspace/test_R/prod/multilevel_brms0b_summary.txt")
text_summ
#tidy(fit1)
sink()
#sjPlot::tab_model(fit1)
fit0b
p_summary <- posterior_summary(fit0b)
write.csv(p_summary, '/home/mrios/workspace/test_R/prod/multilevel_brms0b_psummary.csv')

plot(fit0b)
pp_check(fit0b)

fixef(fit0b)
#coef(fit0)  #coef(fit1)$condition TODO
#coef(fit1)$text[conditions]

conditional_effects(fit0b)


loo0b <- loo(fit0b) 
loo0b

tidy(fit0b)

loo_compare(loo0, loo0b)

tidy(fit0b)
describe_posterior(fit0b)
bayestestR::hdi(fit0b)
##
#multilevel MFD_ST ~ 1 + condition + text + (1 + condition | text)
##

fit1 <- brm(formula = productivity_Matecat~ 1 + condition + text + (1 + condition | text), #(1 + condition | participant)
            data = eyetracking,
            warmup = 1000, iter = 10000, chains = 4, cores = 6,
            control=list(adapt_delta=0.99), seed = bayes_seed,
            backend = "cmdstanr"
            )
fit1
text_summ <-summary(fit1)
sink("/home/mrios/workspace/test_R/prod/multilevel_brms1_summary.txt")
text_summ
#tidy(fit1)
sink()
#sjPlot::tab_model(fit1)
fit1
p_summary <- posterior_summary(fit1)
write.csv(p_summary, '/home/mrios/workspace/test_R/prod/multilevel_brms1_psummary.csv')
randeff <- ranef(fit1)
randeff
write.csv(randeff, '/home/mrios/workspace/test_R/prod/multilevel_brms1_randeff.csv')

plot(fit1)
pp_check(fit1)

fixef(fit1)
coef(fit1)  #coef(fit1)$condition TODO
#coef(fit1)$text[conditions]
extract_random_effects(fit1)


conditional_effects(fit1)



loo1 <- loo(fit1) 
loo1

loo_compare(loo0, loo1)

(model_fit <- eyetracking %>%
    add_predicted_draws(fit1) %>%  # adding the posterior distribution
    ggplot(aes(x = text, y = productivity_Matecat)) +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = eyetracking, colour = "darkseagreen4", size = 3) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ylab("productivity_Matecat\n") +  # latin name for red knot
    xlab("\ntext") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85)))

(model_fit <- eyetracking %>%
    add_predicted_draws(fit1) %>%  # adding the posterior distribution
    ggplot(aes(x = condition, y = productivity_Matecat)) +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = eyetracking, colour = "darkseagreen4", size = 3) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ylab("productivity_Matecat\n") +  # latin name for red knot
    xlab("\ncondition") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85)))


tidy(fit1)
#text level random effects
#https://www.andrewheiss.com/blog/2021/12/01/multilevel-models-panel-data-guide/#intercepts-and-slopes-for-each-country
condition_text_offsets <- ranef(fit1)$text %>%
  as_tibble(rownames = "text") %>% 
  filter(text %in% eyetracking$text) %>% 
  select(text, starts_with("Estimate"))
condition_text_offsets
write.csv(condition_text_offsets, '/home/mrios/workspace/test_R/prod/multilevel_brms1_condoffsets.csv')
#fixed effect + random offset for text-specific intercepts and slopes.
coef(fit1)$text %>%
  as_tibble(rownames = "text") %>% 
  filter(text %in% eyetracking$text) %>% 
  select(text, starts_with("Estimate"))

forest(fit1, pars='conditions')

tidy(fit1)
describe_posterior(fit1)
bayestestR::hdi(fit1)

condition_draws <- fit1 %>%
  spread_draws(b_conditions)

ggplot(condition_draws, aes(x = b_conditions)) +
  stat_halfeye(aes(fill_ramp = stat(x >= 0.5 | x <= -0.5)), fill = "#CCCCCC") +
  scale_fill_ramp_discrete(from = "darkgrey", guide = "none") +
  annotate(geom = "rect", xmin = -0.5, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "darkred", alpha = 0.3) +
  annotate(geom = "label", x = 0, y = 0.75, label = "ROPE") +
  labs(caption = "Point shows median value;\nthick black bar shows 66% credible interval;\nthin black bar shows 95% credible interval")

####
#model participant
#MFD_ST ~ 1 + condition + text + (1 + condition | participant)
####
fit2 <- brm(formula = productivity_Matecat ~ 1 + condition + text + (1 + condition | participant), 
            data = eyetracking,
            warmup = 1000, iter = 10000, chains = 4, cores = 6,
            control=list(adapt_delta=0.9), seed = bayes_seed,
            backend = "cmdstanr"
            )
fit2
tab_model(fit2)
text_summ <-summary(fit2)
sink("/home/mrios/workspace/test_R/prod/multilevel_brms2_summary.txt")
text_summ
sink()
#print(fit2)
#sjPlot::tab_model(fit2)
fit2
p_summary <- posterior_summary(fit2)
p_summary
write.csv(p_summary, '/home/mrios/workspace/test_R/prod/multilevel_brms2_psummary.csv')
#write.csv(text_summ, '/home/mrios/workspace/test_R/quality/multilevel_brms1_summary.txt')
randeff <- ranef(fit2)
randeff
write.csv(randeff, '/home/mrios/workspace/test_R/prod/multilevel_brms2_randeff.csv')
plot(fit2)
pp_check(fit2)

fixef(fit2)
coef(fit2) 
extract_random_effects(fit2)
print(extract_random_effects(fit2), n=42)

conditional_effects(fit2)


loo2 <- loo(fit2)

loo_cpm <- loo_compare(loo1, loo2)
loo_cpm
write.csv(loo_cpm, '/home/mrios/workspace/test_R/prod/multilevel_brms_loo_cpm1-2.csv')
#loo_compare(loo1, loo2)
#loo_compare(loo0, loo0b, loo1, loo2)


(model_fit <- eyetracking %>%
    add_predicted_draws(fit2) %>%  # adding the posterior distribution
    ggplot(aes(x = text, y = MFD_ST)) +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = eyetracking, colour = "darkseagreen4", size = 3) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ylab("cog_load\n") +  # latin name for red knot
    xlab("\ntext") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85)))

(model_fit <- eyetracking %>%
    add_predicted_draws(fit2) %>%  # adding the posterior distribution
    ggplot(aes(x = condition, y = productivity_Matecat)) +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = eyetracking, colour = "darkseagreen4", size = 3) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ylab("productivity_Matecat\n") +  # latin name for red knot
    xlab("\ncondition") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85)))

#participant level random effects
condition_participant_offsets <- ranef(fit2)$participant %>%
  as_tibble(rownames = "participant") %>% 
  filter(participant %in% eyetracking$participant) %>% 
  select(participant, starts_with("Estimate"))
condition_participant_offsets
write.csv(condition_participant_offsets, '/home/mrios/workspace/test_R/prod/multilevel_brms2_condoffsets.csv')
#fixed effect + random offset for text-specific intercepts and slopes.
coef(fit2)$participant %>%
  as_tibble(rownames = "participant") %>% 
  filter(participant %in% eyetracking$participant) %>% 
  select(participant, starts_with("Estimate"))

forest(fit2, pars='conditions')

tidy(fit2)
describe_posterior(fit2)
bayestestR::hdi(fit2)

condition_draws <- fit2 %>%
  spread_draws(b_conditions)
condition_draws


condition_cred_int <- condition_draws %>% 
  median_hdi()
condition_cred_int$.lower


condition_draws %>% 
  # Only look inside the credible interval
  filter(b_conditions >= condition_cred_int$.lower & b_conditions <= condition_cred_int$.upper) %>% 
  summarize(prop_outside_rope = 1 - sum(b_conditions >= -564.52 & b_conditions <= 564.52) / n())
#-564.52, 564.52
ggplot(condition_draws, aes(x = b_conditions)) +
  stat_halfeye(aes(fill_ramp = stat(x >= 564.52 | x <= -564.52)), fill = "#CCCCCC") +
  scale_fill_ramp_discrete(from = "darkgrey", guide = "none") +
  annotate(geom = "rect", xmin = -564.52, xmax = 564.52, ymin = -Inf, ymax = Inf, fill = "darkred", alpha = 0.3) +
  annotate(geom = "label", x = 0, y = 0.75, label = "ROPE") +
  labs(caption = "Point shows median value;\nthick black bar shows 66% credible interval;\nthin black bar shows 95% credible interval")

