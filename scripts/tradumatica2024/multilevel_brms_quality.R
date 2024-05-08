library(brms)
#library(cmdstanr)
library(tidyverse)
library(sjPlot)
#install.packages("tidybayes")
#install.packages('devtools')
#devtools::install_github("mvuorre/brmstools")
#install.packages("bayestestR")
#install.packages("ggdist")
library(ggdist)
library("bayestestR")

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


set.seed(1234) #reproducibility
bayes_seed <- 1234
#set_cmdstan_path(path = '/home/mrios/.cmdstan/cmdstan-2.32.2')
eyetracking <- read.csv('/home/mrios/workspace/imminent_R/prod_cogload_quality_results.csv', header = TRUE, sep = ",")
eyetracking
eyetracking <- eyetracking[!(is.na(eyetracking$no_of_searches_in_external_resources)), ]

eyetracking %>% sample_n_by(condition, text, size = 1)
colnames(eyetracking)
#summary stats
stats <- eyetracking %>%
  group_by(condition, text) %>%
  get_summary_stats(quality_score, type = "mean_sd")
stats
write.csv(stats,'/home/mrios/workspace/imminent_R/qualityTradumatica2024/summary_stats.csv')


stats <- eyetracking %>%
  group_by(condition, bins_translation_experience) %>%
  get_summary_stats(quality_score, type = "mean_sd")
stats
write.csv(stats,'/home/mrios/workspace/imminent_R/qualityTradumatica2024/summary_stats-bins-trexp.csv')
stats <- eyetracking %>%
  group_by(condition, bins_PEMT_experience) %>%
  get_summary_stats(quality_score, type = "mean_sd")
stats
write.csv(stats,'/home/mrios/workspace/imminent_R/qualityTradumatica2024/summary_stats-bins-pemtexp.csv')


eyetracking$bins_PEMT_experience <- factor(eyetracking$bins_PEMT_experience , levels=c("y0", "y1", "y1_2", "y2_5", "y5_8", "y8_11"))
bxp <- ggboxplot(
  eyetracking, x = "bins_PEMT_experience", y = "quality_score",
  facet.by = "condition")  + theme( axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text = element_text(size = 16))
bxp
ggsave("/home/mrios/workspace/imminent_R/qualityTradumatica2024/boxplot_pemtspeed.png")

eyetracking$bins_translation_experience <- factor(eyetracking$bins_translation_experience , levels=c("y3_5", "y5_8", "y8_11", "yp11"))
bxp <- ggboxplot(
  eyetracking, x = "bins_translation_experience", y = "quality_score",
  facet.by = "condition")  + theme( axis.title.x = element_text(size = 20),
                                    axis.title.y = element_text(size = 20),
                                    axis.text = element_text(size = 16))
bxp
ggsave("/home/mrios/workspace/imminent_R/qualityTradumatica2024/boxplot_tr_pemtspeed.png")

#######################
#quality_score ~ 1 + condition + text + no_of_searches_in_external_resources + bins_translation_experience + bins_PEMT_experience + (1 + condition | participant)
# family = lognormal(),
#######################
fit0 <- brm(formula = ter_corpus_score ~  1 + condition + text + no_of_searches_in_external_resources + bins_translation_experience + bins_PEMT_experience + (1 + condition + text + no_of_searches_in_external_resources + bins_translation_experience + bins_PEMT_experience | participant), #(1 + condition | participant)
             data = eyetracking,
             warmup = 1000, iter = 10000, chains = 4, cores = 6,
             seed = bayes_seed
)
fit0
conditional_effects(fit0)
post_fit0 <- describe_posterior(fit0)
post_fit0


#######################
fit0 <- brm(formula = ter_corpus_score ~ 1 + bins_translation_experience + bins_PEMT_experience, #(1 + condition | participant)
             data = eyetracking,
             warmup = 1000, iter = 10000, chains = 4, cores = 6,
             seed = bayes_seed,
             backend = "cmdstanr"
)
fit0
conditional_effects(fit0)
post_fit0 <- describe_posterior(fit0)
post_fit0

#############################
fit0 <- brm(formula = pemt_speed ~ 1 + bins_translation_experience + bins_PEMT_experience, #(1 + condition | participant)
            data = eyetracking,
            warmup = 1000, iter = 10000, chains = 4, cores = 6,
            seed = bayes_seed,
            backend = "cmdstanr"
)
fit0
conditional_effects(fit0)
post_fit0 <- describe_posterior(fit0)
post_fit0

#######################
# quality_score ~ 1 + pemt_speed + (1 + pemt_speed | condition)
##########################
fit0 <- brm(formula = ter_corpus_score ~ 1 + pemt_speed + bins_translation_experience + bins_PEMT_experience + (1 + pemt_speed + bins_translation_experience + bins_PEMT_experience | condition), #(1 + condition | participant)
            data = eyetracking,
            warmup = 1000, iter = 10000, chains = 4, cores = 6,
            seed = bayes_seed,
            backend = "cmdstanr"
)
fit0
conditional_effects(fit0)
post_fit0 <- describe_posterior(fit0)
post_fit0
conditional_effects(fit0, effects="pemt_speed:condition", re_formula = NULL)
conditional_effects(fit0, effects="bins_translation_experience:condition", re_formula = NULL)
conditional_effects(fit0, effects="bins_PEMT_experience:condition", re_formula = NULL)

########################

fit0 <- brm(formula = ter_corpus_score ~ 1 + condition + bins_translation_experience + bins_PEMT_experience, # + (1 + bins_translation_experience + bins_PEMT_experience | condition), 
            data = eyetracking,
            warmup = 1000, iter = 10000, chains = 4, cores = 6,
            seed = bayes_seed,
            backend = "cmdstanr"
)
fit0
conditional_effects(fit0)
post_fit0 <- describe_posterior(fit0)
post_fit0

#conditional_effects(fit0, effects="bins_translation_experience:condition", re_formula = NULL)
#conditional_effects(fit0, effects="bins_PEMT_experience:condition", re_formula = NULL)


(cond_fit <- eyetracking %>%
    group_by(condition) %>%
    add_predicted_draws(fit0) %>%
    ggplot(aes(x = bins_PEMT_experience, y = ter_corpus_score, color = ordered(condition), fill = ordered(condition))) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95), alpha = 1/4) +
    geom_point(data = eyetracking) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    ylab("TER score\n") +
    xlab("\nbins_PEMT_experience") +
    theme_bw() +
    theme(legend.title = element_blank()))

(cond_fit <- eyetracking %>%
    group_by(condition) %>%
    add_predicted_draws(fit0) %>%
    ggplot(aes(x = bins_translation_experience, y = ter_corpus_score, color = ordered(condition), fill = ordered(condition))) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95), alpha = 1/4) +
    geom_point(data = eyetracking) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    ylab("TER score\n") +
    xlab("\nbins_translation_experience") +
    theme_bw() +
    theme(legend.title = element_blank()))
###############################################

fit0 <- brm(formula = ter_corpus_score ~ 1 + condition + text  + (1 + condition   | bins_PEMT_experience), 
            data = eyetracking,
            warmup = 1000, iter = 10000, chains = 4, cores = 6,
            seed = bayes_seed,
            backend = "cmdstanr"
)
fit0
conditional_effects(fit0)
conditional_effects(fit0, effects="condition:bins_PEMT_experience", re_formula = NULL)
forest(fit0,  pars = "conditions")
#forest(fit0,  pars = "textt2")

fit0 <- brm(formula = ter_corpus_score ~ 1 + condition + text  + (1 + condition + text  | bins_translation_experience), 
            data = eyetracking,
            warmup = 1000, iter = 10000, chains = 4, cores = 6,
            seed = bayes_seed,
            backend = "cmdstanr"
)
fit0
conditional_effects(fit0)
conditional_effects(fit0, effects="condition:bins_translation_experience", re_formula = NULL)

##########################


fit0a <- brm(formula = quality_score ~ 1 + bins_translation_experience + bins_PEMT_experience, #(1 + condition | participant)
             data = eyetracking,
             warmup = 1000, iter = 10000, chains = 4, cores = 6,
             seed = bayes_seed,
             backend = "cmdstanr"
)
fit0a
conditional_effects(fit0a)
post_fit0a <- describe_posterior(fit0a)
post_fit0a


fit0a <- brm(formula = MFD_ST ~ 1 + bins_translation_experience + bins_PEMT_experience, #(1 + condition | participant)
             data = eyetracking,
             warmup = 1000, iter = 10000, chains = 4, cores = 6,
             seed = bayes_seed,
             backend = "cmdstanr"
)
fit0a
conditional_effects(fit0a)
#plot(fit0a)

fit0a <- brm(formula = MFD_TT ~ 1 + bins_translation_experience + bins_PEMT_experience, #(1 + condition | participant)
             data = eyetracking,
             warmup = 1000, iter = 10000, chains = 4, cores = 6,
             seed = bayes_seed,
             backend = "cmdstanr"
)
fit0a
conditional_effects(fit0a)
#plot(fit0a)


fit0b <- brm(formula = quality_score ~ 1 + condition + text + no_of_searches_in_external_resources + bins_translation_experience + bins_PEMT_experience, #(1 + condition | participant)
            data = eyetracking,
            warmup = 1000, iter = 10000, chains = 4, cores = 6,
            seed = bayes_seed,
            backend = "cmdstanr"
)
fit0b
text_summ <-summary(fit0b)
sink("/home/mrios/workspace/imminent_R/qualityTradumatica2024/multilevel_brms0b_summary.txt")
text_summ
#tidy(fit1)
sink()
#sjPlot::tab_model(fit1)
fit0b
p_summary <- posterior_summary(fit0b)
write.csv(p_summary, '/home/mrios/workspace/imminent_R/qualityTradumatica2024/multilevel_brms0b_psummary.csv')

plot(fit0b)
pp_check(fit0b)

fixef(fit0b)
#coef(fit0)  #coef(fit1)$condition TODO
#coef(fit1)$text[conditions]

conditional_effects(fit0b)




tidy(fit0b)



post_fit0b <- describe_posterior(fit0b)
post_fit0b
write.csv(post_fit0b, '/home/mrios/workspace/imminent_R/qualityTradumatica2024/multilevel_brms0b_describepost.csv')
bayestestR::hdi(fit0b)

####
#model participant
#quality_score ~ 1 + condition + text + (1 + condition | participant)
####
fit2 <- brm(formula = quality_score ~ 1 + condition + text + no_of_searches_in_external_resources + bins_translation_experience + bins_PEMT_experience + (1 + condition | participant), 
            data = eyetracking,
            warmup = 1000, iter = 10000, chains = 4, cores = 3,
            control=list(adapt_delta=0.9), seed = bayes_seed,
            )
fit2
#tab_model(fit2)
text_summ <-summary(fit2)
sink("/home/mrios/workspace/imminent_R/qualityTradumatica2024/multilevel_brms2_summary.txt")
text_summ
sink()
#print(fit2)
#sjPlot::tab_model(fit2)
fit2
p_summary <- posterior_summary(fit2)
p_summary
write.csv(p_summary, '/home/mrios/workspace/imminent_R/qualityTradumatica2024/multilevel_brms2_psummary.csv')
#write.csv(text_summ, '/home/mrios/workspace/imminent_R/quality/multilevel_brms1_summary.txt')
randeff <- ranef(fit2)
randeff
write.csv(randeff, '/home/mrios/workspace/imminent_R/qualityTradumatica2024/multilevel_brms2_randeff.csv')
plot(fit2)
pp_check(fit2)

fixef(fit2)
coef(fit2) 
extract_random_effects(fit2)


#conditional_effects(fit2)
cond <- conditional_effects(fit2)
plot(cond, theme =  theme( axis.title.x = element_text(size = 20),
                           axis.title.y = element_text(size = 20),
                           axis.text = element_text(size = 16)))

loo2 <- loo(fit2)

#loo_cpm <- loo_compare(loo1, loo2)
#loo_cpm
#write.csv(loo_cpm, '/home/mrios/workspace/imminent_R/qualityTradumatica2024/multilevel_brms_loo_cpm1-2.csv')
#loo_compare(loo0, loo0b, loo1, loo2)




#participant level random effects
condition_participant_offsets <- ranef(fit2)$participant %>%
  as_tibble(rownames = "participant") %>% 
  filter(participant %in% eyetracking$participant) %>% 
  select(participant, starts_with("Estimate"))
condition_participant_offsets
write.csv(condition_participant_offsets, '/home/mrios/workspace/imminent_R/qualityTradumatica2024/multilevel_brms2_condoffsets.csv')
#fixed effect + random offset for text-specific intercepts and slopes.
condition_participant_raneff <- coef(fit2)$participant %>%
  as_tibble(rownames = "participant") %>% 
  filter(participant %in% eyetracking$participant) %>% 
  select(participant, starts_with("Estimate"))
print(condition_participant_raneff, n=21)
forest(fit2, pars='conditions')


post_fit2 <- describe_posterior(fit2)
post_fit2
write.csv(post_fit2, '/home/mrios/workspace/imminent_R/qualityTradumatica2024/multilevel_brms2_describepost.csv')
bayestestR::hdi(fit2)

#TODO??? rstan not cmdstan save_pars = save_pars(all = TRUE)
#bayesfactor_models(fit1, fit2, denominator = 1, verbose = FALSE) 

#ROPE
tidy(fit2)
condition_draws <- fit2 %>%
  spread_draws(b_conditions)

condition_cred_int <- condition_draws %>% 
  median_hdi()
condition_cred_int$.lower


condition_draws %>% 
  # Only look inside the credible interval
  filter(b_conditions >= condition_cred_int$.lower & b_conditions <= condition_cred_int$.upper) %>% 
  summarize(prop_outside_rope = 1 - sum(b_conditions >= -0.66 & b_conditions <= 0.66) / n())

ggplot(condition_draws, aes(x = b_conditions)) +
  stat_halfeye(aes(fill_ramp = stat(x >= 0.66 | x <= -0.66)), fill = "#CCCCCC") +
  scale_fill_ramp_discrete(from = "darkgrey", guide = "none") +
  annotate(geom = "rect", xmin = -0.66, xmax = 0.66, ymin = -Inf, ymax = Inf, fill = "darkred", alpha = 0.3) +
  annotate(geom = "label", x = 0, y = 0.75, label = "ROPE") +
  labs(caption = "Point shows median value;\nthick black bar shows 66% credible interval;\nthin black bar shows 95% credible interval")
#p direction
condition_draws %>% 
  summarize(prop_greater_0 = sum(b_conditions > 0) / n())

(cond_fit <- eyetracking %>%
    group_by(condition) %>%
    add_predicted_draws(fit2) %>%
    ggplot(aes(x = text, y = quality_score, color = ordered(condition), fill = ordered(condition))) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95), alpha = 1/4) +
    geom_point(data = eyetracking) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    ylab("quality_score\n") +
    xlab("\ntext") +
    theme_bw() +
    theme(legend.title = element_blank()))


sessionInfo()
