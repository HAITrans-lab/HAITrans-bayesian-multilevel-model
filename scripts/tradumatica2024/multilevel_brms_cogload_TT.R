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
eyetracking <- read.csv('/home/mrios/workspace/imminent_R/prod_cogload_quality_results.csv', header = TRUE, sep = ",")
eyetracking
eyetracking %>% sample_n_by(condition, text, size = 1)
#delete NA
#eyetracking <- na.omit(eyetracking)
eyetracking <- eyetracking[!(is.na(eyetracking$MFD_TT)), ]
eyetracking <- eyetracking[!(is.na(eyetracking$PEMT_experience)), ]
eyetracking <- eyetracking[!(is.na(eyetracking$no_of_searches_in_external_resources)), ]
colnames(eyetracking)
#summary stats
stats <- eyetracking %>%
  group_by(condition, text) %>%
  get_summary_stats(MFD_TT, type = "mean_sd")
stats
write.csv(stats,'/home/mrios/workspace/imminent_R/cog_load-TTEAMT024/summary_stats_tt.csv')

bxp <- ggboxplot(
  eyetracking, x = "text", y = "MFD_TT",
  facet.by = "condition", short.panel.labs = FALSE
)
bxp


###
#MFD_ST ~ 1 + condition + text 
###

fit0b <- brm(formula = MFD_TT ~ 1 + condition + text + no_of_searches_in_external_resources + PEMT_experience, #(1 + condition | participant)
            data = eyetracking,
            warmup = 1000, iter = 5000, chains = 4, cores = 6,
            seed = bayes_seed,
            backend = "cmdstanr"
)
fit0b
text_summ <-summary(fit0b)
sink("/home/mrios/workspace/imminent_R/cog_load-TTEAMT024/multilevel_brms0b_summary.txt")
text_summ
#tidy(fit1)
sink()
#sjPlot::tab_model(fit1)
fit0b
p_summary <- posterior_summary(fit0b)
write.csv(p_summary, '/home/mrios/workspace/imminent_R/cog_load-TTEAMT024/multilevel_brms0b_psummary.csv')


pp_check(fit0b)

fixef(fit0b)
#coef(fit0)  #coef(fit1)$condition TODO
#coef(fit1)$text[conditions]

conditional_effects(fit0b)



describe_posterior(fit0b)
post_fit0b <- describe_posterior(fit0b)
post_fit0b
write.csv(post_fit0b, '/home/mrios/workspace/imminent_R/cog_load-TTEAMT024/multilevel_brms0b_describepost.csv')
bayestestR::hdi(fit0b)

####
#model participant
#MFD_ST ~ 1 + condition + text + (1 + condition | participant)
####
fit2 <- brm(formula = MFD_TT ~ 1 + condition + text  + no_of_searches_in_external_resources + PEMT_experience + (1 + condition | participant), 
            data = eyetracking,
            warmup = 1000, iter = 10000, chains = 4, cores = 6,
            control=list(adapt_delta=0.9), seed = bayes_seed,
            backend = "cmdstanr"
            )

fit2

text_summ <-summary(fit2)
sink("/home/mrios/workspace/imminent_R/cog_load-TTEAMT024/multilevel_brms2_summary.txt")
text_summ
sink()
#print(fit2)
#sjPlot::tab_model(fit2)
fit2
p_summary <- posterior_summary(fit2)
p_summary
write.csv(p_summary, '/home/mrios/workspace/imminent_R/cog_load-TTEAMT024/multilevel_brms2_psummary.csv')
#write.csv(text_summ, '/home/mrios/workspace/imminent_R/quality/multilevel_brms1_summary.txt')
randeff <- ranef(fit2)
randeff
write.csv(randeff, '/home/mrios/workspace/imminent_R/cog_load-TTEAMT024/multilevel_brms2_randeff.csv')
#plot(fit2)
pp_check(fit2)

fixef(fit2)
coef(fit2) 
extract_random_effects(fit2)
print(extract_random_effects(fit2), n=36)

conditional_effects(fit2)




#participant level random effects
condition_participant_offsets <- ranef(fit2)$participant %>%
  as_tibble(rownames = "participant") %>% 
  filter(participant %in% eyetracking$participant) %>% 
  select(participant, starts_with("Estimate"))
condition_participant_offsets
write.csv(condition_participant_offsets, '/home/mrios/workspace/imminent_R/cog_load-TTEAMT024/multilevel_brms2_condoffsets.csv')
#fixed effect + random offset for text-specific intercepts and slopes.
coef(fit2)$participant %>%
  as_tibble(rownames = "participant") %>% 
  filter(participant %in% eyetracking$participant) %>% 
  select(participant, starts_with("Estimate"))

forest(fit2, pars='conditions')

tidy(fit2)
describe_posterior(fit2)
post_fit2 <- describe_posterior(fit2)
post_fit2
write.csv(post_fit2, '/home/mrios/workspace/imminent_R/cog_load-TTEAMT024/multilevel_brms2_describepost.csv')
bayestestR::hdi(fit2)
#-6.28, 6.28
condition_draws <- fit2 %>%
  spread_draws(b_conditions)

condition_cred_int <- condition_draws %>% 
  median_hdi()
condition_cred_int$.lower


condition_draws %>% 
  # Only look inside the credible interval
  filter(b_conditions >= condition_cred_int$.lower & b_conditions <= condition_cred_int$.upper) %>% 
  summarize(prop_outside_rope = 1 - sum(b_conditions >= -6.14 & b_conditions <= 6.14) / n())

ggplot(condition_draws, aes(x = b_conditions)) +
  stat_halfeye(aes(fill_ramp = stat(x >= 6.14 | x <= -6.14)), fill = "#CCCCCC") +
  scale_fill_ramp_discrete(from = "darkgrey", guide = "none") +
  annotate(geom = "rect", xmin = -6.14, xmax = 6.14, ymin = -Inf, ymax = Inf, fill = "darkred", alpha = 0.3) +
  annotate(geom = "label", x = 0, y = 0.75, label = "ROPE") +
  labs(caption = "Point shows median value;\nthick black bar shows 66% credible interval;\nthin black bar shows 95% credible interval")


(cond_fit <- eyetracking %>%
    group_by(condition) %>%
    add_predicted_draws(fit2) %>%
    ggplot(aes(x = text, y = MFD_TT, color = ordered(condition), fill = ordered(condition))) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95), alpha = 1/4) +
    geom_point(data = eyetracking) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    ylab("MFD_TT\n") +
    xlab("\ntext") +
    theme_bw() +
    theme(legend.title = element_blank()))
