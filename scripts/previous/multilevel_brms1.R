library(brms)
library(tidyverse)
library(sjPlot)
#install.packages("tidybayes")
library(tidybayes)
library(sjstats)
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
#library(mixedup)

set.seed(1234) #reproducibility
bayes_seed <- 1234
eyetracking <- read.csv('/home/mrios/workspace/test_R/ANOVA_table3.xlsx - Sheet1.csv', header = TRUE, sep = ",")
eyetracking
eyetracking %>% sample_n_by(condition, text, size = 1)

#summary stats
stats <- eyetracking %>%
  group_by(condition, text) %>%
  get_summary_stats(quality_score, type = "mean_sd")
stats
write.csv(stats,'/home/mrios/workspace/test_R/quality/summary_stats.csv')

####
#model condition
####

fit0 <- brm(formula = quality_score ~ 1 + condition, #(1 + condition | participant)
            data = eyetracking,
            warmup = 1000, iter = 5000, chains = 4, cores = 6,
            seed = bayes_seed
)
fit0
text_summ <-summary(fit0)
sink("/home/mrios/workspace/test_R/quality/multilevel_brms0_summary.txt")
text_summ
#tidy(fit0)
sink()
#sjPlot::tab_model(fit1)
fit0
p_summary <- posterior_summary(fit0)
write.csv(p_summary, '/home/mrios/workspace/test_R/quality/multilevel_brms0_psummary.csv')

plot(fit0)
pp_check(fit0)

fixef(fit0)
#coef(fit0)  #coef(fit1)$condition TODO
png(file="/home/mrios/workspace/test_R/quality/fit0.png",
    width=600, height=350)
conditional_effects(fit0)
dev.off()


loo0 <- loo(fit0) 
loo0

#(model_fit <- eyetracking %>%
#    add_predicted_draws(fit0) %>%  # adding the posterior distribution
#    ggplot(aes(x = text, y = quality_score)) +  
#    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
#                    alpha = 0.5, colour = "black") +
 #   geom_point(data = eyetracking, colour = "darkseagreen4", size = 3) +   # raw data
#    scale_fill_brewer(palette = "Greys") +
 #   ylab("quality_score\n") +  # latin name for red knot
  #  xlab("\ntext") +
  #  theme_bw() +
   # theme(legend.title = element_blank(),
  #        legend.position = c(0.15, 0.85)))


####

fit1 <- brm(formula = quality_score ~ 1 + condition + text + (1 + condition | text), #(1 + condition | participant)
            data = eyetracking,
            warmup = 1000, iter = 5000, chains = 4, cores = 6,
            control=list(adapt_delta=0.99), seed = bayes_seed
            )
fit1
text_summ <-summary(fit1)
sink("/home/mrios/workspace/test_R/quality/multilevel_brms1_summary.txt")
text_summ
#tidy(fit1)
sink()
#sjPlot::tab_model(fit1)
fit1
p_summary <- posterior_summary(fit1)
write.csv(p_summary, '/home/mrios/workspace/test_R/quality/multilevel_brms1_psummary.csv')
randeff <- ranef(fit1)
randeff
write.csv(randeff, '/home/mrios/workspace/test_R/quality/multilevel_brms1_randeff.csv')

plot(fit1)
pp_check(fit1)

fixef(fit1)
coef(fit1)  #coef(fit1)$condition TODO
extract_random_effects(fit1)

png(file="/home/mrios/workspace/test_R/quality/fit1.png",
    width=600, height=350)
conditional_effects(fit1)
dev.off()


loo1 <- loo(fit1) 
loo1

(model_fit <- eyetracking %>%
    add_predicted_draws(fit1) %>%  # adding the posterior distribution
    ggplot(aes(x = text, y = quality_score)) +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = eyetracking, colour = "darkseagreen4", size = 3) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ylab("quality_score\n") +  # latin name for red knot
    xlab("\ntext") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85)))

(model_fit <- eyetracking %>%
    add_predicted_draws(fit1) %>%  # adding the posterior distribution
    ggplot(aes(x = condition, y = quality_score)) +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = eyetracking, colour = "darkseagreen4", size = 3) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ylab("quality_score\n") +  # latin name for red knot
    xlab("\ncondition") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85)))

####
#model participant
####
fit2 <- brm(formula = quality_score ~ 1 + condition + text + (1 + condition | participant), 
            data = eyetracking,
            warmup = 1000, iter = 5000, chains = 4, cores = 6,
            control=list(adapt_delta=0.99), seed = bayes_seed
)
fit2
text_summ <-summary(fit2)
sink("/home/mrios/workspace/test_R/quality/multilevel_brms2_summary.txt")
text_summ
sink()
#print(fit2)
#sjPlot::tab_model(fit2)
fit2
p_summary <- posterior_summary(fit2)
p_summary
write.csv(p_summary, '/home/mrios/workspace/test_R/quality/multilevel_brms2_psummary.csv')
#write.csv(text_summ, '/home/mrios/workspace/test_R/quality/multilevel_brms1_summary.txt')
randeff <- ranef(fit2)
randeff
write.csv(randeff, '/home/mrios/workspace/test_R/quality/multilevel_brms2_randeff.csv')
plot(fit2)
pp_check(fit2)

fixef(fit2)
coef(fit2) 
extract_random_effects(fit2)
conditional_effects(fit2)
loo2 <- loo(fit2)


loo_compare(loo0, loo1, loo2)


(model_fit <- eyetracking %>%
    add_predicted_draws(fit2) %>%  # adding the posterior distribution
    ggplot(aes(x = text, y = quality_score)) +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = eyetracking, colour = "darkseagreen4", size = 3) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ylab("quality_score\n") +  # latin name for red knot
    xlab("\ntext") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85)))

(model_fit <- eyetracking %>%
    add_predicted_draws(fit2) %>%  # adding the posterior distribution
    ggplot(aes(x = condition, y = quality_score)) +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = eyetracking, colour = "darkseagreen4", size = 3) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ylab("quality_score\n") +  # latin name for red knot
    xlab("\ncondition") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85)))

#####
#condtion * text
#####

#fit3 <- brm(formula = quality_score ~ 1 + condition + text + condition*text, 
#            data = eyetracking,
#            warmup = 1000, iter = 10000, chains = 4, cores = 10,
#            control=list(adapt_delta = 0.999), seed = bayes_seed
#            )
#fit3
#conditional_effects(fit3)

