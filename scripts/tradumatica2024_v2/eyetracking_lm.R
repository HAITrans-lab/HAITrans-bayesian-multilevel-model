library(tidyverse)
#install.packages('ggpubr')
#install.packages('rstatix')
#install.packages('datarium')
#install.packages('report')
#install.packages("lme4")
remove.packages("Matrix")
remove.packages("lme4")
install.packages("lme4", type = "source")
#install.packages("sjPlot")
#install.packages("Matrix")
library(sjPlot)
library(lme4)
library(report)
library(ggpubr)
library(rstatix)

# long format
set.seed(123) #reproducibility
eyetracking <- read.csv('/home/mrios/workspace/imminent_R/prod_cogload_quality_results.csv', header = TRUE, sep = ",")
eyetracking
eyetracking %>% sample_n_by(condition, text, size = 1)

eyetracking <- eyetracking[!(is.na(eyetracking$pemt_speed)), ]
eyetracking <- eyetracking[!(is.na(eyetracking$PEMT_experience)), ]
eyetracking <- eyetracking[!(is.na(eyetracking$no_of_searches_in_external_resources)), ]
#summary stats
eyetracking %>%
  group_by(condition, text) %>%
  get_summary_stats(quality_score, type = "mean_sd")

#viz
bxp <- ggboxplot(
  eyetracking, x = "text", y = "quality_score",
  facet.by = "condition", short.panel.labs = FALSE
)
bxp
#ggsave("3wayaov_boxplot_eyetracking_quality.pdf")

# linear model
results <- lm(quality_score ~ 1 + condition + text + no_of_searches_in_external_resources + PEMT_experience, data = eyetracking)
summary(results)
tidy(results)
sjPlot::tab_model(results)




results <- lm(MFD_ST ~ 1 + condition + text + no_of_searches_in_external_resources + PEMT_experience, data = eyetracking)
summary(results)

tidy(results)
sjPlot::tab_model(results)

results <- lm(MFD_TT ~ 1 + condition + text + no_of_searches_in_external_resources + PEMT_experience, data = eyetracking)
summary(results)

tidy(results)
sjPlot::tab_model(results)

results <- lm(pemt_speed ~ 1 + condition + text + no_of_searches_in_external_resources + PEMT_experience, data = eyetracking)
summary(results)
sjPlot::tab_model(results)

#interactions
results <- lm(quality_score ~ 1 + condition * text, data = eyetracking)
summary(results)
tidy(results)
sjPlot::tab_model(results)

results <- lm(pemt_speed ~ 1 + condition * text, data = eyetracking)
summary(results)
tidy(results)
sjPlot::tab_model(results)

results <- lm(MFD_ST ~ 1 + condition * text, data = eyetracking)
summary(results)
tidy(results)
sjPlot::tab_model(results)

results <- lm(MFD_TT ~ 1 + condition * text, data = eyetracking)
summary(results)
tidy(results)
sjPlot::tab_model(results)

#linear mixed models

results <- lmer(quality_score ~ 1 + condition + text + no_of_searches_in_external_resources + PEMT_experience + (1 + condition | participant), data = eyetracking)
summary(results)
sjPlot::tab_model(results)



results <- lmer(quality_score ~ 1 + condition + text + no_of_searches_in_external_resources + PEMT_experience + (1 + condition  | participant), data = eyetracking)
summary(results)
sjPlot::tab_model(results)

results <- lmer(pemt_speed ~ 1 + condition + text + no_of_searches_in_external_resources + PEMT_experience + (1 + condition | participant), data = eyetracking)
summary(results)
sjPlot::tab_model(results)

results <- lmer(MFD_ST ~ 1 + condition + text + no_of_searches_in_external_resources + PEMT_experience + (1 + condition | participant), data = eyetracking)
summary(results)
sjPlot::tab_model(results)

results <- lmer(MFD_TT ~ 1 + condition + text + no_of_searches_in_external_resources + PEMT_experience + (1 + condition | participant), data = eyetracking)
summary(results)
sjPlot::tab_model(results)


