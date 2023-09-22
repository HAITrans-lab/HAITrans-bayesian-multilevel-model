library(tidyverse)
#install.packages('ggpubr')
#install.packages('rstatix')
#install.packages('datarium')
#install.packages('report')
library(report)
library(ggpubr)
library(rstatix)
library(kableExtra)

# long format
set.seed(123) #reproducibility
eyetracking <- read.csv('/home/mrios/workspace/test_R/ANOVA_table3.xlsx - Sheet1.csv', header = TRUE, sep = ",")
eyetracking
eyetracking %>% sample_n_by(condition, text, size = 1)


#summary stats
stats <- eyetracking %>%
  group_by(condition, text) %>%
  get_summary_stats(quality_score, type = "mean_sd")
stats
write.csv(stats,'/home/mrios/workspace/test_R/quality/summary_stats.csv')
#viz
bxp <- ggboxplot(
  eyetracking, x = "text", y = "quality_score",
  facet.by = "condition", short.panel.labs = FALSE
)
bxp
ggsave("/home/mrios/workspace/test_R/quality/2wayaov_boxplot_quality.pdf")

results <- aov(quality_score ~ condition * text , data = eyetracking)
anova(results)
report(results)
report_aov <- report_text(results)
write.csv(report_aov,'/home/mrios/workspace/test_R/quality/report_aov.csv')


#twoway.model <- lm(quality_score ~ condition*text, data = eyetracking)
#summary(twoway.model)

pwc <- eyetracking %>%
  pairwise_t_test(
    quality_score ~ condition, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc
write.csv(pwc,'/home/mrios/workspace/test_R/quality/pairwise_ttest_condition.csv')

pwc <- eyetracking %>%
  pairwise_t_test(
    quality_score ~ text, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc
write.csv(pwc,'/home/mrios/workspace/test_R/quality/pairwise_ttest_text.csv')

linear_model <- lm(quality_score ~ text*condition, data=eyetracking) 
summary(linear_model)

pwc <- eyetracking %>%
  group_by(text) %>%
  emmeans_test(quality_score ~ condition, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
summary(pwc)
pwc %>%
  kbl(caption = "Pairwise comparisons") %>% 
  kable_material_dark()




pwc <- eyetracking %>%
  group_by(condition) %>%
  emmeans_test(quality_score ~ text, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
summary(pwc)
pwc
write.csv(pwc,'/home/mrios/workspace/test_R/quality/pairwise_ttest_conditiontext1.csv')
pwc %>%
  kbl(caption = "Pairwise comparisons") %>% 
  kable_material_dark()

#Pairwise comparisons
#Compare the different treatments by gender and risk variables:

pwc <- eyetracking %>%
  group_by(condition) %>%
  pairwise_t_test(
    quality_score ~ text, 
    p.adjust.method = "bonferroni"
  )
pwc
write.csv(pwc,'/home/mrios/workspace/test_R/quality/pairwise_ttest_conditiontext2.csv')

eyetracking
one.way <- eyetracking %>%
  group_by(condition) %>%
  anova_test(dv = quality_score, wid = participant, within = text) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way


#model  <- lm(quality_score ~ condition*text, data = eyetracking)
#treatment.effect <- eyetracking %>%
#  group_by(condition) %>%
#  anova_test(quality_score ~ text, error = model)

#treatment.effect %>% kbl(caption = "Treatment effect") %>% 
#  kable_material_dark()
