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
eyetracking <- read.csv('/home/mrios/workspace/test_R/big_data/Big_data_en-de.csv', header = TRUE, sep = ",")
eyetracking
eyetracking %>% sample_n_by(condition, text, size = 1)

colnames(eyetracking)
#summary stats
stats <- eyetracking %>%
  group_by(condition, text) %>%
  get_summary_stats(productivity_HAITrans_PE_speed, type = "mean_sd")
stats
write.csv(stats,'/home/mrios/workspace/test_R/big_data/summary_stats.csv')
#viz
bxp <- ggboxplot(
  eyetracking, x = "text", y = "productivity_HAITrans_PE_speed",
  facet.by = "condition", short.panel.labs = FALSE
)
bxp
ggsave("/home/mrios/workspace/test_R/big_data/2wayaov_boxplot_mfdtt_en-de.pdf")

results <- aov(productivity_HAITrans_PE_speed ~ condition * text , data = eyetracking)
anova(results)
report(results)
report_aov <- report_text(results)
write.csv(report_aov,'/home/mrios/workspace/test_R/big_data/report_aov.csv')


#twoway.model <- lm(quality_score ~ condition*text, data = eyetracking)
#summary(twoway.model)

pwc <- eyetracking %>%
  pairwise_t_test(
    productivity_HAITrans_PE_speed ~ condition, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc
write.csv(pwc,'/home/mrios/workspace/test_R/big_data/pairwise_ttest_condition.csv')

pwc <- eyetracking %>%
  pairwise_t_test(
    productivity_HAITrans_PE_speed ~ text, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc
write.csv(pwc,'/home/mrios/workspace/test_R/big_data/pairwise_ttest_text.csv')

linear_model <- lm(productivity_HAITrans_PE_speed ~ text*condition, data=eyetracking) 
summary(linear_model)

pwc <- eyetracking %>%
  group_by(text) %>%
  emmeans_test(productivity_HAITrans_PE_speed ~ condition, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
summary(pwc)
pwc %>%
  kbl(caption = "Pairwise comparisons") %>% 
  kable_material_dark()




pwc <- eyetracking %>%
  group_by(condition) %>%
  emmeans_test(productivity_HAITrans_PE_speed ~ text, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
summary(pwc)
pwc
write.csv(pwc,'/home/mrios/workspace/test_R/big_data/pairwise_ttest_conditiontext1.csv')
pwc %>%
  kbl(caption = "Pairwise comparisons") %>% 
  kable_material_dark()

#Pairwise comparisons
#Compare the different treatments by gender and risk variables:

pwc <- eyetracking %>%
  group_by(condition) %>%
  pairwise_t_test(
    productivity_HAITrans_PE_speed ~ text, 
    p.adjust.method = "bonferroni"
  )
pwc
write.csv(pwc,'/home/mrios/workspace/test_R/big_data/pairwise_ttest_conditiontext2.csv')

eyetracking
one.way <- eyetracking %>%
  group_by(condition) %>%
  anova_test(dv = productivity_HAITrans_PE_speed, wid = participant, within = text) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way


#model  <- lm(quality_score ~ condition*text, data = eyetracking)
#treatment.effect <- eyetracking %>%
#  group_by(condition) %>%
#  anova_test(quality_score ~ text, error = model)

#treatment.effect %>% kbl(caption = "Treatment effect") %>% 
#  kable_material_dark()
