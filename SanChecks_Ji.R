###########################################################################
# Sanity checks: Evaluate relationship Ji, in-strength, ROI-size
# 
# Written by Hannelore Aerts 
# (UGent, Faculty of Psychology, Department of Data Analysis)
# Date last modification: 05/02/2018
###########################################################################

###----------------------------------------- Ji ~ in-strength checks -------------------------------------------#

### Ji corrected for in-strength 

Anova(aov(J_res_all$J_res ~ results$group), type=3) 
Anova(aov(J_res_all$JiNT_res ~ results$group), type=3) 
Anova(aov(J_res_all$JiT_res ~ results$group), type=3)
TukeyHSD(aov(J_res_all$JiT_res ~ results$group))
plot(J_res_all$JiT_res ~ results$group)

#---------------------------------------------------------------------------------------------------------------#

### In-strength determined by ROI size --> do group differences persist after regressing out 
# in-strength & ROI size (res2)?

Anova(aov(J_res_all2$J_res2 ~ results$group), type=3) 
Anova(aov(J_res_all2$JiNT_res2 ~ results$group), type=3) 
Anova(aov(J_res_all2$JiT_res2 ~ results$group), type=3)
TukeyHSD(aov(J_res_all2$JiT_res2 ~ results$group))
boxplot(J_res_all2$JiT_res2 ~ group)

# --> results similar as regression of in-strength, without taking into account ROI size

#---------------------------------------------------------------------------------------------------------------#

### Because of high correlation between in-strength and ROI size (r=0.80), better only regress out ROI size

Anova(aov(J_res_all3$J_res3 ~ results$group), type=3) 
Anova(aov(J_res_all3$JiNT_res3 ~ results$group), type=3) 
Anova(aov(J_res_all3$JiT_res3 ~ results$group), type=3)
TukeyHSD(aov(J_res_all3$JiT_res3 ~ results$group))
boxplot(J_res_all3$JiT_res3 ~ group)


