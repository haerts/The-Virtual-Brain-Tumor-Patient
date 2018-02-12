###########################################################################
# Sanity checks: Evaluate model fit
# 
# Written by Hannelore Aerts 
# (UGent, Faculty of Psychology, Department of Data Analysis)
# Date last modification: 12/02/2018
###########################################################################

ModelFit_SanChecks$Group = factor(ModelFit_SanChecks$Group)
attach(ModelFit_SanChecks)

myGroups=levels(Group)
myProps=summary(Group)/nrow(ModelFit_SanChecks)
jpeg('ModelFit_SanChecks.jpg', width=800, height=600)
par(mar=c(5,6,4,1)+.1)
plot(Cor ~ Group, data=ModelFit_SanChecks, col=palette, xlab="", outline=FALSE, 
     ylim=c(0.09, 0.55), ylab=expression("Pearson correlation FC" [sim]* "- FC" [emp]),
     cex.axis=1.5, cex.lab=1.8)
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=ModelFit_SanChecks[ModelFit_SanChecks$Group==thisGroup, "Cor"]
  myjitter=jitter(rep(i, dim(thisValues)[1]), amount=myProps[i]/2)
  points(myjitter, thisValues$Cor, pch=16, col=rgb(0,0,0,0.4), cex=1.8)
}
dev.off()

# Test significance
Anova(aov(Cor~Group), type=3) #F(3,140)=36.3399, p=0.00046
TukeyHSD(aov(Cor~Group))

#  Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = Cor ~ Group)
#
#$Group
#                                                 diff          lwr        upr     p adj
#CON avg SC-CON avg params               -0.002771855 -0.054093253 0.04854954 0.9990074
#CON avg SC + PSE-CON avg params          0.055694204  0.004372806 0.10701560 0.0277227
#Individual params & SC-CON avg params    0.062791218  0.011469820 0.11411262 0.0096652
#CON avg SC + PSE-CON avg SC              0.058466059  0.007144660 0.10978746 0.0186260
#Individual params & SC-CON avg SC        0.065563073  0.014241674 0.11688447 0.0062102
#Individual params & SC-CON avg SC + PSE  0.007097014 -0.044224384 0.05841841 0.9840204
