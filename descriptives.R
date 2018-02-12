###################################################################################################################
###                                                                                                             ###
###                 MODELING BRAIN DYNAMICS IN BRAIN TUMOR PATIENTS USING THE VIRTUAL BRAIN                     ###
###               ==============================================================================                ###
###                                       PART 1: DESCRIPTIVE ANALYES                                           ###
###                                                                                                             ###
### Created by Hannelore Aerts                                                                                  ###
### Date last update: 09/02/2018                                                                                ###
### --> use Ji corrected for ROI size only, not ROI size + in-strength                                          ###
###################################################################################################################

library(corrplot)

### Read in data and prepare for analyses ------------------------------------------------------------------------#

setwd("/home/hannelore/Documents/ANALYSES/TVB_global2")
results=read.table(file="RESULTS_ALL.csv", header=TRUE, sep=",")

str(results)
results$subID = as.character(results$subID)
results$sexbin = c(NA*length(results))
results$sexbin[results$sex == "M"] = 0
results$sexbin[results$sex == "F"] = 1
results$group = factor(results$group, ordered=TRUE, levels=c('CON', 'MEN', 'GLI'))
results$lesion_lat = factor(results$lesion_lat, ordered=TRUE, levels=c('None', 'Left', 'Right', 'Bilateral'))
results$fmri_TR = as.factor(results$fmri_TR)
attach(results)


###################################################################################################################
###                                   PART 1: DEMOGRAPHICS DESCRIPTIVES                                         ###
###################################################################################################################

setwd("/home/hannelore/Documents/ANALYSES/TVB_global2/results_descriptives/")
myGroups=levels(group)
myProps=summary(group)/nrow(results)

### Age 
jpeg('age.jpg')
plot(age ~ group, data=results, col=viridis(5), xlab="")
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "age"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=1.2)
}
dev.off()
aggregate(age, list(group), mean)
vars=aggregate(age, list(group), var)
sqrt(vars$x)


### Sex
jpeg('sex.jpg')
plot(sex ~ group, data=results, xlab="",col=viridis(3))
dev.off()
table(sex, group)

### STAI
jpeg('stai.jpg')
plot(STAI ~ group, data=results, xlab="", ylab="state anxiety", col=viridis(5),
     outline=FALSE, ylim=c(20, 80))
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "STAI"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=1.2)
}
dev.off()

### Handedness
table(hand, group)


results_tumor = results[12:36,] #select only patients
results_tumor$group = factor(results_tumor$group, ordered=TRUE, levels=c('MEN', 'GLI'))
results_tumor$lesion_lat = factor(results_tumor$lesion_lat, ordered=TRUE, levels=c('Left', 'Right', 'Bilateral'))
myTumorGroups=levels(results_tumor$group)
myTumorProps=summary(results_tumor$group)/nrow(results_tumor)

### Tumor volume
jpeg('tumor_volume.jpg')
boxplot(results_tumor$lesion_vol_cm3 ~ results_tumor$group,
        col=viridis(5), ylab="tumor volume (cm³)", outline=FALSE, ylim=c(0,90))
for (i in 1:length(myTumorGroups)){
  thisGroup=myTumorGroups[i]
  thisValues=results_tumor[results_tumor$group==thisGroup, "lesion_vol_cm3"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myTumorProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=1.2)
}
dev.off()

### Tumor laterality
jpeg('tumor_laterality.jpg')
plot(results_tumor$lesion_lat ~ results_tumor$group, xlab='', 
     ylab='lesion laterality', col=viridis(5))
dev.off()
rm(list=c('results_tumor', 'myTumorGroups', 'myTumorProps'))




###################################################################################################################
###                                         PART 2: CANTAB DESCRIPTIVES                                         ###
###################################################################################################################

### Cantab distributions -----------------------------------------------------------------------------------------#

setwd("/home/hannelore/Documents/ANALYSES/TVB_global2/results_cantab/")

### MOT
hist(MOT_latency_mean, col="gray", breaks=10)
plot(MOT_latency_mean ~ group, col="gray")
#--> use inverse as "motivation"

results$motivation = (1 / results$MOT_latency_mean) *1000
jpeg('motivation.jpg')
plot(results$motivation ~ group, col=viridis(5), xlab="", ylab="motivation")
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "motivation"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=1.2)
}
dev.off()

### RTI
jpeg('RTI_fiveRT.jpg')
plot(RTI_fiveRT ~ group, col="gray", xlab="")
points(group, RTI_fiveRT)
dev.off()
plot(RTI_fiveMT, RTI_fiveRT)  

### RVP
jpeg('RVP_A.jpg')
plot(RVP_A ~ group, col="gray", xlab="")
points(group, RVP_A)
dev.off()

### SOC
jpeg('SOC_prob_minmoves.jpg')
plot(SOC_prob_minmoves ~ group, col="gray", xlab="")
points(group, SOC_prob_minmoves)
dev.off()

### SSP
jpeg('SSP_spanlength.jpg')
plot(SSP_spanlength ~ group, col="gray", xlab="")
points(group, SSP_spanlength)
dev.off()


### Correlations Cantab - motivation -----------------------------------------------------------------------------#

attach(results)

pairs(motivation ~ RTI_fiveRT + RVP_A + SOC_prob_minmoves + SSP_spanlength)
jpeg('motivation_corrplot.jpg', width=600, height=480)
corrplot(cor(results[,c(77,36,41,46,51)], use="pairwise"), type='lower', method="color", diag=F)
dev.off()


### Correlations Cantab - STAI ----------------------------------------------------------------------------------#

pairs(STAI ~ RTI_fiveRT + RVP_A + SOC_prob_minmoves + SSP_spanlength)
jpeg('STAI_corrplot.jpg', width=600, height=480)
corrplot(cor(results[,c(52,36,41,46,51)], use="pairwise"), type='lower', method="color", diag=F)
dev.off()


### Correlation motivation - STAI -------------------------------------------------------------------------------#

jpeg('STAI_motivation.jpg')
plot(motivation ~ STAI, pch=19)
legend(70, 0.7, round(cor(motivation, STAI, use = "pairwise"),2))
dev.off()
#--> 2 seperate things, both related to performance on cantab!!


### Correlation cantab - covariates -----------------------------------------------------------------------------#
corrplot(cor(results[,c(42,47,52,57,2,4,7:9,58,83)], use="pairwise.complete.obs"), 
         type='lower', diag=F, method='number')
  #--> important continuous covariates: age, lesion volume, STAI & motivation



###################################################################################################################
###                                     PART 3: GTA RESULTS DESCRIPTIVES                                        ###
###################################################################################################################

setwd("/home/hannelore/Documents/ANALYSES/TVB_global2/results_GTA/")

summary(results[,59:69])


### ABSOLUTE THRESHOLD

# Pairwise correlations
pairs(~ SC_thrA_density + SC_thrA_clust + SC_thrA_Eloc + SC_thrA_Q + SC_thrA_Eglob + SC_thrA_comm 
      + SC_thrA_degree + SC_thrA_strength + SC_thrA_EBC + SC_thrA_BC + SC_thrA_PC)

jpeg('SC_thrA_corplot.jpg', width = 600, height = 480)
corrplot(cor(results[,59:69]), type='lower', diag=FALSE, method='number')
dev.off()
cor(results[,59:69])

# Remove variables with correlations > 0.80:
#Density & degree highly related --> remove degree
#Clust & Eloc highly related --> remove clust
#EBC & BC highly related --> remove EBC
#Eglob & comm highly related --> remove comm
#Eglob & strength highly related --> remove strength
#Eloc & Eglob highly related --> remove Eloc

corrplot(cor(results[,c(59,62,63,68,69)]), type='lower', diag=FALSE, method='color')


# PCA
palette(c("blue4", "firebrick"))
results_SC_thrA = results[,c(59:69)]
colnames(results_SC_thrA)=
  c("density", "clustering", "Eloc", "Q", "Eglob", "communicability", "degree", "strength", "EBC", "BC", "PC")
SC_thrA_pca = prcomp(results_SC_thrA, center=TRUE, scale=TRUE)
jpeg('SC_thrA_pca1.jpg')
biplot(SC_thrA_pca, xlabs=rep("*", nrow(results)), cex=1.3)
dev.off()

results_SC_thrA = results[,c(62,63,69)]
colnames(results_SC_thrA)=c("Q", "Eglob", "PC")
SC_thrA_pca = prcomp(results_SC_thrA, center=TRUE, scale=TRUE)
jpeg('SC_thrA_pca2.jpg')
biplot(SC_thrA_pca, xlabs=rep("*", nrow(results)), cex=1.3)
dev.off()


#--> FINAL: Q, Eglob, PC (most segregated in PCA space)


# Check group differences
jpeg('SC_thrA_Q.jpg')
boxplot(SC_thrA_Q ~ group, col="gray", ylab='Maximized modularity (Q)')
points(group, SC_thrA_Q)
dev.off()
hist(SC_thrA_Q)
Q_groupdif = aov(SC_thrA_Q ~ group)
summary(Q_groupdif)

jpeg('SC_thrA_Eglob.jpg')
boxplot(SC_thrA_Eglob ~ group, col="gray", ylab='Global efficiency')
points(group, SC_thrA_Eglob)
dev.off()
hist(SC_thrA_Eglob)
Eglob_groupdif = aov(SC_thrA_Eglob ~ group)
summary(Eglob_groupdif) 

jpeg('SC_thrA_PC.jpg')
boxplot(SC_thrA_PC ~ group, col="gray", ylab="Participation coefficient")
points(group, SC_thrA_PC)
dev.off()
hist(SC_thrA_PC)
PC_groupdif = aov(SC_thrA_PC ~ group)
summary(PC_groupdif) 

rm(list=c('Eglob_groupdif', 'PC_groupdif', 'Q_groupdif', 'SC_thrA_pca', 'results_SC_thrA'))


# Important confounders to correct for:
corrplot(cor(results[,c(62,63,69,2,4,7:9,58)], use="pairwise.complete.obs"), 
         type='lower', diag=F, method='number')
#--> important continuous covariates: age, STAI 



###################################################################################################################
###                                  PART 4: MODELING PARAMETERS DESCRIPTIVES                                   ###
###################################################################################################################

### Modeling parameters distributions ----------------------------------------------------------------------------#

setwd("/home/hannelore/Documents/ANALYSES/TVB_global2/results_TVBii/")

myGroups=levels(group)
myProps=summary(group)/nrow(results)

# G - J relationship
palette(viridis(4))
par(mfrow=c(1,1))

jpeg('G-Jmd_dist_thrA.jpg')
plot(J_thrA_md ~ G_thrA, col=group, type="p", pch=19)
legend(1.1, 1.50, levels(group), col=1:length(group), pch=19)
dev.off()

jpeg('G-JiTmd_dist_thrA.jpg')
plot(Ji_tumor_thrA_md ~ G_thrA, col=group, type="p", pch=19)
legend(1.05,2.2, levels(group), col=1:length(group), pch=19)
dev.off()


### Average model parameters per group --------------------------------------------------------------------#

aggregate(results[, c('G_thrA', 'G_thrR', 'J_thrA', 'J_thrR')], list(results$group), mean)
aggregate(results[, c('G_thrA', 'G_thrR', 'J_thrA', 'J_thrR')], list(results$group), median)
aggregate(results[, c('Ji_tumor_thrA', 'Ji_tumor_thrR', 'Ji_DM_thrA', 'Ji_DM_thrR')], list(results$group), mean)
aggregate(results[, c('Ji_tumor_thrA', 'Ji_tumor_thrR', 'Ji_DM_thrA', 'Ji_DM_thrR')], list(results$group), median)


### Linear regression of modeling parameters by group, no covariates ---------------------------------------------#

# G
aov_G_thrA=aov(G_thrA ~ group)
summary(aov_G_thrA) 
shapiro.test(aov_G_thrA$residuals) 

# J (median whole-brain)
aov_Jmd_thrA=aov(J_thrA_md ~ group)
summary(aov_Jmd_thrA) 
shapiro.test(aov_Jmd_thrA$residuals) 

# Ji - non-tumor (median)
results$Ji_nontumor_thrA_md[1:11]=results$J_thrA_md[1:11]
attach(results)
aov_JiNT_thrA=aov(Ji_nontumor_thrA_md ~ group)
summary(aov_JiNT_thrA) 
shapiro.test(aov_JiNT_thrA$residuals) 
leveneTest(Ji_nontumor_thrA_md, group) 

# Ji - tumor (median)
# --> without sham tumor regions, just using whole-brain median Ji for controls
results$Ji_tumor_thrA_md[1:11]=results$J_thrA_md[1:11]
attach(results)
aov_Jimd_tumor_thrA=aov(Ji_tumor_thrA_md ~ group)
summary(aov_Jimd_tumor_thrA) 
shapiro.test(aov_Jimd_tumor_thrA$residuals) 
TukeyHSD(aov(Ji_tumor_thrA_md~group))
leveneTest(Ji_tumor_thrA_md, group)


## Now using residuals (v3: only regressing out ROI size, not in-strength)
results$J_res = J_res_all3$JiBrain
results$JiNT_res = J_res_all3$JiNT
results$JiT_res = J_res_all3$JiT
results$JiNT_res[1:11] = results$J_res[1:11]
results$JiT_res[1:11] = results$J_res[1:11]
attach(results)


# JiNT_res (median): use KW because residuals non-normal
kruskal.test(JiNT_res ~ group)
leveneTest(JiNT_res, group) 
pairwise.wilcox.test(JiNT_res, group, p.adjust.method = "BH")

# JiT_res (median)
kruskal.test(JiT_res ~ group) 
leveneTest(JiT_res, group) 
pairwise.wilcox.test(JiT_res, group, p.adjust.method = "BH")




### Link-wise correlations ~ group ------------------------------------------------------------------------------#

#FCsim-FCemp (PSE; thrA)
summary(Fcsim_emp_thrA)
sqrt(var(Fcsim_emp_thrA))
aov_PSE1 = aov(Fcsim_emp_thrA ~ group)
summary(aov_PSE1) 
shapiro.test(aov_PSE1$residuals) 
etaSquared(aov_PSE1, type=3, anova=FALSE) 

aov_PSE1x = aov(Fcsim_emp_thrA ~ lesion_vol_cm3)
summary(aov_PSE1x) 
shapiro.test(aov_PSE1x$residuals) 
etaSquared(aov_PSE1x, type=3, anova=FALSE) 

#SC-FCemp
summary(SC_FC)
aov_PSE2 = aov(SC_FC ~ group)
summary(aov_PSE2) 
shapiro.test(aov_PSE2$residuals) 
etaSquared(aov_PSE2, type=3, anova=FALSE)

#FCsim-FCemp <=> SC_FC
cor.test(Fcsim_emp_thrA, SC_FC)
hist(Fcsim_emp_thrA)
hist(SC_FC)



###################################################################################################################
###                                           PART 5: CLEAN DATASET                                             ###
###################################################################################################################


# Regress out important confounding effects of GTA and CANTAB scores

lm_RTI = lm(RTI_fiveRT ~ motivation + STAI + sexbin + age + lesion_vol)
summary(lm_RTI) 
RTI_res=lm_RTI$residuals
shapiro.test(RTI_res) 
results$RTI_res = c(RTI_res[1:27], NA , RTI_res[28:35])

lm_RVP = lm(RVP_A ~ motivation + STAI + sexbin + age + lesion_vol)
summary(lm_RVP) 
RVP_res=lm_RVP$residuals
shapiro.test(RVP_res) 
results$RVP_res = c(RVP_res[1:12], NA , RVP_res[13:26], NA, RVP_res[27:34])

lm_SOC = lm(SOC_prob_minmoves ~ motivation + STAI + sexbin + age + lesion_vol)
summary(lm_SOC) 
SOC_res=lm_SOC$residuals
shapiro.test(SOC_res) 
results$SOC_res = c(SOC_res[1:27], NA , SOC_res[28:35])

lm_SSP = lm(SSP_spanlength ~ motivation + STAI + sexbin + age + lesion_vol)
summary(lm_SSP) 
SSP_res=lm_SSP$residuals
shapiro.test(SSP_res) 
results$SSP_res = c(SSP_res[1:27], NA , SSP_res[28:35])


lm_Eglob = lm(SC_thrA_Eglob ~ age + sexbin + hand + lesion_vol + intnorm_scaling + fmri_MD + fmri_TR + STAI)
summary(lm_Eglob)
Eglob_res = lm_Eglob$residuals
shapiro.test(Eglob_res) 
hist(Eglob_res)
results$Eglob_res=Eglob_res

lm_Q = lm(SC_thrA_Q ~ age + sexbin + hand + lesion_vol + intnorm_scaling + fmri_MD + fmri_TR + STAI)
summary(lm_Q)
Q_res = lm_Q$residuals
shapiro.test(Q_res) 
hist(Q_res)
results$Q_res=Q_res

lm_PC = lm(SC_thrA_PC ~ age + sexbin + hand + lesion_vol + intnorm_scaling + fmri_MD + fmri_TR + STAI)
summary(lm_PC)
PC_res = lm_PC$residuals
shapiro.test(PC_res) 
hist(PC_res)
results$PC_res=PC_res

rm(list=c('RTI_res', 'RVP_res', 'SOC_res', 'SSP_res', 'PC_res', 'Eglob_res', 'Q_res', 
          'lm_RTI', 'lm_RVP', 'lm_SOC', 'lm_SSP', 'lm_Eglob', 'lm_Q', 'lm_PC'))

attach(results)

# Test for group differences:
Anova(aov(RTI_res~group), type=3) 
Anova(aov(RVP_res~group), type=3) 
Anova(aov(SOC_res~group), type=3) 
Anova(aov(SSP_res~group), type=3) 
Anova(aov(Eglob_res~group), type=3) 
Anova(aov(Q_res~group), type=3)
Anova(aov(PC_res~group), type=3) 
TukeyHSD(aov(PC_res~group))

# Quasi normalize using mean/sd of controls 
results$RTI_resnorm = (results$RTI_res - mean(results$RTI_res[c(1:11)])) / sqrt(var(results$RTI_res[c(1:11)]))
results$RVP_resnorm = (results$RVP_res - mean(results$RVP_res[c(1:11)])) / sqrt(var(results$RVP_res[c(1:11)]))
results$SOC_resnorm = (results$SOC_res - mean(results$SOC_res[c(1:11)])) / sqrt(var(results$SOC_res[c(1:11)]))
results$SSP_resnorm = (results$SSP_res - mean(results$SSP_res[c(1:11)])) / sqrt(var(results$SSP_res[c(1:11)]))
results$Eglob_resnorm = (results$Eglob_res - mean(results$Eglob_res[c(1:11)])) / sqrt(var(results$Eglob_res[c(1:11)]))
results$Q_resnorm = (results$Q_res - mean(results$Q_res[c(1:11)])) / sqrt(var(results$Q_res[c(1:11)]))
results$PC_resnorm = (results$PC_res - mean(results$PC_res[c(1:11)])) / sqrt(var(results$PC_res[c(1:11)]))


# Save dataset with ALL variables
write.table(x=results, file="RESULTS_ALL3.csv", quote=TRUE, sep=';', dec='.', row.names=FALSE)

# Final dataset
cbind(c(1:96), colnames(results))
vars_final = c(1,2,3,81,4:13,15,17,58,82,90:96)
results_cleaned = results[,vars_final]
str(results_cleaned)

write.table(x=results_cleaned, file="C:/Users/Hannelore/Dropbox/PhD/Papers/TVB_global/Results/RESULTS_cleaned2.csv",
            quote=TRUE, sep=';', dec='.', row.names=FALSE)

rm(list=ls())


###################################################################################################################
###                                           PART 6: VISUALIZATION OF RESULTS                                  ###
###################################################################################################################

setwd("/home/hannelore/Documents/ANALYSES/TVB_global2/")
results=read.table(file="RESULTS_cleaned3.csv", header=TRUE, sep=";")
results$group = factor(results$group, ordered=TRUE, levels=c('CON', 'MEN', 'GLI'))
attach(results)
library(car)
library(lsr)


#-----------------------------------------------------------------------------------------------------------------#

# Plot CANTAB results by group

myGroups=levels(group)
myProps=summary(group)/nrow(results)
colors=palette(c("slateblue", "cyan", "magenta"))
attach(results)

setwd('/home/hannelore/Documents/ANALYSES/TVB_global2/results_cantab')

jpeg('CANTAB_all_resnorm4.jpg', width=1200, height = 1000)
par(mfrow=c(2,2), mar=c(5,8,4,4)+.1, mgp = c(5, 2, 0))
plot(RTI_resnorm ~ group, col=colors, xlab="", ylab="Reaction time",
     ylim=c(-2.1,3.3), outline=FALSE, cex.lab=3.5, cex.axis=2.8, frame=FALSE)
grid(nx=NA, ny=NULL, col="darkgray", lty="longdash")
plot(RTI_resnorm ~ group, col=colors, xlab="", ylab="Reaction time",
     ylim=c(-2.1,3.3), outline=FALSE, cex.lab=3.5, cex.axis=2.8, frame=FALSE, add=TRUE)
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "RTI_resnorm"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=3)
  rm(myjitter)
}
plot(RVP_resnorm ~ group, col=colors, xlab="", ylab="Sustained attention",
     outline=FALSE, ylim=c(-2.8,3), cex.lab=3.5, cex.axis=2.8, frame=FALSE)
grid(nx=NA, ny=NULL, col="darkgray", lty="longdash")
plot(RVP_resnorm ~ group, col=colors, xlab="", ylab="Sustained attention",
     outline=FALSE, ylim=c(-2.8,3), cex.lab=3.5, cex.axis=2.8, frame=FALSE, add=TRUE)
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "RVP_resnorm"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=3)
  rm(myjitter)
}
plot(SOC_resnorm ~ group, col=colors, xlab="", ylab="Planning accuracy",
     outline=FALSE, ylim=c(-2.2,2), cex.lab=3.5, cex.axis=2.8, frame=FALSE)
grid(nx=NA, ny=NULL, col="darkgray", lty="longdash")
plot(SOC_resnorm ~ group, col=colors, xlab="", ylab="Planning accuracy",
     outline=FALSE, ylim=c(-2.2,2), cex.lab=3.5, cex.axis=2.8, frame=FALSE, add=TRUE)
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "SOC_resnorm"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=3)
  rm(myjitter)
}
plot(SSP_resnorm ~ group, col=colors, xlab="", ylab="Working memory capacity",
     outline=FALSE, ylim=c(-2.1,2), cex.lab=3.5, cex.axis=2.8, frame=FALSE)
grid(nx=NA, ny=NULL, col="darkgray", lty="longdash")
plot(SSP_resnorm ~ group, col=colors, xlab="", ylab="Working memory capacity",
     outline=FALSE, ylim=c(-2.1,2), cex.lab=3.5, cex.axis=2.8, frame=FALSE, add=TRUE)
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "SSP_resnorm"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=3)
  rm(myjitter)
}
dev.off()


#-----------------------------------------------------------------------------------------------------------------#

### Plot SC results by group

setwd('/home/hannelore/Documents/ANALYSES/TVB_global2/results_GTA/')

jpeg('GTA_thrA_resnorm3.jpg', width=1200, height=800)
par(mfrow=c(1,3), mar=c(5,8,4,4)+.1, mgp = c(5, 2, 0))
plot(Eglob_resnorm ~ group, col=colors, xlab="", ylab="Global efficiency", outline=FALSE, 
     ylim=c(-2,3), cex.lab=4.5, cex.axis=3.5, frame=FALSE)
grid(nx=NA, ny=NULL, col="darkgray", lty="longdash")
plot(Eglob_resnorm ~ group, col=colors, xlab="", ylab="Global efficiency", outline=FALSE, 
     ylim=c(-2,3), cex.lab=4.5, cex.axis=3.5, frame=FALSE, add=TRUE)
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "Eglob_resnorm"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=3)
  rm(myjitter)
}
plot(Q_resnorm ~ group, col=colors, xlab="", ylab="Modularity", outline=FALSE, 
     ylim=c(-2,2), cex.lab=4.5, cex.axis=3.5, frame=FALSE) #yaxt=c(-1,2,3),
grid(nx=NA, ny=NULL, col="darkgray", lty="longdash")
plot(Q_resnorm ~ group, col=colors, xlab="", ylab="Modularity", frame=FALSE, outline=FALSE, 
     ylim=c(-2,2), cex.lab=4.5, cex.axis=3.5, yaxt="n", yaxt=c(-1,2,3), add=TRUE)
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "Q_resnorm"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=3)
  rm(myjitter)
}
plot(PC_resnorm ~ group, col=colors, xlab="", ylab="Participation coefficient", outline=FALSE, 
     ylim=c(-2,3), cex.lab=4.5, cex.axis=3.5, frame=FALSE)
grid(nx=NA, ny=NULL, col="darkgray", lty="longdash")
plot(PC_resnorm ~ group, col=colors, xlab="", ylab="Participation coefficient", outline=FALSE, 
     ylim=c(-2,3), cex.lab=4.5, cex.axis=3.5, frame=FALSE, add=TRUE)
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "PC_resnorm"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=3)
  rm(myjitter)
}
dev.off()


#-----------------------------------------------------------------------------------------------------------------#

### Plot modeling parameters by group

setwd("/home/hannelore/Documents/ANALYSES/TVB_global2/results_TVBii/")

myGroups=levels(group)
myProps=summary(group)/nrow(results)


### Plot J raw ~ Group (TUMOR) (A)

jpeg('J_groupT_v2.jpg', width=1000, height=1000)
par(mfrow=c(1,1), mar=c(5,8,4,4)+.1, mgp = c(5, 2, 0))
plot(Ji_tumor_thrA_md ~ group, col=colors, xlab="", ylab="", outline=FALSE, 
     cex.axis=3.5, cex.lab=4, ylim=c(1,2.4), frame=FALSE);
grid(nx=NA, ny=NULL, col="darkgray", lty="longdash")
plot(Ji_tumor_thrA_md ~ group, col=colors, xlab="", ylab="", outline=FALSE, 
     cex.axis=3.5, cex.lab=4, ylim=c(1,2.4), frame=FALSE, add=TRUE);
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "Ji_tumor_thrA_md"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=3)
  rm(myjitter)
}
title(expression("J"[brain]), cex.main=5, adj=0.15)
title(expression("J"[tumor]), cex.main=5) 
title(expression("J"[tumor]), cex.main=5, adj=0.85) 
dev.off()


### Plot J raw (HEALTHY) ~ Group (B)

jpeg('J_groupH_v2.jpg', width=1000, height=1000)
par(mfrow=c(1,1), mar=c(5,8,4,4)+.1, mgp = c(5, 2, 0))
plot(Ji_nontumor_thrA_md ~ group, col=colors, xlab="", ylab="", outline=FALSE, 
     cex.axis=3.5, cex.lab=4, ylim=c(1.26,1.55), frame=FALSE);
grid(nx=NA, ny=NULL, col="darkgray", lty="longdash")
plot(Ji_nontumor_thrA_md ~ group, col=colors, xlab="", ylab="", outline=FALSE, 
     cex.axis=3.5, cex.lab=4, ylim=c(1.26,1.55), frame=FALSE, add=TRUE);
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "Ji_nontumor_thrA_md"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=3)
  rm(myjitter)
}
title(expression("J"[brain]), cex.main=5, adj=0.15)
title(expression("J"[non-tumor]), cex.main=5) 
title(expression("J"[non-tumor]), cex.main=5, adj=0.95) 
dev.off()


### Plot J res ~ Group (TUMOR) (C)

jpeg('Jres_groupT_v2.jpg', width=1000, height=1000)
par(mfrow=c(1,1), mar=c(5,8,4,4)+.1, mgp = c(5, 2, 0))
plot(JiT_res ~ group, col=colors, xlab="", ylab="", outline=FALSE, 
     cex.axis=3.5, cex.lab=4, ylim=c(-1.9,0.9), frame=FALSE);
grid(nx=NA, ny=NULL, col="darkgray", lty="longdash")
plot(JiT_res ~ group, col=colors, xlab="", ylab="", outline=FALSE, 
     cex.axis=3.5, cex.lab=4, ylim=c(-1.9,0.9), frame=FALSE, add=TRUE);
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "JiT_res"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=3)
  rm(myjitter)
}
title(expression("J"[brain]), cex.main=5, adj=0.15)
title(expression("J"[tumor]), cex.main=5) 
title(expression("J"[tumor]), cex.main=5, adj=0.85) 
dev.off()


### Plot J res ~ Group (HEALTHY) (D)

jpeg('Jres_groupH_v2.jpg', width=1000, height=1000)
par(mfrow=c(1,1), mar=c(5,8,4,4)+.1, mgp = c(5, 2, 0))
plot(JiNT_res ~ group, col=colors, xlab="", ylab="", outline=FALSE, 
     cex.axis=3.5, cex.lab=4, ylim=c(0.45,0.65), frame=FALSE);
grid(nx=NA, ny=NULL, col="darkgray", lty="longdash")
plot(JiNT_res ~ group, col=colors, xlab="", ylab="", outline=FALSE, 
     cex.axis=3.5, cex.lab=4, ylim=c(0.45,0.65), frame=FALSE, add=TRUE);
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "JiNT_res"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=3)
  rm(myjitter)
}
title(expression("J"[brain]), cex.main=5, adj=0.15)
title(expression("J"[non-tumor]), cex.main=5) 
title(expression("J"[non-tumor]), cex.main=5, adj=0.95) 
dev.off()


### Plot G ~ Group (E)

palette=c("slateblue", "cyan", "magenta")

jpeg('G-Group5.jpg', height=1200, width=400)
par(mar=c(5,6,4,1)+.1)
plot(G_thrA ~ group, col=palette, ylab="", xlab="", cex.axis=2.1, frame=F); 
grid(nx=NA, ny=NULL, col="darkgray", lty="longdash")
plot(G_thrA ~ group, col=palette, ylab="", xlab="", cex.axis=2.1, frame=F, add=TRUE); 
title("G", cex.main=3, font.main=1)
for (i in 1:length(myGroups)){
  thisGroup=myGroups[i]
  thisValues=results[results$group==thisGroup, "G_thrA"]
  myjitter=jitter(rep(i, length(thisValues)), amount=myProps[i]/2)
  points(myjitter, thisValues, pch=16, col=rgb(0,0,0,0.4), cex=2.5)
}
dev.off()
