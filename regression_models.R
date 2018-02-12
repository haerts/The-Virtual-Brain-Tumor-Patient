###################################################################################################################
###                                                                                                             ###
###                   MODELING BRAIN DYNAMICS IN BRAIN TUMOR PATIENTS USING THE VIRTUAL BRAIN                   ###
###                ==============================================================================               ###
###    				 Use regression model to check effect of GTA/cognition on modeling params    			    ###
###                                                                                                             ###
### Created by Hannelore Aerts                                                                                  ###
### Date last update: 06/02/2018                                                                                ###
### --> use Ji residuals, corrected for ROI size only (not for in-strength)
###################################################################################################################

### Read in data and prepare for analyses ------------------------------------------------------------------------#

setwd("/home/hannelore/Documents/ANALYSES/TVB_global2/")
results=read.table(file="RESULTS_cleaned3.csv", header=TRUE, sep=";")
results$group = factor(results$group, ordered=TRUE, levels=c('CON', 'MEN', 'GLI'))
attach(results)

library(car)
library(lsr)


### Test for group differences in CANTAB and GTA metrics ---------------------------------------------------------#
Anova(aov(RTI_resnorm~group), type=3) 
Anova(aov(RVP_resnorm~group), type=3) 
Anova(aov(SOC_resnorm~group), type=3)
Anova(aov(SSP_resnorm~group), type=3) 
Anova(aov(Eglob_resnorm~group), type=3) 
Anova(aov(Q_resnorm~group), type=3) 
Anova(aov(PC_resnorm~group), type=3) 
TukeyHSD(aov(PC_resnorm~group))



### Model relationship with most simple model possible -----------------------------------------------------------#

### G ###

### G - CANTAB
lm_G_cantab=lm(G_thrA ~ RTI_resnorm + RVP_resnorm + SOC_resnorm + SSP_resnorm)
summary(lm_G_cantab)
etaSquared(lm_G_cantab, anova=FALSE, type = 3)

### G - GTA
lm_G_GTA=lm(G_thrA ~ Eglob_resnorm + Q_resnorm + PC_resnorm)
summary(lm_G_GTA)
etaSquared(lm_G_GTA, type=3, anova=FALSE)


### whole-brain median J ###

### J - CANTAB
lm_J_cantab=lm(J_thrA_md ~ RTI_resnorm + RVP_resnorm + SOC_resnorm + SSP_resnorm)
summary(lm_J_cantab)
etaSquared(lm_J_cantab, anova=FALSE, type = 3)

### J - GTA
lm_J_GTA=lm(J_thrA_md ~ Eglob_resnorm + Q_resnorm + PC_resnorm)
summary(lm_J_GTA)
etaSquared(lm_J_GTA, type=3, anova=FALSE)


### whole-brain median J -- FOR CONTROLS ONLY -- ###

### J - CANTAB
lm_J_cantab_CON=lm(J_thrA_md[1:11] ~ RTI_resnorm[1:11] + RVP_resnorm[1:11] + SOC_resnorm[1:11] + SSP_resnorm[1:11])
summary(lm_J_cantab_CON)
etaSquared(lm_J_cantab_CON, anova=FALSE, type = 3)

### J - GTA
lm_J_GTA_CON=lm(J_thrA_md[1:11] ~ Eglob_resnorm[1:11] + Q_resnorm[1:11] + PC_resnorm[1:11])
summary(lm_J_GTA_CON)
etaSquared(lm_J_GTA_CON, type=3, anova=FALSE)


### Non-tumor median J -- IN PATIENTS ONLY -- ###

### JiNT - CANTAB
lm_JiNT_cantab=lm(Ji_nontumor_thrA_md[12:36] ~ RTI_resnorm[12:36] + RVP_resnorm[12:36] + 
                    SOC_resnorm[12:36] + SSP_resnorm[12:36])
summary(lm_JiNT_cantab)
etaSquared(lm_JiNT_cantab, anova=FALSE, type = 3)

### JiNT - GTA
lm_JiNT_GTA=lm(Ji_nontumor_thrA_md[12:36] ~ Eglob_resnorm[12:36] + Q_resnorm[12:36] + PC_resnorm[12:36])
summary(lm_JiNT_GTA)
etaSquared(lm_JiNT_GTA, anova=FALSE, type = 3)


### Tumor median J -- IN PATIENTS ONLY -- ###

### JiT - CANTAB 
lm_JiT_cantab=lm(Ji_tumor_thrA_md[12:36] ~ RTI_resnorm[12:36] + RVP_resnorm[12:36] + 
                   SOC_resnorm[12:36] + SSP_resnorm[12:36])
summary(lm_JiT_cantab)
etaSquared(lm_JiT_cantab, anova=FALSE, type = 3)

### JiT - GTA
lm_JiT_GTA=lm(Ji_tumor_thrA_md[12:36] ~ Eglob_resnorm[12:36] + Q_resnorm[12:36] + PC_resnorm[12:36])
summary(lm_JiT_GTA)
etaSquared(lm_JiT_GTA, anova=FALSE, type = 3)



### whole-brain median J -- RES ###

### Jres - CANTAB
lm_Jres_cantab=lm(J_res ~ RTI_resnorm + RVP_resnorm + SOC_resnorm + SSP_resnorm)
summary(lm_Jres_cantab)
etaSquared(lm_Jres_cantab, anova=FALSE, type = 3)

### Jres - GTA
lm_Jres_GTA=lm(J_res ~ Eglob_resnorm + Q_resnorm + PC_resnorm)
summary(lm_Jres_GTA)
etaSquared(lm_J_GTA, type=3, anova=FALSE)


### whole-brain median Jres -- FOR CONTROLS ONLY -- ###

### Jres - CANTAB
lm_Jres_cantab_CON=lm(J_res[1:11] ~ RTI_resnorm[1:11] + RVP_resnorm[1:11] + SOC_resnorm[1:11] + SSP_resnorm[1:11])
summary(lm_Jres_cantab_CON)
etaSquared(lm_Jres_cantab_CON, anova=FALSE, type = 3)

### Jres - GTA
lm_Jres_GTA_CON=lm(J_res[1:11] ~ Eglob_resnorm[1:11] + Q_resnorm[1:11] + PC_resnorm[1:11])
summary(lm_Jres_GTA_CON)
etaSquared(lm_Jres_GTA_CON, type=3, anova=FALSE)


### Non-tumor median Jres -- IN PATIENTS ONLY -- ###

### JiNTres - CANTAB
lm_JiNTres_cantab=lm(JiNT_res[12:36] ~ RTI_resnorm[12:36] + RVP_resnorm[12:36] + 
                    SOC_resnorm[12:36] + SSP_resnorm[12:36])
summary(lm_JiNTres_cantab)
etaSquared(lm_JiNTres_cantab, anova=FALSE, type = 3)

### JiNTres - GTA

lm_JiNTres_GTA=lm(JiNT_res[12:36] ~ Eglob_resnorm[12:36] + Q_resnorm[12:36] + PC_resnorm[12:36])
summary(lm_JiNTres_GTA)
etaSquared(lm_JiNTres_GTA, anova=FALSE, type = 3)


### Tumor median Jres -- IN PATIENTS ONLY -- ###

### JiTres - CANTAB 
lm_JiTres_cantab=lm(JiT_res[12:36] ~ RTI_resnorm[12:36] + RVP_resnorm[12:36] + 
                   SOC_resnorm[12:36] + SSP_resnorm[12:36])
summary(lm_JiTres_cantab)
etaSquared(lm_JiTres_cantab, anova=FALSE, type = 3)

### JiT - GTA
lm_JiTres_GTA=lm(JiT_res[12:36] ~ Eglob_resnorm[12:36] + Q_resnorm[12:36] + PC_resnorm[12:36])
summary(lm_JiTres_GTA)
etaSquared(lm_JiTres_GTA, anova=FALSE, type = 3)


### Plot significant associations --------------------------------------------------------------------------------#

palette(c("slateblue", "cyan", "magenta"))
setwd("/home/hannelore/Documents/ANALYSES/TVB_global2/results_regression")


### G & JiNT - Eglob together

jpeg('G_JiNT-Eglob.jpg', width=1200, height=600)
par(mar=c(5,6,4,1)+.1, mfrow=c(1,2), oma=c(4,1,1,1))
plot(G_thrA ~ Eglob_resnorm, pch=21, bg=group, cex=2, cex.lab=2, cex.axis=2,
     xlab="Global efficiency", ylab="G")
legend(1.2,2.2, legend=paste('r =',round(cor(Eglob_resnorm, G_thrA),2)), box.lty=0, cex=2)
model=(lm(G_thrA ~ Eglob_resnorm))
myPredict=predict(model, interval="confidence")
ix <- sort(Eglob_resnorm,index.return=T)$ix
lines(Eglob_resnorm[ix], myPredict[ix , 1], col="darkgray", lwd=2 )  
polygon(c(rev(Eglob_resnorm[ix]), Eglob_resnorm[ix]), c(rev(myPredict[ix,3]), myPredict[ix,2]), 
        col = rgb(0.7,0.7,0.7,0.4) , border = NA)
rm(list=c("model", "myPredict", "ix"))

plot(JiNT_res[12:36] ~ Eglob_resnorm_PAT, pch=21, bg=group_PAT, cex=2, cex.lab=2, cex.axis=2,
     xlab="Global efficiency", ylab=expression("J"[non-tumor]), xaxt='n')
axis(side = 1, at = c(-1,0,1,2), cex.axis=2)
legend(1,0.65, legend=paste('r =',round(cor(Eglob_resnorm_PAT, JiNT_res[12:36]),2)), box.lty=0, cex=2)
model=(lm(JiNT_res[12:36] ~ Eglob_resnorm_PAT))
myPredict=predict(model, interval="confidence")
ix <- sort(Eglob_resnorm_PAT,index.return=T)$ix
lines(Eglob_resnorm_PAT[ix], myPredict[ix , 1], col="darkgray", lwd=2 )  
polygon(c(rev(Eglob_resnorm_PAT[ix]), Eglob_resnorm_PAT[ix]), c(rev(myPredict[ix,3]), myPredict[ix,2]), 
        col = rgb(0.7,0.7,0.7,0.4) , border = NA)
rm(list=c("model", "myPredict", "ix"))

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("topright", levels(group), xpd=TRUE, pt.bg=1:length(group), box.lty=0, pch=21, horiz=TRUE, cex=1.8)

dev.off()


# J - RTI & RVP together 

jpeg('J-RTI_RVP_resnorm5.jpg', width=1200, height=1200)

par(mar=c(5,6,4,1)+.1, mfrow=c(2,2))

RTI_resnorm_PAT = RTI_resnorm[c(12:27,29:35)]
plot(JiT_res[c(12:27,29:35)] ~ RTI_resnorm_PAT, pch=21, bg='blue', cex=2.5, cex.lab=2.5, cex.axis=2.5,
     xlab="Reaction time", ylab=expression("J"[tumor]))
RTI_resnorm_PAT = RTI_resnorm[c(12:27,29:35)]
legend(-1.7,-1.1, legend=paste('r =',round(cor(RTI_resnorm_PAT, JiT_res[c(12:27,29:35)]),2)), box.lty=0, cex=2.5)
model=(lm(JiT_res[c(12:27,29:35)] ~ RTI_resnorm_PAT))
myPredict=predict(model, interval="confidence")
ix <- sort(RTI_resnorm_PAT,index.return=T)$ix
lines(RTI_resnorm_PAT[ix], myPredict[ix , 1], col="darkgray", lwd=2 )  
polygon(c(rev(RTI_resnorm_PAT[ix]), RTI_resnorm_PAT[ix]), c(rev(myPredict[ix,3]), myPredict[ix,2]), 
        col = rgb(0.7,0.7,0.7,0.4) , border = NA)
rm(list=c("model", "myPredict", "ix"))

RTI_resnorm_PAT = RTI_resnorm[c(12:27,29:35)]
plot(JiNT_res[c(12:27,29:35)] ~ RTI_resnorm_PAT, pch=21, bg='blue', cex=2.5, cex.lab=2.5, cex.axis=2.5,
     xlab="Reaction time", ylab=expression("J"[non-tumor]))
RTI_resnorm_PAT = RTI_resnorm[c(12:27,29:35)]
legend(-1.75,0.51, legend=paste('r =',round(cor(RTI_resnorm_PAT, JiNT_res[c(12:27,29:35)], use="pairwise.complete.obs"),2)), box.lty=0, cex=2.5)
model=(lm(JiNT_res[c(12:27,29:35)] ~ RTI_resnorm_PAT))
myPredict=predict(model, interval="confidence")
ix <- sort(RTI_resnorm_PAT,index.return=T)$ix
lines(RTI_resnorm_PAT[ix], myPredict[ix , 1], col="darkgray", lwd=2 )  
polygon(c(rev(RTI_resnorm_PAT[ix]), RTI_resnorm_PAT[ix]), c(rev(myPredict[ix,3]), myPredict[ix,2]), 
        col = rgb(0.7,0.7,0.7,0.4) , border = NA)
rm(list=c("model", "myPredict", "ix"))
#left out PAT31 because of outlier... now r=0.32
#with outlier r=-0.19


RVP_resnorm_PAT = RVP_resnorm[c(12,14:27,29:36)]

plot(JiT_res[c(12,14:27,29:36)] ~ RVP_resnorm_PAT, pch=21, bg='blue', cex=2.5, cex.lab=2.5, cex.axis=2.5,
     xlab="Sustained attention", ylab=expression("J"[tumor]))
legend(-3.2,-1.35, legend=paste('r =',round(cor(RVP_resnorm_PAT, JiT_res[c(12,14:27,29:36)]),2)), box.lty=0, cex=2.5)
model=(lm(JiT_res[c(12,14:27,29:36)] ~ RVP_resnorm_PAT))
myPredict=predict(model, interval="confidence")
ix <- sort(RVP_resnorm_PAT,index.return=T)$ix
lines(RVP_resnorm_PAT[ix], myPredict[ix , 1], col="darkgray", lwd=2 )  
polygon(c(rev(RVP_resnorm_PAT[ix]), RVP_resnorm_PAT[ix]), c(rev(myPredict[ix,3]), myPredict[ix,2]), 
        col = rgb(0.7,0.7,0.7,0.4) , border = NA)
rm(list=c("model", "myPredict", "ix"))

plot(JiNT_res[c(12,14:27,29:36)] ~ RVP_resnorm_PAT, pch=21, bg='blue', cex=2.5, cex.lab=2.5, cex.axis=2.5,
     xlab="Sustained attention", ylab=expression("J"[non-tumor]))
legend(-3.2,0.51, legend=paste('r = ',round(cor(RVP_resnorm_PAT, JiNT_res[c(12,14:27,29:36)]),2)), box.lty=0, cex=2.5)
model=(lm(JiNT_res[c(12,14:27,29:36)] ~ RVP_resnorm_PAT))
myPredict=predict(model, interval="confidence")
ix <- sort(RVP_resnorm_PAT,index.return=T)$ix
lines(RVP_resnorm_PAT[ix], myPredict[ix , 1], col="darkgray", lwd=2 )  
polygon(c(rev(RVP_resnorm_PAT[ix]), RVP_resnorm_PAT[ix]), c(rev(myPredict[ix,3]), myPredict[ix,2]), 
        col = rgb(0.7,0.7,0.7,0.4) , border = NA)
rm(list=c("model", "myPredict", "ix"))

dev.off()




