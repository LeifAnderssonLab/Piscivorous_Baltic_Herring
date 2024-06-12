################################################################################
################################################################################
#
#    R code for Goodall et al. - Herring otolith chemistry
#
#     Yvette Heimbrand SLU Aqua 2024-06-10
#     Email: yvette.heimbrand@slu.se
#
################################################################################
################################################################################
#
rm(list=ls()) # Empty data frame
options("install.lock"=FALSE) # Run in case packages will not install in library

data <- read.csv("Goodall_et_al__Herring_otolith_chemistry_data.csv", sep = ";") # Read data file

# Load libraries
library(viridis)
library(car)
library(lme4)
library(effects)
library(sjPlot)
library(r2glmm)
library(ggplot2)

str(data) # Check data

# Transform data
data$Dist<-as.numeric(data$Distance..µm)
data$Batch<-as.factor(data$Batch)
data$ID_nr<-as.factor(data$ID)
str(data) # Check data

#Plot data of individual Sr:Ca profiles
data$Batch = factor(data$Batch, levels=c('Slåttersill','Autumn-spawning','Spring-spawning'))
plot1<-ggplot(data=data, aes(x=Dist, y=Sr.Ca_Avg, group=ID_nr)) +
        geom_line(aes(color=Batch),size=1.3)+xlab("")+ylab("Sr:Ca")+
        facet_wrap(.~ID_nr)+theme_bw()+
        theme(legend.text = element_text(size=18),legend.title.align = 0,legend.position = "bottom", 
        legend.title = element_text(size=14, vjust = .5, hjust = .3))
plot1 + scale_colour_manual(values = c("red", "orange", "dodgerblue"))


### Plot lifelong Sr:Ca profiles from the otolith core to the edge per spawning type.
data$Batch = factor(data$Batch, levels=c('Slåttersill','Autumn-spawning','Spring-spawning'))
plot2<-ggplot(data=data, aes(x=Dist, y=Sr.Ca_Avg, group=ID_nr,colour=Batch)) +
  geom_line(aes(),size=1,alpha=0.4)+theme_classic(base_size = 25) + 
  theme(plot.title = element_text(hjust = 0.5))+#scale_color_brewer(palette="Set1")+
  geom_smooth(aes(group=Batch), method = lm, formula = y ~ splines::bs(x, 25),colour="Black")+
  ylab("Sr:Ca")+xlab("Distance from the otolith core (µm)")+facet_wrap(.~Batch)+ theme(legend.position='none')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot2 + scale_colour_manual(values = c("red", "orange", "dodgerblue")) # Figure for paper

#*******************************************************************************
#
### Select data corresponding to adult herring -> Distance from the core > = 15000

df <- subset(data, Dist >=1500)

# Check distribution of Sr:Ca
summary(df$Sr.Ca_Avg)
hist(df$Sr.Ca_Avg)
qqPlot(df$Sr.Ca_Avg)

# Log transform Sr:Ca
df$Sr_log<-log(df$Sr.Ca_Avg)
hist(df$Sr_log) # Better
qqPlot(df$Sr_log) # Better

# Plot data
df$Batch_f = factor(df$Batch, levels=c('Slåttersill','Autumn-spawning','Spring-spawning'))
plot<-ggplot(data=df, aes(x=Dist, y=Sr_log, group=ID_nr,colour=Batch_f)) +
  geom_line(aes(),size=1,alpha=0.4)+theme_classic(base_size = 25) + 
  theme(plot.title = element_text(hjust = 0.5))+#scale_color_brewer(palette="Set1")+
  geom_smooth(aes(group=Batch_f), method = lm,colour="Black")+
  ylab("Sr:Ca")+xlab("Distance")+facet_wrap(.~Batch_f)+ theme(legend.position='none')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot+ scale_colour_manual(values = c("red", "orange", "dodgerblue")) 

### Repeated measurements mixed model 
#
#Compare models
g3.mixed<-lmer(Sr.Ca_Avg~ Dist*Batch +(1|ID),data=df,na.action=na.exclude)
g4.mixed<-lmer(Sr.Ca_Avg~ Dist+Batch +(1|ID),data=df,na.action=na.exclude)
g5.mixed<-lmer(log(Sr.Ca_Avg)~ Dist*Batch +(1|ID),data=df,na.action=na.exclude)
g6.mixed<-lmer(log(Sr.Ca_Avg)~ Dist+Batch +(1|ID),data=df,na.action=na.exclude)
g7.mixed<-lmer(log(Sr.Ca_Avg)~ Dist+Batch +Dist*Batch+(1|ID),data=df,na.action=na.exclude)
AIC(g3.mixed,g4.mixed,g5.mixed,g6.mixed,g7.mixed) # selecting g7.mixed model

Anova(g7.mixed)
#Check normality
plot(g7.mixed)
qqnorm(resid(g7.mixed))

# Pairwise comparison of model slopes for adult herring spawning types
library(emmeans)
cs <- emtrends(g7.mixed,specs=pairwise~Batch,var="Dist")
css <- summary(cs,infer=TRUE)
css$contrasts # Table SXXX for paper
css$emtrends

# Plot effect
p<-emmip(g7.mixed, Batch ~ Dist, cov.reduce = range)
p

df$Batch = factor(df$Batch, levels=c('Slåttersill','Autumn-spawning','Spring-spawning'))
  
q<-p+theme_classic()+geom_line(aes(colour=Batch), linewidth=2)+
  theme_classic(base_size = 23)+
  ylab("Linear prediction of log Sr:Ca")+
  xlab("Distance from the otolith core (µm)")+
  theme(legend.position='top')+
  theme(legend.title=element_blank())

q + scale_colour_manual(values = c("red", "orange", "dodgerblue")) #  Figure for paper


#Check Multicollinearity
vif(g7.mixed) #All vif < 5 = No multicolinearity

# Check effects
sjPlot::plot_model(g7.mixed)

plot(allEffects(g7.mixed))

# Check r2
r2 <- r2beta(model=g7.mixed,partial=TRUE,method='sgv')
print(r2)
plot(x=r2)

################################################################################
################################################################################