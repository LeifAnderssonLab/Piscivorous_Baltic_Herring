# Import the dataset
Dioxine_Data_Table <- read.delim2("/pathtodata/Dioxine_Data_Table.txt")
Dioxine_Data_Table <- Dioxine_Data_Table[ -c(2,3) ,]

custom_colors <- c("limegreen", "indianred3")

############################
# Figure 2C script.       ##
############################

# Fat content
pdf("/pathtodata/Fat_Content.pdf", height = 5, width = 5)
boxplot(Dioxine_Data_Table$Fat_content... ~ Dioxine_Data_Table$Population, ylab="% Fat Content", xlab = "", ylim = c(0,12), col = custom_colors, axes = F)
box(bty="l")
axis(2)
axis(1, at = 1:2, labels = c("Slåttersill", "Strömming"))
dev.off()

############################
# Figure 2D script.       ##
############################

# Dioxine content comparison
pdf("/pathtodata/Dioxine_Only.pdf", height = 5, width = 5)
boxplot(Dioxine_Data_Table$X.PCDD.F_.pg.TEQ.g.vv. ~ Dioxine_Data_Table$Population, ylab="dioxine", xlab = "", col = custom_colors, axes = F)
box(bty="l")
axis(2)
axis(1, at = 1:2, labels = c("Slåttersill", "Strömming"))
dev.off()

############################
# Figure 2E script.       ##
############################
rm(list=ls()) # Empty data frame
options("install.lock"=FALSE) # Run in case packages will not install in library

data <- read.delim("/pathtodata/Raw_Data_Edit.txt")

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
data$Dist<-as.numeric(data$Distance_um)
data$Batch<-as.factor(data$Batch)
data$ID_nr<-as.factor(data$ID)
str(data) # Check data

#Plot data of individual Sr:Ca profiles
ggplot(data=data, aes(x=Dist, y=Sr_Ca_Avg, group=ID_nr)) +
  geom_line(aes(color=Batch),size=1.3)+xlab("")+ylab("Sr:Ca")+
  facet_wrap(.~ID_nr)+theme_bw()+
  theme(
    legend.text = element_text(size=18),
    legend.title.align = 0,
    legend.position = "bottom", 
    legend.title = element_text(size=14, vjust = .5, hjust = .3))

# Suset data Dist > = 0
B <- subset(data, Dist>="0")

### Plot Sr:Ca profiles fromm the otolith core to the anterior rostrum edge per spawning type.
B$Batch_f = factor(B$Batch, levels=c('Slåttersill','Autumn-spawning','Spring-spawning'))
plot<-ggplot(data=B, aes(x=Dist, y=Sr_Ca_Avg, group=ID_nr,colour=Batch_f)) +
  geom_line(aes(),size=1,alpha=0.4)+theme_classic(base_size = 25) + 
  theme(plot.title = element_text(hjust = 0.5))+#scale_color_brewer(palette="Set1")+
  geom_smooth(aes(group=Type), method = lm, formula = y ~ splines::bs(x, 25),colour="Black")+
  ylab("Sr:Ca")+xlab("Distance")+facet_wrap(.~Batch_f)+ theme(legend.position='none')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot+ scale_colour_manual(values = c("dodgerblue", "orange", "red")) # Figure 2F for paper

############################
# Figure 2F script.       ##
############################

library(ggplot2)

Isotope_Data <- read.delim2("~/pathtodata/Isotope_Data.txt")
Isotope_Data <- Isotope_Data[ -c(1:2) ,]

# Filter the dataset for Stromming and Slattersill groups
stromming_data <- subset(Isotope_Data, FID == "Stromming")
slattersill_data <- subset(Isotope_Data, FID == "Slattersill")

pdf("/pathtoexport/Dotplot_Figure.pdf", height = 5, width = 5)

# Create dot plot with ggplot
dotplot <- ggplot(data = Isotope_Data, aes(x = d13C, y = d15N, color = FID)) +
  xlim(-27, -20) + ylim(9, 12) +
  geom_point(size = 3) + # Adjust size of points here
  labs(x = "d13C", y = "d15N") +
  scale_color_manual(values = c("Stromming" = "indianred3", "Slattersill" = "limegreen"), name = "") + # Add color legend with species names
  theme_classic() +
  theme(
    legend.position = c(0.45, 0.001), # Adjust the position of the legend within the plot
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "transparent") # Make legend background transparent
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE) # Arrange legend into one row with two columns
  )

# Add ellipse
dotplot_with_ellipse <- dotplot +
  stat_ellipse(type = "norm", level = 0.95, geom = "polygon", alpha = 0) 

# Print the plot
print(dotplot_with_ellipse)

dev.off()