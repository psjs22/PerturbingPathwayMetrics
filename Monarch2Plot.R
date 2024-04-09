# set working directory and load necessary packages
setwd("") # same location as the Matlab file MonarchExamplePPCM2.m
library(R.matlab)
library(tidyverse)

# load in monarch data from Matlab
Monarch2data <- read.csv('Monarch2data.csv',header = FALSE)
colnames(Monarch2data) <- c('Pathway','Value','Pert','C_metric','xValue')

# create plot 
ggplot() + 
  geom_point(data=Monarch2data, aes(x=xValue, y=Value, shape=Pert, fill=C_metric, color=C_metric))+
  scale_shape_manual(values = c(25,24), breaks=c('T','C'),name='Perturbation Applied',labels=c('Threats','Threats and Conservation')) +
  scale_color_manual(values = c('black','blue','red'), breaks=c('H','S','M'),name='Contribution Metric',label=c('Habitat','Subpopulation','Metapopulation')) +
  scale_fill_manual(values = c('black','blue','red'), breaks=c('H','S','M')) +
  guides(fill = 'none') + # remove fill legend
  labs(x='Pathway', y='Change in Contribution Metric') +
  theme_minimal() + # remove axes
  geom_hline(aes(yintercept=0)) + # add y axis
  scale_x_continuous(breaks = seq(0,29,2), minor_breaks = seq(0,29,1)) +
  theme(legend.position = 'bottom', legend.box = 'vertical', legend.margin = margin()) +
  theme(text=element_text(size=12), axis.title = element_text(size=14)) # change font size

# save high resolution image
library(gridExtra)
#save .png image
png(".../Monarch2plot.png", width=18, height=12, units='cm', res=300) #set size and res #add in folder of where to save file
grid.arrange (Fig) # make plot
dev.off()
#save .tiff image
tiff(".../Monarch2plot.tiff", width=18, height=12, units='cm', res=300) #set size and res #add in folder of where to save file
grid.arrange (Fig) # make plot
dev.off()
#save .eps image
ggsave(".../Monarch2plot.eps", width=18, height=12, units='cm', dpi=300) #set size and res #add in folder of where to save file

