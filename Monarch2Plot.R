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
  theme(legend.position = 'bottom')
