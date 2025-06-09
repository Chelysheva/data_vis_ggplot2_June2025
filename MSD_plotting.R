###Data visualisation with ggplot2
##Irina & Rao
##10.06.2025
#MSD cytokine dataset

library(tidyverse)
####Loading data####
cyto<-read.csv("MSD_data.csv")
#Remove negative values
#cyto <- cyto %>%
#  filter(if_all(c(IFNy, IL1, IL2, IL10), ~ .x > 0 | is.na(.x)))
#cyto <- cyto %>%
#  mutate(across(c(IFNy, IL1, IL2, IL10), ~ ifelse(.x < 0, 0, .x)))
cyto <- cyto %>%
  mutate(across(
    c(IFNy, IL1, IL2, IL10),
    ~ ifelse(.x < 0,
             0.5 * min(.x[.x > 0], na.rm = TRUE),
             .x)))

#Convert groups into factors
cyto$Randomisation <- as.factor(cyto$Randomisation)
levels(cyto$Randomisation)
c(1,2) #note the difference
levels(cyto$Randomisation)<-c("ChAdOx1", "Control")

#Wide to long format
#Use gather from tidyr
cyto_melted <- pivot_longer(cyto, cols = -c(ID, Randomisation), names_to = "cytokine", values_to = "value")
cyto_melted<-na.omit(cyto_melted)

####Plotting####

#####Barplot - we plot the medians####
ggplot(cyto_melted, aes(fill=Randomisation, y=value, x=cytokine)) + 
  geom_bar(position="dodge", stat = "summary", fun = "median") +
  theme_bw()

#Stacked barplot
ggplot(cyto_melted, aes(fill=Randomisation, y=value, x=cytokine)) + 
  geom_bar(position="stack", stat = "summary",) +
  theme_bw()

#Percent stacked barplot
ggplot(cyto_melted, aes(fill=Randomisation, y=value, x=cytokine)) + 
  geom_bar(position="fill", stat="summary", fun = "mean") +
  theme_bw()

#Reverse the groupping - check the relative aboundance of each cytokine by randomisation group
ggplot(cyto_melted, aes(fill=cytokine, y=value, x=Randomisation)) + 
  geom_bar(position="fill", stat="identity", width = 0.5) +
  theme_bw()

#Add error bars
ggplot(cyto_melted, aes(fill=Randomisation, y=(value), x=cytokine)) + 
  geom_bar(position="dodge", stat = "summary", fun = "mean") + 
  stat_summary(geom = "errorbar", fun.data = mean_se, 
               position = position_dodge(width = 0.90), width = 0.5) +
  theme_bw()

#Bar plots is not a good way to represent our data!

#####Boxplots####
ggplot(cyto_melted, aes(fill=Randomisation, y=value, x=cytokine)) + 
  geom_boxplot() +
  theme_bw()

#Scale - log transform
ggplot(cyto_melted, aes(fill=Randomisation, y=value, x=cytokine)) + 
  scale_y_log10() +
  geom_boxplot() +
  theme_bw()

#Facets
ggplot(cyto_melted, aes(fill=Randomisation, y=value, x=cytokine)) + 
  scale_y_log10() +
  geom_boxplot() +
  theme_bw() +
  #theme(legend.position = "none")+ #remove the legend as it is redundant
  facet_wrap(~Randomisation)

#Try: geom_dotplot() - why this is not a best approach for my data?

####Jitter plots####
ggplot(cyto_melted, aes(y=value, x=cytokine)) +
  scale_y_log10() +
  geom_jitter(aes(shape=Randomisation, color=Randomisation), #specifying aes only for jitter
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              size = 1.2) +
  theme_bw()

#Adding stat
ggplot(cyto_melted, aes(group=Randomisation, y=log(value), x=cytokine)) + #group is needed for stat_summary function
  geom_jitter(aes(shape=Randomisation, color=Randomisation), #specifying aes only for jitter
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              size = 1.2) + 
  theme_bw() +
  stat_summary(
               fun.data="mean_sdl",  fun.args = list(mult=1), 
               geom = "pointrange",  size = 0.4,
               position = position_dodge(0.6), show.legend = F)

##Finalising graph format
ggplot(cyto_melted, aes(group=Randomisation, y=value, x=cytokine)) +
  scale_y_log10() +
  geom_jitter(aes(color=Randomisation),
              position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5),
              size = 1.2, alpha=0.5) + #alpha - transparency
  theme_bw() +
  scale_color_manual(values=c("red","blue")) +
  xlab("") +
  ylab("Cytokine concentration") +
  stat_summary(
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.2,
    position = position_dodge(0.5), show.legend = FALSE, alpha=0.5)


####Performing statistics and adding significance levels####
#install.packages("ggpubr")
library(ggpubr)
#Summary
test<-compare_means(data = cyto_melted, formula = value ~ Randomisation, group.by = "cytokine")
#Adding to the ggplot
ggplot(cyto_melted, aes(group=Randomisation, y=value, x=cytokine)) +
  scale_y_log10() +
  geom_jitter(aes(color=Randomisation),
              position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5),
              size = 1.2, alpha=0.5) + 
  theme_pubr() +
  scale_color_manual(values=c("red","blue")) +
  xlab("") +
  ylab("Cytokine concentration") +
  stat_summary(
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.2,
    position = position_dodge(0.5), show.legend = FALSE, alpha=0.5)+
  stat_compare_means(aes(group = Randomisation), method = "wilcox", hide.ns = T, paired = F,
                     label = "p.signif",bracket.size = 0.3)
#alternatively - label = "p.signif"

library(rstatix)
#Use p-adjusted values instead
stat.test <- cyto_melted %>%
  group_by(cytokine) %>%
  wilcox_test(value ~ Randomisation) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
View(stat.test)
stat.test <- stat.test %>% add_xy_position(x = "cytokine")
#Finalising graph format + stats
ggplot(cyto_melted, aes(group=Randomisation, y=value, x=cytokine)) +
  scale_y_log10() +
  geom_jitter(aes(color=Randomisation),
              position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5),
              size = 1.2, alpha=0.5) + 
  theme_pubr() +
  scale_color_manual(values=c("red","blue")) +
  xlab("") +
  ylab("Cytokine concentration") +
  stat_summary(
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4, alpha = 0.5,
    position = position_dodge(0.5), show.legend = FALSE) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0,hide.ns = TRUE, inherit.aes = FALSE, y.position = 3)

#Try: reorder the cytokines on this graph like by alphabet (check out how to reorder factor)

#Create interactive plot
# install.packages("plotly")
library(plotly)
p<-ggplot(cyto_melted, aes(group=Randomisation, y=value, x=cytokine)) +
  scale_y_log10() +
  geom_jitter(aes(color=Randomisation),
              position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5),
              size = 1.2, alpha = 0.5) + 
  theme_pubr() +
  scale_color_manual(values=c("red","blue")) +
  xlab("") +
  ylab("Cytokine concentration") +
  stat_summary(
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4, alpha = 0.5,
    position = position_dodge(0.5), show.legend = FALSE) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0,hide.ns = TRUE, inherit.aes = FALSE, y.position = 3)
p
ggplotly(p)

