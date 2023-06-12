rm(list = ls())

library(readxl)
library(tidyverse)
library(dbplyr)
library(dplyr)
library(patchwork)
library(ggplot2)
library(stringr)
library(ggpubr)
library(reshape2)
library(rio)
library(xtable)
library(gridExtra)
library(gtable)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(tidyr)
library(here)


setwd(here('/Volumes/research/metabogut/BMI_2021/enrichment_analysis/new_analysis_2208'))
here::i_am("data_preparation/significant_kos_gmms.R")

#read results
res_gmm <- read_excel(here('results/results_gmm/main_analyses', 'all_results_main_v221104.xlsx'))

#read the files where I get GMM annotations
abund_gmm <- read.csv(here('results/output_files/gmm/abundance.csv'))

#overlap them to take the description
merge1 <- merge(res_gmm,abund_gmm, by.x = 'gmm', by.y = 'Module', all = TRUE)
merge1_subset <- merge1[,c(1:13)]

#there are columns with :maleTRUE or :femaleTRUE. Split and merge again
#subset the rows where : is present
gmm_subset_colon <- merge1_subset[grepl(":", merge1_subset[["gmm"]]),]
#split to two columns 
gmm_split_colon_1 <- unlist(lapply(strsplit(as.character(gmm_subset_colon$gmm), ":"), '[[', 1))
gmm_split_colon_2 <- unlist(lapply(strsplit(as.character(gmm_subset_colon$gmm), ":"), '[[', 2))
gmm_subset_colon <- gmm_subset_colon %>% mutate(split1 = gmm_split_colon_1, split2 = gmm_split_colon_2)


#merge the dataframe with res_gmm to get the description
transient_gmm_with_desc <- merge(gmm_subset_colon, abund_gmm, by.x = 'split1', by.y = 'Module', all = FALSE)
transient_gmm_with_desc <- transient_gmm_with_desc[,c(2,16)]
#make description and GMMs unique
transient_gmm_with_desc <- transient_gmm_with_desc[!duplicated(transient_gmm_with_desc[,c('gmm','Description.y')]),]

#merge back with the original results dataframe
gmm_with_desc <- merge(res_gmm, transient_gmm_with_desc, by = 'gmm', all = TRUE)
gmm_with_desc$Description.y <-  ifelse(is.na(gmm_with_desc$Description.y), merge1_subset$Description, gmm_with_desc$Description.y)

#fix final structure
colnames(gmm_with_desc)[13] <- 'Description'
res_gmm <- gmm_with_desc[,c(1,13, 2:12)]
res_gmm <- res_gmm[order(res_gmm$p_fdr_per_outcome),]

#save results
export(res_gmm, file = here('results/output_files/gmm', 'all_res_annotated.xlsx'))

#check significant results
sig_res_gmm <- subset(res_gmm, p_fdr_per_outcome <= 0.05) #860 for all models

#check sig BMI
length(which(sig_res_gmm$model_name == 'full_M1_BMI' & sig_res_gmm$cohort == 'all_cohorts' & sig_res_gmm$analysis == 'main')) #68

length(which(sig_res_gmm$model_name == 'full_M3_waist_adjBMI' & sig_res_gmm$cohort == 'all_cohorts' & sig_res_gmm$analysis == 'main')) #22

length(which(sig_res_gmm$model_name == 'full_M2_WHR_adjBMI' & sig_res_gmm$cohort == 'all_cohorts' & sig_res_gmm$analysis == 'main')) #15


#save to variables
sig_res_gmm_all_main_full_BMI <- sig_res_gmm[which(sig_res_gmm$model_name == 'full_M1_BMI' & sig_res_gmm$cohort == 'all_cohorts' & sig_res_gmm$analysis == 'main'),]
sig_res_gmm_all_main_full_waist <- sig_res_gmm[which(sig_res_gmm$model_name == 'full_M3_waist_adjBMI' & sig_res_gmm$cohort == 'all_cohorts' & sig_res_gmm$analysis == 'main'),]
sig_res_gmm_all_main_full_whr <- sig_res_gmm[which(sig_res_gmm$model_name == 'full_M2_WHR_adjBMI' & sig_res_gmm$cohort == 'all_cohorts' & sig_res_gmm$analysis == 'main'),]

all_sig_pheno_all_main_full <- rbind(sig_res_gmm_all_main_full_BMI,sig_res_gmm_all_main_full_waist, sig_res_gmm_all_main_full_whr)

all_sig_pheno_all_main_full <- all_sig_pheno_all_main_full %>% mutate(phenotype = ifelse(model_name == 'full_M1_BMI', 'BMI', 
                                                                                         ifelse(model_name == 'full_M3_waist_adjBMI', 'Waist_adjBMI','WHR_adj_BMI'))) #74 significant modules (some overlap among the 3 phenotypes)

#save results
export(all_sig_pheno_all_main_full, file = here('results/output_files/gmm', 'all_sig_res_annotated.xlsx'), overwrite = TRUE)


P <- all_sig_pheno_all_main_full %>%  group_by (phenotype) %>% ggplot(aes(x=fct_reorder(Description,est), y=est, ymin=ci_low, ymax=ci_high)) +
  geom_pointrange(position = position_dodge(width = 1), size=0.1,aes(color=phenotype), shape = 20) + 
  scale_color_manual(values = c("black", "red", "blue"), labels = c("BMI", "Waist_adjBMI", "WHR_adjBMI")) + 
  xlab("GMM") + 
  ylab('Beta') + 
  coord_flip() + theme_minimal() + 
  theme(axis.text.y = element_text(size= 4.5)) 
P 


###############################################################################
###############################################################################
###############################################################################
# Marlena suggested that I add estimates of GMMs that are significant in one phenotype but not the other
#first take res_gmm which has everything and choose only the models of our interest

res_gmm_all_main_full_BMI <- res_gmm[which(res_gmm$model_name == 'full_M1_BMI' & res_gmm$cohort == 'all_cohorts' & res_gmm$analysis == 'main'),]
res_gmm_all_main_full_waist <- res_gmm[which(res_gmm$model_name == 'full_M3_waist_adjBMI' & res_gmm$cohort == 'all_cohorts' & res_gmm$analysis == 'main'),]
res_gmm_all_main_full_whr <- res_gmm[which(res_gmm$model_name == 'full_M2_WHR_adjBMI' & res_gmm$cohort == 'all_cohorts' & res_gmm$analysis == 'main'),]

##### it should be 103 modules for all (sig and nonsig)

#add them all on top of each other
all_pheno_all_main_full <- rbind(res_gmm_all_main_full_BMI,res_gmm_all_main_full_waist, res_gmm_all_main_full_whr)
all_pheno_all_main_full <- all_pheno_all_main_full %>% mutate(phenotype = ifelse(model_name == 'full_M1_BMI', 'BMI', 
                                                                                 ifelse(model_name == 'full_M3_waist_adjBMI', 'Waist_adjBMI','WHR_adj_BMI')))

#merge with all_sig_pheno_all_main_full to get the estimates of GMMs that are for instance significant in BMI but not in whradjbmi

merge_all_with_sig <- merge(all_sig_pheno_all_main_full, all_pheno_all_main_full , by = 'gmm', all = FALSE) #74 unique GMMs

# subset the needed data
merge_all_with_sig_subset <- merge_all_with_sig %>% 
  transmute(gmm = gmm, Description = Description.y, analysis = analysis.y, outcome = outcome.y,
            cohort = cohort.y, model_name = model_name.y, 
            n_obs = n_obs.y, est = est.y,
            se = se.y,ci_low = ci_low.y, ci_high = ci_high.y, p = p.y, 
            p_fdr_per_outcome = p_fdr_per_outcome.y,
            phenotype = phenotype.y)

merge_all_with_sig_subset <- merge_all_with_sig_subset %>% 
  mutate(phenotype = ifelse(phenotype == 'WHR_adj_BMI', 'WHRadjBMI', 
                            ifelse(phenotype == 'Waist_adjBMI', 'WaistadjBMI', 'BMI')))


#remove duplicated values 
merge_all_with_sig_subset <- unique(merge_all_with_sig_subset)
#save to excel to edit the CIs
export(merge_all_with_sig_subset, file = here('results/output_files/gmm', 'merge_all_with_sig_subset.xlsx'),overwrite = TRUE)


#reread the file to use in the plot
#merge_all_with_sig_subset <- read_excel('/Users/ma8244mi/Desktop/BMI-paper/enrichment_analysis/new_analysis_2208/results/output_files/gmm/merge_all_with_sig_subset.xlsx')

to_table <- merge_all_with_sig_subset[,c(1,2,8,13,14)]

#to_table <- to_table  %>% arrange(est)

#x <- merge_all_with_sig_subset %>% 
# group_by (phenotype) %>% 
#arrange(Description, est)


#use spread from tidyr to transform column values to colnames #there are 222 unique estimates (74 for each phenotype)

transform_to_table <- to_table %>% spread(phenotype, p_fdr_per_outcome)
##########I have a list of 74 gmms (MFxxxx) where each one is repeated 3 times (T = 221.. 1 module for each phenotype 
######I want to replace the NA values with the p_fdr_per_outcome values so that they are all on the same line for each gmm (see example_output in excelsheet). 
transform_to_table<- transform_to_table %>%
  group_by(gmm) %>%
  summarise( Description = Description,
             BMI = na.omit(BMI),
             WaistadjBMI = na.omit(WaistadjBMI),
             WHRadjBMI = na.omit(WHRadjBMI))

transform_to_table <- unique(transform_to_table) # 74 unique Modules

transform_to_table <- transform_to_table %>% 
  mutate(BMI = signif(BMI,2),
         WaistadjBMI = signif(WaistadjBMI,2),
         WHRadjBMI = signif(WHRadjBMI,2))

######### START PLOTTING EACH PHENOTYPE SEPARATELY

data_bmi_plot <- merge_all_with_sig_subset[merge_all_with_sig_subset$phenotype == 'BMI',]

P_bmi <- data_bmi_plot %>% ggplot(aes(x=reorder(factor(Description),est), y=est, ymin=ci_low, ymax=ci_high)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_pointrange(size=0.1,aes(color=phenotype), shape = 20) + 
  scale_color_manual(values = 'grey30', labels = 'BMI') + 
  xlab("") + 
  ylab('est (95% CI)') + 
  ggtitle("GMM") + 
  # coord_cartesian(ylim = c(-5, 30))+ 
  coord_flip() + theme_minimal() +
  theme(axis.text.y = element_text(size= 8.5),
        axis.text.x = element_text(size = 6.5),
        axis.title.x = element_text(size = 9),
        plot.title = element_text(size = 9.2,hjust = -0.3),
        legend.position = 'none',
        #legend.position = c(0.9,1.5),
        #legend.position = 'bottom',
        # legend.direction = 'horizontal',
        plot.margin=grid::unit(c(2,0,0,0), "mm")) 
P_bmi

#get y axis labels in the order they were plotted
desc_in_P_bmi <- ggplot_build(P_bmi)$layout$panel_params[[1]]$y$get_labels()
#flip bottom to top
desc_in_P_bmi<- rev(desc_in_P_bmi)


###### before plotting the heatmap, make sure that the order of the table Description matches the ones in the plot (P)

# Extract the levels of the text_column  in df1 in the same order as the unique values in df2
levels_df1_bmi <- levels(as.factor(transform_to_table$Description))[match(desc_in_P_bmi, levels(as.factor(transform_to_table$Description)))]

# Find the indices of the key values in the text_column column of the data frame
indices_bmi <- match(levels_df1_bmi, transform_to_table$Description)

# Reorder the values in the text_column column based on the indices
transform_to_table <- transform_to_table[indices_bmi,]


#heatmap
################################### BMI

data_bmi_hp <- transform_to_table[,3]
data_bmi_hp[nrow(data_bmi_hp) +1, ] <- 0 
data_bmi_hp[nrow(data_bmi_hp) +1, ] <- 1 
data_bmi_hp <- as.matrix(data_bmi_hp)
#colnames(data) <- NULL

#create this foo variable just to add numbers that are not numeric to the table because some numeric numbers are being converted to e notation and others are not
#foo_bmi_hp <- as.data.frame(transform_to_table %>% 
#                      mutate(BMI = formatC(BMI, format = 'e', digits = 2)))
#names(foo) <- NULL

#foo_bmi_hp[nrow(foo_bmi_hp) +1, ] <- 0 
#foo_bmi_hp[nrow(foo_bmi_hp) +1, ] <- 1 
#setting colors
clrsp <- colorRampPalette(c("red", "white"))   
clrs <- clrsp(nrow(data_bmi_hp)) 
#set breaks from 0 to 0.05
breaks1 <- seq(-log10(0+1e-50), -log10(1), length.out = nrow(data_bmi_hp))

data_bmi_hp[75,1] <- data_bmi_hp[75,1] + 1e-50


p_bmi_hp <- grid.grabExpr(draw(ComplexHeatmap::pheatmap(-log10(data_bmi_hp),cluster_rows = F, cluster_cols = F, 
                                                        cellwidth=42, 
                                                        cellheight = 10.2,
                                                        display_numbers = data_bmi_hp, 
                                                        column_names_side = c("top"),
                                                        name = "p value",
                                                        angle_col = c("0"),  
                                                        number_color = "black", 
                                                        border_color = NA,
                                                        fontsize_col = 9.3, 
                                                        fontsize = 10.5,
                                                        legend = FALSE,
                                                        col =clrs, breaks = breaks1)))
all_hp_bmi <- grid.arrange( P_bmi,p_bmi_hp, widths = c(2, 1))

#ggsave('/Users/ma8244mi/Desktop/GMM_hp_bmi.pdf', all_hp_bmi, width = 30,height = 50, units = 'cm')

########################## waist # match the description in waist plot to that of the heatmap so that they are aligned

data_waist_plot <- merge_all_with_sig_subset[merge_all_with_sig_subset$phenotype == 'WaistadjBMI',]

# First match the description in waist with that in bmi (because all plots should be aligned with the same descriptions to the left)
levels_df1_waist <- levels(as.factor(data_waist_plot$Description))[match(desc_in_P_bmi, levels(as.factor(data_waist_plot$Description)))]

# Find the indices of the key values in the text_column column of the data frame to fix the heatmap
#I am using transform_to_table dataframe because it contains the exponentiated pvalues which are to be added in heatmap cells
indices_waist_hp <- match(levels_df1_waist, transform_to_table$Description)

# Reorder the values in the text_column column based on the indices
transform_to_table_waist_hp <- transform_to_table[indices_waist_hp,] #this is for the heatmap

#----- #the task here is to match the description in waist dataframe to that in BMI so that when they are plotted in same plot, the descriptions are aligned in both phenotypes

# pr = this is used for point range plot
indices_waist_pr <- match(levels_df1_waist, data_waist_plot$Description)

# Reorder the values in the text_column column based on the indices
data_waist_plot <- data_waist_plot[indices_waist_pr,]

#this is to force plot the y axis labels the order they are in the dataframe
data_waist_plot$Description <- factor(data_waist_plot$Description, levels=data_waist_plot$Description)



P_waist <- data_waist_plot %>% ggplot(aes(x=Description, y=est, ymin=ci_low, ymax=ci_high)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_pointrange(size=0.1,aes(color=phenotype), shape = 20) + 
  scale_color_manual(values = 'mediumorchid', labels = 'WaistadjBMI') + scale_x_discrete(limits=rev)+
  xlab("") + 
  ylab('est (95% CI)') + 
  ggtitle("") + 
  # coord_cartesian(ylim = c(-5, 30))+ 
  coord_flip() + theme_minimal() +
  theme(#axis.text.y = element_text(size= 6.5),
    axis.text.x = element_text(size = 6.5),
    axis.title.x = element_text(size = 9),
    plot.title = element_text(size = 9.2,hjust = -0.09),
    #plot.title = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    legend.position = 'none',
    #legend.position = c(0.9,1.5),
    #legend.position = 'bottom',
    # legend.direction = 'horizontal',
    plot.margin=grid::unit(c(2,0,0,0), "mm")) 
P_waist


# heatmap
data_waist_hp <- transform_to_table_waist_hp[,4]
data_waist_hp[nrow(data_waist_hp) +1, ] <- 0 
data_waist_hp[nrow(data_waist_hp) +1, ] <- 1 


data_waist_hp <- as.matrix(data_waist_hp)
#colnames(data) <- NULL

#create this foo variable just to add numbers that are not numeric to the table because some numeric numbers are being converted to e notation and others are not
#foo_waist_hp <- as.data.frame(transform_to_table %>% 
#                           mutate(WaistadjBMI = formatC(WaistadjBMI, format = 'e', digits = 2)))
#names(foo) <- NULL

#foo_waist_hp[nrow(foo_waist_hp) +1, ] <- 0 
#foo_waist_hp[nrow(foo_waist_hp) +1, ] <- 1 

data_waist_hp[75,1] <- data_waist_hp[75,1] + 1e-50

p_waist_hp <- grid.grabExpr(draw(ComplexHeatmap::pheatmap(-log10(data_waist_hp),cluster_rows = F, cluster_cols = F, 
                                                          cellwidth=42, 
                                                          cellheight = 10.2,
                                                          display_numbers = data_waist_hp, 
                                                          column_names_side = c("top"),
                                                          name = "p value",
                                                          angle_col = c("0"),  
                                                          number_color = "black", 
                                                          border_color = NA,
                                                          fontsize_col = 9.3, 
                                                          fontsize = 10.5,
                                                          legend = FALSE,
                                                          col =clrs, breaks = breaks1)))

grid.arrange( P_bmi,p_bmi_hp, P_waist, p_waist_hp, widths = c(2,1,1,1))

###########################WHR match the description in waist plot to that of the heatmap so that they are aligned

data_whr_plot <- merge_all_with_sig_subset[merge_all_with_sig_subset$phenotype == 'WHRadjBMI',]

# First match the description in whr with that in bmi (because all plots should be aligned with the same descriptions to the left)
levels_df1_whr <- levels(as.factor(data_whr_plot$Description))[match(desc_in_P_bmi, levels(as.factor(data_whr_plot$Description)))]

# Find the indices of the key values in the text_column column of the data frame to fix the heatmap
#I am using transform_to_table dataframe because it contains the exponentiated pvalues which are to be added in heatmap cells
indices_whr_hp <- match(levels_df1_whr, transform_to_table$Description)

# Reorder the values in the text_column column based on the indices
transform_to_table_whr_hp <- transform_to_table[indices_whr_hp,] #this is for the heatmap

#----- #the task here is to match the description in whr dataframe to that in BMI so that when they are plotted in same plot, the descriptions are aligned in both phenotypes

# pr = this is used for point range plot
indices_whr_pr <- match(levels_df1_whr, data_whr_plot$Description)

# Reorder the values in the text_column column based on the indices
data_whr_plot <- data_whr_plot[indices_whr_pr,]

data_whr_plot$Description <- factor(data_whr_plot$Description, levels=data_whr_plot$Description)


P_whr <- data_whr_plot %>% ggplot(aes(x=Description, y=est, ymin=ci_low, ymax=ci_high)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_pointrange(size=0.1,aes(color=phenotype), shape = 20) + 
  scale_color_manual(values = 'orange', labels = 'WHRadjBMI') + scale_x_discrete(limits=rev)+
  xlab("") + 
  ylab('est (95% CI)') + 
  ggtitle("") + 
  # coord_cartesian(ylim = c(-5, 30))+ 
  coord_flip() + theme_minimal() +
  theme(#axis.text.y = element_text(size= 6.5),
    axis.text.x = element_text(size = 6.5),
    axis.title.x = element_text(size = 9),
    plot.title = element_text(size = 9.2,hjust = -0.09),
    #plot.title = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    legend.position = 'none',
    #legend.position = c(0.9,1.5),
    #legend.position = 'bottom',
    # legend.direction = 'horizontal',
    plot.margin=grid::unit(c(2,0,0,0), "mm")) 
P_whr



data_whr_hp <- transform_to_table_whr_hp[,5]
data_whr_hp[nrow(data_whr_hp) +1, ] <- 0 
data_whr_hp[nrow(data_whr_hp) +1, ] <- 1 

data_whr_hp <- as.matrix(data_whr_hp)
#colnames(data) <- NULL

#create this foo variable just to add numbers that are not numeric to the table because some numeric numbers are being converted to e notation and others are not
#foo_whr_hp <- as.data.frame(transform_to_table %>% 
#                               mutate(WHRadjBMI = formatC(WHRadjBMI, format = 'e', digits = 2)))
#names(foo) <- NULL

#foo_whr_hp[nrow(foo_whr_hp) +1, ] <- 0 
#foo_whr_hp[nrow(foo_whr_hp) +1, ] <- 1 

data_whr_hp[75,1] <- data_whr_hp[75,1] + 1e-50

p_whr_hp <- grid.grabExpr(draw(ComplexHeatmap::pheatmap(-log10(data_whr_hp),cluster_rows = F, cluster_cols = F, 
                                                        cellwidth=42, 
                                                        cellheight = 10.2,
                                                        display_numbers = data_whr_hp, 
                                                        column_names_side = c("top"),
                                                        name = "p value",
                                                        angle_col = c("0"),  
                                                        number_color = "black", 
                                                        fontsize_col = 9.3, 
                                                        fontsize = 10.5,
                                                        border_color = NA,
                                                        col =clrs, breaks = breaks1)))

#add white rectangle
rect_whr <- rectGrob(gp=gpar(fill='white', col = 'white'), x=unit(-1.4,"npc"), 
                     y=unit(0.007,"npc"),
                     width=unit(1.3,"npc"), 
                     height=unit(0.048,"npc"))

rect_waist <- rectGrob(gp=gpar(fill='white', col = 'white'), x=unit(-8,"npc"), 
                       y=unit(0.007,"npc"),
                       width=unit(1.5,"npc"), 
                       height=unit(0.048,"npc"))

rect_bmi <- rectGrob(gp=gpar(fill='white', col = 'white'), x=unit(-15,"npc"), 
                     y=unit(0.007,"npc"),
                     width=unit(1.5,"npc"), 
                     height=unit(0.048,"npc"))


all_plots <- grid.arrange( P_bmi,p_bmi_hp,P_waist, p_waist_hp, P_whr, 
                           p_whr_hp,rect_whr, rect_waist, rect_bmi, widths = c(3.5,1,2,1,2,1,0.5,0.5,0.5))

ggsave('/Users/ma8244mi/Desktop/BMI-paper/plots/GMM-heatmap.tiff', all_plots, width = 50,height = 28, units = 'cm')

