# install packages 
install.packages("ggplot2")
install.packages("showtext")
install.packages("extrafont")

# load libraries
library(ggplot2)
library(readr)
library(showtext)
library(extrafont)
library(dplyr)


# Efruxi_CHX_GSEA Barplot -------------------------------------------------


# import data - this case two files one for upregulated pathways at desired condition, one for down 
data_up <- read.delim("C:\\Users\\oberl\\R_projects_07-02-25\\cac_fgf21_ingWATRNAseq\\GSEA\\Efruxi_vs_CHX_mh_all.Gsea.1760973137337\\gsea_report_for_CHX_1760973137337.tsv", sep = "\t")
data_down <- read.delim("C:\\Users\\oberl\\R_projects_07-02-25\\cac_fgf21_ingWATRNAseq\\GSEA\\Efruxi_vs_CHX_mh_all.Gsea.1760973137337\\gsea_report_for_Efruxi_1760973137337.tsv", sep = "\t")

# combine both dataframes
data_combined <- rbind(data_up, data_down)

# set classes for filtering conditions - Normalized enrichment score and false discovery rate
data_combined$NES<- as.numeric(data_combined$NES)
data_combined$FDR.q.val <- as.numeric(data_combined$FDR.q.val)
str(data_combined)

# set cutoff values for filtering conditions
NES_cutoff <- 1.5 # enrichment score should be min. 1.5
FDR_cutoff <- 0.25 # FDR should not exceed 0.25 == 25%

# filter dataframe accordingly - NES should be bigger or equal than cutoff; FDR should be smaller than cutoff 
significant_pathways <- data_combined %>%
  filter(abs(NES) >= NES_cutoff & FDR.q.val < FDR_cutoff)

# remove the hallmark gsea annotation
significant_pathways$NAME <- gsub("^HALLMARK_", "", significant_pathways$NAME)

# order data beginnning with the biggest NES value in descending order
significant_pathways$NAME <- factor(significant_pathways$NAME,
                                    levels = significant_pathways$NAME[order(significant_pathways$NES, decreasing = FALSE)])


### Barplot ###

ggp <-ggplot(significant_pathways, aes(x = NAME, y = NES, fill = NES))+ 
  geom_bar(stat = "identity", 
           width = 0.7, # bar width
           alpha = 0.8) + # transparency
  coord_flip() + # flip coordinates to be plotted horizontally
  
  scale_y_continuous(
    breaks = seq(-3, 3, by = 1),
    limits = c(min(significant_pathways$NES, -3), 
               max(significant_pathways$NES, 3))) +
  
  scale_fill_gradient2(
    low = "gray65",
    mid = "white",
    high = "lightslateblue",
    midpoint = 0,
    name = "NES") + # set colors
  
  labs(title = "ingWAT, Hallmark GSEA Pathways, Efruxi_CHX", x = "", y = "NES") + # set labels for x, y axis and title
  theme_minimal(base_family = "Arial") + # define theme
  theme(panel.grid = element_blank(), # no grid
        axis.text.y = element_blank(), # remove y labels, we wanted them to be integrated to save space
        axis.title.x = element_text(size = 9, hjust =0.5), # define x axis title
        axis.ticks.x = element_line(linewidth = 0.7), # tick lines settings
        axis.line.x.bottom = element_line(color = "black", size = 0.7), # extra line for definition
        legend.text = element_text(size = 10), # legend settings
        plot.title = element_text(size = 10, hjust = 0.5)) + # plot title settings
  
  geom_text(aes(y = ifelse(NES < 0, 0.03, -0.03), label = NAME), # define distance to middle axis
            hjust = ifelse(significant_pathways$NES < 0, 0, 1), # display pathway labels on opposing site 
            size = 2.6) 


# plot
ggp


# save plot
ggsave("GSEA_Plot_Efruxi.tiff", ggp, width = 6, height = 5, dpi = 600)



# AAV_CHX_GSEA Barplot ----------------------------------------------------


# import data - this case two files one for upregulated pathways at desired condition, one for down 
data_up <- read.delim("C:\\Users\\oberl\\R_projects_07-02-25\\cac_fgf21_ingWATRNAseq\\GSEA\\AAV_vs_CHX_mh_all.Gsea.1760973747843\\gsea_report_for_AAV_1760973747843.tsv", sep = "\t")
data_down <- read.delim("C:\\Users\\oberl\\R_projects_07-02-25\\cac_fgf21_ingWATRNAseq\\GSEA\\AAV_vs_CHX_mh_all.Gsea.1760973747843\\gsea_report_for_CHX_1760973747843.tsv", sep = "\t")

# combine both dataframes
data_combined <- rbind(data_up, data_down)

# set classes for filtering conditions - Normalized enrichment score and false discovery rate
data_combined$NES<- as.numeric(data_combined$NES)
data_combined$FDR.q.val <- as.numeric(data_combined$FDR.q.val)
str(data_combined)

# set cutoff values for filtering conditions
NES_cutoff <- 1.5 # enrichment score should be min. 1.5
FDR_cutoff <- 0.25 # FDR should not exceed 0.25 == 25%

# filter dataframe accordingly - NES should be bigger or equal than cutoff; FDR should be smaller than cutoff 
significant_pathways <- data_combined %>%
  filter(abs(NES) >= NES_cutoff & FDR.q.val < FDR_cutoff)

# remove the hallmark gsea annotation
significant_pathways$NAME <- gsub("^HALLMARK_", "", significant_pathways$NAME)

# order data beginnning with the biggest NES value in descending order
significant_pathways$NAME <- factor(significant_pathways$NAME,
                                    levels = significant_pathways$NAME[order(significant_pathways$NES, decreasing = FALSE)])


### Barplot ###

ggp <-ggplot(significant_pathways, aes(x = NAME, y = NES, fill = NES))+ 
  geom_bar(stat = "identity", 
           width = 0.7, # bar width
           alpha = 0.8) + # transparency
  coord_flip() + # flip coordinates to be plotted horizontally
  
  scale_y_continuous(
    breaks = seq(-3, 3, by = 1),
    limits = c(min(significant_pathways$NES, -3), 
               max(significant_pathways$NES, 3))) +
  
  scale_fill_gradient2(
    low = "gray65",
    mid = "white",
    high = "maroon3",
    midpoint = 0,
    name = "NES") + # set colors
  
  labs(title = "WAT, Hallmark GSEA Pathways, AAV_CHX", x = "", y = "NES") + # set labels for x, y axis and title
  theme_minimal(base_family = "Arial") + # define theme
  theme(panel.grid = element_blank(), # no grid
        axis.text.y = element_blank(), # remove y labels, we wanted them to be integrated to save space
        axis.title.x = element_text(size = 9, hjust =0.5), # define x axis title
        axis.ticks.x = element_line(linewidth = 0.7), # tick lines settings
        axis.line.x.bottom = element_line(color = "black", size = 0.7), # extra line for definition
        legend.text = element_text(size = 10), # legend settings
        plot.title = element_text(size = 10, hjust = 0.5)) + # plot title settings
  
  geom_text(aes(y = ifelse(NES < 0, 0.03, -0.03), label = NAME), # define distance to middle axis
            hjust = ifelse(significant_pathways$NES < 0, 0, 1), # display pathway labels on opposing site 
            size = 2.6) 


# plot
ggp


# save plot
ggsave("GSEA_Plot_AAV.tiff", ggp, width = 5.5, height = 5, dpi = 600)


# CHX_HC GSEA Barplot -----------------------------------------------------


# import data - this case two files one for upregulated pathways at desired condition, one for down 
data_up <- read.delim("C:\\Users\\oberl\\R_projects_07-02-25\\cac_fgf21_ingWATRNAseq\\GSEA\\CHX_vs_HC_mh_all.Gsea.1760973880028\\gsea_report_for_CHX_1760973880028.tsv", sep = "\t")
data_down <- read.delim("C:\\Users\\oberl\\R_projects_07-02-25\\cac_fgf21_ingWATRNAseq\\GSEA\\CHX_vs_HC_mh_all.Gsea.1760973880028\\gsea_report_for_HC_1760973880028.tsv", sep = "\t")

# combine both dataframes
data_combined <- rbind(data_up, data_down)

# set classes for filtering conditions - Normalized enrichment score and false discovery rate
data_combined$NES<- as.numeric(data_combined$NES)
data_combined$FDR.q.val <- as.numeric(data_combined$FDR.q.val)
str(data_combined)

# set cutoff values for filtering conditions
NES_cutoff <- 1.5 # enrichment score should be min. 1.5
FDR_cutoff <- 0.25 # FDR should not exceed 0.25 == 25%

# filter dataframe accordingly - NES should be bigger or equal than cutoff; FDR should be smaller than cutoff 
significant_pathways <- data_combined %>%
  filter(abs(NES) >= NES_cutoff & FDR.q.val < FDR_cutoff)

# remove the hallmark gsea annotation
significant_pathways$NAME <- gsub("^HALLMARK_", "", significant_pathways$NAME)

# order data beginnning with the biggest NES value in descending order
significant_pathways$NAME <- factor(significant_pathways$NAME,
                                    levels = significant_pathways$NAME[order(significant_pathways$NES, decreasing = FALSE)])


### Barplot ###

ggp <-ggplot(significant_pathways, aes(x = NAME, y = NES, fill = NES))+ 
  geom_bar(stat = "identity", 
           width = 0.7, # bar width
           alpha = 0.8) + # transparency
  coord_flip() + # flip coordinates to be plotted horizontally
  
  scale_y_continuous(
    breaks = seq(-3, 3, by = 1),
    limits = c(min(significant_pathways$NES, -3), 
               max(significant_pathways$NES, 3))) +
  
  scale_fill_gradient2(
    low = "navyblue",
    mid = "white",
    high = "red3",
    midpoint = 0,
    name = "NES") + # set colors
  
  labs(title = "WAT, Hallmark GSEA Pathways, CHX_HC", x = "", y = "NES") + # set labels for x, y axis and title
  theme_minimal(base_family = "Arial") + # define theme
  theme(panel.grid = element_blank(), # no grid
        axis.text.y = element_blank(), # remove y labels, we wanted them to be integrated to save space
        axis.title.x = element_text(size = 9, hjust =0.5), # define x axis title
        axis.ticks.x = element_line(linewidth = 0.7), # tick lines settings
        axis.line.x.bottom = element_line(color = "black", size = 0.7), # extra line for definition
        legend.text = element_text(size = 10), # legend settings
        plot.title = element_text(size = 10, hjust = 0.5)) + # plot title settings
  
  geom_text(aes(y = ifelse(NES < 0, 0.03, -0.03), label = NAME), # define distance to middle axis
            hjust = ifelse(significant_pathways$NES < 0, 0, 1), # display pathway labels on opposing site 
            size = 2.6) 


# plot
ggp


# save plot
ggsave("GSEA_Plot_CHX.tiff", ggp, width = 5.5, height = 5, dpi = 600)





