################################################################################
##Create folder for project organization

if(!dir.exists("input_data")){dir.create("input_data")}
if(!dir.exists("output_data")){dir.create("output_data")}
if(!dir.exists("figures")){dir.create("figures")}
if(!dir.exists("scripts")){dir.create("scripts")}

################################################################################
##load packages 

source(file = "scripts/install_load_packages.r")

################################################################################
##Load data

#Microbiome Data
load("input_data/tse_function.RData")

#Read Removal Efficiency file 
RE = read_excel("input_data/Removal_Efficiency.xlsx")

################################################################################
#plot removal efficiency 
RE$week_year <- paste(RE$Week_Number, RE$Year, sep = "_")

RE$week_year <- factor(RE$week_year, levels = c("43_2014","46_2014","51_2014","2_2015","4_2015","10_2015","14_2015","19_2015","22_2015","26_2015",
                                                "30_2015","34_2015","39_2015","42_2015","46_2015","49_2015",
                                                "52_2015","3_2016","6_2016","12_2016","16_2016","17_2016","22_2016","24_2016",
                                                "27_2016","30_2016","34_2016","39_2016","42_2016","47_2016","52_2016","6_2017") )

RE_long <- pivot_longer(RE, 
                        cols = c("N_removal","BOD_removal","P_removal","Averaged_removal_efficiency"), 
                        names_to = "quality_type", 
                        values_to = "value")

#Create the plot with colorblind-friendly colors
png(filename="figures/removal_efficiency.png" ,units = 'in',width=9, height=6, res=1000)
ggplot(RE_long, aes(x = week_year, y = value, color = quality_type, group = quality_type)) +
  geom_line() +                  # Add lines
  geom_point() +                 # Add dots
  theme_minimal() +              # Use a minimal theme
  labs(x = "Week-Year", y = "Removal Efficiency [%]", color = "Legend") + # Labels
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey90"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) + # Tilt x-axis labels
  scale_color_brewer(palette = "Set1") # Colorblind-friendly palette
dev.off()
################################################################################
##Add removal efficiency to colData in tse_function 

colData(tse_function) <- DataFrame(cbind(colData(tse_function), RE))
################################################################################
##Calculate alpha diversity measures 
##Pre-process data by correlation analysis for each alpha diversity measure

################################################################################
#Richness  
tse_function <- mia::estimateRichness(tse_function, 
                                      assay.type = "counts", 
                                      index =   c("ace", "chao1", "hill", "observed"), 
                                      name=  c("ace", "chao1", "hill", "observed"))

Richness = as.data.frame(colData(tse_function)[ , c("SampleFileName", "ace", "chao1", "hill", "observed")])

Richness_corr <- round(cor(Richness[ ,2:5]), 1)
testRes_rich = cor.mtest(Richness[,2:5],conf.level = 0.95)

#png("figures/corr_richness.png", units="in", width=5, height=5, res=1000)
ggcorrplot(Richness_corr, hc.order = TRUE, type = "upper", lab = TRUE, 
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"), p.mat = testRes_rich$p)
#dev.off()
################################################################################
#Diversity
tse_function <- mia::estimateDiversity(tse_function, 
                                       assay.type = "counts",
                                       index =  c("coverage", "fisher", "gini_simpson", "inverse_simpson",
                                                  "shannon"), 
                                       name =  c("coverage", "fisher", "gini_simpson", "inverse_simpson",
                                                 "shannon"))

Diversity = as.data.frame(colData(tse_function)[ , c("SampleFileName", "coverage", "fisher", "gini_simpson", "inverse_simpson",
                                                     "shannon")])

Diversity_corr <- round(cor(Diversity[,2:6]), 2)
testRes_diversity = cor.mtest(Diversity[,2:6],conf.level = 0.95)

#png("figures/corr_diversity.png", units="in", width=5, height=5, res=1000)
ggcorrplot(Diversity_corr, hc.order = TRUE, type= "upper", lab = TRUE, 
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"), p.mat= testRes_diversity$p)
#dev.off()
################################################################################
#Evenness 
tse_function <- estimateEvenness(tse_function, 
                                 assay.type = "counts", 
                                 index=c("camargo", "pielou", "simpson_evenness", "evar", "bulla"),
                                 name = c("camargo", "pielou", "simpson_evenness", "evar", "bulla"))


Evenness = as.data.frame(colData(tse_function)[ , c("SampleFileName", "camargo", "pielou", "simpson_evenness", "evar", "bulla")])


Evenness_corr <- round(cor(Evenness[,2:6]), 1)
testRes_even = cor.mtest(Evenness[,2:6],conf.level = 0.95)

#png("figures/corr_evenness.png", units="in", width=5, height=5, res=1000)
ggcorrplot(Evenness_corr, hc.order = TRUE, type= "upper", lab = TRUE, 
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"), p.mat = testRes_even$p)
#dev.off()

################################################################################
#Dominance
tse_function <- estimateDominance(tse_function, 
                                  assay.type = "counts", 
                                  index=c("absolute", "dbp", "core_abundance", "gini", "dmn", "relative",
                                          "simpson_lambda"), 
                                  name = c("absolute", "dbp", "core_abundance", "gini", "dmn", "relative",
                                           "simpson_lambda"))

Dominance = as.data.frame(colData(tse_function)[ , c("SampleFileName", "absolute","dbp", "core_abundance",
                                                     "gini", "dmn", "relative",
                                                     "simpson_lambda")])


Dominance_corr <- round(cor(Dominance[,2:8]), 1)
testRes_dominance = cor.mtest(Dominance[,2:8],conf.level = 0.95)

#png("figures/corr_dominance.png", units="in", width=5, height=5, res=1000)
ggcorrplot(Dominance_corr, hc.order = TRUE, type= "upper", lab = TRUE, 
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"), p.mat=testRes_dominance$p)
#dev.off()

################################################################################
##Rarity 
tse_function <- mia::estimateDiversity(tse_function, 
                                       assay.type = "counts",
                                       index = "log_modulo_skewness")

#Another calculation for rarity 
df = data.frame(assay(tse_function,2))
#check
#dim(df)
#sapply(df, sum)
rarity_0.01 = as.data.frame(colSums(df * (df <0.01)))
rownames(rarity_0.01) <- NULL
colnames(rarity_0.01) = "rarity_0.01"
#Same as above based on counts/reads and not relative abundance 
#rarity_0.01_n = colSums(df != 0)

#Get rare species
#x = df * (df <0.01)
#y = cbind(rowData, x)
#You can loop to get all rare phylum 
#y %>%
#  filter( X4717_003 > 0) %>%
#  select(Domain) %>%
#  unique()
Rarity = as.data.frame(colData(tse_function)[ , c("SampleFileName", "log_modulo_skewness")])
Rarity = cbind(Rarity, rarity_0.01)


Rarity_corr <- round(cor(Rarity[,2:3]), 2)
testRes_rarity = cor.mtest(Rarity[,2:3],conf.level = 0.95)

#png("figures/corr_rarity.png", units="in", width=5, height=5, res=1000)
ggcorrplot(Rarity_corr, hc.order = TRUE, type= "upper", lab = TRUE, 
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"), p.mat = testRes_rarity$p)
#dev.off()
################################################################################
#all alpha diversity measures
alpha = as.data.frame(colData(tse_function)[ ,17:40])
alpha = cbind(alpha,Rarity[ , c("rarity_0.01")])
colnames(alpha)[colnames(alpha) == "Rarity[, c(\"rarity_0.01\")]"] <- "rarity_0.01"
################################################################################
##Divergence 

tse_function <- addDivergence(tse_function)
#colnames(colData(tse))

Divergence  = as.data.frame(colData(tse_function)[ , c("SampleFileName", "divergence")])


alpha = cbind(alpha,Divergence[2])

################################################################################
#Checking correlation in the resulting alpha diversity measure to check that 
# other types of alpha diversity don't correlate to each other 

alpha = subset(alpha, select = -c(ace_se, chao1_se) )

alpha_corr <- round(cor(alpha), 1)
testRes_alpha = cor.mtest(alpha,conf.level = 0.95)

#png("figures/corr_alpha_all.png", units="in", width=5, height=5, res=1000)
ggcorrplot(alpha_corr,hc.order = TRUE,lab = TRUE,lab_size = 2.5,
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"),
           tl.cex = 10, p.mat = testRes_alpha$p)
#dev.off()

################################################################################
#cross correlation analysis 

total_corr = round(cor(alpha), 1)

png(filename="figures/corr_alpha_all.png" ,units = 'in',width=8, height=6, res=1000)
corrplot(total_corr, order ='hclust', hclust.method ='single' ,addrect = 3, rect.col = "green",
         rect.lwd = 3 ,method = 'circle', addCoef.col = 'white',
         number.digits = 3,number.cex = 0.3,tl.pos ='lt', tl.srt=45, tl.cex = 0.5)
dev.off()

alpha_clean = alpha[ ,c("observed","log_modulo_skewness","gini_simpson", "shannon")]
################################################################################
plot_list <- list()
for (i in 1:length(alpha_clean)){ 
  #Add name of the alpha diversity used in CCF
  n=colnames(alpha_clean)
  text = paste("cross correlation function for", n[i] , "and removal efficiency")
  #CCF
  plot <- ggCcf(alpha_clean[ ,i],
                RE$Averaged_removal_efficiency,
                type = "correlation",
                na.action = na.contiguous) +
    theme_minimal() + 
    scale_x_continuous(limits = c(-7, 0), breaks = seq(-7,0,1)) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.2)) +
    labs(title=text,
         x ="Lag", y = "Correlation")+
    theme(plot.title = element_text(size = 9))
  plot_list[[i]] <- plot
  #ggplotly(plot)
  #Save plot 
  cff_name = paste('figures/','ccf_', n[i],".png")
  png(cff_name ,units = 'in',width=5, height=5, res=1000)
  print(plot)
  dev.off()
  
}


# Combine all plots into one figure
combined_plot <- grid.arrange(grobs = plot_list, ncol = 2)  # Adjust ncol for desired layout

# Save the combined plot
png("figures/combined_ccf_plots.png", units = 'in', width = 10, height = 6, res = 1000)
print(grid.arrange(grobs = plot_list, ncol = 2))
dev.off()


write.csv(as.data.frame(colData(tse_function)[ , c("Sample_ID","SampleFileName","SRA_accession","Biosample_accession", "gini_simpson")]), "output_data/alpha_diversity.csv", row.names=FALSE)










