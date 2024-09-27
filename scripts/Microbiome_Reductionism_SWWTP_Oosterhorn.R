################################################################################
##Create folder for project organization
if(!dir.exists("input_data")){dir.create("input_data")}
if(!dir.exists("output_data")){dir.create("output_data")}
if(!dir.exists("figures")){dir.create("figures")}
if(!dir.exists("scripts")){dir.create("scripts")}

################################################################################
##load packages 
suppressMessages({
  suppressWarnings({
    source(file = "scripts/install_load_packages.r")
  })
})
################################################################################
##Load data

#Microbiome Data
load("input_data/tse_function.RData")

#Read Removal Efficiency file 
RE = read_excel("input_data/Removal_Efficiency.xlsx")

################################################################################
#Number of unique functional genes existing in the 32 samples  
length(unique(rowData(tse_function)$Order))

#get top 10 functional genes with their respective taxonomy mapping 
tse_genes <- agglomerateByRank(tse_function, rank ="Order")
top_genes <- getTopFeatures(tse_genes, top = 10, assay.type = "relabundance")
# Create an empty data frame to store results
tax_map_table <- data.frame(Kingdom = character(),
                            Phylum = character(),
                            Class = character(),
                            Order = character(),
                            stringsAsFactors = FALSE)
for (i in top_genes){
  i <- sub("^Order:", "", i)
  tax_map = mapTaxonomy(tse_genes, taxa = i)
  tax_map_table  <- rbind(tax_map_table , data.frame(
    Kingdom = as.data.frame(tax_map)[ ,1],
    Phylum = as.data.frame(tax_map)[ ,2],
    Class = as.data.frame(tax_map)[ ,3],
    Order = as.data.frame(tax_map)[ ,4],
    stringsAsFactors = FALSE
  ))
}

write_xlsx(tax_map_table, "output_data/taxonomy_top10_genes.xlsx")  

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
##DIVE: Diversity-Informed Valuation of Ecosystem functionality

alpha = dive(tse_function,as.data.frame(RE$Averaged_removal_efficiency))

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

#no title plots for paper
suppressMessages({
  suppressWarnings({
    for (i in 1:length(alpha)){ 
      # CCF
      plot <- ggCcf(alpha[ ,i],
                    RE$Averaged_removal_efficiency,
                    type = "correlation",
                    na.action = na.contiguous) +
        theme_minimal() + 
        scale_x_continuous(limits = c(-4, 0), breaks = seq(-4, 0, 1)) +
        scale_y_continuous(limits = c(-0.45, 0.45), breaks = seq(-0.45, 0.45, 0.1)) +
        labs(x = "Lag", y = "Correlation Coefficient")  +  # Set title to NULL
        theme(plot.title = element_blank())
      plot_list[[i]] <- plot
    }
    
    
    combined_plot <- plot_grid(
        plotlist = plot_list[c(7,11,13,21)], 
        labels = c("a", "b", "c", "d","e","f","g"),   # Add labels
        ncol = 2                          # Specify the number of columns
    )
    
    # Save the combined plot
    png("figures/combined_ccf_plots.png", units = 'in', width = 10, height = 6, res = 1000)
    print(plot_grid(plotlist = plot_list[c(7,11,13,21)], labels = c("a", "b", "c", "d","e","f","g"),ncol = 2))
    dev.off()
  })
})

write.csv(as.data.frame(colData(tse_function)[ , c("Sample_ID","SampleFileName","SRA_accession","Biosample_accession", "gini_simpson")]), "output_data/alpha_diversity.csv", row.names=FALSE)
################################################################################
