# File: functions.R


## Boxplot Functions

### FUNCT: cil_v_sec_boxplot_points()



cil_v_sec_boxplot_points <- function(gene = "FOXJ1", point_size = 3, dataset = all_q_norm){
  set.seed(42)
  vector <- as.vector(assayDataElement(dataset[gene, ], elt = "q_norm"))
  plot <- ggplot(pData(dataset),
                 aes(x = segment,
                     y = vector,
                     fill = segment
                 )
  )+
    geom_boxplot(outlier.shape = NA)+
    labs(title = TeX(paste0("\\textit{", gene, "}")),
         y = TeX(paste0("\\textit{", gene, "} Expression")),
         x = "Cell Type")+
    
    geom_jitter(size = point_size, width = 0.25, stroke = 0.3)+
    scale_fill_manual(values =c("green", "red"))+
    
    
    theme_bw()+
    
    stat_compare_means(
      comparisons = list(c("Ciliated", "Secretory")),
      method = "t.test", 
      bracket.nudge.y = -0.5,
      size = point_size+5)+
    
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))
  
  
  return(plot)
  
}





### FUNCT: make_boxplot_legend()





make_boxplot_legend <- function(dataset = NA, width = 0.5){
  
  plot <- Cil_v_Sec_Anat_Boxplot_Points("ESR1", dataset = dataset)+
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") +
    theme(
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.spacing = unit(20, "pt")
    )
  
  
  # guides(shape = guide_legend(override.aes = list(size = 6)))+
  # theme(legend.key.height= unit(0.8, 'cm'),
  #    legend.key.width= unit(width, 'cm'))
  
  legend <- get_legend(plot)
  
  return(legend)
  
}






### FUNCT: Cil_v_Sec_Anat_Boxplot_Points



Cil_v_Sec_Anat_Boxplot_Points <- function(gene = "FOXJ1", point_size = 3, dataset = all_q_norm){
  
  # %in% checks to see if a name is in a list. 
  # I define not in (%ni%) to do the opposite (check if NOT in the list)
  `%ni%` <- Negate(`%in%`)
  gene_not_found <- c()
  if(gene %ni% rownames(assayDataElement(dataset, elt = "q_norm"))){
    return(paste0(gene, " not in list"))
  }
  vector <- as.vector(assayDataElement(dataset[gene, ], elt = "q_norm"))
  plot <- ggplot(
    pData(dataset), 
    aes(x = region,
        y = vector,
        fill = region
    )
  )+
    geom_boxplot(outlier.shape = NA) +
    labs(title = TeX(paste0("\\textit{", gene, "} Expression")), 
         y = TeX(paste0("\\textit{", gene, "} Expression")),
         x = "Anatomical Region")+
    facet_wrap(~segment)+
    
    geom_jitter(size = point_size, width = 0.25, stroke = 1.3, aes(shape = Patient))+
    
    scale_shape_manual(values = c(3, 4, 13, 21, 22, 23, 24))+
    
    theme_bw()+
    
    # scale_x_discrete(labels=c("Fim","Inf","Amp", "Isth"))+
    scale_fill_manual(values= c("#7CAE00", "#00C1C6", "#F8766D", "#C77CFF"))+
    
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          legend.text = element_text(size =12),
          legend.title = element_text(size = 12),
          strip.text = element_text(size = 15))
  
  
  return(plot)
}




adj_text_size <- function(plot, pnt = 20){
  plot <- plot + theme(
    text = element_text(size = pnt),
    axis.text.x = element_text(size = pnt),
    axis.text.y = element_text(size = pnt),
    axis.title.y = element_text(size = pnt),
    axis.title.x = element_text(size = pnt),
    legend.text = element_text(size = pnt),
    legend.title = element_text(size = pnt),
    strip.text = element_text(size = pnt)
  )
  
  return(plot)
}



### FUNCT: add_sig



my_comparisons = list(
  c("Amp", "Isth"),
  c("Inf", "Amp"),
  c("Fim", "Inf"), 
  c("Fim", "Amp"),
  c("Fim", "Isth")
)


add_sig <- function(plot, comparisons_wanted = my_comparisons, pnt = 5){
  plot <- plot + 
    stat_compare_means(
      comparisons = comparisons_wanted, 
      method = "t.test", 
      bracket.nudge.y = -3,
      size = pnt)+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))
  return (plot)
}


# add_sig_pwc <- function(plot, pnt = 5){
#   plot <- plot + geom_pwc(aes(group = region), method = "t_test", size = pnt)+
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
#   return(plot)
# }






### FUNCT: get_legend replacement function

# Get legend from the cowplot function broke with the latest ggplot update. Here is a temporary fix from 
# 
# https://github.com/wilkelab/cowplot/issues/202. 
# 
# Which may eventually be fixed. Until then, I will try this function. 


get_legend <- function(plot, legend = NULL) {
  
  gt <- ggplotGrob(plot)
  
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  
  indices <- grep(pattern, gt$layout$name)
  
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}




### Functions for creating volcano plots Fig 3



#### FUNCT: volcano_plot

# This function takes a data frame containing log_fold, p_value, and minuslog10_Pvalue columns and plots a volcano plot.
# 
# Optionally, you can label selected genes in purple by providing the plot with a list of genes names in `highlight_genes`.
# 
# Change `flip_log_fold` to `TRUE` to flip which side of the volcano plot is upregulated. 
# 
# To label genes, use the function add_labels below. 


volcano_plot <- function(
    df = NULL,
    highlight_genes = c(), #empty list
    title_name = NULL,
    title_size = 8,
    upreg = NULL,
    downreg = NULL,
    xmin = -4,
    xmax = 4,
    lab_size = 3,
    FC_cutoff = 0.5,
    flip_log_fold = F
){
  if(is.null(upreg) & is.null(downreg)){
    FC_name = TeX("$log_2(FC)$")
  }
  else{
    FC_name = TeX(paste0("$log_2(\\frac{\\mu_{", upreg, "}}{\\mu_{", downreg, "}})$"))
  }
  if(is.null(title_name)){
    title_name = deparse(substitute(df))
  }
  if(flip_log_fold){
    df$log_fold <- -df$log_fold
  }
  df$color <- "NS or FC < 0.5"
  
  df$color[df$log_fold > FC_cutoff & df$p_value < 0.05] <- "Upregulated"
  df$color[df$log_fold < -FC_cutoff & df$p_value < 0.05] <- "Downregulated"
  df$color[df$Marker.name %in% highlight_genes] <- "Markers of Mature Ciliated Cells"
  df$color <- factor(df$color,
                     levels = c("NS or FC < 0.5", "Upregulated", "Downregulated", "Markers of Mature Ciliated Cells"))
  
  genes_to_label <- filter(df, df$color != "NS or FC < 0.5")
  
  
  ggplot(df, aes(x = log_fold, y = minuslog10_Pvalue))+
    geom_point(aes(color = color)) +
    labs(title = title_name, x = FC_name, y = TeX("$-log_{10}(p_{value})$"))+
    xlim(xmin, xmax)+
    
    # coord_cartesian(ylim= c(0, 10))+
    
    scale_color_manual(values = c(`Downregulated` = "blue",
                                  `Upregulated` = "red",
                                  `NS or FC < 0.5` = "gray",
                                  `Markers of Mature Ciliated Cells` = "forestgreen"),
                       guide = guide_legend(override.aes = list(size = 6)))+
    
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "blue")+
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "blue")+
    
    # add lines at p-value cut-offs
    geom_hline(yintercept = 1.3, linetype = "dashed", color = "red")+
    
    
    
    # highlight desired genes
    # must be called after the original geom_point to ensure highlighted genes are visible
    geom_point(data = df[df$Marker.name %in% highlight_genes, ],
               aes(x = log_fold, y = minuslog10_Pvalue),
               color = "forestgreen"
    )+
    
    
    theme_light(base_size = 24, base_family = "sans")+
    
    # remove dumb grid-lines
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(color = "Color")+
    
    theme(plot.title = element_text(size = title_size))
}




#### FUNCT: add_labels

# 
# labels genes in the volcano plots.


add_labels <- function(plot = NA, df = NA, n_genes = 10, method = "manh", highlight_genes = c(), lab_size = 3, flip_log_fold = F){
  if(flip_log_fold){
    df$log_fold <- -df$log_fold
  }
  
  
  # find the top n genes using the selected method
  # return names
  top_n <- top_genes(df = df, sort_method = method, n_genes = n_genes)$Marker.name
  
  # merge the name list, removing duplicate values
  gene_list <- unique(c(top_n, highlight_genes))
  
  # filter dataset to get only genes to label
  
  genes_to_label <- filter(df, df$Marker.name %in% gene_list)
  
  plot_out <- plot +
    
    # Now add labels to both the up and down regulated genes
    
    # for log_fold > 0, plot between 05 and 3
    
    geom_label_repel(data = subset(genes_to_label, log_fold > 0),
                     aes(label = Marker.name,
                         x = log_fold,
                         y = minuslog10_Pvalue),
                     xlim = c(0.5, 3),
                     # xlim = c(1, NA)
                     # # nudge_x = c(0.5, -0.5),
                     # nudge_y = c(0.2, -0.2),
                     # box.padding = unit(0.5, "lines"),
                     # point.padding = unit(.3+3*0.1, "lines"),
                     # max.time = 100,
                     # max.iter = 10^6,
                     size = lab_size,
                     max.overlaps = 40
    )+
    
    # for log_fold < 0, plot -0.5 to -3
    geom_label_repel(data = subset(genes_to_label, log_fold < 0),
                     aes(label = Marker.name,
                         x = log_fold,
                         y = minuslog10_Pvalue),
                     xlim = c(-3, -0.5),
                     # xlim = c(1, NA)
                     # nudge_x = c(0.5, -0.5),
                     # nudge_y = c(0.2, -0.2),
                     # box.padding = unit(0.5, "lines"),
                     # point.padding = unit(.3+3*0.1, "lines"),
                     # max.time = 100,
                     # max.iter = 10^6,
                     size = lab_size,
                     max.overlaps = 20
    )
  
  
  return(plot_out)
  
}





#### FUNCT: top_genes

# Identifies top genes for the volcano plot according to the selected method. 



top_genes <- function(df = NA, sort_method = "manh", n_genes = 10, use_cutoff = T, log_fold_cutoff = 0.5, p_value_cutoff = 0.05){
  
  # remove all rows that do not meet the threshold requirement
  if(use_cutoff){
    df <- df[abs(df$log_fold) >= log_fold_cutoff, ]
    df <- df[df$p_value <= p_value_cutoff, ]
  }
  
  if(sort_method == "manh"){
    # Make a new column called man with the manhattan distance
    df <- df |> mutate(man = abs(log_fold) + abs(minuslog10_Pvalue)) |> arrange(desc(man))}
  else if(sort_method == "euc"){
    # Make a new column called man with the Euclidian distance
    df <- df |> mutate(euc = sqrt(log_fold^2 + minuslog10_Pvalue^2)) |> arrange(desc(euc))
  }
  else if (sort_method == "fc"){
    df <- df |> arrange(desc(abs(log_fold)))
  }
  else if(sort_method == "sig"){
    df <- df |> arrange(desc(minuslog10_Pvalue))
  }
  else{
    message("Not a valid sort method. Select `mnh`, `euc`, `fc`, or `sig`.")
  }
  
  df_out <- head(df, n = n_genes)
  
  return(df_out)
}




### FUNCT: remove_legend()

# This function removes unwanted legends from a ggplot object.



remove_legend <- function(plot){
  plot + theme(legend.position = "none")
}






## 3. Volcano Plot Function 


#### FUNCT: volcano_plot

# This function takes a data frame containing log_fold, p_value, and minuslog10_Pvalue columns and plots a volcano plot.
# 
# Optionally, you can label selected genes in purple by providing the plot with a list of genes names in `highlight_genes`.
# 
# Change `flip_log_fold` to `TRUE` to flip which side of the volcano plot is upregulated. 
# 
# To label genes, use the function add_labels below. 


volcano_plot <- function(
    df = NULL,
    highlight_genes = c(), #empty list
    title_name = NULL,
    title_size = 8,
    upreg = NULL,
    downreg = NULL,
    xmin = -3.5,
    xmax = 3.5,
    lab_size = 3,
    FC_cutoff = 0.5,
    flip_log_fold = F
){
  if(is.null(upreg) & is.null(downreg)){
    FC_name = TeX("$log_2(FC)$")
  }
  else{
    FC_name = TeX(paste0("$log_2(\\frac{\\mu_{", upreg, "}}{\\mu_{", downreg, "}})$"))
  }
  if(is.null(title_name)){
    title_name = deparse(substitute(df))
  }
  if(flip_log_fold){
    df$log_fold <- -df$log_fold
  }
  df$color <- "NS or FC < 0.5"

  df$color[df$log_fold > FC_cutoff & df$p_value < 0.05] <- "Upregulated"
  df$color[df$log_fold < -FC_cutoff & df$p_value < 0.05] <- "Downregulated"
  df$color[df$Marker.name %in% highlight_genes] <- "Highlighted"
  df$color <- factor(df$color,
                        levels = c("NS or FC < 0.5", "Upregulated", "Downregulated", "Highlighted"))
  
  genes_to_label <- filter(df, df$color != "NS or FC < 0.5")
  

  ggplot(df, aes(x = log_fold, y = minuslog10_Pvalue))+
    geom_point(aes(color = color)) +
    labs(title = title_name, x = FC_name, y = TeX("$-log_{10}(p_{value})$"))+
    xlim(xmin, xmax)+

    # coord_cartesian(ylim= c(0, 10))+

    scale_color_manual(values = c(`Downregulated` = "blue",
                                  `Upregulated` = "red",
                                  `NS or FC < 0.5` = "gray",
                                  `Highlighted` = "darkorchid"),
                       guide = guide_legend(override.aes = list(size = 6)))+

    geom_vline(xintercept = 0.5, linetype = "dashed", color = "blue")+
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "blue")+

    # add lines at p-value cut-offs
    geom_hline(yintercept = 1.3, linetype = "dashed", color = "red")+
    

    
    # highlight desired genes
    # must be called after the original geom_point to ensure highlighted genes are visible
    geom_point(data = df[df$Marker.name %in% highlight_genes, ],
               aes(x = log_fold, y = minuslog10_Pvalue),
               color = "darkorchid"
               )+


    theme_light(base_size = 16, base_family = "sans")+

    # remove dumb grid-lines
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(color = "Color")+
      
    theme(plot.title = element_text(size = title_size))
}




### FUNCT: add_labels

add_labels <- function(plot = NA, df = NA, n_genes = 10, method = "manh", highlight_genes = c(), lab_size = 3, flip_log_fold = F){
    if(flip_log_fold){
      df$log_fold <- -df$log_fold
      }
    
  
    # find the top n genes using the selected method
    # return names
    top_n <- top_genes(df = df, sort_method = method, n_genes = n_genes)$Marker.name
    
    # merge the name list, removing duplicate values
    gene_list <- unique(c(top_n, highlight_genes))
    
    # filter dataset to get only genes to label
  
    genes_to_label <- filter(df, df$Marker.name %in% gene_list)
  
    plot_out <- plot +

    # Now add labels to both the up and down regulated genes
      
    # for log_fold > 0, plot between 05 and 3

    geom_label_repel(data = subset(genes_to_label, log_fold > 0),
                     aes(label = Marker.name,
                         x = log_fold,
                         y = minuslog10_Pvalue),
                     xlim = c(0.5, 4),
                     # xlim = c(1, NA)
                    # # nudge_x = c(0.5, -0.5),
                    # nudge_y = c(0.2, -0.2),
                    # box.padding = unit(0.5, "lines"),
                    # point.padding = unit(.3+3*0.1, "lines"),
                    # max.time = 100,
                    # max.iter = 10^6,
                    size = lab_size,
                    max.overlaps = 40
                    )+
      
      # for log_fold < 0, plot -0.5 to -3
      geom_label_repel(data = subset(genes_to_label, log_fold < 0),
                     aes(label = Marker.name,
                         x = log_fold,
                         y = minuslog10_Pvalue),
                     xlim = c(-4, -0.5),
                     # xlim = c(1, NA)
                    # nudge_x = c(0.5, -0.5),
                    # nudge_y = c(0.2, -0.2),
                    # box.padding = unit(0.5, "lines"),
                    # point.padding = unit(.3+3*0.1, "lines"),
                    # max.time = 100,
                    # max.iter = 10^6,
                    size = lab_size,
                    max.overlaps = 20
                    )

    
    return(plot_out)
  
    }




### FUNCT: top_genes


top_genes <- function(df = NA, sort_method = "manh", n_genes = 10, use_cutoff = T, log_fold_cutoff = 0.5, p_value_cutoff = 0.05){
  
  # remove all rows that do not meet the threshold requirement
  if(use_cutoff){
    df <- df[abs(df$log_fold) >= log_fold_cutoff, ]
    df <- df[df$p_value <= p_value_cutoff, ]
  }
  
  if(sort_method == "manh"){
  # Make a new column called man with the manhattan distance
  df <- df |> mutate(man = abs(log_fold) + abs(minuslog10_Pvalue)) |> arrange(desc(man))}
  else if(sort_method == "euc"){
  # Make a new column called man with the Euclidian distance
  df <- df |> mutate(euc = sqrt(log_fold^2 + minuslog10_Pvalue^2)) |> arrange(desc(euc))
  }
  else if (sort_method == "fc"){
    df <- df |> arrange(desc(abs(log_fold)))
  }
  else if(sort_method == "sig"){
    df <- df |> arrange(desc(minuslog10_Pvalue))
  }
  else{
    message("Not a valid sort method. Select `mnh`, `euc`, `fc`, or `sig`.")
  }
  
  df_out <- head(df, n = n_genes)
  
  return(df_out)
}



### FUNCT: get_common_genes()



get_common_genes <- function(df_1 = NA, df_2 = NA, type = "upreg", FC_cutoff = 0.5, p_value_cutoff = 0.05){

  # look for the upreg values (log_fold >= 0.5)
  if(type == "upreg"){
  list_1 <- df_1[(df_1$log_fold >= FC_cutoff & df_1$p_value <= p_value_cutoff), ]$Marker.name
  list_2 <- df_2[(df_2$log_fold >= FC_cutoff & df_2$p_value <= p_value_cutoff), ]$Marker.name
  gene_list <- intersect(list_1, list_2)
  return(gene_list)
  }
    # look for the downreg values (log_fold <= -0.5)
  if(type == "downreg"){
  list_1 <- df_1[(df_1$log_fold <= -FC_cutoff & df_1$p_value <= p_value_cutoff), ]$Marker.name
  list_2 <- df_2[(df_2$log_fold <= -FC_cutoff & df_2$p_value <= p_value_cutoff), ]$Marker.name
  gene_list <- intersect(list_1, list_2)
  return(gene_list)
  }
  
      # look for the downreg values (log_fold <= -0.5)
  if(type == "both"){
  list_1 <- df_1[(df_1$log_fold >= FC_cutoff & df_1$p_value <= p_value_cutoff), ]$Marker.name
  list_2 <- df_2[(df_2$log_fold >= FC_cutoff & df_2$p_value <= p_value_cutoff), ]$Marker.name
  list_up <- intersect(list_1, list_2)
    
  list_3 <- df_1[(df_1$log_fold <= -FC_cutoff & df_1$p_value <= p_value_cutoff), ]$Marker.name
  list_4 <- df_2[(df_2$log_fold <= -FC_cutoff & df_2$p_value <= p_value_cutoff), ]$Marker.name
  
  list_down <- intersect(list_3, list_4)
  
  gene_list <- append(list_up, list_down)
  
  return(gene_list)
    
  }
  
  
  if(!type %in% c("upreg", "downreg", "both")){
    return("Error, type must be `upreg`, `downreg`, or `both`.")
  }
  
}


#### FUNCT: remove_legend



remove_legend <- function(plot){
  plot + theme(legend.position = "none")
}

#### FUNCT: get_legend()


get_legend <- function(plot, legend = NULL) {
  
  gt <- ggplotGrob(plot)
  
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  
  indices <- grep(pattern, gt$layout$name)

  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}






