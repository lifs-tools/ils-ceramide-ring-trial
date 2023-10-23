
# Plotting all samples (with replicate StDev as error bars) of all labs. 
plot_labscatter <- function(d_samples, d_assays, variable, sample_type, labid_col, normalisation_sample = NULL, 
                            excluded_labs = "", save_plot = TRUE, C16_max = NA, C18_max = NA, C24_max = NA, C241_max = NA, filename_prefix){
  
  var_meanreplicates <- paste0(variable, "_mean")
  
  d_samples <- d_samples |> filter(!(LabId %in% excluded_labs))
  d_assays <- d_assays |> filter(!(LabId %in% excluded_labs))
  
  if(!is.null(normalisation_sample)) {
    
    if(variable == "C_Adj"){
      d_stat_filt <- get_study_stats(flag_outlier(d_assays) |> filter(!Outlier_mp))
    } else if(variable == "C_SinglePoint"){
      d_stat_filt <- get_study_stats(flag_outlier(d_assays) |> filter(!Outlier_sp))
    } else {stop("Variable not present or supported!")}
    
    d_assays <- d_assays |>
      left_join(d_stat_filt |> filter(SampleType == "SRM") |> select(ceramideName, C_SinglePoint_interlabMEAN, C_Adj_interlabMEAN)) |> 
      group_by(ceramideName, LabId) |>
      mutate(!!ensym(var_meanreplicates) := !!ensym(var_meanreplicates)/median((!!ensym(var_meanreplicates))[SampleType == normalisation_sample]) * (!!sym(paste0(variable, "_interlabMEAN")))) 
    
    
    d_samples <- d_samples |> 
      left_join(d_stat_filt |> filter(SampleType == "SRM") |> select(ceramideName, C_SinglePoint_interlabMEAN, C_Adj_interlabMEAN)) |> 
      group_by(ceramideName, LabId, Sample) |> 
      mutate(lab_mean_conc = median((!!ensym(var_meanreplicates))[SampleType == normalisation_sample])) |> 
      group_by(ceramideName, LabId) |>  
      mutate(!!ensym(var_meanreplicates) := !!ensym(var_meanreplicates)/mean(lab_mean_conc[SampleType == normalisation_sample]) * (!!sym(paste0(variable, "_interlabMEAN")))) 
    
    View(d_samples)
  }
  
  if(variable == "C_Adj"){
    d_assays_filt <- flag_outlier(d_assays) |> filter(!Outlier_mp)
    d_stat_all <- get_study_stats(d_assays |> filter(SampleType == sample_type))
    d_stat_filt <- get_study_stats(d_assays_filt |> filter(SampleType == sample_type))
    conc_text <- "calibration curve"
  } else if(variable == "C_SinglePoint"){
    d_assays_filt <- flag_outlier(d_assays) |> filter(!Outlier_sp)
    d_stat_all <- get_study_stats(d_assays |> filter(SampleType == sample_type))
    d_stat_filt <- get_study_stats(d_assays_filt |> filter(SampleType == sample_type))
    conc_text <- "ISTD"
  } else {
    stop("Variable not present or supported!")
  }
  
  
  d_samples_subset <- d_samples |> filter(SampleType == sample_type) 
  set.seed(1334)
  pos_points <- position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2 )
  
  x_labels <- if(labid_col == "LabNo") x_labels <- d_samples_subset$LabNo |> unique() else x_labels <- d_samples_subset$LabNum |> unique() |> str_pad( 2, pad = "0")
  
  plt2 <- ggplot(data = d_samples_subset , mapping = aes(x = LabNum+0.5, y = !!ensym(var_meanreplicates), group = MethodNo, color = Protocol, fill = Protocol)) +
    geom_rect(data = d_stat_filt, 
              mapping = aes(ymin =  !!sym(paste0(variable, "_interlabMEAN")) - 2 *  !!sym(paste0(variable, "_interlabSD")), 
                            ymax =  !!sym(paste0(variable, "_interlabMEAN")) + 2 *  !!sym(paste0(variable, "_interlabSD"))), 
              xmin = 0, xmax = 38, inherit.aes = FALSE, fill = "grey70", alpha = 0.046) +
    geom_hline(data = d_stat_filt, 
               mapping = aes(yintercept = !!sym(paste0(variable, "_interlabMEAN"))), 
               linewidth = 0.4, color = "grey70") +  
    geom_hline(data = d_stat_filt, 
               mapping = aes(yintercept = !!sym(paste0(variable, "_interlabMEAN")) + 2 * !!sym(paste0(variable, "_interlabSD"))), 
               linewidth = 0.4, color = "#67878f", linetype = "dashed") +  
    geom_hline(data = d_stat_filt, 
               mapping = aes(yintercept = !!sym(paste0(variable, "_interlabMEAN")) - 2 * !!sym(paste0(variable, "_interlabSD"))), 
               linewidth = 0.4, color = "#67878f", linetype = "dashed") +  
    geom_hline(data = d_stat_all, 
               mapping = aes(yintercept = !!sym(paste0(variable, "_limitQ1"))), 
               linewidth = 0.4, color = "#34bf44", linetype = "dotted") +  
    geom_hline(data = d_stat_all, 
               mapping = aes(yintercept = !!sym(paste0(variable, "_limitQ3"))), 
               linewidth = 0.4, color = "#34bf44", linetype = "dotted") +  
    geom_linerange(aes(x = LabNum+0.5, 
                       ymin = !!sym(paste0(variable, "_mean")) - !!sym(paste0(variable, "_SD")), 
                       ymax = !!sym(paste0(variable, "_mean")) + !!sym(paste0(variable, "_SD"))), 
                    position = pos_points, size = .3,  alpha = 0.25) + 
    stat_summary(position = pos_points,
                 fun = "mean",
                 geom = "crossbar",inherit.aes = TRUE,
                 width = .7,
                 #color = "darkred",
                 linewidth = 0.2, alpha = 0.7)+
    geom_point(size = 1.5, shape = 21, position = pos_points, alpha = 0.5) + 
    facet_wrap(vars(ceramideName), ncol = 1, scales = "free") +
    scale_y_continuous(limits = c(0,NA))+
    scale_x_continuous(limits = c(.5,NA), 
                       breaks = c(1:(length(d_assays$LabNo |> unique()))), 
                       labels = x_labels,
                       expand = expansion(mult = c(0,0.01)))+
    scale_fill_manual(values = c("SOP" = "#faafaf", "OTHER" = "#b4c5fa")) +
    scale_color_manual(values = c("SOP" = "#e30202", "OTHER" = "darkblue"))+
    labs(x = "Laboratory", y = glue("\U003BCmol/L ({conc_text})")) + 
    labs(color=NULL, fill = NULL)+
    theme_classic(base_size = 8) +  
    ggh4x::facetted_pos_scales(
      y = list(
        ceramideName == "Cer 18:1;O2/16:0" ~ scale_y_continuous(limits = c(0, C16_max)),
        ceramideName == "Cer 18:1;O2/18:0" ~ scale_y_continuous(limits = c(0, C18_max)),
        ceramideName == "Cer 18:1;O2/24:0" ~ scale_y_continuous(limits = c(0, C24_max)),
        ceramideName == "Cer 18:1;O2/24:1" ~ scale_y_continuous(limits = c(0, C241_max))
      )) +
    theme(
      axis.text.x = element_text(size = 7, angle = 0, vjust = 1.5, hjust = -0.25), 
      panel.grid.major.x = element_line(colour = "grey70", linewidth = 0.2, linetype = "dotted"),
      strip.background =element_rect(fill="grey99")
    )  +
    theme(legend.position = c(.94,.97), 
          legend.key.size = unit(0.2, "cm"),
          legend.spacing.y = unit(.001, 'cm'),
          legend.text=element_text(size=8)) 
  
  if(save_plot) ggsave(plot = plt2, filename = here("manuscript/output", paste0(filename_prefix, "LabScatterPlot_", sample_type, "_", variable, if_else(is.null(normalisation_sample), "", paste0("_norm_by_",sample_type)),".png")), device = "png",
                       dpi = 600, width = 160, height =210, units = "mm")
  return(list(plt = plt2, stats = d_stat_filt))
}



# Plotting all samples (with replicate StDev as error bars) of all labs. 
plot_comparison_norm <- function(d_assays, variable, sample_type, normalisation_sample = NULL, 
                            excluded_labs = "", save_plot = TRUE, filename_prefix){
  
  var_meanreplicates <- paste0(variable, "_mean")
  
  d_assays <- d_assays |> filter(!(LabId %in% excluded_labs))
  
  if(variable == "C_Adj"){
  d_stat_filt <- get_study_stats(flag_outlier(d_assays) |> filter(!Outlier_mp))
  } else if(variable == "C_SinglePoint"){
  d_stat_filt <- get_study_stats(flag_outlier(d_assays) |> filter(!Outlier_sp))
  } else {stop("Variable not present or supported!")}
  
  d_assays <- d_assays |>
  left_join(d_stat_filt |> filter(SampleType == "SRM") |> select(ceramideName, C_SinglePoint_interlabMEAN, C_Adj_interlabMEAN)) |> 
  group_by(ceramideName, LabId) |>
  mutate(!!sym(paste0(var_meanreplicates, "_norm")) := !!sym(var_meanreplicates)/median((!!sym(var_meanreplicates))[SampleType == normalisation_sample]) * (!!sym(paste0(variable, "_interlabMEAN")))) 
  

  d_plot <- d_assays |> 
    filter(SampleType %in% sample_type) |>
    select(SampleType, ceramideName, LabId, Protocol, C_Adj_mean, C_Adj_mean_norm) |> 
    pivot_longer(cols = C_Adj_mean:C_Adj_mean_norm, names_to = "norm_type",values_to = "Concentration")

  d_plot_stats <- d_plot |> 
    group_by(SampleType, ceramideName) |> 
    summarise(
      MEAN = mean(Concentration[norm_type=="C_Adj_mean"]),
      MEAN_norm = mean(Concentration[norm_type=="C_Adj_mean_norm"]),
      CV = sd(Concentration[norm_type=="C_Adj_mean"])/mean(Concentration[norm_type=="C_Adj_mean"])*100,
      CV_norm = sd(Concentration[norm_type=="C_Adj_mean_norm"])/mean(Concentration[norm_type=="C_Adj_mean_norm"])*100
    )
  
  write_csv(x = d_plot_stats, file = here("manuscript/output/Stats_norm.csv"))
  
  
    
  set.seed(1334)
  pos_points <- position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2 )

  plt <- ggplot(d_plot |> arrange(norm_type),
                   aes(x = norm_type, y = Concentration, color = Protocol, fill = Protocol, group = LabId))+
    #geom_boxplot(aes(x = ISTD, y = Concentration),width = .3, linewidth = 0.2, color = "red", outlier.shape = NA, alpha =0.5) +
    geom_point(position = pos, size = 1.5, shape = 21, stroke = 0.2,alpha = 0.9)  +

    #ggbeeswarm::geom_beeswarm(cex = 2.5,size = 1.2, stroke = .2, method = "swarm", priority = "random", alpha = .7, shape = 21) +
    stat_summary(aes(x = norm_type, y = Concentration), geom = "errorbar",inherit.aes = FALSE,
                 width = .3,
                 color = "darkred",
                 linewidth = 0.25,
                 alpha = 0.6,
                 fun.min = \(x) mean(x) - sd(x),
                 fun.max = \(x) mean(x) + sd(x)) +

    stat_summary(aes(x = norm_type, y = Concentration),
                 fun = "mean",
                 geom = "crossbar",inherit.aes = FALSE,
                 width = .4,
                 color = "darkred",
                 linewidth = 0.2, alpha = 0.7)+


    facet_wrap(vars(ceramideName), scales = "free", nrow = 1) +
    geom_line(aes(group = LabId), position = pos, alpha = 0.2, linewidth = 0.2, linetype = "longdash") +
    scale_fill_manual(values = c("SOP" = "#ffc2c2", "OTHER" = "#d1dcff"), labels=c("SOP", "OTHER")) +
    scale_color_manual(values = c("SOP" = "#e30202", "OTHER" = "darkblue"), labels=c("SOP", "OTHER")) +
    expand_limits(y=0) +
    
    scale_x_discrete(labels = c("C_Adj_mean" = "None",
                                "C_Adj_mean_norm" = "SRM1950")) +

    #scale_linetype_manual(values = c("SOP" = "solid", "OTHER" = "dotted"))+
    labs(color=NULL, fill = NULL)+
    theme_classic(base_size = 7) +
    labs(y="\U003BCmol/L (calibration curve)", x="Recalibration")+
    #scale_x_discrete(drop = FALSE, expand = expansion(mult = c(.03, .05))) +
    theme(
      axis.text.x = element_text(size = 6, face = "plain", angle = 0),
      axis.text.y = element_text(size = 6, face = "plain", angle = 0),
      axis.title.x = element_text(size = 7, face = "bold", angle = 0, margin = margin(t = 8, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size = 6, face = "bold", angle = 90, margin = margin(t = 0, r = 7, b = 0, l = 8)),
      #panel.grid.major.x = element_line(colour = "grey70", linewidth = 0.2, linetype = "dotted"),
      strip.background = element_rect(fill="grey99")) +
     theme(legend.position = c(.96,.08),
          legend.key.size = unit(0.2, "cm"),
          legend.text =  element_text(size = 5, face = "plain", angle = 0),
          legend.spacing.y = unit(.001, "cm"))
  
  plt
}


