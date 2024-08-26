get_study_stats <- function(d, limit_iqr = 1.5){ 
  d_res <- d %>%
    group_by(ceramideName, SampleType) |> 
    dplyr::summarize(
      n = n(),
      n_SOP = length(RatioLipidToIS_mean[Protocol == "SOP"]),
      RatioLipidToIS_interlabMEAN = mean(RatioLipidToIS_mean, na.rm = FALSE),
      RatioLipidToIS_interlabSD = sd(RatioLipidToIS_mean, na.rm = FALSE), 
      RatioLipidToIS_interlabCV = RatioLipidToIS_interlabSD/RatioLipidToIS_interlabMEAN *100,
      RatioLipidToIS_median_intralabCV = median(RatioLipidToIS_intralabCV, na.rm = FALSE),
      RatioLipidToIS_interlabMEDIAN = median(RatioLipidToIS_mean, na.rm = FALSE),      
      RatioLipidToIS_IQR = IQR(RatioLipidToIS_mean, na.rm = FALSE),
      RatioLipidToIS_limitQ1 = quantile(RatioLipidToIS_mean, 0.25, na.rm = FALSE) -  (RatioLipidToIS_IQR * limit_iqr),
      RatioLipidToIS_limitQ3 = quantile(RatioLipidToIS_mean, 0.75, na.rm = FALSE) +  (RatioLipidToIS_IQR * limit_iqr),
      C_SinglePoint_interlabMEAN = mean(C_SinglePoint_mean, na.rm = TRUE),
      C_SinglePoint_interlabSEM = sd(C_SinglePoint_mean, na.rm = TRUE)/sqrt(n), 
      C_SinglePoint_interlabSD = sd(C_SinglePoint_mean, na.rm = TRUE), 
      C_SinglePoint_interlabCV = C_SinglePoint_interlabSD/C_SinglePoint_interlabMEAN *100,
      C_SinglePoint_median_intralabCV = median(C_SinglePoint_intralabCV, na.rm = FALSE),
      C_SinglePoint_IQR = IQR(C_SinglePoint_mean, na.rm = TRUE),
      C_SinglePoint_limitQ1 = quantile(C_SinglePoint_mean, 0.25, na.rm = TRUE) -  (C_SinglePoint_IQR * limit_iqr),
      C_SinglePoint_limitQ3 = quantile(C_SinglePoint_mean, 0.75, na.rm = TRUE) +  (C_SinglePoint_IQR * limit_iqr),
      C_Adj_interlabMEAN = mean(C_Adj_mean, na.rm = TRUE),
      C_Adj_interlabMEDIAN = median(C_Adj_mean, na.rm = TRUE),
      C_Adj_interlabSD = sd(C_Adj_mean, na.rm = TRUE), 
      C_Adj_interlabMAD = mad(C_Adj_mean, na.rm = TRUE), 
      C_Adj_interlabSEM = sd(C_Adj_mean, na.rm = TRUE)/sqrt(n), 
      C_Adj_interlabCV = C_Adj_interlabSD/C_Adj_interlabMEAN *100,
      C_Adj_interlabRCV = C_Adj_interlabMAD/C_Adj_interlabMEDIAN *100,
      C_Adj_C16_interlabMEAN = mean(C_Adj_C16_mean, na.rm = TRUE),
      C_Adj_C16_interlabSD = sd(C_Adj_C16_mean, na.rm = TRUE), 
      C_Adj_C16_interlabCV = C_Adj_C16_interlabSD/C_Adj_C16_interlabMEAN *100,
      C_Adj_C18_interlabMEAN = mean(C_Adj_C18_mean, na.rm = TRUE),
      C_Adj_C18_interlabSD = sd(C_Adj_C18_mean, na.rm = TRUE), 
      C_Adj_C18_interlabCV = C_Adj_C18_interlabSD/C_Adj_C18_interlabMEAN *100,
      C_Adj_C24_interlabMEAN = mean(C_Adj_C24_mean, na.rm = TRUE),
      C_Adj_C24_interlabSD = sd(C_Adj_C24_mean, na.rm = TRUE), 
      C_Adj_C24_interlabCV = C_Adj_C24_interlabSD/C_Adj_C24_interlabMEAN *100,
      C_Adj_C241_interlabMEAN = mean(C_Adj_C241_mean, na.rm = TRUE),
      C_Adj_C241_interlabSD = sd(C_Adj_C241_mean, na.rm = TRUE), 
      C_Adj_C241_interlabCV = C_Adj_C241_interlabSD/C_Adj_C241_interlabMEAN *100,
      C_Adj_IQR = IQR(C_Adj_mean, na.rm = TRUE),
      C_Adj_limitQ1 = quantile(C_Adj_mean, 0.25, na.rm = TRUE) -  (C_Adj_IQR * limit_iqr),
      C_Adj_limitQ3 = quantile(C_Adj_mean, 0.75, na.rm = TRUE) +  (C_Adj_IQR * limit_iqr),
      C_Adj_median_intralabCV = median(C_Adj_intralabCV, na.rm = TRUE)) |> 
    ungroup() |> 
    arrange(SampleType)
  d_res
}

# Flag outliers, based on Tukey’s 3xIQR
flag_outlier <- function(data, include_calibdata, limit_iqr = 1.5){
  data <- data |>
    group_by(ceramideName, SampleType) %>% 
    mutate(
      IQR_sp = IQR(C_SinglePoint_mean, na.rm = TRUE),
      Q1_sp = quantile(C_SinglePoint_mean, 0.25, na.rm = TRUE),
      Q3_sp = quantile(C_SinglePoint_mean, 0.75, na.rm = TRUE), 
      IQR_mp = IQR(C_Adj_mean, na.rm = TRUE),
      Q1_mp = quantile(C_Adj_mean, 0.25, na.rm = TRUE),
      Q3_mp = quantile(C_Adj_mean, 0.75, na.rm = TRUE), 
      Outlier_sp = !between(C_SinglePoint_mean,(Q1_sp - limit_iqr*IQR_sp),(Q3_sp + limit_iqr*IQR_sp)), 
      Outlier_mp = !between(C_Adj_mean,(Q1_mp - limit_iqr*IQR_mp),(Q3_mp + limit_iqr*IQR_mp))
    ) |> 
    ungroup()
  data
}

sig_annot = function(x){
  #if(x < 1e-9) {"p < 1e-9"}
  if(x < 1e-3) {paste0("p = ", formatC(x, format = "e", digits = 2))}
  else if(x <= 1) {paste0("p = ", formatC(x, digits = 3))}      
  else{NA}
}

# formula standard uncertainty (Ghorasaini et al, Analytical Chemistry, 2021)
std_uncert <-
  function(x)
    sqrt(pi / (2 * length(x))) * mad(x, constant = 1.4826, na.rm = TRUE)

labScatterPlot <- function(sample_type, d_sample_mean, d_srm_mean, variable = "C_Adj", normalisation_sample = "SRM", name_prefix) {
  print(paste("Creating plots of", sample_type, "normalized on SRM 1950"))
  C16_max <- NA
  C18_max <- NA
  C24_max <- NA
  C241_max <- NA
  if (variable=="C_Adj") {
    C16_max <- d_sample_mean |> filter(ceramideName=="Cer 18:1;O2/16:0") |> summarise(C16_max=max(na.omit(C_Adj_mean))) |> unlist()
    C18_max <- d_sample_mean |> filter(ceramideName=="Cer 18:1;O2/18:0") |> summarise(C18_max=max(na.omit(C_Adj_mean))) |> unlist()
    C24_max <- d_sample_mean |> filter(ceramideName=="Cer 18:1;O2/24:0") |> summarise(C24_max=max(na.omit(C_Adj_mean))) |> unlist()
    C241_max <- d_sample_mean |> filter(ceramideName=="Cer 18:1;O2/24:1") |> summarise(C241_max=max(na.omit(C_Adj_mean))) |> unlist()
  } else {
    C16_max <- d_sample_mean |> filter(ceramideName=="Cer 18:1;O2/16:0") |> summarise(C16_max=max(na.omit(C_SinglePoint_mean))) |> unlist()
    C18_max <- d_sample_mean |> filter(ceramideName=="Cer 18:1;O2/18:0") |> summarise(C18_max=max(na.omit(C_SinglePoint_mean))) |> unlist()
    C24_max <- d_sample_mean |> filter(ceramideName=="Cer 18:1;O2/24:0") |> summarise(C24_max=max(na.omit(C_SinglePoint_mean))) |> unlist()
    C241_max <- d_sample_mean |> filter(ceramideName=="Cer 18:1;O2/24:1") |> summarise(C241_max=max(na.omit(C_SinglePoint_mean))) |> unlist()
  }
  
  base_name <- glue("{name_prefix}_{variable}_{sample_type}_combined_beforeafterNorm_to_{normalisation_sample}")
  p1 <- plot_labscatter(d_sample_mean = d_sample_mean, d_srm_mean = d_srm_mean, variable = variable, labid_col = "LabNum", sample_type = sample_type, normalisation_sample = NULL, excluded_labs = "", save_plot = FALSE, C16_max = C16_max, C18_max = C18_max, C24_max = C24_max, C241_max = C241_max, base_name)
  p2 <- plot_labscatter(d_sample_mean = d_sample_mean, d_srm_mean = d_srm_mean, variable = variable, labid_col = "LabNum", sample_type = sample_type, normalisation_sample = normalisation_sample, excluded_labs = "", save_plot = FALSE, C16_max = C16_max, C18_max = C18_max, C24_max = C24_max, C241_max = C241_max, base_name)
  
  p_comb <- p1$plt + p2$plt
  
  ggsave(plot = p_comb, filename = here("manuscript/output", glue("{base_name}.png")), device = "png", dpi = 600, width = 290, height =210, units = "mm")
  list(p_comb=p_comb, p1=p1, p2=p2)
}

# Plotting all samples (with replicate StDev as error bars) of all labs. 
plot_labscatter <- function(d_sample_mean, d_srm_mean, variable, sample_type, labid_col, normalisation_sample = NULL, 
                            excluded_labs = "", save_plot = TRUE, C16_max = NA, C18_max = NA, C24_max = NA, C241_max = NA, filename_prefix, as_pdf = FALSE){
  
  var_meanreplicates <- paste0(variable, "_mean")
  
  d_sample_mean <- d_sample_mean |> filter(!(LabId %in% excluded_labs))
  d_srm_mean <- d_srm_mean |> filter(!(LabId %in% excluded_labs))
  
  if(!is.null(normalisation_sample)) {
    
    if(variable == "C_Adj"){
      d_stat_filt <- get_study_stats(flag_outlier(d_srm_mean) |> filter(!Outlier_mp))
    } else if(variable == "C_SinglePoint"){
      d_stat_filt <- get_study_stats(flag_outlier(d_srm_mean) |> filter(!Outlier_sp))
    } else {stop("Variable not present or supported!")}
    
    d_srm_mean <- d_srm_mean |>
      left_join(d_stat_filt |> filter(SampleType == "SRM") |> select(ceramideName, C_SinglePoint_interlabMEAN, C_Adj_interlabMEAN)) |> 
      group_by(ceramideName, LabId) |>
      mutate(!!ensym(var_meanreplicates) := !!ensym(var_meanreplicates)/median((!!ensym(var_meanreplicates))[SampleType == normalisation_sample]) * (!!sym(paste0(variable, "_interlabMEAN")))) 
    
    
    d_sample_mean <- d_sample_mean |> 
      left_join(d_stat_filt |> filter(SampleType == "SRM") |> select(ceramideName, C_SinglePoint_interlabMEAN, C_Adj_interlabMEAN)) |> 
      group_by(ceramideName, LabId, Sample) |> 
      mutate(lab_mean_conc = median((!!ensym(var_meanreplicates))[SampleType == normalisation_sample])) |> 
      group_by(ceramideName, LabId) |>  
      mutate(!!ensym(var_meanreplicates) := !!ensym(var_meanreplicates)/mean(lab_mean_conc[SampleType == normalisation_sample]) * (!!sym(paste0(variable, "_interlabMEAN")))) 
  }
  
  if(variable == "C_Adj"){
    d_srm_mean_filt <- flag_outlier(d_srm_mean) |> filter(!Outlier_mp)
    d_stat_all <- get_study_stats(d_srm_mean |> filter(SampleType == sample_type))
    d_stat_filt <- get_study_stats(d_srm_mean_filt |> filter(SampleType == sample_type))
    conc_text <- "\U003BCmol/L (multi-point calibration)"
  } else if(variable == "C_SinglePoint"){
    d_srm_mean_filt <- flag_outlier(d_srm_mean) |> filter(!Outlier_sp)
    d_stat_all <- get_study_stats(d_srm_mean |> filter(SampleType == sample_type))
    d_stat_filt <- get_study_stats(d_srm_mean_filt |> filter(SampleType == sample_type))
    conc_text <- "\U003BCmol/L (single-point calibration)"
  } else {
    stop("Variable not present or supported!")
  }
  
  d_sample_mean_subset <- d_sample_mean |> filter(SampleType == sample_type) 
  
  d_sample_mean_subset <- d_sample_mean_subset |> 
    select(
      SampleType,	
      ceramideName,
      Sample,
      LabId,
      LabNum,
      MethodNo,
      Protocol,
      C_SinglePoint_mean,
      C_SinglePoint_SD,
      C_Adj_mean,
      C_Adj_SD)
  
  
  d_stat_filt_csv <- d_stat_filt |> select(SampleType,
                                           ceramideName, 
                                           n, 
                                           n_SOP,
                                           MEAN = C_Adj_interlabMEAN, 
                                           SD=C_Adj_interlabSD, 
                                           CV_interlab = C_Adj_interlabCV, 
                                           `Q1-1.5×IQR`  = C_Adj_limitQ1,	
                                           `Q3+1.5×IQR`  = C_Adj_limitQ3)
  
  readr::write_excel_csv(d_sample_mean_subset, file = here(paste0("manuscript/output/", filename_prefix, "_sourcedata_points.csv")))
  readr::write_excel_csv(d_stat_filt_csv, file = here(paste0("manuscript/output/", filename_prefix, "_sourcedata_stats.csv")))
  
  
  set.seed(1334)
  pos_points <- position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2 )
  
  x_labels <- if(labid_col == "LabId_group") x_labels <- d_sample_mean_subset$LabId_group |> unique() else x_labels <- d_sample_mean_subset$LabNum |> unique() |> str_pad( 2, pad = "0")
  
  plt2 <- ggplot(data = d_sample_mean_subset , mapping = aes(x = LabNum+0.5, y = !!ensym(var_meanreplicates), group = MethodNo, color = Protocol, fill = Protocol)) +
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
                    position = pos_points, linewidth = .3,  alpha = 0.25) + 
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
                       breaks = c(1:(length(d_srm_mean$LabId_group |> unique()))), 
                       labels = x_labels,
                       expand = expansion(mult = c(0,0.01)))+
    scale_fill_manual(values = c("SOP" = "#faafaf", "OTHER" = "#b4c5fa")) +
    scale_color_manual(values = c("SOP" = "#e30202", "OTHER" = "darkblue"))+
    labs(x = "Laboratory", y = glue("{conc_text}")) + 
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
    theme(legend.position.inside = c(.94,.97), 
          legend.key.size = unit(0.2, "cm"),
          legend.spacing.y = unit(.001, 'cm'),
          legend.text=element_text(size=8)) 
  
  
  if(save_plot) {
    if (as_pdf) 
      ggsave(plot = plt2, filename = here("manuscript/output", paste0(filename_prefix, "_LabScatterPlot_", sample_type, "_", variable, "_", if_else(is.null(normalisation_sample), "", paste0("_norm_by_",sample_type)),".pdf")), device=cairo_pdf,
                       dpi = 600, width = 160, height =210, units = "mm")
    else
      ggsave(plot = plt2, filename = here("manuscript/output", paste0(filename_prefix, "_LabScatterPlot_", sample_type, "_", variable, "_", if_else(is.null(normalisation_sample), "", paste0("_norm_by_",sample_type)),".png")), device = "png",
                         dpi = 600, width = 160, height =210, units = "mm")
  }
  return(list(plt = plt2, stats = d_stat_filt))
}



# Plotting all samples (with replicate StDev as error bars) of all labs. 
plot_comparison_norm <- function(d_srm_mean, variable, sample_type, normalisation_sample = NULL, 
                            excluded_labs = "", save_plot = TRUE, filename_prefix = "SRM-normalization",
                            selected_ceramides = c("Cer 18:1;O2/16:0", "Cer 18:1;O2/18:0", "Cer 18:1;O2/24:0", "Cer 18:1;O2/24:1"),
                            pos = ggbeeswarm::position_beeswarm(dodge.width = 0.07, method = "center"), show_n = TRUE){
  
  var_meanreplicates <- paste0(variable, "_mean")
  var_meanreplicates_norm <- paste0(var_meanreplicates, "_norm")
  
  d_srm_mean <- d_srm_mean |> filter(!(LabId %in% excluded_labs))
  
  if(variable == "C_Adj"){
  d_stat_filt <- get_study_stats(flag_outlier(d_srm_mean) |> filter(!Outlier_mp))
    y_label <- "\U003BCmol/L (multi-point calibration)"
  } else if(variable == "C_SinglePoint"){
  d_stat_filt <- get_study_stats(flag_outlier(d_srm_mean) |> filter(!Outlier_sp))
    y_label <- "\U003BCmol/L (single-point calibration)"
  } else {stop("Variable not present or supported!")}
  
  d_srm_mean <- d_srm_mean |>
  left_join(d_stat_filt |> filter(SampleType == "SRM") |> select(ceramideName, C_SinglePoint_interlabMEAN, C_Adj_interlabMEAN)) |> 
  group_by(ceramideName, LabId) |>
  mutate(!!sym(var_meanreplicates_norm) := !!sym(var_meanreplicates)/median((!!sym(var_meanreplicates))[SampleType == normalisation_sample]) * (!!sym(paste0(variable, "_interlabMEAN")))) 

  d_plot <- d_srm_mean |> 
    filter(SampleType %in% sample_type, ceramideName %in% selected_ceramides) |>
    select(SampleType, ceramideName, LabId, Protocol, !!sym(var_meanreplicates), !!sym(var_meanreplicates_norm)) |> 
    pivot_longer(cols = !!sym(var_meanreplicates):!!sym(var_meanreplicates_norm), names_to = "norm_type",values_to = "Concentration")

  d_plot_stats <- d_plot |> 
    group_by(SampleType, ceramideName) |> 
    summarise(
      MEAN = mean(Concentration[norm_type==var_meanreplicates]),
      MEAN_norm = mean(Concentration[norm_type==var_meanreplicates_norm]),
      SD = sd(Concentration[norm_type==var_meanreplicates]),
      SD_norm = sd(Concentration[norm_type==var_meanreplicates_norm]),
      SEM = sd(Concentration[norm_type==var_meanreplicates])/length(Concentration[norm_type==var_meanreplicates]),,
      SEM_norm = sd(Concentration[norm_type==var_meanreplicates_norm])/length(Concentration[norm_type==var_meanreplicates_norm]),
      CV = sd(Concentration[norm_type==var_meanreplicates])/mean(Concentration[norm_type==var_meanreplicates])*100,
      CV_norm = sd(Concentration[norm_type==var_meanreplicates_norm])/mean(Concentration[norm_type==var_meanreplicates_norm])*100,
      p_value = t.test(y = Concentration[norm_type==var_meanreplicates], x = Concentration[norm_type==var_meanreplicates_norm], paired = T)$p.value
    )
  

  readr::write_excel_csv(d_plot, file = here(paste0("manuscript/output/", filename_prefix, "",  sample_type, "_sourcedata_points.csv")))
  readr::write_excel_csv(d_plot_stats, file = here(paste0("manuscript/output/", filename_prefix, "", sample_type, "_sourcedata_statistics.csv")))
  
  
  d_sum <- d_plot |> 
    ungroup() |> 
    mutate(y_max = max(Concentration), .by = ceramideName) |> 
    summarise(
      n_labs = paste0("n=",n()), 
      conc_min = min(Concentration) - mean(y_max)/20,
      .by = c(norm_type, ceramideName)
    )
  
  #write_csv(x = d_plot_stats, file = here(glue("manuscript/output/Stats_norm_{sample_type}.csv")))
    
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
    
    scale_x_discrete(labels = c("None", "SRM1950"))
  
  if(show_n){
    plt <- plt + geom_text(
      data = d_sum,
      aes(label = n_labs, x = norm_type , y = conc_min),
      size = 1.3, nudge_x = 0.0,inherit.aes = FALSE, color = "grey70", fontface = "italic"
    )
  }
  
  plt <- plt + 
    #scale_linetype_manual(values = c("SOP" = "solid", "OTHER" = "dotted"))+
    labs(color=NULL, fill = NULL)+
    theme_classic(base_size = 7) +
    labs(y=y_label, x="Recalibration") +
    #scale_x_discrete(drop = FALSE, expand = expansion(mult = c(.03, .05))) +
    theme(
      axis.text.x = element_text(size = 6, face = "plain", angle = 0),
      axis.text.y = element_text(size = 6, face = "plain", angle = 0),
      axis.title.x = element_text(size = 7, face = "bold", angle = 0, margin = margin(t = 8, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size = 6, face = "bold", angle = 90, margin = margin(t = 0, r = 7, b = 0, l = 8)),
      #panel.grid.major.x = element_line(colour = "grey70", linewidth = 0.2, linetype = "dotted"),
      strip.background = element_rect(fill="grey99")) +
     theme(legend.position.inside = c(.96,.08),
          legend.key.size = unit(0.2, "cm"),
          legend.text =  element_text(size = 5, face = "plain", angle = 0),
          legend.spacing.y = unit(.001, "cm"))
  
  plt
}


