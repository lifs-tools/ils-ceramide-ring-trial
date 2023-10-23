#overrides for ggpubr::stat_regline_equation which does not allow to set the figures of merit for R etc.
stat_regline_equation <- function(
  mapping = NULL, data = NULL, formula = y~x,
  label.x.npc = "left", label.y.npc = "top",
  label.x = NULL, label.y = NULL, output.type = "expression",
  geom = "text", position = "identity",  na.rm = FALSE, show.legend = NA,
  inherit.aes = TRUE, ...
)
{
  
  parse <- ifelse(output.type == "expression", TRUE, FALSE)
  
  layer(
    stat = StatReglineEquation, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(formula = formula, label.x.npc  = label.x.npc , label.y.npc  = label.y.npc,
                  label.x = label.x, label.y = label.y,
                  output.type = output.type, parse = parse, na.rm = na.rm, ...)
  )
}

# use options("digits") to get the number of digits to display using signif
StatReglineEquation<- ggproto("StatReglineEquation", Stat,
    required_aes = c("x", "y"),
    default_aes = aes(label = after_stat(eq.label), hjust = ..hjust.., vjust = ..vjust..),
    compute_group = function(data, scales, formula, label.x.npc, label.y.npc,
                             label.x, label.y, output.type)
    {
      
      force(data)
      if (length(unique(data$x)) < 2) {
        return(data.frame()) # Not enough data to perform test
      }
      .test <- .stat_lm(formula, data, output.type = output.type, signif = getOption("digits"))
      # Returns a data frame with label: x, y, hjust, vjust
      .label.pms <- .label_params(data = data, scales = scales,
                                  label.x.npc = label.x.npc, label.y.npc = label.y.npc,
                                  label.x = label.x, label.y = label.y ) |>
        mutate(hjust = 0)
      cbind(.test, .label.pms)
    }
)

# Compute regression line equation
.stat_lm <- function(formula, data, output.type = "expression", signif = 2){
  
  res.lm <- stats::lm(formula, data)
  coefs <- stats::coef(res.lm)
  
  formula.rhs.chr <- as.character(formula)[3]
  if (grepl("-1", formula.rhs.chr) || grepl("- 1", formula.rhs.chr)) {
    coefs <- c(0, coefs)
  }
  
  rr <- summary(res.lm)$r.squared |> signif(signif)
  adj.rr <- summary(res.lm)$adj.r.squared |> signif(signif)
  AIC <- stats::AIC(res.lm) |> signif(signif)
  BIC <- stats::BIC(res.lm) |> signif(signif)
  
  # Build model equation
  eq.char <- as.character(signif(polynom::as.polynomial(coefs), signif))
  eq.char <- gsub("e([+-]?[0-9]*)", "%*%10^\\1", eq.char)
  if (output.type %in% c("latex", "tex", "tikz")) {
    eq.char <- gsub("*", " ", eq.char, fixed = TRUE)
  }
  # Add y
  if (output.type == "expression") {
    lhs <- "italic(y)~`=`~"
  } else if (output.type %in% c("latex", "tex", "tikz", "text")) {
    lhs <- "y = "
  }
  eq.char <- paste(lhs, eq.char, sep = "")
  
  # Build data frame with the output
  if (output.type == "expression") {
    eq.x.rhs = "~italic(x)"
  } else {
    eq.x.rhs = " x"
  }
  
  if (output.type == "expression") {
    z <- data.frame(eq.label = gsub("x", eq.x.rhs, eq.char, fixed = TRUE),
                    rr.label = paste("italic(R)^2", rr, sep = "~`=`~"),
                    adj.rr.label = paste("italic(R)[adj]^2",
                                         adj.rr, sep = "~`=`~"),
                    AIC.label = paste("AIC", AIC, sep = "~`=`~"),
                    BIC.label = paste("BIC", BIC, sep = "~`=`~"))
  } else if (output.type %in% c("latex", "tex", "text")) {
    z <- data.frame(eq.label = gsub("x", eq.x.rhs, eq.char, fixed = TRUE),
                    rr.label = paste("R^2", rr, sep = " = "),
                    adj.rr.label = paste("R_{adj}^2",adj.rr, sep = " = "),
                    AIC.label = paste("AIC", AIC, sep = " = "),
                    BIC.label = paste("BIC", BIC, sep = " = "))
  }
  
  z <- z |>
    mutate(rr = rr, adj.rr = adj.rr, AIC = AIC, BIC = BIC) |>
    dplyr::select(rr, adj.rr, AIC, BIC, everything())
  
  z
}

plotDatasetSummary <- function(
  datasetSummary, 
  standardReports, 
  preferredReports,
  theme = mytheme,
  outputDirectory
) {
  allReports <- data.frame("LabId"=c(standardReports, preferredReports))
  filteredSummary <- datasetSummary %>% dplyr::semi_join(allReports, by = c("LabId"))
  readr::write_csv(filteredSummary, file=file.path(outputDirectory,"filteredDataSummary.csv"))
  # filteredSummary <- readr::read_csv(file="figures/filteredDataSummary.csv")
  filteredSummaryLong <- filteredSummary %>% 
    dplyr::select(-`Processing Status`) %>%
    dplyr::group_by(LabId) %>%
    tidyr::pivot_longer(
      cols = 2:last_col(),
      names_to = "Category",
      values_to = "Reported"
    ) |> 
    dplyr::mutate(
      Reported=forcats::fct_relevel(Reported, c("Standard", "Preferred", "No", "Yes (with remarks)", "Yes")),
      Category=forcats::fct_relevel(Category, c("Blank", "Calibration Line 1", "Calibration Line 2", "LLQC","LQC","MQC","HQC","HLQC", "NIST SRM", "NIST hTAG", "NIST T1D", "NIST Young AA"))
    )
  heatmap <- ggplot2::ggplot(
    filteredSummaryLong, ggplot2::aes(LabId, Category, fill = Reported)
    ) + 
    ggplot2::geom_tile(color = "black") + 
    ggplot2::scale_fill_viridis_d() + 
    theme + 
    ggplot2::coord_fixed()
  ggplot2::ggsave(filename = file.path(outputDirectory, "dataSummaryHeatmap.png"), plot = heatmap, width = 8, height = 4)
  heatmap
}

meanSdCalibrationLinesPlots <- function(calibLineDataLmPred, theme=mytheme, avgLipidToISRatios, suffix="", cunit="", d_theoratios, outputDirectory) {
  library(ggrepel)
  # Plot "average" calibration line
  calibLineDataLmPredPlot <- calibLineDataLmPred |> 
    group_by(
      ceramideId,
      ceramideName,
      Concentration
    ) |> 
    mutate(.fitted = .fitted) |>
    summarise(
      ymin=min(.fitted),
      ymax=max(.fitted),
      ymed=median(.fitted),
      yavg=mean(.fitted),
      ysd=sd(.fitted),
      yavgmin=max(ymin, yavg-ysd),
      yavgmax=yavg+ysd,
      ymedmin=max(ymin, ymed-ysd),
      ymedmax=ymed+ysd,
    ) |> mutate(
      labelavg=paste(format(yavg, digits=3, nsmall=3),"+/-", format(ysd, digits=2, nsmall=2)),
      labelmed=paste(format(ymed, digits=3, nsmall=3),"+/-", format(ysd, digits=2, nsmall=2)),
      labelconc=paste(format(Concentration, digits=3, nsmall=3)),
      labelsd=paste("\u03C3 +/-", format(ysd, digits=2, nsmall=2))
    )

  write_csv(calibLineDataLmPredPlot, here::here(file.path(outputDirectory,"temp.csv")))
  
  calibLineDataLmPredPlot <- calibLineDataLmPredPlot %>% 
    left_join(d_theoratios, by=c("ceramideId", "ceramideName", "Concentration"))
  
  if(!is.null(avgLipidToISRatios)) {
    calibLineDataLmPredPlot <- calibLineDataLmPredPlot |> left_join(avgLipidToISRatios, by=c("ceramideId", "ceramideName"))
  }
  
  avgMedCalibrationLinesLog10 <- ggplot(data=calibLineDataLmPredPlot, mapping=aes(x=Concentration, y=yavg)) +
    geom_ribbon(aes(x=Concentration, ymin=yavgmin, ymax=yavgmax), fill="blue", alpha=0.1, show.legend = FALSE) + 
    geom_line(aes(x=Concentration, y=TheoreticalRatio), color="green", linewidth = 1.3, show.legend = FALSE) + 
    geom_line(aes(x=Concentration, y=yavg), color="blue", show.legend = FALSE) + 
    geom_line(aes(x=Concentration, y=ymed), color="black", linetype=2, show.legend = FALSE, alpha=0.2) +
    geom_point(aes(x=Concentration, y=yavg), color="blue", alpha=0.4) +
    geom_point(aes(x=Concentration, y=ymed), color="black", alpha=0.4) +
    geom_errorbar(aes(x=Concentration, ymin=yavgmin, ymax=yavgmax), color="black", alpha=0.5, linetype=1) +
    geom_vline(aes(xintercept=Concentration), linetype=1, color="black", alpha=0.3) +
    geom_text(aes(x=Concentration, y=max(yavg+(ymax-ymin)*0.05), angle=90, label=labelconc), size=3, vjust = -0.5, hjust=0.5, check_overlap = TRUE, color="black") +
    geom_text_repel(aes(x=Concentration, y=yavgmin, label=labelsd), min.segment.length = 1, size=3, nudge_y = -0.1, point.padding=4, direction="y", color="black")
  if(!is.null(avgLipidToISRatios)) {
    avgMedCalibrationLinesLog10 <- avgMedCalibrationLinesLog10 +
      geom_ribbon(aes(x=Concentration, y=meanAreaRatios,ymin=meanAreaRatios-sdAreaRatios,ymax=meanAreaRatios+sdAreaRatios), fill="red", alpha=0.1, show.legend = FALSE) +
      geom_hline(aes(yintercept=meanAreaRatios), color="red", show.legend = FALSE, linetype=3)
  }
  avgMedCalibrationLinesLog10 <- avgMedCalibrationLinesLog10 +
    scale_color_discrete(breaks=c("blue","black"), labels=c("Mean","Median")) +
      facet_wrap(.~ceramideName, scales = "free_x") + 
    xlab(paste0("Concentration ", cunit)) + ylab("Ratio Light/Heavy") +
    scale_x_log10() + scale_y_log10() + 
    theme
  ggsave(file.path(outputDirectory, paste0("CalibrationLinesAvgMedRangeLog10",suffix,".png")), avgMedCalibrationLinesLog10, width = 9, height = 9)
  
  avgMedCalibrationLinesLinear <- ggplot(data=calibLineDataLmPredPlot, mapping=aes(x=Concentration, y=yavg)) +
    geom_ribbon(aes(x=Concentration, ymin=yavgmin, ymax=yavgmax), fill="blue", alpha=0.1, show.legend = FALSE) + 
    geom_line(aes(x=Concentration, y=TheoreticalRatio), color="green", linewidth = 1.3, show.legend = FALSE) + 
    geom_line(aes(x=Concentration, y=yavg), color="blue", show.legend = FALSE) + 
    geom_line(aes(x=Concentration, y=ymed), color="black", linetype=2, show.legend = FALSE, alpha=0.2) +
    geom_point(aes(x=Concentration, y=yavg), color="blue", alpha=0.4) +
    geom_point(aes(x=Concentration, y=ymed), color="black", alpha=0.4) +
    geom_errorbar(aes(x=Concentration, ymin=yavgmin, ymax=yavgmax), color="black", alpha=0.5, linetype=1) +
    geom_vline(aes(xintercept=Concentration), linetype=1, color="black", alpha=0.3) +
    geom_text(aes(x=Concentration, y=max(yavg+(ymax-ymin)*0.05), angle=90, label=labelconc), size=3, vjust = -0.5, hjust=0.5, check_overlap = TRUE, color="black") +
    geom_text_repel(aes(x=Concentration, y=yavgmin, label=labelsd), min.segment.length = 1, size=3, nudge_y = -0.1, point.padding=4, direction="y", color="black")
  if(!is.null(avgLipidToISRatios)) {
    avgMedCalibrationLinesLinear <- avgMedCalibrationLinesLinear + 
      geom_ribbon(aes(x=Concentration, y=meanAreaRatios,ymin=meanAreaRatios-sdAreaRatios,ymax=meanAreaRatios+sdAreaRatios), fill="red", alpha=0.1, show.legend = FALSE) +
      geom_hline(aes(yintercept=meanAreaRatios), color="red", show.legend = FALSE, linetype=3)
  }
  avgMedCalibrationLinesLinear <- avgMedCalibrationLinesLinear + 
    scale_color_discrete(breaks=c("blue","black"), labels=c("Mean","Median")) +
    facet_wrap(.~ceramideName, scales = "free_x") +
    xlab(paste0("Concentration ", cunit)) + ylab("Ratio Light/Heavy") +
    theme
  ggsave(file.path(outputDirectory, paste0("CalibrationLinesAvgMedRangeLinear",suffix,".png")), avgMedCalibrationLinesLinear)
  return(list(avgMedCalibrationLinesLog10, avgMedCalibrationLinesLinear))
}

calibrationLineSurveyPlot <- function(intraAssayTable, theme=mytheme, colorscale=mycolorscale, outlierShape=19, outputDirectory) {
  calibrationLinesPlot <-
    ggplot(
      data = intraAssayTable,
      mapping = aes(y = area, x = ceramideName, color = isotope)
    ) + 
    geom_boxplot(outlier.shape=outlierShape) + 
    geom_point(position = position_jitterdodge()) + 
    facet_grid(SampleType ~ LabId) + 
    theme + colorscale
  ggsave(file.path(outputDirectory, "CalibrationLinesSurvey.png"), calibrationLinesPlot, width=16)
  calibrationLinesPlot
}

qcSurveyPlot <- function(intraAssayQC, theme=mytheme, colorscale=mycolorscale, outlierShape=19, outputDirectory) {
  QCPlot <-
    ggplot(
      data = intraAssayQC,
      mapping = aes(y = area, x = ceramideName, color = isotope)
    ) + 
    geom_boxplot(outlier.shape=outlierShape) + 
    geom_point(position = position_jitterdodge()) + 
    facet_grid(SampleType ~ LabId) + 
    theme + colorscale
  ggsave(file.path(outputDirectory, "QCSamplesSurvey.png"), QCPlot, width=16, height=10)
  QCPlot
}

nistSurveyPlot <- function(IntraAssayNIST, theme=mytheme, colorscale=mycolorscale, outlierShape=19, outputDirectory) {
  NISTSurveyPlot <-
    ggplot(
      data = IntraAssayNIST,
      mapping = aes(y = area, x = ceramideName, color = isotope)
    ) + 
    geom_boxplot(outlier.shape=outlierShape) + 
    geom_point(position = position_jitterdodge()) + 
    facet_grid(SampleType ~ LabId) + 
    theme + colorscale
  ggsave(file.path(outputDirectory, "NISTSamplesSurvey.png"), NISTSurveyPlot, width=16)
  NISTSurveyPlot
}

scatterRatioPlots <- function(linRegFormula, data, ceramideIds, outputDirectory) {
  ceramideIds %>% 
    map(
      ~ scatterRatioPlot(
        linRegFormula=linRegFormula,
        data=data,
        selectedCeramideId=.x,
        outputDirectory=outputDirectory
      )
    )
}

scatterRatioPlot <- function(linRegFormula, data, selectedCeramideId, outputDirectory) {
  dataCer <- data %>% dplyr::filter(ceramideId == selectedCeramideId)
  stopifnot(length(unique(dataCer$ceramideId))==1)
  scatterRatioPlotCer <- ggplot(
    data = dataCer, 
    aes(x = RatioLipidToIS, y = Concentration, color=LabId),
    alpha = 0.1
  ) + 
    geom_point() + 
    facet_grid(SampleType~ceramideName, scales = "free") + 
    stat_smooth(method="lm", formula = linRegFormula, mapping = aes(weight = 1/(Concentration^2))) +
    ylab(paste0("Actual Concentration ", unique(data$Unit))) +
    coord_flip()
  ggsave(file.path(outputDirectory, paste0("CalibrationLinesScatterRatioCer", selectedCeramideId,".png")), scatterRatioPlotCer)
  scatterRatioPlotCer
}

nistAreaRatioPlots <- function(data, sampleTypes, theme=mytheme, fillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
  sampleTypes %>% 
    map(
      ~ nistAreaRatioPlot(
        data=data,
        selectedSampleType=.x,
        theme=theme,
        fillscale=fillscale,
        outlierShape=outlierShape,
        outputDirectory=outputDirectory
      )
    )
}

nistAreaRatioPlot <- function(data, selectedSampleType, theme=mytheme, fillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
  dataSampleType <- data %>% dplyr::filter(SampleType == selectedSampleType)
  stopifnot(length(unique(dataSampleType$SampleType))==1)
  plot <- ggplot(data=dataSampleType, 
         mapping = aes(x=LabId, y=RatioLipidToIS, fill=ceramideName)) +
    geom_boxplot(outlier.shape=outlierShape) + ylab("Cer / ISCer (area ratio)") +
    facet_grid(ceramideName~SampleType, scales="free_y") + theme + fillscale
  ggsave(file.path(outputDirectory, paste0("NIST_Light_Heavy_AreasRatioPlot-", selectedSampleType,".png")), plot)
  plot
}

calibrationLinePlot <- function(linRegFormula, lmPred, adjAnalyteConcentrations, theme=mytheme, labId, digits = 4, cunit = "", d_theoratios, outputDirectory) {
  oldDigits <- getOption("digits")
  options(digits = digits)
  filteredLmPred <- lmPred %>% filter(LabId==labId)
  filteredAdjAnalyteConcentrations <- adjAnalyteConcentrations %>% filter(LabId==labId)
  
  filteredAdjAnalyteConcentrations <- filteredAdjAnalyteConcentrations %>% 
    left_join(d_theoratios, by=c("ceramideId", "ceramideName", "Concentration"))
  
  Cl4Plot <- ggplot(
    data=filteredAdjAnalyteConcentrations, 
    mapping = aes(x=Concentration , y=RatioLipidToIS)
  ) + 
    geom_point(aes(color=LabId, shape=ceramideName)) +
    # 1.96 is for 5% and 95% CIs
    geom_ribbon(data = filteredLmPred, aes(x= Concentration, y = .fitted, ymin=.fitted-1.96*.se.fit, ymax=.fitted+1.96*.se.fit), alpha=0.2) +
    geom_line(aes(x=Concentration, y=TheoreticalRatio), color="green", linewidth = 1.3, show.legend = FALSE) + 
    geom_line(data = filteredLmPred, aes(x = Concentration ,y = .fitted, color = LabId), linewidth = 1) +
    stat_regline_equation(formula = linRegFormula, label.y = max(filteredLmPred$.fitted)+(max(filteredLmPred$.fitted)-min(filteredLmPred$.fitted))*0.4, aes(label =  paste(after_stat(eq.label), sep = "~~"))) + 
    stat_regline_equation(formula = linRegFormula, label.y = max(filteredLmPred$.fitted)+(max(filteredLmPred$.fitted)-min(filteredLmPred$.fitted))*0.2, aes(label =  paste(after_stat(adj.rr.label), sep = "~~"))) +
    facet_wrap(LabId+SampleType~ceramideName, scales = "free", ncol=4) + 
    theme + 
    labs(subtitle = paste("Results for Lab", labId)) +
    xlab(paste0("Actual Concentration ", cunit)) + ylab("Ratio Light/Heavy")
  ggsave(file.path(outputDirectory, paste0("CalibrationLineAnalyteConcentrationsLab",labId, ".png")), plot=Cl4Plot, width=10)
  options(digits = oldDigits)
  filteredAdjAnalyteConcentrations
}

calibrationLineVsSinglePointPlot <- function(analyteConcentrationsFromCalibLines, theme=mytheme, mycolorscale, outputDirectory) {
  plot <- ggplot(
    aes(x=Concentration, y=C_Adj/C_SinglePoint, color=ceramideName, label=LabId, group=LabId), 
    data=analyteConcentrationsFromCalibLines) + 
    geom_smooth(method = "lm") +
    facet_wrap(~ceramideName, scale="free") + 
    theme +
    mycolorscale
  ggsave(file.path(outputDirectory, paste0("figures/CalibrationLineVsSinglePoint.png")), plot=plot, width=8)
  plot
}

calibrationAnalyteConcentrationsPlot <- function(analyteConcentrationsFromCalibLines, calibLineDataLmPred, theme=mytheme, cunit, linRegFormula, outputDirectory) {
ClPlot <- ggplot(
  data=analyteConcentrationsFromCalibLines, 
  mapping = aes(x = Concentration, y=RatioLipidToIS)
) + 
  geom_point(aes(color=LabId, shape=ceramideName)) +
  # 1.96 is for 5% and 95% CIs
  geom_ribbon(data = calibLineDataLmPred, aes(x= Concentration, y = .fitted, ymin=.fitted-1.96*.se.fit, ymax=.fitted+1.96*.se.fit), alpha=0.2) +
  geom_line(data = calibLineDataLmPred, aes(x = Concentration ,y = .fitted, color = LabId), linewidth = 1) +
  stat_regline_equation(formula = linRegFormula, label.y = max(calibLineDataLmPred$.fitted)+(max(calibLineDataLmPred$.fitted)-min(calibLineDataLmPred$.fitted))*0.2, aes(label =  paste(after_stat(eq.label), sep = "~~"))) + 
  stat_regline_equation(formula = linRegFormula, label.y = max(calibLineDataLmPred$.fitted)+(max(calibLineDataLmPred$.fitted)-min(calibLineDataLmPred$.fitted))*0.1, aes(label =  paste(after_stat(adj.rr.label), sep = "~~"))) +
  facet_wrap(LabId+SampleType~ceramideName, scales = "free", ncol=8) + 
  labs(subtitle = paste("Results for all labs")) +
  xlab(paste0("Concentration ", cunit)) + ylab("Ratio Light/Heavy") +
  theme
  ggsave(file.path(outputDirectory, "CalibrationLineAnalyteConcentrations.png"), plot=ClPlot, width=12, height=60, limitsize = FALSE)
  return(ClPlot)
}

################################################################################
# QC concentrations by sample type Plot
################################################################################
qcConcentrationsPlot <- function(qcAveragedConcentrationsLong, theme=mytheme, fillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
  qcConcentrationsPlot <- ggplot(aes(x=SampleType,y=Adj_Conc,fill=Calibration), data=qcAveragedConcentrationsLong) + 
  geom_boxplot(outlier.shape = outlierShape) + 
  scale_x_discrete(name="Sample Type", labels=c("LLQC","LQC","MQC","HQC","HLQC")) +
  facet_grid(ceramideName~facetLabel, scales="free_y") + 
  ylab(paste0("Adj. Concentration ",unique(qcAveragedConcentrationsLong$Unit))) +
  theme + fillscale
  ggsave(file.path(outputDirectory, "QCSamplesConcentrations.png"), qcConcentrationsPlot)
  qcConcentrationsPlot
}

################################################################################
# QC concentrations by lab id Plot
################################################################################
qcConcentrationsByLabIdPlot <- function(qcAveragedConcentrations, theme=mytheme, fillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
  qcConcentrationsByLabIdPlot <- ggplot(aes(x=LabId,y=Avg_C_Adj,fill=ceramideName), data=qcAveragedConcentrations) + 
    geom_boxplot() + 
    facet_grid(ceramideName~SampleType, scales="free_y") + 
    ylab(paste0("Adj. Concentration ",unique(qcAveragedConcentrations$Unit))) +
    theme + fillscale
  ggsave(file.path(outputDirectory, "QCSamplesConcentrationsByLabId.png"), qcConcentrationsByLabIdPlot, width=35)
  qcConcentrationsByLabIdPlot
}

################################################################################
# QC Sample Stats Plot
################################################################################
qcSampleStatsPlot <- function(qcSampleStats, theme=mytheme, fillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
  #"LabId,","SampleType","ceramideId","ceramideName","Unit","Protocol","Instrument","MassAnalyzerType"
  qcSampleStatsLong <- qcSampleStats %>% 
    ungroup() %>% 
    select(
      -SdAvg_C_Adj, 
      -CV_C_Adj,
      -SdAvg_C_SinglePoint,
      -CV_C_SinglePoint) %>% pivot_longer(
        cols = c("AvgAvg_C_Adj", "Avg_C_SinglePoint"),
        names_to = c("CalibrationType"),
        values_to = c("MeanAdjustedConcentration")
      ) %>%
    mutate(
      CalibrationType = forcats::fct_recode(
        CalibrationType,
        "Calibration Curve" = "AvgAvg_C_Adj",
        "Single Point" = "Avg_C_SinglePoint"
      )#,
      # LabId = fct_reorder(LabId, MassAnalyzerType)
    )
  qcSampleStatsPlot <- ggplot(aes(x=MeanAdjustedConcentration, y=LabId, shape=MassAnalyzerType, color=MassAnalyzerType), data=qcSampleStatsLong) +
    geom_point() + 
    xlab(paste0("Mean Adj. Conc. ",unique(qcSampleStatsLong$Unit))) +
    facet_grid(SampleType+CalibrationType~ceramideName, scales="free_x") + 
    theme + fillscale
  ggsave(file.path(outputDirectory, "QCSamplesStatsByLabId.png"), qcSampleStatsPlot, height=35)
  qcSampleStatsPlot  
}

################################################################################
# QC concentrations density by lab id Plot
################################################################################
qcConcentrationsByLabIdDensityPlot <- function(qcAveragedConcentrations, theme=mytheme, fillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
  qcConcentrationsByLabIdDensityPlot <- ggplot(aes(x=Avg_C_Adj,fill=ceramideName), data=qcAveragedConcentrations) + 
    geom_density(trim=TRUE, alpha=0.9) + 
    geom_rug(alpha=0.9) +
    geom_density(aes(x=C_SinglePoint), trim=TRUE, alpha=0.25) + 
    geom_rug(aes(x=C_SinglePoint), alpha=0.25) +
    facet_grid(ceramideName~SampleType, scales="free_y") + 
    xlab(paste0("Adj. Concentration ",unique(qcAveragedConcentrations$Unit))) +
    scale_x_log10() +
    theme + fillscale
  ggsave(file.path(outputDirectory, "QCSamplesConcentrationsByLabIdDensity.png"), qcConcentrationsByLabIdDensityPlot, width=10, height=6)
  qcConcentrationsByLabIdDensityPlot
}

################################################################################
# Stds Box Plot
################################################################################
stdsBoxplot <- function(averagedConcentrations, selectCeramideId, theme=mythemeXRot, fillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
  dta <- averagedConcentrations %>% filter(ceramideId==selectCeramideId)
  loDf <- data.frame(
    LoD=dta$Avg_C_LoD_C, 
    LoQ=dta$Avg_C_LoQ_C, 
    LabId=dta$LabId, 
    Sample=dta$Sample,
    label=paste(dta$Sample, " ", dta$Concentration, " ", dta$Unit)
  ) %>% filter(
    Sample=="STD 6"
  ) %>% distinct()
  loDfLong <- loDf %>% pivot_longer(cols = c("LoD", "LoQ"), names_to = "FoM", values_to = "y")
  # boxplot of averaged adjusted concentrations
  boxPlot <- ggplot(aes(x = LabId, y = Avg_C_Adj, fill = Sample), data=dta) + 
    geom_boxplot(show.legend = FALSE, outlier.shape = outlierShape) + 
    geom_point(aes(x = LabId, y = Avg_C_Adj), shape = 20, position=position_jitterdodge(), show.legend = FALSE) + 
    facet_wrap(c("Sample","ceramideName"), scales = "free") +
    geom_point(aes(x = LabId, y = y, shape = FoM), data=loDfLong, show.legend = TRUE) +
    geom_hline(aes(yintercept = Concentration), color="red", linetype = 2, show.legend = FALSE) +
    scale_shape_discrete(solid = FALSE) + # use non-solid characters for FoMs LoQ and LoD
    guides(fill = "none") + # do not show points in legend for Sample
    ylab(paste0("Adj. Concentration ",unique(dta$Unit))) +
    labs(caption = paste0("Avg. LoD from calib. lines for STD 6, where LoD=3.3*sd(a_l/a_h)/slope\n",
           "Avg. LoQ from calib. lines for STD 6, where LoQ=10*sd(a_l/a_h)/slope")) +
    theme + 
    fillscale
  ggsave(file.path(outputDirectory, paste0("CalibrationLineSTDsCer",selectCeramideId,".png")), plot=boxPlot, width=20)
  boxPlot
}

################################################################################
# Stds Variance Histogram Plot
################################################################################
stdsVarHistoPlot <- function(averagedConcentrations, selectCeramideId, theme=mythemeXRot, fillscale=scale_fill_nejm(), outputDirectory) {
  dta <- averagedConcentrations %>% 
    filter(ceramideId==selectCeramideId)
  # boxplot of averaged adjusted concentrations
  boxPlot <- ggplot(aes(x = Avg_C_Adj, fill = Sample), data=dta) + 
    geom_density(alpha=0.9) +
    geom_density(aes(x = Avg_C_SinglePoint, fill = Sample), data=dta, alpha=0.25) +
    geom_rug(alpha=0.9) +
    geom_rug(aes(x = Avg_C_SinglePoint), data=dta, alpha=0.25) +
    geom_vline(aes(xintercept = Concentration), color="red", linetype = 2, show.legend = FALSE) +
    facet_wrap(c("label","ceramideName"), scales = "free") +
    guides(fill = "none") + # do not show points in legend for Sample
    xlab(paste0("Adj. Concentration ",unique(dta$Unit))) +
    theme + 
    fillscale
  ggsave(file.path(outputDirectory, paste0("CalibrationLineSTDsVarHistoCer",selectCeramideId,".png")), plot=boxPlot, width=20)
  boxPlot
}

################################################################################
# Stds Recovery Percent Histogram Plot
################################################################################
stdsRecoveryPercentHistoPlot <- function(averagedConcentrations, selectCeramideId, theme=mythemeXRot, fillscale=scale_fill_nejm(), outputDirectory) {
  dta <- averagedConcentrations %>% 
    filter(ceramideId==selectCeramideId) %>% 
    pivot_longer(cols=c(Avg_C_Adj_Rec_Perc, Avg_C_SinglePoint_Rec_Perc), names_to="Calibration_Type", values_to="Adj_C_Recovery_Perc")
  readr::write_csv(file=file.path(outputDirectory, "stdsRecoveryPercentHistoData.csv"), dta)
  # boxplot of averaged adjusted concentrations
  boxPlot <- ggplot(aes(y = Adj_C_Recovery_Perc, x = LabId, fill = Sample), data=dta) + 
    geom_boxplot(aes(y = Adj_C_Recovery_Perc, fill = Calibration_Type), data=dta, alpha=0.9) +
    geom_hline(aes(yintercept = 100), color="red", linetype = 2, show.legend = FALSE) +
    geom_hline(aes(yintercept = 80), color="red", linetype = 2, show.legend = FALSE, alpha=0.5) +
    geom_hline(aes(yintercept = 120), color="red", linetype = 2, show.legend = FALSE, alpha=0.5) +
    facet_wrap(c("label","ceramideName"), scales = "free") +
    guides(fill = "none") + # do not show points in legend for Sample
    ylab(paste0("Adj. Concentration Recovery % (log10)")) +
    scale_y_log10() +
    theme + 
    fillscale
  ggsave(file.path(outputDirectory, paste0("CalibrationLineSTDsRecPercHistoCer",selectCeramideId,".png")), plot=boxPlot, width=20)
  boxPlot
}

################################################################################
# NIST Concentrations Plot
################################################################################
nistConcentrationsPlot <- function(nistAveragedConcentrationsLong, theme=mythemeXRot45, fillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
  nistConcentrationsPlot <- ggplot(aes(x=SampleType, y=Adj_Conc, fill=Calibration), data=nistAveragedConcentrationsLong) + 
    geom_boxplot(outlier.shape=outlierShape) + 
    facet_grid(ceramideName~facetLabel, scales="free_y") + 
    ylab(paste0("Adj. Concentration ",unique(nistAveragedConcentrationsLong$Unit))) +
    theme + fillscale
  ggsave(file.path(outputDirectory, "NISTSamplesConcentrations.png"), nistConcentrationsPlot, width=6)
  nistConcentrationsPlot
}

################################################################################
# NIST CV Plot
################################################################################
nistCvPlot <- function(nistConcentrationsStatsLong, theme=mythemeXRot45, fillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
  # TODO: filter should be removed once the calibration line for lab 28 does not produce a negative concentration
  # for the NIST YAA sample
  nistCvPlot <- ggplot(aes(x=SampleType, y=CV, fill=Calibration), data=nistConcentrationsStatsLong %>%
                         filter(AvgAvg_C_Adj>0)) + 
    geom_boxplot(outlier.shape=outlierShape) + 
    facet_grid(ceramideName~facetLabel, scales="free_y") + 
    ylab("CV") + 
    theme + fillscale
  ggsave(file.path(outputDirectory, "NISTSamplesConcentrationsCV.png"), nistCvPlot, width=6)
  nistCvPlot
}

################################################################################
# NIST concentrations by sample type Plot
################################################################################
nistConcentrationsPlotBySampleType <- function(data, namesuffix, width=6, unit, theme=mythemeXRot, xfillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
    plot <- ggplot(aes(x=LabId,y=Adj_Conc,fill=Calibration), data=data) + 
      geom_boxplot(outlier.shape=outlierShape) + 
      facet_grid(ceramideName~., scales="free_y") +
      ylab(paste0("Adj. Concentration ", unit)) + 
      labs(title=unique(data$SampleType)) + 
      theme + xfillscale
    ggsave(file.path(outputDirectory, paste0("figures/NISTSamplesConcentrationsBySampleType-", namesuffix, ".png")), plot, width=width)
    plot
}

################################################################################
# Concentrations by Sample type, LabId and Mass Analyzer Plot
################################################################################
nistConcentrationsPlotByMa <- function(data, namesuffix, width=6, unit, theme=mythemeXRot, xfillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
  plot <- ggplot(aes(x=LabId,y=Adj_Conc,fill=Calibration), data=data) + 
    geom_boxplot(outlier.shape=outlierShape) + 
    facet_grid(ceramideName~MassAnalyzerType, scales="free") +
    ylab(paste0("Adj. Concentration ", unit)) + 
    labs(title=unique(data$SampleType)) + 
    theme + xfillscale
  ggsave(file.path(outputDirectory, paste0("NISTSamplesConcentrationsBySampleAndMassAnalyzerType-", namesuffix, ".png")), plot, width=2*width)
  plot
}
################################################################################
# Concentrations by LabId Plot
################################################################################
nistConcentrationsByLabIdPlot <- function(nistAveragedConcentrationsLong, theme=mythemeXRot, xfillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
  nistConcentrationsByLabIdPlot <- ggplot(aes(x=LabId,y=Adj_Conc,fill=Calibration), data=nistAveragedConcentrationsLong) + 
  geom_boxplot(outlier.shape=outlierShape) + 
  facet_grid(ceramideName~SampleType, scales="free_y") + 
  ylab(paste0("Adj. Concentration ",unique(nistAveragedConcentrationsLong$Unit))) + 
  mythemeXRot + scale_fill_nejm()
  ggsave(file.path(outputDirectory, "NISTSamplesConcentrationsByLabId.png"), nistConcentrationsByLabIdPlot, width=35)
  nistConcentrationsByLabIdPlot
}

################################################################################
# CV by LabId Plot
################################################################################
nistCvByLabIdPlot <- function(nistConcentrationsStatsLong, theme=mythemeXRot, xfillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
  nistCvByLabIdPlot <- ggplot(aes(x=LabId,y=CV*100,fill=Calibration), data=nistConcentrationsStatsLong) + 
  geom_boxplot(outlier.shape = outlierShape) + 
  scale_y_log10() +
  facet_grid(ceramideName~SampleType, scales="free_y") + 
  ylab(paste0("CV")) + 
  mythemeXRot + scale_fill_nejm()
  ggsave(file.path(outputDirectory, "NISTSamplesCvByLabId.png"), nistCvByLabIdPlot, width=35)
  nistCvByLabIdPlot
}

################################################################################
# Zscores by LabId Plot
################################################################################
nistZscoreByLabIdPlot <- function(nistConcentrationsZScoreStatsLong, theme=mythemeXRot, xfillscale=scale_fill_nejm(), outlierShape=19, outputDirectory) {
  nistZscoreByLabIdPlot <- ggplot(aes(x=LabId,y=ZScore,fill=Calibration), data=nistConcentrationsZScoreStatsLong) + 
  geom_boxplot(outlier.shape = outlierShape) + 
  facet_grid(ceramideName~SampleType, scales="free_y") + 
  ylab(paste0("ZScore")) + 
  mythemeXRot + scale_fill_nejm()
  ggsave(file.path(outputDirectory, "NISTSamplesZScoreByLabId.png"), nistZscoreByLabIdPlot, width=35)
  nistZscoreByLabIdPlot
}

################################################################################
# Ring Trial comparison plot
################################################################################
ringTrialsComparisonPlot <- function(
  ringTrialsComparisonTbl,
  theme=mythemeXRot,
  outlierShape=outlierShape,
  outputDirectory
) {
  # formula standard uncertainty (Ghorasaini et al, Analytical Chemistry, 2021)
  std_uncert <- function(x) sqrt(pi/(2*length(x))) * mad(x, constant = 1.4826, na.rm = TRUE)
  
  ringTrialsComparisonTblFiltered <- ringTrialsComparisonTbl |>
    filter(Study!="Giera 2021 (SP)")
  
  d_lipidyzer <- ringTrialsComparisonTbl |>
    filter(Study=="Giera 2021 (SP)") |> # added because we use the row bound values
    select(Study, Compound, everything()) |>
    mutate(MADM = MAD/MEDM *100)
  
  d_sum_1 <- ringTrialsComparisonTblFiltered |> 
    group_by(Compound, Study) |> 
    summarise(
      n_labs = n(),
      MEAN = mean(Conc_mean, na.rm = TRUE),
      SD = sd(Conc_mean, na.rm = TRUE),
      CV = SD/MEAN *100,
      SEM = SD/sqrt(n_labs),
      CV_SEM = SEM/MEAN *100,
      MAD = mad(Conc_mean, constant = 1.4826, na.rm = TRUE),
      MEDM = median(Conc_mean, na.rm = TRUE),
      MADM = MAD/MEDM * 100, 
      Uncertainty = std_uncert(Conc_mean), 
      COD = Uncertainty/MEDM *100, 
      Conc_max = max(Conc_mean, na.rm = TRUE))
  
  
  write_csv(x = d_sum_1, file = file.path(outputDirectory, "ILS-Ceramides_Conc_Summary_20220808.csv"))
  
  d_sum <- d_sum_1 |>
    bind_rows(d_lipidyzer) |>
    mutate(Study = factor(Study, levels = c("Quehenberger 2010", "Bowden 2017 (MP)",
                                            "Moseley 2019 (SP)", "Giera 2021 (SP)",
                                            "ILS Ceramides (MP)")))|>
    arrange(Study) |>
    # Conc_max used for label positioning only
    mutate(Conc_max = if_else(is.na(Conc_max), MEDM + MAD, Conc_max)) |>
    group_by(Compound) |>
    mutate(Conc_max = Conc_max + max(Conc_max)*0.05)
  p_Fig1 <- ggplot(ringTrialsComparisonTblFiltered, aes(x = Study, y = Conc_mean)) +
    
    ggbeeswarm::geom_beeswarm(
      shape = 21, size = 2, stroke  = 0.5, dodge.width = 1, cex = 1.5, color = "grey50") +
    geom_crossbar(
      data = d_sum,
      aes(y = MEDM, ymin = MEDM - Uncertainty, ymax = MEDM + Uncertainty),
      inherit.aes = TRUE, 
      fatten = .1, size = 0.2, fill = "grey95", color = "red", width = .6, alpha = 0.1) +
    geom_crossbar(
      data = d_sum,
      aes(y = MEDM, ymin = MEDM, ymax = MEDM),
      inherit.aes = TRUE, fatten = 0, size = 1.1, fill = "grey95", color = "red", width = 0.7) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(.03, .1))) +
    geom_text(
      data = d_sum,
      aes(label = n_labs, y = Conc_max),
      size = 3, nudge_x = 0.0, inherit.aes = TRUE, color = "#d40202", fontface = "italic") +
    facet_wrap(vars(Compound), scales = "free", nrow = 1) +
    scale_x_discrete(drop = FALSE, expand = expansion(mult = c(.03, .25))) +
    theme_classic() +
    labs(x = "", y = "\U003BCmol/L") +
    theme(
      axis.text.x = element_text(size = 7, angle = 46, vjust = 1, hjust = 1)
    )
  
  ggsave(
    plot = p_Fig1,
    filename = file.path(outputDirectory, "ils-to-other-ring-trials-comparison-nist-srm1950.png"),
    dpi = 300, width = 220, height = 110, units = "mm"
  )
  p_Fig1
}
