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
  filteredSummary <- datasetSummary |> dplyr::semi_join(allReports, by = c("LabId"))
  readr::write_csv(filteredSummary, file=file.path(outputDirectory,"filteredDataSummary.csv"))
  filteredSummaryLong <- filteredSummary |> 
    dplyr::select(-`Processing Status`) |>
    dplyr::group_by(LabId) |>
    tidyr::pivot_longer(
      cols = 2:last_col(),
      names_to = "Category",
      values_to = "Reported"
    ) |> 
    dplyr::mutate(
      Reported=as.factor(Reported),
      Category=as.factor(Category),
    ) |>
    dplyr::mutate(
      Reported=forcats::fct_relevel(Reported, c("Standard", "Preferred", "Yes", "Yes (with remarks)", "No")),
      Category=forcats::fct_relevel(Category, c("Type", "Blank", "Calibration Line 1", "Calibration Line 2", "LLQC","LQC","MQC","HQC","HLQC", "NIST SRM", "NIST hTAG", "NIST T1D", "NIST Young AA"))
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
      labelsd=paste("\u03BC +/-", format(ysd, digits=2, nsmall=2))
    )

  #write_csv(calibLineDataLmPredPlot, here::here(file.path(outputDirectory,"temp.csv")))
  
  calibLineDataLmPredPlot <- calibLineDataLmPredPlot |> 
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
    geom_errorbar(aes(x=Concentration, ymin=yavgmin, ymax=yavgmax), color="blue", alpha=0.5, linetype=1) +
    geom_vline(aes(xintercept=Concentration), linetype=1, color="black", alpha=0.3) +
    geom_text(aes(x=Concentration, y=max(yavg+(ymax-ymin)*0.05), angle=90, label=labelconc), size=3, vjust = -0.5, hjust=0.5, check_overlap = TRUE, color="black") +
    ggrepel::geom_text_repel(aes(x=Concentration, y=yavgmin, label=labelsd), min.segment.length = 1, size=3, nudge_y = -0.1, point.padding=4, direction="y", color="black")
  if(!is.null(avgLipidToISRatios)) {
    avgMedCalibrationLinesLog10 <- avgMedCalibrationLinesLog10 +
      geom_ribbon(aes(x=Concentration, y=meanAreaRatios,ymin=meanAreaRatios-sdAreaRatios,ymax=meanAreaRatios+sdAreaRatios), fill="red", alpha=0.1, show.legend = FALSE) +
      geom_hline(aes(yintercept=meanAreaRatios), color="red", show.legend = FALSE, linetype=3)
  }
  avgMedCalibrationLinesLog10 <- avgMedCalibrationLinesLog10 +
    scale_color_discrete(breaks=c("blue","black"), labels=c("Mean","Median")) +
      facet_wrap(.~ceramideName, scales = "free_x") + 
    xlab(paste0("Concentration ", cunit)) + ylab("Ratio Non-labeled/Labeled") +
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
    geom_errorbar(aes(x=Concentration, ymin=yavgmin, ymax=yavgmax), color="blue", alpha=0.5, linetype=1) +
    geom_vline(aes(xintercept=Concentration), linetype=1, color="black", alpha=0.3) +
    geom_text(aes(x=Concentration, y=max(yavg+(ymax-ymin)*0.05), angle=90, label=labelconc), size=3, vjust = -0.5, hjust=0.5, check_overlap = TRUE, color="black") +
    ggrepel::geom_text_repel(aes(x=Concentration, y=yavgmin, label=labelsd), min.segment.length = 1, size=3, nudge_y = -0.1, point.padding=4, direction="y", color="black")
  if(!is.null(avgLipidToISRatios)) {
    avgMedCalibrationLinesLinear <- avgMedCalibrationLinesLinear + 
      geom_ribbon(aes(x=Concentration, y=meanAreaRatios,ymin=meanAreaRatios-sdAreaRatios,ymax=meanAreaRatios+sdAreaRatios), fill="red", alpha=0.1, show.legend = FALSE) +
      geom_hline(aes(yintercept=meanAreaRatios), color="red", show.legend = FALSE, linetype=3)
  }
  avgMedCalibrationLinesLinear <- avgMedCalibrationLinesLinear + 
    scale_color_discrete(breaks=c("blue","black"), labels=c("Mean","Median")) +
    facet_wrap(.~ceramideName, scales = "free_x") +
    xlab(paste0("Concentration ", cunit)) + ylab("Ratio Non-labeled/Labeled") +
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

scatterRatioPlot <- function(linRegFormula, data, selectedCeramideId, outputDirectory) {
  dataCer <- data |> dplyr::filter(ceramideId == selectedCeramideId)
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
  sampleTypes |> 
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
  print(paste("Creating area ratio plot for", selectedSampleType))
  dataSampleType <- data |> dplyr::filter(SampleType == selectedSampleType)
  stopifnot(length(unique(dataSampleType$SampleType))==1)
  plot <- ggplot(data=dataSampleType, 
         mapping = aes(x=LabId, y=RatioLipidToIS, fill=ceramideName)) +
    geom_boxplot(outlier.shape=outlierShape) + ylab("Cer / ISCer (area ratio)") +
    facet_grid(ceramideName~SampleType, scales="free_y") + theme + fillscale
  ggsave(file.path(outputDirectory, paste0("NIST_Non-labeled_Labeled_AreasRatioPlot-", selectedSampleType,".png")), plot)
  plot
}

calibrationLinePlot <- function(linRegFormula, lmPred, adjAnalyteConcentrations, theme=mytheme, labId, digits = 4, cunit = "", d_theoratios, outputDirectory) {
  print(paste("Creating calibration line plot for lab", labId))
  oldDigits <- getOption("digits")
  options(digits = digits)
  filteredLmPred <- lmPred |> filter(LabId==labId)
  filteredAdjAnalyteConcentrations <- adjAnalyteConcentrations |> filter(LabId==labId)
  
  filteredAdjAnalyteConcentrations <- filteredAdjAnalyteConcentrations |> 
    left_join(d_theoratios, by=c("ceramideId", "ceramideName", "Concentration"))
  
  plot <- ggplot(
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
    xlab(paste0("Actual Concentration ", cunit)) + ylab("Ratio Non-labeled/Labeled")
  ggsave(file.path(outputDirectory, paste0("CalibrationLineAnalyteConcentrationsLab",labId, ".png")), plot=plot, width=10)
  options(digits = oldDigits)
  plot
}

calibrationLineResidualsPlot <- function(linRegFormula, lmPred, adjAnalyteConcentrations, theme=mytheme, labId, digits = 4, cunit = "", d_theoratios, outputDirectory) {
  print(paste("Creating calibration line residuals plot for lab", labId))
  oldDigits <- getOption("digits")
  options(digits = digits)
  filteredLmPred <- lmPred |> filter(LabId==labId)
  plot <- ggplot(
    data=filteredLmPred, 
    mapping = aes(x=RatioLipidToIS , y=.resid)
  ) + 
    #geom_point(aes(color=LabId, shape=ceramideName)) +
    geom_point(data = filteredLmPred, aes(x= Concentration, y = .std.resid, shape = ceramideName)) +
    facet_wrap(LabId+SampleType~ceramideName, scales = "free", ncol=4) + 
    theme + 
    labs(subtitle = paste("Results for Lab", labId)) +
    xlab(paste0("Concentration ", cunit)) + ylab("Std. Residuals")
  ggsave(file.path(outputDirectory, paste0("CalibrationLineConcentrationsVsStdResidualsLab",labId, ".png")), plot=plot, width=10)
  options(digits = oldDigits)
  plot
}

calibrationLineVsSinglePointPlot <- function(analyteConcentrationsFromCalibLines, selectedCeramideId, suffix="", theme=mytheme, mycolorscale, outputDirectory) {
  dta <- analyteConcentrationsFromCalibLines |> filter(ceramideId==selectedCeramideId)
  plot <- ggplot(
      aes(x=C_SinglePoint, y=C_Adj, color=LabId, fill=LabId, label=LabId, group=LabId), 
      data=dta
    ) + 
    geom_abline(slope = 1, intercept = 0, color = "#6c9694", linewidth=.3, linetype = "longdash") +
    geom_point(shape = 21, stroke = 0.33, size = 1.2, alpha = .6) + 
    geom_hline(aes(yintercept = Concentration), color = "lightgray", linewidth=.3, linetype = "dashed", alpha = .5) +
    geom_vline(aes(xintercept = Concentration), color = "lightgray", linewidth=.3, linetype = "dashed", alpha = .5) +
    stat_smooth(method = "lm", aes( x = C_SinglePoint, y = C_Adj, color=LabId, fill=LabId), linewidth=.15, linetype = "dotted", inherit.aes = FALSE, se = FALSE, ) +
    facet_grid(rows = vars(ceramideName), cols = vars(Protocol, MassAnalyzerType), scales = "free") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "\U003BCmol/L (single-point calibration)", y = "\U003BCmol/L (multi-point calibration)") +
    labs(color=NULL, fill = NULL)+
    theme_classic( base_size =  9) +  
    theme(
      axis.text = element_text(size = 7), 
      axis.title = element_text(size = 8, face = "bold"), 
      strip.background =element_rect(fill="grey99")
    ) + theme(aspect.ratio = 1) +
    theme(legend.position.inside="bottom",
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.1, "cm"),
          legend.spacing.y = unit(.001, 'cm')
    )
  ggsave(file.path(outputDirectory, paste0("CalibrationLineVsSinglePoint-", selectedCeramideId, suffix, ".png")), width = 8, height = 4, plot=plot)
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
  xlab(paste0("Concentration ", cunit)) + ylab("Ratio Non-labeled/Labeled") +
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
  qcSampleStatsLong <- qcSampleStats |> 
    ungroup() |> 
    select(
      -SdAvg_C_Adj, 
      -CV_C_Adj,
      -SdAvg_C_SinglePoint,
      -CV_C_SinglePoint) |> pivot_longer(
        cols = c("AvgAvg_C_Adj", "Avg_C_SinglePoint"),
        names_to = c("CalibrationType"),
        values_to = c("MeanAdjustedConcentration")
      ) |>
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
  dta <- averagedConcentrations |> filter(ceramideId==selectCeramideId)
  loDf <- data.frame(
    LoD=dta$Avg_C_LoD_C, 
    LoQ=dta$Avg_C_LoQ_C, 
    LabId=dta$LabId, 
    Sample=dta$Sample,
    label=paste(dta$Sample, " ", dta$Concentration, " ", dta$Unit)
  ) |> filter(
    Sample=="STD 6"
  ) |> distinct()
  loDfLong <- loDf |> pivot_longer(cols = c("LoD", "LoQ"), names_to = "FoM", values_to = "y")
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
    labs(caption = paste0("Avg. LoD from calib. lines for STD 6, where LoD=3.3*SE(intercept)/slope\n",
           "Avg. LoQ from calib. lines for STD 6, where LoQ=10*SE(intercept)/slope")) +
    theme + 
    fillscale
  ggsave(file.path(outputDirectory, paste0("CalibrationLineSTDsCer",selectCeramideId,".png")), plot=boxPlot, width=20)
  boxPlot
}

################################################################################
# Stds Variance Histogram Plot
################################################################################
stdsVarHistoPlot <- function(averagedConcentrations, selectCeramideId, theme=mythemeXRot, fillscale=scale_fill_nejm(), outputDirectory) {
  dta <- averagedConcentrations |> 
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
  dta <- averagedConcentrations |> 
    filter(ceramideId==selectCeramideId) |> 
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
  nistCvPlot <- ggplot(aes(x=SampleType, y=CV, fill=Calibration), data=nistConcentrationsStatsLong |>
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
    ggsave(file.path(outputDirectory, paste0("NISTSamplesConcentrationsBySampleType-", namesuffix, ".png")), plot, width=width)
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
  geom_hline(aes(yintercept = 20), color="red", linetype = 2, show.legend = FALSE) +
  geom_hline(aes(yintercept = 10), color="yellow", linetype = 2, show.legend = FALSE) +
  geom_hline(aes(yintercept = 5), color="green", linetype = 2, show.legend = FALSE) +
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
