plot_survival <- function(data, factors, cluster_name, sample_filter = NULL) {
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(survival)
  library(survminer)
  library(ggpubr)

  if (!"sampInfo" %in% names(data)) stop("data must contain 'sampInfo'")

  sampInfo <- data$sampInfo

  # Finding 3 fix: Use explicit sample_filter parameter instead of undefined 'whitelist'
 # If sample_filter is provided as a logical vector, use it; otherwise keep all samples
  if (!is.null(sample_filter)) {
    if (is.logical(sample_filter) && length(sample_filter) == nrow(sampInfo)) {
      sampInfo <- sampInfo[sample_filter, , drop = FALSE]
    } else if (is.character(sample_filter) || is.numeric(sample_filter)) {
      sampInfo <- sampInfo[sample_filter, , drop = FALSE]
    }
  }

  # Remove samples with NA cluster assignment (e.g., those filtered during clustering)
  sampInfo <- sampInfo[!is.na(sampInfo[[cluster_name]]), , drop = FALSE]
  
  sampInfo$cluster <- sampInfo[[cluster_name]]
  
  if(length(table(sampInfo$cluster))>1){
    
    nclass <- length(unique(sampInfo$cluster))
    cols <- determine_colors(factors, nclass)
    col.colors <- cols$col.colors
    
    # Fit survival curves
    fit <- survfit(Surv(time, event) ~ cluster, data = sampInfo)
    fit2 <- coxph(Surv(time, event) ~ as.factor(cluster), data = sampInfo)
    bic <- BIC(fit2)
    hr <- summary(fit2)$coefficients[, 2]
    
    # Survival plot and risk table
    splot <- ggsurvplot(fit, data = sampInfo, risk.table = TRUE, pval = TRUE)
    
    risk.table <- ggtexttable(
      t(
        splot$data.survtable %>%
          select(strata, time, n.risk) %>%
          pivot_wider(names_from = strata, values_from = n.risk, id_cols = time)
      ),
      theme = ttheme(tbody.style = tbody_style(
        color = c("black", col.colors),
        fill = "white",
        face = "bold"
      ))
    )
    
    # Build annotation label
    hr_labels <- paste0("HR", 1:length(hr), " = ", round(hr, 2), collapse = "\n")
    # title_text <- paste0(str_extract(cluster_name, "alpha=[^\\.]+"), 
    #                      " factor ", gsub(".*factor(\\d+).*", "\\1", cluster_name), 
    #                      ", ", labels[as.numeric(gsub(".*factor(\\d+).*", "\\1", cluster_name))])
    
    surv.plot <- splot$plot +
      geom_text(x = max(splot$data.survtable$time) - 5, y = 0.85, size = 5,
                label = paste0("BIC = ", round(bic, 2), "\n", hr_labels)) +
      scale_color_manual(labels = as.character(1:nclass), values = col.colors)
    
    final_plot <- ggarrange(surv.plot, risk.table, nrow = 2, heights = c(3, 1))

    print(final_plot)

  }
}
