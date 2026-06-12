# Functions to interpret sensitivities

## plotDAO
# Plot profiles of DIC, Alkalinity, and TOC in separate panels, with a co-variable
plotDAO <- function(runs,
                    unitplot,
                    breaks,
                    labels,
                    varia){
  ### Plot - DIC, Alkalinity - Organic C
  plot_df <- map_dfr(runs, function(run) {
    
    tibble(
      depthcm = run$output$depth * 100,        # Depth to cm
      DIC     = run$output$y[, "DIC"] / 1000,  # mmol m-3 -> mM
      ALK     = run$output$y[, "ALK"] / 1000,  # mmol m-3 -> mM
      TOC     = run$output$OCpercent        ,  # TOC as percent
      w       = run[[1]]
    )
    
  }) %>%
    filter(depthcm <= 25) %>%
    pivot_longer(
      cols = c(DIC, ALK, TOC),
      names_to = "variable",
      values_to = "value"
    )
  
  plot_df$variable <- factor(
    plot_df$variable,
    levels = c("DIC", "ALK", "TOC"),
    labels = c(
      "DIC (mM)",
      "Alkalinity (mM)",
      "Organic Carbon (%)"
    )
  )
  
  dfout <- plot_df
  
  plotout <- ggplot(
    plot_df,
    aes(
      x = value,
      y = depthcm,
      group = interaction(w, variable),
      colour = w
    )
  ) +
    geom_path(linewidth = 1) +
    facet_wrap2(~variable, scales = "free_x") +
    facetted_pos_scales(
      x = list(
        NULL,              # DIC
        NULL,              # ALK
        scale_x_log10()    # OC
      )
    ) +
    scale_y_reverse(limits = c(25, 0),
                    expand = c(0,0)) +
    
    scale_colour_viridis_c(
      name = parse(text = paste0(varia, unitplot))[[1]],
      option = "plasma",
      breaks = breaks,
      labels = labels,
      trans = "log10"
    ) +
    labs(
      x = NULL,
      y = "Depth (cm)",
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.6),
      axis.ticks = element_line(colour = "black", linewidth = 0.6),
      axis.text = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      panel.spacing = unit(1.2, "lines")
    )
  return(list(dfout, plotout))
}


## PlotCbudg
# Plot DIC production by different OM mineralization processes along a covariable
# as a barplot.
plotCbudg <- function(runs,
                      unitplot,
                      varia){
  
  ### C mineralization per mineralization process
  budget_df <- map_dfr(runs, function(run){
    
    bud <- FEMSDIA_get_budgetC(run$output)
    
    tibble(
      w = run[[1]],
      OxicMineralisation = bud$Rates[["OxicMineralisation"]],
      Denitrification    = bud$Rates[["Denitrification"]],
      MnReduction        = bud$Rates[["MnReduction"]],
      FeReduction        = bud$Rates[["IronReduction"]],
      SulphateReduction  = bud$Rates[["SulphateReduction"]],
      Methanogenesis     = bud$Rates[["Methanogenesis"]]
    )
  })
  
  
  
  budget_long <- budget_df %>%
    tidyr::pivot_longer(
      cols = -w,
      names_to = "Process",
      values_to = "Rate"
    )
  
  budget_long$Process <- factor(
    budget_long$Process,
    levels = rev(c(
      "OxicMineralisation",
      "Denitrification",
      "MnReduction",
      "FeReduction",
      "SulphateReduction",
      "Methanogenesis"
    ))
  )
  
  dfout <- budget_df
  
  plotout <- ggplot(budget_long, aes(x = factor(w), y = Rate, fill = Process)) +
    geom_col(width = 0.85, colour = "black", width = 0.2) +
    
    scale_fill_viridis_d(
      option = "plasma",
      begin = 0.,
      end = 1,
      alpha = 0.8
    )  +
    
    labs(
      x = parse(text = paste0(varia, unitplot))[[1]],
      y = expression("Mineralisation rate (mmol  " * m^{-2} * d^{-1} * ")")
    ) +
    theme_classic(base_size = 13) +
    theme(axis.line = element_line(colour = "black",
                                   linewidth = 0.6),
          axis.ticks = element_line(colour = "black",
                                    linewidth = 0.6),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          legend.position = "right"
    )
  return(list(dfout, plotout))
}

## PlotDICinorg
# Plot DIC consumption, dissolution, and burial by processes related to PIC dynamics
# Contrasts net dissolution (DIC production - positive), with PIC burial (negative).
plotDICinorg <- function(runs,
                         unitplot,
                         varia){
  
  ### Carbonate budget
  budget_carbonate_df <- map_dfr(runs, function(run){
    
    tibble(
      w = run[[1]],
      CarbonateDissolution   = run$output$TotCALCdiss,
      CarbonatePrecipitation = run$output$TotCALCprod,
      AragoniteDissolution   = run$output$TotARAGdiss,
      CarbonateBurial        = -run$output$CALCdeepflux,
      AragoniteBurial        = -run$output$ARAGdeepflux,
      CarbonateNetDIC        = run$output$TotCALCdiss - run$output$TotCALCprod,
      AragoniteNetDIC        = run$output$TotARAGdiss
    )
  })
  
  
  budget_carbonate_long <- budget_carbonate_df %>%
    tidyr::pivot_longer(
      cols = -w,
      names_to = "Process",
      values_to = "Rate"
    )
  
  dfout <- budget_carbonate_df
  
  plotout <- ggplot(subset(budget_carbonate_long, Process %in% c("CarbonateNetDIC", "AragoniteNetDIC", "CarbonateBurial", "AragoniteBurial")),
                    aes(x = factor(w), y = Rate, fill = Process)) +
    geom_col(width = 0.85, colour = "black", width = 0.2) +
    
    scale_fill_viridis_d(
      option = "plasma",
      begin = 0.,
      end = 1,
      alpha = 0.8
    )  +
    
    labs(
      x = parse(text = paste0(varia, unitplot))[[1]],
      y = expression("PIC dissolution or burial (mmol " * m^{-2} * d^{-1} * ")")
    ) +
    theme_classic(base_size = 13) +
    theme(axis.line = element_line(colour = "black",
                                   linewidth = 0.6),
          axis.ticks = element_line(colour = "black",
                                    linewidth = 0.6),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          legend.position = "right"
    )
  return(list(dfout, plotout))
}

## Get JAlkmin Krumins netprocesses
# Get net Alkalinity generation by mineralization processes

get_JAlk_mineralisation <- function(out) {
  
  tibble::tibble(
    Process = c(
      "Net ammonification",
      "Denitrification",
      "Net iron reduction",
      "Net sulfate reduction",
      "Net methanogenesis"
    ),
    
    JAlk = c(
      out$TotNH3prodMin - out$TotNitri1 - out$TotAnammox,  # net NH4+ accumulation/production
      0.8 * out$TotDenit,
      8 * out$TotFered - 2 * out$TotFeoxid - out$TotFeMnoxid,
      out$TotBSR - 2 * out$TotH2Soxid -
        2 * out$TotFeSoxid -
        4 * out$TotFeS2oxid +
        4 * out$TotH2SFeoxid +
        2 * out$TotH2SMnoxid,
      0
    )
  )
}

## PlotAlkProd netprocesses
# Plot net alkalinity generation by mineralization processes

plotAlkProd <- function(runs,
                        unitplot,
                        varia){
  
  JAlk_min_df <- purrr::map_dfr(runs, function(run) {
    get_JAlk_mineralisation(run$output) |>
      dplyr::mutate(w = run[[1]])
  })
  
  JAlk_min_wide <- JAlk_min_df %>%
    pivot_wider(
      names_from = names(JAlk_min_df)[3],
      values_from = JAlk
    )
  
  dfout <- JAlk_min_wide
  
  plotout <- ggplot(JAlk_min_df, aes(x = factor(w), y = JAlk, fill = Process)) +
    geom_col(width = 0.85, colour = "black", linewidth = 0.2) +
    scale_fill_viridis_d(option = "plasma",
                         begin = 0.,
                         end = 1,
                         alpha = 0.8) +
    labs(
      x = parse(text = paste0(varia, unitplot))[[1]],
      y = expression(J[Alk]~"(mmol eq m"^{-2}~d^{-1}*")"),
      fill = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(axis.line = element_line(colour = "black",
                                   linewidth = 0.6),
          axis.ticks = element_line(colour = "black",
                                    linewidth = 0.6),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          legend.position = "right"
    )
  return(list(dfout, plotout))
}

## Alkbudget- KRUMINS
# Get alkalinity budget following Krumins 2013 (Jalk)

get_alk_krumins_style <- function(out) {
  
  # Net process terms
  net_ammonification <- out$TotNH3prodMin - out$TotNitri1 -out$TotAnammox
  
  net_fe_reduction <- 8 * out$TotFered - 2 * out$TotFeoxid - 1 * out$TotFeMnoxid
  net_sulfate_reduction <- out$TotBSR - 2 * out$TotH2Soxid - 2 * out$TotFeSoxid -
    4 * out$TotFeS2oxid + 4 * out$TotH2SFeoxid + 2 * out$TotH2SMnoxid
  
  # Carbonate dissolution
  calcite_diss <- 2 * out$TotCALCdiss
  aragonite_diss <- 2 * out$TotARAGdiss
  
  # Burial terms
  alk_burial <- -2 * (out$CALCdeepflux + out$ARAGdeepflux)
  
  reduced_burial <- 2 * out$FeSdeepflux + 4 * out$FeS2deepflux # is not in output + 2 * out$FeCO3deepflux
  
  # Water-column reoxidation of escaping sulfide (see Krumins)
  sulfide_wc_oxid <- -2 * abs(out$H2Sflux)
  
  tibble(
    Process = c(
      "Net sulfate reduction",
      "Net ammonification",
      "Denitrification",
      "Net iron reduction",
      "Calcite dissolution",
      "Aragonite dissolution",
      "Alk burial",
      "Reduced S / Fe burial",
      "Sulfides oxidized in water column"
    ),
    
    JAlk = c(
      net_sulfate_reduction,
      net_ammonification,
      0.8 * out$TotDenit,
      net_fe_reduction,
      calcite_diss,
      aragonite_diss,
      abs(alk_burial),
      reduced_burial,
      abs(sulfide_wc_oxid)
    )
  )
}

## PlotAlkbudget KRUMINS
plot_AlkKrumins <- function(runs,
                            unitplot,
                            varia){
  
  alk_df <- map_dfr(runs, function(run) {
    get_alk_krumins_style(run$output) %>%
      mutate(w = run[[1]])
  })
  
  alk_df_wide <- alk_df %>%
    pivot_wider(
      names_from = names(alk_df)[3],
      values_from = JAlk
    )
  
  dfout <- alk_df_wide
  
  # left bar: exclude water-column sulfide oxidation
  # right bar: include water-column sulfide oxidation
  alk_plot_df <- alk_df %>%
    mutate(
      Bar = case_when(
        Process %in% c(
          "Net sulfate reduction",
          "Net ammonification",
          "Denitrification",
          "Net iron reduction",
          "Calcite dissolution",
          "Aragonite dissolution"
        ) ~ "JAlk",
        
        Process %in% c(
          "Alk burial",
          "Reduced S / Fe burial",
          "Sulfides oxidized in water column"
        ) ~ "JAlk*"
      )
    ) %>%
    filter(!is.na(Bar)) %>%
    mutate(
      Bar = factor(Bar, levels = c("JAlk", "JAlk*")),
      Process = factor(
        Process,
        levels = c(
          "Net sulfate reduction",
          "Net ammonification",
          "Denitrification",
          "Net iron reduction",
          "Calcite dissolution",
          "Aragonite dissolution",
          "Alk burial",
          "Reduced S / Fe burial",
          "Sulfides oxidized in water column"
        )
      )
    )
  
  alk_plot_df <- alk_plot_df %>%
    mutate(
      w_factor = factor(w, levels = sort(unique(w))),
      x_base   = as.numeric(w_factor),
      x_pos    = ifelse(Bar == "JAlk", x_base - 0.17, x_base + 0.17)
    )
  
  plotout <- ggplot(
    alk_plot_df,
    aes(
      x = x_pos,
      y = JAlk,
      fill = Process
    )
  ) +
    geom_col(
      width = 0.32,
      colour = "black",
      linewidth = 0.2
    ) +
    scale_x_continuous(
      breaks = unique(alk_plot_df$x_base),
      labels = levels(alk_plot_df$w_factor)
    )+
    scale_fill_viridis_d(option = "plasma",
                         begin = 0.,
                         end = 1,
                         alpha = 0.8) +
    labs(
      x = parse(text = paste0(varia, unitplot))[[1]],
      y = expression("Alkalinity generation and fate (mmol eq " * m^{-2} * d^{-1} * ")"),
      fill = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(axis.line = element_line(colour = "black",
                                   linewidth = 0.6),
          axis.ticks = element_line(colour = "black",
                                    linewidth = 0.6),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          legend.position = "right"
    )
  return(list(dfout, plotout))
}
