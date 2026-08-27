
# Functions to interpret sensitivities                                       #
# Last checked 24/8/2026                                                     #
# author: Emil De Borger                                                     #
# Functions to extract data and plots from sensitivity runs and case studies.#


# Functions for the sensitivity analyses ####
## 1. Function to extract net alkalinity balance from only carbonate dynamics ####
# and pyrite burial.

get_netalklongterm <- function(out) {
  
  # Carbonate dissolution
  carbonate_diss <- 2 * (out$TotCALCdiss + out$TotARAGdiss)
  carbonate_prod <- - 2 * (out$TotCALCprod)
  
  # Burial terms
  # alk_burial <- -2 * (out$CALCdeepflux + out$ARAGdeepflux)
  pyrite_burial <- 4 * out$FeS2deepflux
  
  tibble(
    Process = c(
      "Carbonate dissolution",
      "Carbonate precipitation",
      "Pyrite burial"
    ),
    
    JAlk = c(
      carbonate_diss,
      carbonate_prod,
      pyrite_burial
    )
  )
}

## 2. Functions that use function above to plot net long term alkalinity flux.####

### 2.1. Plot net alkalinity budget for MAR (unit is converted from w to MAR). ####
plot_netalklongtermMAR <- function(runs,
                                unitplot,
                                varia){
  
  JAlk_lt_df <- purrr::map_dfr(runs, function(run) {
    get_netalklongterm(run$output) |>
      dplyr::mutate(w = run[[1]])
  })
  
  JAlk_lt_wide <- JAlk_lt_df %>%
    pivot_wider(
      names_from = names(JAlk_lt_df)[3],
      values_from = JAlk
    )

  
  dfout <- JAlk_lt_wide
    
  net_df <- JAlk_lt_df %>%
    dplyr::group_by(w) %>%
    dplyr::summarise(net_JAlk = sum(JAlk), .groups = "drop")
  
  plotout <- ggplot(JAlk_lt_df, aes(x = factor(w), y = JAlk, fill = Process)) +
    geom_col(width = 0.85, colour = "black", linewidth = 0.2) +
    geom_point(
      data = net_df,
      aes(x = factor(w), y = net_JAlk),
      inherit.aes = FALSE,
      colour = "black",
      size = 3
    ) +
    scale_fill_viridis_d(option = "viridis", #"plasma"
                         begin = 0.,
                         end = 1,
                         alpha = 0.7) +
    scale_x_discrete(
      labels = function(x) {
        MAR_g <- as.numeric(as.character(x)) * 0.5 * 10 * 1000
        format(round(MAR_g, 1), big.mark = ",", scientific = FALSE)
      }
    ) +
    labs(
      x = expression(
        "Mass accumulation rate (g m"^{-2}~"yr"^{-1}*")"
      )
    ) +
    labs(
      x = parse(text = paste0(varia, unitplot))[[1]],
      y = expression("Net benthic AT flux"~~~"(mmol m"^{-2}~d^{-1}*")"),
      fill = "Process"
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

### 2.2. Plot alkalinity budget for other variables. ####

plot_netalklongterm <- function(runs,
                                   unitplot,
                                   varia){
  
  JAlk_lt_df <- purrr::map_dfr(runs, function(run) {
    get_netalklongterm(run$output) |>
      dplyr::mutate(w = run[[1]])
  })
  
  JAlk_lt_wide <- JAlk_lt_df %>%
    pivot_wider(
      names_from = names(JAlk_lt_df)[3],
      values_from = JAlk
    )
  
  dfout <- JAlk_lt_wide
  
  net_df <- JAlk_lt_df %>%
    dplyr::group_by(w) %>%
    dplyr::summarise(net_JAlk = sum(JAlk), .groups = "drop")
  
  plotout <- ggplot(JAlk_lt_df, aes(x = factor(w), y = JAlk, fill = Process)) +
    geom_col(width = 0.85, colour = "black", linewidth = 0.2) +
    geom_point(
      data = net_df,
      aes(x = factor(w), y = net_JAlk),
      inherit.aes = FALSE,
      colour = "black",
      size = 3
    ) +
    scale_fill_viridis_d(option = "viridis", #"plasma"
                         begin = 0.,
                         end = 1,
                         alpha = 0.7) + 
    labs(
      x = parse(text = paste0(varia, unitplot))[[1]],
      y = expression("Net benthic " * A[T] ~ "flux"~~~"(mmol m"^{-2}~d^{-1}*")"),
      fill = "Process"
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

# Functions for the simulations ####
## 1. Extract data to plot during fitting ####
## Not used in chapter itself.

prepare_profile_comparison <- function(
    model_output,
    measured_file,
    variables,
    sheet = "Profiles",
    measured_filters = list(),
    measured_depth_col = "MidDepth",
    error_columns = NULL,
    max_depth_cm = 20
) {
  
  model_long <- as.data.frame(model_output$y) %>%
    dplyr::mutate(depth_cm = model_output$depth * 100) %>%
    dplyr::filter(depth_cm <= max_depth_cm) %>%
    dplyr::select(depth_cm, dplyr::all_of(variables)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(variables),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Error = NA_real_,
      Source = "Model"
    )
  
  measured <- readxl::read_excel(
    measured_file,
    sheet = sheet
  )
  
  if (length(measured_filters) > 0) {
    for (filter_name in names(measured_filters)) {
      measured <- measured %>%
        dplyr::filter(
          .data[[filter_name]] %in%
            measured_filters[[filter_name]]
        )
    }
  }
  
  measured_long <- purrr::map_dfr(
    variables,
    function(variable_name) {
      
      error_name <- if (is.null(error_columns)) {
        paste0(variable_name, "_ERR")
      } else {
        error_columns[[variable_name]]
      }
      
      error_values <- if (
        !is.null(error_name) &&
        error_name %in% names(measured)
      ) {
        as.numeric(measured[[error_name]])
      } else {
        rep(NA_real_, nrow(measured))
      }
      
      tibble::tibble(
        depth_cm = as.numeric(
          measured[[measured_depth_col]]
        ),
        Variable = variable_name,
        Value = as.numeric(
          measured[[variable_name]]
        ),
        Error = error_values,
        Source = "Measured"
      )
    }
  ) %>%
    dplyr::filter(
      !is.na(Value),
      depth_cm <= max_depth_cm
    )
  
  list(
    model = model_long,
    measured = measured_long,
    combined = dplyr::bind_rows(
      model_long,
      measured_long
    )
  )
}

## 2. Plot data to fit ####
plot_profile_comparison <- function(
    profile_data,
    variable_labels = NULL,
    max_depth_cm = 20,
    x_label = expression(
      "Concentration (mmol m"^{-3}*")"
    ),
    y_label = "Depth (cm)",
    nrow = 1,
    model_linewidth = 1,
    point_size = 2.5,
    error_linewidth = 0.8,
    error_cap_height = 0.15,
    show_errors = TRUE
) {
  
  model_df <- profile_data$model
  measured_df <- profile_data$measured
  
  facet_labeller <- if (is.null(variable_labels)) {
    label_value
  } else {
    as_labeller(variable_labels)
  }
  
  p <- ggplot() +
    
    geom_path(
      data = model_df,
      aes(
        x = Value,
        y = depth_cm,
        group = Variable
      ),
      linewidth = model_linewidth
    )
  
  if (
    show_errors &&
    any(!is.na(measured_df$Error))
  ) {
    
    error_df <- measured_df %>%
      filter(
        !is.na(Error),
        Error > 0
      ) %>%
      mutate(
        xmin = Value - Error,
        xmax = Value + Error
      )
    
    p <- p +
      
      # Horizontal error line
      geom_segment(
        data = error_df,
        aes(
          x = xmin,
          xend = xmax,
          y = depth_cm,
          yend = depth_cm
        ),
        linewidth = error_linewidth
      ) +
      
      # Left-hand cap
      geom_segment(
        data = error_df,
        aes(
          x = xmin,
          xend = xmin,
          y = depth_cm - error_cap_height,
          yend = depth_cm + error_cap_height
        ),
        linewidth = error_linewidth
      ) +
      
      # Right-hand cap
      geom_segment(
        data = error_df,
        aes(
          x = xmax,
          xend = xmax,
          y = depth_cm - error_cap_height,
          yend = depth_cm + error_cap_height
        ),
        linewidth = error_linewidth
      )
  }
  
  p +
    geom_point(
      data = measured_df,
      aes(
        x = Value,
        y = depth_cm
      ),
      shape = 21,
      fill = "white",
      size = point_size,
      stroke = 0.8
    ) +
    facet_wrap(
      ~Variable,
      scales = "free_x",
      nrow = nrow,
      labeller = facet_labeller
    ) +
    scale_y_reverse(
      limits = c(max_depth_cm, 0),
      expand = c(0, 0)
    ) +
    labs(
      x = x_label,
      y = y_label
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = grid::unit(1.2, "lines")
    )
}

# Additional functions not used in chapter itself ####

## 1. A version of the above that includes carbonate burial (not used in the book chapter).  ####

get_netalklongterm2 <- function(out) {
  
  # Carbonate dissolution
  carbonate_diss <- 2 * (out$TotCALCdiss + out$TotARAGdiss)
  carbonate_prod <- - 2 * (out$TotCALCprod)
  carbonate_burial <- - 2 * (out$CALCdeepflux + out$ARAGdeepflux)
  
  # Burial terms
  # alk_burial <- -2 * (out$CALCdeepflux + out$ARAGdeepflux)
  pyrite_burial <- 4 * out$FeS2deepflux
  
  tibble(
    Process = c(
      "Carbonate dissolution",
      "Carbonate precipitation",
      "Carbonate burial",
      "Pyrite burial"
    ),
    
    JAlk = c(
      carbonate_diss,
      carbonate_prod,
      carbonate_burial,
      pyrite_burial
    )
  )
}


## PlotAlkbudget KRUMINS
plot_netalklongterm2 <- function(runs,
                                 unitplot,
                                 varia){
  
  JAlk_lt_df <- purrr::map_dfr(runs, function(run) {
    get_netalklongterm2(run$output) |>
      dplyr::mutate(w = run[[1]])
  })
  
  JAlk_lt_wide <- JAlk_lt_df %>%
    pivot_wider(
      names_from = names(JAlk_lt_df)[3],
      values_from = JAlk
    )
  
  dfout <- JAlk_lt_wide
  
  net_df <- JAlk_lt_df %>%
    dplyr::group_by(w) %>%
    dplyr::summarise(net_JAlk = sum(JAlk), .groups = "drop")
  
  plotout <- ggplot(JAlk_lt_df, aes(x = factor(w), y = JAlk, fill = Process)) +
    geom_col(width = 0.85, colour = "black", linewidth = 0.2) +
    geom_point(
      data = net_df,
      aes(x = factor(w), y = net_JAlk),
      inherit.aes = FALSE,
      colour = "black",
      size = 3
    ) +
    scale_fill_viridis_d(option = "plasma",
                         begin = 0.,
                         end = 1,
                         alpha = 0.8) +
    labs(
      x = parse(text = paste0(varia, unitplot))[[1]],
      y = expression("Alkalinity generation"~"(mmol eq m"^{-2}~d^{-1}*")"),
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

## 2. plot DIC, alkalininty, and organic C profiles of sensitivity runs. ####
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


## 3. Plot DIC production ####
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

## 4. Plot DIC consumption and production ####
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

## 5. Get JAlkmin similar to Krumins et al., 2013, netprocesses ####
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

## 6. Get alkalinty budget ####
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
  
  reduced_burial <- 2 * out$FeSdeepflux + 4 * out$FeS2deepflux 
  # is not in output but theoretically  + 2 * out$FeCO3deepflux
  
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

