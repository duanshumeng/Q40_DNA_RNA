# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Read data（df_quartet.mcr from Fig3-4-WES_combine.R)

df_new <- df_quartet.mcr

# Create visualization function (adapted to your data structure)
plot_mcr_by_rmsk_class <- function(data, 
                                   rmsk_class,
                                   type = "SNV",
                                   title = NULL,
                                   colors = c("steelblue", "lightcoral"),
                                   text_size = 1,
                                   show_stats = TRUE,
                                   stats_position = "bottom",
                                   save_plot = FALSE,
                                   output_dir = ".",
                                   file_format = "png",
                                   width = 10,
                                   height = 7,
                                   dpi = 300) {
  
  # Filter data
  df_filtered <- data %>%
    filter(type == !!type, Rmsk_class == !!rmsk_class)
  
  if (nrow(df_filtered) == 0) {
    warning(paste("No data found for type", type, "and Rmsk_class", rmsk_class))
    return(NULL)
  }
  
  # If MCR column doesn't exist but Mendelian_Concordance_Rate exists, use it
  if (!"MCR" %in% names(df_filtered) && "Mendelian_Concordance_Rate" %in% names(df_filtered)) {
    df_filtered$MCR <- df_filtered$Mendelian_Concordance_Rate
  }
  
  # Calculate averages by Quality group
  grouped <- df_filtered %>%
    dplyr::group_by(Quality) %>%
    dplyr::summarise(
      Total_Variants = mean(Total_Variants),
      Mendelian_Concordant_Variants = mean(Mendelian_Concordant_Variants),
      MCR = mean(MCR),
      .groups = "drop"
    )
  
  # Check if both Q30 and Q40 data are present
  if (!all(c("Q30", "Q40") %in% grouped$Quality)) {
    warning_message <- paste("Data does not contain both Q30 and Q40 quality levels for", rmsk_class, "-", type)
    warning(warning_message)
    show_stats <- FALSE
  }
  
  # Initialize variables
  line_color <- "darkgreen"  # Default line color
  text_str <- NULL
  
  # Calculate percentage change and determine line color (if both Q30 and Q40 exist)
  if (all(c("Q30", "Q40") %in% grouped$Quality) && show_stats) {
    tv_q30 <- grouped %>% filter(Quality == "Q30") %>% pull(Total_Variants)
    tv_q40 <- grouped %>% filter(Quality == "Q40") %>% pull(Total_Variants)
    mcv_q30 <- grouped %>% filter(Quality == "Q30") %>% pull(Mendelian_Concordant_Variants)
    mcv_q40 <- grouped %>% filter(Quality == "Q40") %>% pull(Mendelian_Concordant_Variants)
    mcr_q30 <- grouped %>% filter(Quality == "Q30") %>% pull(MCR)
    mcr_q40 <- grouped %>% filter(Quality == "Q40") %>% pull(MCR)
    
    tv_increase <- ((tv_q40 - tv_q30) / tv_q30) * 100
    mcv_increase <- ((mcv_q40 - mcv_q30) / mcv_q30) * 100
    mcr_change <- mcr_q40 - mcr_q30
    
    # Determine line color based on MCR change direction
    if (mcr_change > 0) {
      line_color <- "#2E8B57"  # Green for increase (SeaGreen)
      change_symbol <- "↑"
      change_direction <- "increased"
    } else if (mcr_change < 0) {
      line_color <- "#DC143C"  # Red for decrease (Crimson)
      change_symbol <- "↓"
      change_direction <- "decreased"
    } else {
      line_color <- "darkgreen"  # Default for no change
      change_symbol <- "→"
      change_direction <- "unchanged"
    }
    
    # Create statistics text with color-coded direction
    text_str <- sprintf(
      "Q40 vs Q30:\n• Total variants change: %.1f%%\n• Concordant variants change: %.1f%%\n• MCR: %.3f (Q30) → %.3f (Q40) [%s%.3f, %s]",
      tv_increase, mcv_increase, mcr_q30, mcr_q40, 
      ifelse(mcr_change > 0, "+", ""), mcr_change, change_direction
    )
  }
  
  # Convert data to long format for plotting
  grouped_long <- grouped %>%
    select(-MCR) %>%
    pivot_longer(
      cols = c(Total_Variants, Mendelian_Concordant_Variants),
      names_to = "Variant_Type",
      values_to = "Count"
    )
  
  # Set factor levels to ensure correct order
  grouped_long$Variant_Type <- factor(grouped_long$Variant_Type,
                                      levels = c("Total_Variants", "Mendelian_Concordant_Variants"))
  
  # Set color mapping
  bar_colors <- setNames(colors, c("Total_Variants", "Mendelian_Concordant_Variants"))
  
  # Calculate scaling factor
  max_count <- max(grouped_long$Count, na.rm = TRUE)
  max_mcr <- max(grouped$MCR, na.rm = TRUE)
  min_mcr <- min(grouped$MCR, na.rm = TRUE)
  
  # Avoid division by zero
  if (max_mcr <= 0) {
    scale_factor <- 1
  } else {
    scale_factor <- max_count / max_mcr
  }
  
  # Set Y-axis range
  y_buffer <- max_count * 0.15  # Space for bar labels
  y_max <- max_count + y_buffer
  
  # Generate title automatically
  if (is.null(title)) {
    title <- paste(rmsk_class, "Region:", type, "- Variant Counts and MCR by Quality")
  }
  
  # Create base plot
  p <- ggplot() +
    # Bar plot layer
    geom_bar(data = grouped_long,
             aes(x = Quality, y = Count, fill = Variant_Type),
             stat = "identity", position = position_dodge(width = 0.7),
             width = 0.6, alpha = 0.8) +
    scale_fill_manual(values = bar_colors,
                      labels = c("Total Variants", "Mendelian Concordant Variants")) +
    # Add bar value labels
    geom_text(data = grouped_long,
              aes(x = Quality, y = Count, label = round(Count, 0), group = Variant_Type),
              position = position_dodge(width = 0.7),
              vjust = -0.5, size = 3.5 * text_size) +
    # Line plot layer (using scaled MCR values)
    # Line color determined by MCR change direction
    geom_line(data = grouped,
              aes(x = Quality, y = MCR * scale_factor, group = 1),
              color = line_color, size = 1.5) +
    # Points remain with default color
    geom_point(data = grouped,
               aes(x = Quality, y = MCR * scale_factor),
               size = 3) +
    # Add line plot value labels (showing original MCR values)
    geom_text(data = grouped,
              aes(x = Quality, y = MCR * scale_factor, label = sprintf("%.3f", MCR)),
              vjust = -1, size = 4 * text_size, fontface = "bold") +
    # Dual axis setup - use line_color for the secondary axis to match line color
    scale_y_continuous(
      name = "Variant Count",
      limits = c(0, y_max),
      sec.axis = sec_axis(~./scale_factor, 
                          name = "MCR (Mendelian Concordant / Total)",
                          breaks = seq(round(min_mcr, 2), round(max_mcr, 2), length.out = 5))
    ) +
    labs(x = "Quality", fill = "Variant Type", title = title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14 * text_size, face = "bold", 
                                margin = margin(b = 20)),
      axis.title.y = element_text(size = 12 * text_size, color = "black"),
      # Use line_color for the secondary axis label to match line color
      axis.title.y.right = element_text(size = 12 * text_size, color = line_color),
      axis.text.y.right = element_text(color = line_color),
      axis.title.x = element_text(size = 12 * text_size),
      legend.position = "top",
      legend.title = element_text(size = 11 * text_size),
      legend.text = element_text(size = 10 * text_size),
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank()
    )
  
  # Add statistics box
  if (show_stats && !is.null(text_str)) {
    # Determine statistics box position
    if (stats_position == "top") {
      y_pos <- y_max * 0.9
      vjust_val <- 1
    } else {
      y_pos <- y_max * 0.1
      vjust_val <- 0
    }
    
    p <- p +
      annotate("text", 
               x = 1.5,  # Centered
               y = y_pos,
               label = text_str, 
               hjust = 0.5, 
               vjust = vjust_val,
               size = 3.5 * text_size, 
               color = "black",
               fontface = "plain",
               alpha = 0.9)
  }
  
  # Save plot
  if (save_plot) {
    # Ensure output directory exists
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Create filename (remove special characters)
    safe_class <- gsub("[^a-zA-Z0-9]", "_", rmsk_class)
    safe_type <- gsub("[^a-zA-Z0-9]", "_", type)
    filename <- paste0(safe_class, "_", safe_type, "_MCR_comparison.", file_format)
    filepath <- file.path(output_dir, filename)
    
    if (file_format == "png") {
      ggsave(filepath, plot = p, width = width, height = height, dpi = dpi)
    } else if (file_format == "pdf") {
      ggsave(filepath, plot = p, width = width, height = height)
    } else {
      warning("Unsupported file_format, using PNG format")
      ggsave(paste0(filepath, ".png"), plot = p, width = width, height = height, dpi = dpi)
    }
    
    message(paste("Plot saved to:", filepath))
  }
  
  # Print data summary
  cat("\n=== Data Summary (", rmsk_class, "-", type, ") ===\n")
  print(grouped)
  if (!is.null(text_str)) {
    cat("\n=== Statistics ===\n")
    cat(text_str, "\n")
    # Also print the line color used
    cat("Line color (MCR change direction):", line_color, "\n")
  }
  cat("\n")
  
  return(p)
}

# Function to create combined plot for all Rmsk_class for a specific type
create_combined_plot_by_type <- function(data, 
                                         type = "SNV",
                                         output_dir = "./combined_plots",
                                         file_format = "pdf",
                                         width = 16,
                                         height = 12) {
  
  # Get all unique Rmsk_class values
  rmsk_classes <- unique(data$Rmsk_class)
  
  cat("Creating combined plot for type:", type, "\n")
  cat("Number of Rmsk_class:", length(rmsk_classes), "\n")
  cat("Rmsk_class:", paste(rmsk_classes, collapse = ", "), "\n\n")
  
  # Create list to store plots
  plot_list <- list()
  
  # Generate individual plots for each Rmsk_class
  for (i in seq_along(rmsk_classes)) {
    class <- rmsk_classes[i]
    cat(sprintf("[%d/%d] Processing: %s\n", i, length(rmsk_classes), class))
    
    # Create individual plot
    p <- plot_mcr_by_rmsk_class(
      data = data,
      rmsk_class = class,
      type = type,
      title = paste(class, "-", type),
      show_stats = TRUE,
      save_plot = FALSE
    )
    
    if (!is.null(p)) {
      # Remove legend from all but the first plot for cleaner layout
      if (i > 1) {
        p <- p + theme(legend.position = "none")
      }
      plot_list[[class]] <- p
      cat("  ✓ Completed\n")
    } else {
      cat("  ⓘ No data, skipping\n")
    }
  }
  
  # Check if we have any plots
  if (length(plot_list) == 0) {
    warning(paste("No valid plots generated for type", type))
    return(NULL)
  }
  
  # Combine plots using patchwork
  # Determine layout (approximately square)
  n_plots <- length(plot_list)
  n_cols <- ceiling(sqrt(n_plots))
  n_rows <- ceiling(n_plots / n_cols)
  
  cat("\nCombining", n_plots, "plots into", n_rows, "rows ×", n_cols, "columns\n")
  
  # Create combined plot
  combined_plot <- wrap_plots(plot_list, ncol = n_cols, nrow = n_rows) +
    plot_annotation(
      title = paste(type, "Analysis: MCR and Variant Counts by Rmsk_class and Quality"),
      subtitle = "Line color indicates MCR change direction: Green = Increase (Q40 > Q30), Red = Decrease (Q40 < Q30)",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40")
      )
    )
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save combined plot
  filename <- paste0("Combined_", type, "_Analysis.", file_format)
  filepath <- file.path(output_dir, filename)
  
  if (file_format == "pdf") {
    ggsave(filepath, plot = combined_plot, width = width, height = height)
  } else {
    ggsave(filepath, plot = combined_plot, width = width, height = height, dpi = 300)
  }
  
  cat("Combined plot saved to:", filepath, "\n\n")
  
  return(combined_plot)
}

# Function to create summary plot with color-coded MCR changes
create_summary_plot <- function(data) {
  # Prepare data with MCR change information
  summary_data <- data %>%
    # First, get mean values for each Rmsk_class, type, and Quality
    dplyr::group_by(Rmsk_class, type, Quality) %>%
    dplyr::summarise(
      Mean_MCR = mean(MCR, na.rm = TRUE),
      Mean_Total_Variants = mean(Total_Variants, na.rm = TRUE),
      Mean_Concordant_Variants = mean(Mendelian_Concordant_Variants, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Reshape to calculate Q40 vs Q30 changes
    pivot_wider(
      id_cols = c(Rmsk_class, type),
      names_from = Quality,
      values_from = c(Mean_MCR, Mean_Total_Variants, Mean_Concordant_Variants)
    ) %>%
    # Calculate changes
    mutate(
      MCR_change = Mean_MCR_Q40 - Mean_MCR_Q30,
      Change_direction = case_when(
        MCR_change > 0 ~ "Increase",
        MCR_change < 0 ~ "Decrease",
        TRUE ~ "No change"
      ),
      Change_color = case_when(
        MCR_change > 0 ~ "#2E8B57",  # Green for increase
        MCR_change < 0 ~ "#DC143C",  # Red for decrease
        TRUE ~ "darkgray"            # Gray for no change
      )
    )
  
  # Create summary plot with color-coded bars
  p <- ggplot(summary_data, aes(x = Rmsk_class, y = Mean_MCR_Q40, fill = Change_direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_point(aes(y = Mean_MCR_Q30), size = 3, color = "black") +
    geom_segment(aes(x = Rmsk_class, xend = Rmsk_class,
                     y = Mean_MCR_Q30, yend = Mean_MCR_Q40),
                 size = 1, color = "black") +
    facet_wrap(~type, scales = "free_x") +
    geom_text(aes(label = sprintf("Q40: %.3f\nQ30: %.3f\nΔ: %+.3f", 
                                  Mean_MCR_Q40, Mean_MCR_Q30, MCR_change)),
              vjust = -0.5, size = 2.8, lineheight = 0.9) +
    labs(title = "Summary: MCR Comparison between Q30 and Q40",
         subtitle = "Bar color indicates MCR change direction: Green = Increase, Red = Decrease",
         x = "Rmsk_class",
         y = "MCR (Q40)",
         fill = "MCR Change Direction") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    ) +
    scale_fill_manual(values = c("Increase" = "#2E8B57", 
                                 "Decrease" = "#DC143C",
                                 "No change" = "darkgray"))
  
  return(p)
}

# Main execution function
main <- function() {
  cat("==========================================\n")
  cat("MCR Data Visualization Script\n")
  cat("==========================================\n\n")
  
  # 1. Check data structure
  cat("1. Checking data structure...\n")
  print(str(df_new))
  cat("\nData preview:\n")
  print(head(df_new))
  
  # 2. Check unique values
  cat("\n2. Unique value statistics:\n")
  cat("Rmsk_class:", paste(unique(df_new$Rmsk_class), collapse = ", "), "\n")
  cat("type:", paste(unique(df_new$type), collapse = ", "), "\n")
  cat("Quality:", paste(unique(df_new$Quality), collapse = ", "), "\n")
  
  # 3. Create output directory
  output_dir <- "./mcr_visualizations"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 4. Create combined plots for SNV and InDel
  cat("\n3. Creating combined plots...\n")
  
  # Combined plot for SNV
  cat("\nCreating combined plot for SNV...\n")
  combined_snv <- create_combined_plot_by_type(
    data = df_new,
    type = "SNV",
    output_dir = output_dir,
    file_format = "pdf",
    width = 16,
    height = 12
  )
  
  # Combined plot for InDel
  cat("\nCreating combined plot for InDel...\n")
  combined_indel <- create_combined_plot_by_type(
    data = df_new,
    type = "InDel",
    output_dir = output_dir,
    file_format = "pdf",
    width = 16,
    height = 12
  )
  
  # Save individual plots as PDF if needed
  cat("\n4. Saving individual plots as PDF...\n")
  
  # Get all Rmsk_class and type combinations
  rmsk_classes <- unique(df_new$Rmsk_class)
  types <- unique(df_new$type)
  
  # Create directory for individual plots
  individual_dir <- file.path(output_dir, "individual_plots")
  if (!dir.exists(individual_dir)) {
    dir.create(individual_dir, recursive = TRUE)
  }
  
  for (class in rmsk_classes) {
    for (t in types) {
      cat(sprintf("Processing: %s - %s\n", class, t))
      
      tryCatch({
        p <- plot_mcr_by_rmsk_class(
          data = df_new,
          rmsk_class = class,
          type = t,
          title = paste(class, "Region:", t, "- Variant Counts and MCR by Quality"),
          show_stats = TRUE,
          save_plot = TRUE,
          output_dir = individual_dir,
          file_format = "pdf",
          width = 10,
          height = 7
        )
        
        if (!is.null(p)) {
          cat("  ✓ Saved\n")
        } else {
          cat("  ⓘ No data, skipping\n")
        }
      }, error = function(e) {
        cat(sprintf("  ✗ Error: %s\n", e$message))
      })
    }
  }
  
  # 5. Create enhanced summary plot
  cat("\n5. Creating enhanced summary plot...\n")
  summary_plot <- create_summary_plot(df_new)
  
  # Save summary plot as PDF
  summary_path <- file.path(output_dir, "summary_plot.pdf")
  ggsave(summary_path, 
         plot = summary_plot, 
         width = 14, 
         height = 8)
  cat("Summary plot saved to:", summary_path, "\n")
  
  # Display a preview of the summary plot
  print(summary_plot)
  
  # 6. Generate detailed report
  cat("\n6. Generating detailed analysis report...\n")
  report_file <- file.path(output_dir, "analysis_report.txt")
  
  sink(report_file)
  cat("MCR Visualization Analysis Report\n")
  cat("Generated at:", Sys.time(), "\n")
  cat("==========================================\n\n")
  
  cat("DATA OVERVIEW:\n")
  cat("Total rows:", nrow(df_new), "\n")
  cat("Number of Rmsk_class:", length(unique(df_new$Rmsk_class)), "\n")
  cat("Number of types:", length(unique(df_new$type)), "\n")
  cat("Quality types:", paste(unique(df_new$Quality), collapse = ", "), "\n\n")
  
  cat("COLOR CODING SCHEME:\n")
  cat("• Green line (#2E8B57): MCR increased from Q30 to Q40\n")
  cat("• Red line (#DC143C): MCR decreased from Q30 to Q40\n")
  cat("• Dark green line (default): MCR unchanged or comparison not available\n\n")
  
  cat("STATISTICS BY RMSK_CLASS:\n")
  for (class in unique(df_new$Rmsk_class)) {
    cat("\n", toupper(class), ":\n")
    cat("=", rep("=", nchar(class)+2), "\n", sep="")
    
    class_data <- df_new[df_new$Rmsk_class == class, ]
    cat("  Sample count:", nrow(class_data), "\n")
    
    for (t in unique(class_data$type)) {
      type_data <- class_data[class_data$type == t, ]
      cat("  ", t, ":\n")
      cat("    • Quality distribution:", paste(table(type_data$Quality), collapse = ", "), "\n")
      
      # Calculate MCR statistics if Q30 and Q40 both exist
      q30_mcr <- type_data[type_data$Quality == "Q30", ]$MCR
      q40_mcr <- type_data[type_data$Quality == "Q40", ]$MCR
      
      if (length(q30_mcr) > 0 && length(q40_mcr) > 0) {
        avg_q30 <- mean(q30_mcr, na.rm = TRUE)
        avg_q40 <- mean(q40_mcr, na.rm = TRUE)
        mcr_change <- avg_q40 - avg_q30
        
        cat("    • MCR Q30:", sprintf("%.3f", avg_q30), "\n")
        cat("    • MCR Q40:", sprintf("%.3f", avg_q40), "\n")
        cat("    • MCR change:", sprintf("%+.3f", mcr_change), 
            ifelse(mcr_change > 0, "(INCREASE)", 
                   ifelse(mcr_change < 0, "(DECREASE)", "(NO CHANGE)")), "\n")
      } else {
        cat("    • MCR range:", sprintf("%.3f", min(type_data$MCR, na.rm = TRUE)), 
            "-", sprintf("%.3f", max(type_data$MCR, na.rm = TRUE)), "\n")
      }
    }
  }
  
  sink()
  
  cat("Report saved to:", report_file, "\n")
  cat("\n==========================================\n")
  cat("Script execution completed!\n")
  cat("Output directory:", normalizePath(output_dir), "\n")
  cat("Files generated:\n")
  cat("  • Combined_SNV_Analysis.pdf\n")
  cat("  • Combined_InDel_Analysis.pdf\n")
  cat("  • summary_plot.pdf\n")
  cat("  • individual_plots/ (directory with individual plots)\n")
  cat("  • analysis_report.txt\n")
  cat("==========================================\n")
  
  # Return list of combined plots for potential further use
  return(list(
    combined_snv = combined_snv,
    combined_indel = combined_indel,
    summary_plot = summary_plot
  ))
}

# Execute main function
results <- main()
