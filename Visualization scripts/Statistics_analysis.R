library(dplyr)
library(tidyr)
#Statistics analysis for All Data Quality---------

calculate_unpaired_cohens_d <- function(q30_values, q40_values) {
  # 计算均值和标准差
  mean_q30 <- mean(q30_values, na.rm = TRUE)
  mean_q40 <- mean(q40_values, na.rm = TRUE)
  
  sd_q30 <- sd(q30_values, na.rm = TRUE)
  sd_q40 <- sd(q40_values, na.rm = TRUE)
  
  n1 <- length(q30_values)
  n2 <- length(q40_values)
  
  # 计算合并标准差（pooled standard deviation）
  pooled_sd <- sqrt(((n1 - 1) * sd_q30^2 + (n2 - 1) * sd_q40^2) / (n1 + n2 - 2))
  
  # Cohen's d 计算公式
  cohens_d <- (mean_q40-mean_q30) / pooled_sd
  
  # 计算置信区间（使用t分布）
  if (n1 >= 2 && n2 >= 2 && pooled_sd > 0) {
    # 标准误差
    se_d <- sqrt(1/n1 + 1/n2)
    
    # t临界值（95% CI）
    t_critical <- qt(0.975, df = n1 + n2 - 2)
    
    # 置信区间
    ci_lower <- cohens_d - t_critical * se_d
    ci_upper <- cohens_d + t_critical * se_d
  } else {
    ci_lower <- NA
    ci_upper <- NA
  }
  
  result_df <- data.frame(
    cohens_d = cohens_d,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    mean_q30 = mean_q30,
    mean_q40 = mean_q40,
    sd_q30 = sd_q30,
    sd_q40 = sd_q40,
    n1 = n1,
    n2 = n2
  )
  
  return(result_df)
}


get_test_type_QC <- function(df) {
  # 定义内部函数：计算配对Cohen's d效应量
  calculate_paired_cohens_d <- function(q30_values, q40_values) {
    # 计算配对差值
    differences <- q40_values - q30_values
    
    # 计算均值差和标准差
    mean_diff <- mean(differences, na.rm = TRUE)
    sd_diff <- sd(differences, na.rm = TRUE)
    
    # Cohen's d（配对版本）
    if (!is.na(sd_diff) && sd_diff > 0) {
      cohens_d <- mean_diff / sd_diff
    } else {
      cohens_d <- NA
    }
    
    # 计算置信区间（使用t分布）
    n <- length(differences)
    if (n >= 2 && !is.na(sd_diff) && sd_diff > 0 && !is.na(cohens_d)) {
      # 标准误差
      se_d <- sqrt((1/n) + (cohens_d^2)/(2*n))
      
      # t临界值（95% CI）
      t_critical <- qt(0.975, df = n-1)
      
      ci_lower <- cohens_d - t_critical * se_d
      ci_upper <- cohens_d + t_critical * se_d
    } else {
      ci_lower <- NA
      ci_upper <- NA
    }
    
    return(list(
      cohens_d = cohens_d,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      mean_diff = mean_diff,
      sd_diff = sd_diff,
      n_pairs = n
    ))
  }
  
  metrics <- unique(df$Metric)
  
  # Create an empty list to store results
  results_list <- list()
  
  # Loop through each Metric
  for (metric in metrics) {
    cat("\n=== Processing indicator:", metric, "===\n")
    
    # Filter data for the current Metric
    metric_data <- df %>% filter(Metric == metric)
    
    # Group by Quality, prepare paired data
    # Assume each sample has paired data under Q30 and Q40
    metric_data_wide <- metric_data %>%
      dplyr::select(Sample, Quality, Value, source) %>%
      dplyr::mutate(pair_id = gsub("_R[12].*|_ELE|_ILM|Q30_|Q40_|Fresh_", "", source)) %>% 
      dplyr::group_by(pair_id, Quality) %>%
      dplyr::summarise(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop') %>%
      pivot_wider(names_from = Quality, values_from = mean_value)
    
    if (!("Q30" %in% names(metric_data_wide) & "Q40" %in% names(metric_data_wide))) {
      cat("  Warning: Current indicator lacks Q30 or Q40 data\n")
      next
    }
    
    # Remove NA values
    metric_data_wide <- metric_data_wide %>%
      filter(!is.na(Q30) & !is.na(Q40))
    
    if (nrow(metric_data_wide) < 2) {
      cat("  Warning: Insufficient valid paired samples\n")
      next
    }
    
    # 计算Q30和Q40各自的均值和标准差
    mean_Q30 <- mean(metric_data_wide$Q30, na.rm = TRUE)
    sd_Q30 <- sd(metric_data_wide$Q30, na.rm = TRUE)
    mean_Q40 <- mean(metric_data_wide$Q40, na.rm = TRUE)
    sd_Q40 <- sd(metric_data_wide$Q40, na.rm = TRUE)
    
    cat("  Q30: Mean =", round(mean_Q30, 3), "SD =", round(sd_Q30, 3), "\n")
    cat("  Q40: Mean =", round(mean_Q40, 3), "SD =", round(sd_Q40, 3), "\n")
    
    # 计算配对Cohen's d效应量
    effect_size <- calculate_paired_cohens_d(metric_data_wide$Q30, metric_data_wide$Q40)
    
    # 计算差值（用于正态性检验）
    differences <- metric_data_wide$Q40 - metric_data_wide$Q30
    
    # 初始化变量
    shapiro_p <- NA
    shapiro_note <- "Not tested"
    shapiro_result <- "Not tested"
    
    # 1. 正态性检验
    if (length(differences) < 3) {
      cat("  Sample size too small for Shapiro-Wilk test (n < 3)\n")
      shapiro_p <- NA
      shapiro_note <- "Not tested (n<3)"
      shapiro_result <- "Not tested (n<3)"
      
      # 样本量太小，直接使用Wilcoxon检验
      cat("  Using Wilcoxon signed-rank test\n")
      test_result <- tryCatch({
        wilcox.test(metric_data_wide$Q40, metric_data_wide$Q30,
                    paired = TRUE, alternative = "two.sided", 
                    exact = FALSE, conf.int = TRUE)
      }, error = function(e) {
        cat("  Error in Wilcoxon test:", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(test_result)) next
      
      test_method <- "Wilcoxon signed-rank test"
      test_statistic <- ifelse(!is.null(test_result$statistic), test_result$statistic, NA)
      df_value <- NA
      p_value <- test_result$p.value
      
    } else {
      # 检查差值是否全部相同
      if (length(unique(differences)) == 1) {
        cat("  All differences are identical, skipping Shapiro-Wilk test\n")
        shapiro_p <- NA
        shapiro_note <- "Not tested (all diffs identical)"
        shapiro_result <- "Not tested (all diffs identical)"
      } else {
        # Shapiro-Wilk正态性检验
        shapiro_test <- tryCatch({
          shapiro.test(differences)
        }, error = function(e) {
          cat("  Error in Shapiro-Wilk test:", e$message, "\n")
          return(NULL)
        })
        
        if (is.null(shapiro_test)) {
          shapiro_p <- NA
          shapiro_note <- "Test failed"
          shapiro_result <- "Test failed"
        } else {
          shapiro_p <- shapiro_test$p.value
          shapiro_note <- ifelse(shapiro_p > 0.05, "Normal", "Non-normal")
          shapiro_result <- shapiro_note
          cat("  Shapiro-Wilk normality test: W =", round(shapiro_test$statistic, 4), 
              ", p =", format.pval(shapiro_test$p.value, digits = 3), "\n")
        }
      }
      
      # 2. 根据正态性检验结果选择检验方法
      if (is.na(shapiro_p) || shapiro_p > 0.01) {
        cat("  Using paired t-test\n")
        test_result <- tryCatch({
          t.test(metric_data_wide$Q40, metric_data_wide$Q30,
                 paired = TRUE, alternative = "two.sided")
        }, error = function(e) {
          cat("  Error in t-test:", e$message, "\n")
          return(NULL)
        })
        
        if (is.null(test_result)) next
        
        test_method <- "Paired t-test"
        test_statistic <- test_result$statistic
        df_value <- test_result$parameter
        p_value <- test_result$p.value
        
      } else {
        cat("  Using Wilcoxon signed-rank test\n")
        test_result <- tryCatch({
          wilcox.test(metric_data_wide$Q40, metric_data_wide$Q30,
                      paired = TRUE, alternative = "two.sided", 
                      exact = FALSE, conf.int = TRUE)
        }, error = function(e) {
          cat("  Error in Wilcoxon test:", e$message, "\n")
          return(NULL)
        })
        
        if (is.null(test_result)) next
        
        test_method <- "Wilcoxon signed-rank test"
        test_statistic <- ifelse(!is.null(test_result$statistic), test_result$statistic, NA)
        df_value <- NA
        p_value <- test_result$p.value
      }
    }
    
    # 确定效应量大小分类
    if (!is.na(effect_size$cohens_d)) {
      effect_size_magnitude <- case_when(
        abs(effect_size$cohens_d) < 0.2 ~ "Negligible",
        abs(effect_size$cohens_d) < 0.5 ~ "Small",
        abs(effect_size$cohens_d) < 0.8 ~ "Medium",
        TRUE ~ "Large"
      )
    } else {
      effect_size_magnitude <- "Not calculated"
    }
    
    # 确定显著性标记
    significance <- case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
    
    # 格式化输出字符串
    mean_Q30_str <- sprintf("%.3f ± %.3f", mean_Q30, sd_Q30)
    mean_Q40_str <- sprintf("%.3f ± %.3f", mean_Q40, sd_Q40)
    mean_diff_CI <- sprintf("%.3f ± %.3f", effect_size$mean_diff, effect_size$sd_diff)
    
    if (!is.na(effect_size$ci_lower) && !is.na(effect_size$ci_upper)) {
      cohens_d_CI <- sprintf("%.3f [%.3f, %.3f]", 
                             effect_size$cohens_d, effect_size$ci_lower, effect_size$ci_upper)
    } else {
      cohens_d_CI <- sprintf("%.3f [NA, NA]", effect_size$cohens_d)
    }
    
    # 存储结果 - 按照你提供的参考格式
    results_list[[metric]] <- data.frame(
      Test = test_method,
      Shapiro_p = shapiro_p,
      Shapiro_note = shapiro_note,
      Test_statistic = test_statistic,
      df = df_value,
      p_value = p_value,
      Mean_Q30 = mean_Q30,
      SD_Q30 = sd_Q30,
      Mean_Q40 = mean_Q40,
      SD_Q40 = sd_Q40,
      Mean_diff = effect_size$mean_diff,
      SD_diff = effect_size$sd_diff,
      Cohens_d = effect_size$cohens_d,
      CI_lower_d = effect_size$ci_lower,
      CI_upper_d = effect_size$ci_upper,
      n_pairs = effect_size$n_pairs,
      Significance = significance,
      Shapiro_result = shapiro_result,
      Effect_size_magnitude = effect_size_magnitude,
      Mean_Q30_str = mean_Q30_str,
      Mean_Q40_str = mean_Q40_str,
      Mean_diff_CI = mean_diff_CI,
      Cohens_d_CI = cohens_d_CI,
      stringsAsFactors = FALSE
    )
    
    # 打印简要结果
    cat("  Test result: p =", format.pval(p_value, digits = 3), 
        " (", significance, ")\n", sep = "")
    if (!is.na(effect_size$cohens_d)) {
      cat("  Cohen's d:", round(effect_size$cohens_d, 3), 
          "[", round(effect_size$ci_lower, 3), ",", 
          round(effect_size$ci_upper, 3), "] (", effect_size_magnitude, ")\n")
    }
    cat("  Number of sample pairs:", effect_size$n_pairs, "\n")
  }
  
  # Combine all results into one dataframe
  if (length(results_list) == 0) {
    cat("\nNo results to return.\n")
    return(NULL)
  }
  
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  # 添加Metric列作为第一列
  results_df <- cbind(Metric = names(results_list), results_df)
  
  # 重新排列列顺序，将字符串列放在后面
  results_df <- results_df %>%
    select(Metric, Test, Shapiro_p, Shapiro_note, Test_statistic, df, p_value,
           Mean_Q30, SD_Q30, Mean_Q40, SD_Q40, Mean_diff, SD_diff, Cohens_d,
           CI_lower_d, CI_upper_d, n_pairs, Significance, Shapiro_result,
           Effect_size_magnitude, Mean_Q30_str, Mean_Q40_str, Mean_diff_CI, Cohens_d_CI)
  
  # Display results
  cat("\n\n=== Complete Results ===\n")
  print(results_df)
  
  # 打印简洁摘要
  cat("\n\n=== Summary Table ===\n")
  print(results_df[, c("Metric", "Test", "Shapiro_result", "p_value", 
                       "Significance", "Effect_size_magnitude", 
                       "Mean_Q30_str", "Mean_Q40_str", "Mean_diff_CI", 
                       "Cohens_d_CI", "n_pairs")])
  
  return(results_df)
}

#Statistics analysis for WES----------
prepare_difference_data <- function(df) {
  df$source_clean <- gsub('_ILM|_ELE', '', df$source)
  
  diff_data <- df %>%
    dplyr::group_by(Depth, source_clean) %>%
    dplyr::filter(all(c("Q30", "Q40") %in% Quality)) %>%
    dplyr::summarise(
      Q30_value = F1.score[Quality == "Q30"],
      Q40_value = F1.score[Quality == "Q40"],
      difference =  Q40_value-Q30_value,
      .groups = 'drop'
    )
  
  return(diff_data)
}


plot_all_qqplots <- function(df, title = "Q-Q Plots of Differences by Depth") {

  diff_data <- prepare_difference_data(df)

  depths <- unique(diff_data$Depth)

  depth_order <- c("10X", "15X", "20X", "25X", "30X", "60X", "90X", "120X")
  depths <- depths[order(match(depths, depth_order))]

  plot_list <- list()
  
  for (depth_i in depths) {
    print(paste0('Processing....', depth_i))
    
    depth_diff <- diff_data %>%
      filter(Depth == depth_i) %>%
      pull(difference)
    
    if (length(depth_diff) < 2) {
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = paste("Insufficient data\nfor", depth_i), 
                 size = 4) +
        theme_void()
      plot_list[[depth_i]] <- ggplotGrob(p)
    } else {

      qq_data <- as.data.frame(qqnorm(depth_diff, plot.it = FALSE))

      theoretical <- qq_data$x
      actual <- qq_data$y

      qq_fit <- lm(actual ~ theoretical)
      intercept <- coef(qq_fit)[1]
      slope <- coef(qq_fit)[2]
      
      shapiro_test <- shapiro.test(depth_diff)
      shapiro_p <- shapiro_test$p.value

      is_normal <- shapiro_p > 0.01

      p <- ggplot(qq_data, aes(x = theoretical, y = actual)) +
        geom_point(size = 3, shape = 21, 
                   fill = ifelse(is_normal, "#4DAF4A", "#E41A1C"),
                   color = "black", alpha = 0.7) +
        geom_abline(intercept = intercept, slope = slope, 
                    color = "#377EB8", size = 0.8, linetype = "solid") +
        geom_abline(intercept = 0, slope = 1, 
                    color = "#999999", size = 0.5, linetype = "dashed") +
        labs(
          title = paste("Depth:", depth_i),
          subtitle = paste("Shapiro-Wilk p =", format.pval(shapiro_p, digits = 2),
                           ifelse(is_normal, "(Normal)", "(Non-normal)")),
          x = "Theoretical Quantiles",
          y = "Sample Quantiles"
        ) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "gray70", size = 0.5),
          plot.margin = margin(5, 5, 5, 5)
        )

      p <- p + annotate("text", 
                        x = min(theoretical), 
                        y = max(actual),
                        label = paste("n =", length(depth_diff)),
                        hjust = 0, vjust = 1, size = 3)

      plot_list[[depth_i]] <- ggplotGrob(p)
    }
  }
  
  combined_plot <- wrap_plots(plot_list, ncol = 4, nrow = 2) +
    plot_annotation(
      title = title,
      subtitle = "Comparison of Q40-Q30 F1.score differences across sequencing depths",
      caption = "Red points: non-normal (p ≤ 0.01); Green points: normal (p > 0.01)\nDashed line: y = x (perfect normality); Solid line: data trend",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        plot.caption = element_text(size = 9, hjust = 0)
      )
    )
  
  return(combined_plot)
}

# Method 1: manually compute paired Cohen's d
calculate_paired_cohens_d <- function(q30_values, q40_values) {
  # compute paired differences
  differences <-  q40_values-q30_values
  
  # compute mean difference and sd
  mean_diff <- mean(differences, na.rm = TRUE)
  sd_diff <- sd(differences, na.rm = TRUE)
  
  # Cohen's d (paired version)
  cohens_d <- mean_diff / sd_diff
  
  # compute CI (t-distribution)
  n <- length(differences)
  if (n >= 2) {
    # standard error
    se_d <- sqrt((1/n) + (cohens_d^2)/(2*n))
    
    # t critical value (95% CI)
    t_critical <- qt(0.975, df = n-1)
    
    ci_lower <- cohens_d - t_critical * se_d
    ci_upper <- cohens_d + t_critical * se_d
  } else {
    ci_lower <- NA
    ci_upper <- NA
  }
  
  return(list(
    cohens_d = cohens_d,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    mean_diff = mean_diff,
    sd_diff = sd_diff,
    n_pairs = n
  ))
}


get_test_type_with_effectsize <- function(df) {
  df$source <- gsub('_ELE|_ILM|ELE_|ILM_','',df$source)
  depths <- unique(df$Depth)
  
  results_list <- list()
  
  for (i in seq_along(depths)) {
    depth_i <- depths[i]
    
    cat("\n=== Processing Depth:", depth_i, "===\n")
    
    # filter data for current Depth
    combo_data <- df %>% filter(Depth == depth_i)
    
    # prepare paired data
    samples_with_both <- combo_data %>%
      dplyr::group_by(source) %>%
      dplyr::filter(all(c("Q30", "Q40") %in% Quality)) %>%
      ungroup() %>%
      distinct(source)
    
    if (nrow(samples_with_both) < 2) {
      cat("  Warning: Insufficient paired samples\n")
      next
    }
    
    # check whether each source has equal numbers of Q30 and Q40
    sample_counts <- combo_data %>%
      dplyr::filter(source %in% samples_with_both$source) %>%
      dplyr::filter(Quality %in% c("Q30", "Q40")) %>%
      dplyr::group_by(source, Quality) %>%
      dplyr::summarise(count = n(), .groups = 'drop') %>%
      pivot_wider(
        names_from = Quality, 
        values_from = count,
        values_fill = 0
      )
    
    # identify sources with unequal Q30/Q40 counts
    unequal_samples <- sample_counts %>%
      filter(Q30 != Q40) %>%
      pull(source)
    
    if (length(unequal_samples) > 0) {
      cat("  Warning: Unequal sample sizes detected for sources:", 
          paste(unequal_samples, collapse = ", "), "\n")
      cat("  Skipping paired test due to unequal sample sizes\n")
      next
    }
    
    metric_data_wide <- combo_data %>%
      filter(source %in% samples_with_both$source) %>%
      select(source, Quality, F1.score) %>%
      pivot_wider(
        names_from = Quality, 
        values_from = F1.score,
        values_fn = mean
      ) %>%
      filter(!is.na(Q30) & !is.na(Q40))
    
    if (nrow(metric_data_wide) < 2) {
      cat("  Warning: Insufficient valid paired samples\n")
      next
    }
    
    # compute mean and sd for Q30 and Q40 separately
    mean_Q30 <- mean(metric_data_wide$Q30, na.rm = TRUE)
    sd_Q30 <- sd(metric_data_wide$Q30, na.rm = TRUE)
    mean_Q40 <- mean(metric_data_wide$Q40, na.rm = TRUE)
    sd_Q40 <- sd(metric_data_wide$Q40, na.rm = TRUE)
    
    # compute differences
    differences <- metric_data_wide$Q40 - metric_data_wide$Q30
    
    # === key fix: check whether differences are constant ===
    if (length(unique(differences)) == 1) {
      cat("  Warning: All paired differences are identical (value =", 
          unique(differences), ")\n")
      cat("  Skipping Shapiro-Wilk test, proceeding directly with Wilcoxon test\n")
      
      # for constant differences normality test is meaningless, use Wilcoxon directly
      # note: Wilcoxon may still warn but will run
      test_result <- suppressWarnings(
        wilcox.test(metric_data_wide$Q40, metric_data_wide$Q30,
                    paired = TRUE, alternative = "two.sided",
                    exact = FALSE, conf.int = TRUE)
      )
      
      test_method <- "Wilcoxon signed-rank test"
      test_statistic <- test_result$statistic
      df_value <- NA
      shapiro_p <- NA
      shapiro_note <- "Not tested (constant differences)"
      
    } else if (nrow(metric_data_wide) < 3) {
      # when sample size < 3, Shapiro-Wilk cannot be performed
      cat("  Warning: Sample size too small for Shapiro-Wilk test (n < 3)\n")
      cat("  Proceeding with Wilcoxon test directly\n")
      
      test_result <- wilcox.test(metric_data_wide$Q40, metric_data_wide$Q30,
                                 paired = TRUE, alternative = "two.sided",
                                 exact = FALSE, conf.int = TRUE)
      
      test_method <- "Wilcoxon signed-rank test"
      test_statistic <- test_result$statistic
      df_value <- NA
      shapiro_p <- NA
      shapiro_note <- "Not tested (n < 3)"
      
    } else {
      # normal case: perform Shapiro-Wilk test
      shapiro_test <- shapiro.test(differences)
      shapiro_p <- shapiro_test$p.value
      shapiro_note <- ifelse(shapiro_p > 0.05, "Normal", "Non-normal")
      
      # choose test method
      if (shapiro_test$p.value > 0.05) {
        cat("  Using paired t-test\n")
        test_result <- t.test(metric_data_wide$Q40, metric_data_wide$Q30, 
                              paired = TRUE, alternative = "two.sided")
        
        test_method <- "Paired t-test"
        test_statistic <- test_result$statistic
        df_value <- test_result$parameter
        
      } else {
        cat("  Using Wilcoxon signed-rank test\n")
        test_result <- wilcox.test(metric_data_wide$Q40, metric_data_wide$Q30, 
                                   paired = TRUE, alternative = "two.sided",
                                   exact = FALSE, conf.int = TRUE)
        
        test_method <- "Wilcoxon signed-rank test"
        test_statistic <- test_result$statistic
        df_value <- NA
      }
    }
    
    # compute effect size
    # make sure function calculate_paired_cohens_d is defined
    effect_size <- calculate_paired_cohens_d(
      metric_data_wide$Q30, 
      metric_data_wide$Q40
    )
    
    # store results
    results_list[[i]] <- data.frame(
      Depth = depth_i,
      Test = test_method,
      Shapiro_p = shapiro_p,
      Shapiro_note = if(exists("shapiro_note")) shapiro_note else NA,
      Test_statistic = test_statistic,
      df = df_value,
      p_value = test_result$p.value,
      # add descriptive stats for Q30 and Q40
      Mean_Q30 = mean_Q30,
      SD_Q30 = sd_Q30,
      Mean_Q40 = mean_Q40,
      SD_Q40 = sd_Q40,
      # descriptive stats for differences
      Mean_diff = effect_size$mean_diff,
      SD_diff = effect_size$sd_diff,
      Cohens_d = effect_size$cohens_d,
      CI_lower_d = effect_size$ci_lower,
      CI_upper_d = effect_size$ci_upper,
      n_pairs = effect_size$n_pairs,
      stringsAsFactors = FALSE
    )
    
    # print results
    cat("  Q30: Mean =", round(mean_Q30, 3), "SD =", round(sd_Q30, 3), "\n")
    cat("  Q40: Mean =", round(mean_Q40, 3), "SD =", round(sd_Q40, 3), "\n")
    cat("  p-value:", format.pval(test_result$p.value, digits = 3), "\n")
    if (!is.null(effect_size$cohens_d) && !is.na(effect_size$cohens_d)) {
      cat("  Cohen's d:", round(effect_size$cohens_d, 3), 
          "[", round(effect_size$ci_lower, 3), ",", 
          round(effect_size$ci_upper, 3), "]\n")
    }
  }
  
  # combine all results
  if (length(results_list) == 0) return(NULL)
  
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  # add significance flags and interpretation
  results_summary <- results_df %>%
    mutate(
      Significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      Shapiro_result = case_when(
        is.na(Shapiro_p) ~ "Not tested (n<3)",
        Shapiro_p > 0.05 ~ "Normal",
        TRUE ~ "Non-normal"
      ),
      Effect_size_magnitude = case_when(
        abs(Cohens_d) < 0.2 ~ "Negligible",
        abs(Cohens_d) < 0.5 ~ "Small",
        abs(Cohens_d) < 0.8 ~ "Medium",
        TRUE ~ "Large"
      ),
      # format output
      Mean_Q30_str = sprintf("%.3f ± %.3f", Mean_Q30, SD_Q30),
      Mean_Q40_str = sprintf("%.3f ± %.3f", Mean_Q40, SD_Q40),
      Mean_diff_CI = sprintf("%.3f [NA, NA]", Mean_diff),
      Cohens_d_CI = sprintf("%.3f [%.3f, %.3f]", 
                            Cohens_d, CI_lower_d, CI_upper_d)
    )
  
  return(results_summary)
}


get_test_type_with_effectsize_RNAseq <- function(df) {
  # remove platform info from source names
  df$source <- gsub('ELE|ILM','', rownames(df))
  
  cat("\n=== Processing All Data (No Depth Classification) ===\n")
  
  # prepare paired data - ensure each sample has both Q30 and Q40
  samples_with_both <- df %>%
    dplyr::group_by(source) %>%
    dplyr::filter(all(c("Q30", "Q40") %in% Quality)) %>%
    ungroup() %>%
    distinct(source)
  
  if (nrow(samples_with_both) < 2) {
    cat("Warning: Insufficient paired samples\n")
    return(NULL)
  }
  
  # check whether each source has equal numbers of Q30 and Q40
  sample_counts <- df %>%
    dplyr::filter(source %in% samples_with_both$source) %>%
    dplyr::filter(Quality %in% c("Q30", "Q40")) %>%
    dplyr::group_by(source, Quality) %>%
    dplyr::summarise(count = n(), .groups = 'drop') %>%
    pivot_wider(
      names_from = Quality, 
      values_from = count,
      values_fill = 0
    )
  
  # identify sources with unequal Q30/Q40 counts
  unequal_samples <- sample_counts %>%
    filter(Q30 != Q40) %>%
    pull(source)
  
  if (length(unequal_samples) > 0) {
    cat("Warning: Unequal sample sizes detected for sources:", 
        paste(unequal_samples, collapse = ", "), "\n")
    cat("Skipping paired test due to unequal sample sizes\n")
    return(NULL)
  }
  
  # create wide-format data - one row per source with Q30 and Q40 SNR values
  # note: if a source has multiple replicates, take the mean
  metric_data_wide <- df %>%
    dplyr::filter(source %in% samples_with_both$source) %>%
    dplyr::select(source, Quality, SNR) %>%
    dplyr::group_by(source, Quality) %>%
    dplyr::summarise(SNR = mean(SNR, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(
      names_from = Quality, 
      values_from = SNR
    ) %>%
    filter(!is.na(Q30) & !is.na(Q40))
  
  if (nrow(metric_data_wide) < 2) {
    cat("Warning: Insufficient valid paired samples\n")
    return(NULL)
  }
  
  # compute mean and sd for Q30 and Q40 separately
  mean_Q30 <- mean(metric_data_wide$Q30, na.rm = TRUE)
  sd_Q30 <- sd(metric_data_wide$Q30, na.rm = TRUE)
  mean_Q40 <- mean(metric_data_wide$Q40, na.rm = TRUE)
  sd_Q40 <- sd(metric_data_wide$Q40, na.rm = TRUE)
  
  # compute differences
  differences <- metric_data_wide$Q40 - metric_data_wide$Q30
  
  # === key fix: check whether differences are constant ===
  if (length(unique(differences)) == 1) {
    cat("Warning: All paired differences are identical (value =", 
        unique(differences), ")\n")
    cat("Skipping Shapiro-Wilk test, proceeding directly with Wilcoxon test\n")
    
    # for constant differences normality test is meaningless, use Wilcoxon directly
    test_result <- suppressWarnings(
      wilcox.test(metric_data_wide$Q40, metric_data_wide$Q30,
                  paired = TRUE, alternative = "two.sided",
                  exact = FALSE, conf.int = TRUE)
    )
    
    test_method <- "Wilcoxon signed-rank test"
    test_statistic <- test_result$statistic
    df_value <- NA
    shapiro_p <- NA
    shapiro_note <- "Not tested (constant differences)"
    
  } else if (nrow(metric_data_wide) < 3) {
    # when sample size < 3, Shapiro-Wilk cannot be performed
    cat("Warning: Sample size too small for Shapiro-Wilk test (n < 3)\n")
    cat("Proceeding with Wilcoxon test directly\n")
    
    test_result <- wilcox.test(metric_data_wide$Q40, metric_data_wide$Q30,
                               paired = TRUE, alternative = "two.sided",
                               exact = FALSE, conf.int = TRUE)
    
    test_method <- "Wilcoxon signed-rank test"
    test_statistic <- test_result$statistic
    df_value <- NA
    shapiro_p <- NA
    shapiro_note <- "Not tested (n < 3)"
    
  } else {
    # normal case: perform Shapiro-Wilk test
    shapiro_test <- shapiro.test(differences)
    shapiro_p <- shapiro_test$p.value
    shapiro_note <- ifelse(shapiro_p > 0.05, "Normal", "Non-normal")
    
    # choose test method
    if (shapiro_test$p.value > 0.05) {
      cat("Using paired t-test\n")
      test_result <- t.test(metric_data_wide$Q40, metric_data_wide$Q30, 
                            paired = TRUE, alternative = "two.sided")
      
      test_method <- "Paired t-test"
      test_statistic <- test_result$statistic
      df_value <- test_result$parameter
      
    } else {
      cat("Using Wilcoxon signed-rank test\n")
      test_result <- wilcox.test(metric_data_wide$Q40, metric_data_wide$Q30, 
                                 paired = TRUE, alternative = "two.sided",
                                 exact = FALSE, conf.int = TRUE)
      
      test_method <- "Wilcoxon signed-rank test"
      test_statistic <- test_result$statistic
      df_value <- NA
    }
  }
  
  # compute effect size
  # make sure function calculate_paired_cohens_d is defined
  effect_size <- calculate_paired_cohens_d(
    metric_data_wide$Q30, 
    metric_data_wide$Q40
  )
  
  # create results data frame
  results_df <- data.frame(
    Test = test_method,
    Shapiro_p = shapiro_p,
    Shapiro_note = if(exists("shapiro_note")) shapiro_note else NA,
    Test_statistic = test_statistic,
    df = df_value,
    p_value = test_result$p.value,
    # add descriptive stats for Q30 and Q40
    Mean_Q30 = mean_Q30,
    SD_Q30 = sd_Q30,
    Mean_Q40 = mean_Q40,
    SD_Q40 = sd_Q40,
    # descriptive stats for differences
    Mean_diff = effect_size$mean_diff,
    SD_diff = effect_size$sd_diff,
    Cohens_d = effect_size$cohens_d,
    CI_lower_d = effect_size$ci_lower,
    CI_upper_d = effect_size$ci_upper,
    n_pairs = effect_size$n_pairs,
    stringsAsFactors = FALSE
  )
  
  # print results
  cat("\n=== Results ===\n")
  cat("Q30: Mean =", round(mean_Q30, 3), "SD =", round(sd_Q30, 3), "\n")
  cat("Q40: Mean =", round(mean_Q40, 3), "SD =", round(sd_Q40, 3), "\n")
  cat("p-value:", format.pval(test_result$p.value, digits = 3), "\n")
  if (!is.null(effect_size$cohens_d) && !is.na(effect_size$cohens_d)) {
    cat("Cohen's d:", round(effect_size$cohens_d, 3), 
        "[", round(effect_size$ci_lower, 3), ",", 
        round(effect_size$ci_upper, 3), "]\n")
  }
  
  # add significance flags and interpretation
  results_summary <- results_df %>%
    mutate(
      Significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      Shapiro_result = case_when(
        is.na(Shapiro_p) ~ "Not tested (n<3 or constant differences)",
        Shapiro_p > 0.05 ~ "Normal",
        TRUE ~ "Non-normal"
      ),
      Effect_size_magnitude = case_when(
        abs(Cohens_d) < 0.2 ~ "Negligible",
        abs(Cohens_d) < 0.5 ~ "Small",
        abs(Cohens_d) < 0.8 ~ "Medium",
        TRUE ~ "Large"
      ),
      # format output
      Mean_Q30_str = sprintf("%.3f ± %.3f", Mean_Q30, SD_Q30),
      Mean_Q40_str = sprintf("%.3f ± %.3f", Mean_Q40, SD_Q40),
      Mean_diff_CI = sprintf("%.3f ± %.3f", Mean_diff, SD_diff),
      Cohens_d_CI = sprintf("%.3f [%.3f, %.3f]", 
                            Cohens_d, CI_lower_d, CI_upper_d)
    )
  
  return(results_summary)
}


get_test_type_with_effectsize_fpkm = function(df_clean){
  result_list <- list()
  
  for(group_val in unique(df_clean$group)) {
    
    # select data for current group
    df_group <- df_clean %>% filter(group == group_val)
    
    # extract independent-sample data
    q30_values <- df_group$cv_value[df_group$Quality == "Q30"]
    q40_values <- df_group$cv_value[df_group$Quality == "Q40"]
    
    # sample sizes
    n_q30 <- length(q30_values)
    n_q40 <- length(q40_values)
    
    # use Kolmogorov-Smirnov test instead of Shapiro-Wilk (supports large samples)
    ks_q30 <- ks.test(q30_values, "pnorm", 
                      mean = mean(q30_values), 
                      sd = sd(q30_values))
    ks_q40 <- ks.test(q40_values, "pnorm", 
                      mean = mean(q40_values), 
                      sd = sd(q40_values))
    
    # determine normality test result (if either p < 0.05, mark as non-normal)
    shapiro_p <- min(ks_q30$p.value, ks_q40$p.value)
    shapiro_note <- ifelse(shapiro_p < 0.05, "Non-normal", "Normal")
    shapiro_result <- ifelse(shapiro_p < 0.05, "Non-normal", "Normal")
    
    # note: KS test was used
    shapiro_note <- ifelse(shapiro_p < 0.05, 
                           "Non-normal (KS test)", 
                           "Normal (KS test)")
    
    # Wilcoxon independent-sample test (Mann-Whitney U test)
    test_result <- wilcox.test(cv_value ~ Quality, data = df_group, 
                               paired = FALSE, exact = FALSE)
    
    # use custom function to compute Cohen's d (independent samples)
    cohens_d_result <- calculate_unpaired_cohens_d(q30_values, q40_values)
    
    # compute means and sds
    mean_q30 <- mean(q30_values, na.rm = TRUE)
    sd_q30 <- sd(q30_values, na.rm = TRUE)
    mean_q40 <- mean(q40_values, na.rm = TRUE)
    sd_q40 <- sd(q40_values, na.rm = TRUE)
    mean_diff <- mean_q30 - mean_q40
    
    # compute pooled sd for SD_diff
    n1 <- n_q30
    n2 <- n_q40
    pooled_sd <- sqrt(((n1-1)*sd_q30^2 + (n2-1)*sd_q40^2) / (n1+n2-2))
    sd_diff <- pooled_sd  # use pooled sd as SD_diff
    
    # compute CI for mean difference (independent-sample t-test)
    t_test <- t.test(q30_values, q40_values, 
                     paired = FALSE, 
                     var.equal = FALSE)  # unequal variances
    
    mean_diff_ci_lower <- t_test$conf.int[1]
    mean_diff_ci_upper <- t_test$conf.int[2]
    mean_diff_ci_radius <- (mean_diff_ci_upper - mean_diff)  # CI radius
    
    # format strings
    mean_q30_str <- sprintf("%.3f ± %.3f", mean_q30, sd_q30)
    mean_q40_str <- sprintf("%.3f ± %.3f", mean_q40, sd_q40)
    mean_diff_ci_str <- sprintf("%.3f ± %.3f", mean_diff, mean_diff_ci_radius)
    cohens_d_ci_str <- sprintf("%.3f [%.3f, %.3f]", 
                               cohens_d_result$cohens_d,
                               cohens_d_result$ci_lower,
                               cohens_d_result$ci_upper)
    
    # classify effect size magnitude (based on absolute Cohen's d)
    d_abs <- abs(cohens_d_result$cohens_d)
    effect_size_magnitude <- ifelse(d_abs >= 0.8, "Large",
                                    ifelse(d_abs >= 0.5, "Medium",
                                           ifelse(d_abs >= 0.2, "Small", "Very small")))
    
    # significance flag
    significance <- ifelse(test_result$p.value < 0.001, "***",
                           ifelse(test_result$p.value < 0.01, "**",
                                  ifelse(test_result$p.value < 0.05, "*", "ns")))
    
    # extract relevant statistics
    result_list[[group_val]] <- data.frame(
      Group=group_val,
      Test = "Wilcoxon rank-sum test",
      Shapiro_p = shapiro_p,
      Shapiro_note = shapiro_note,
      Test_statistic = as.numeric(test_result$statistic),
      df = NA,
      p_value = test_result$p.value,
      Mean_Q30 = mean_q30,
      SD_Q30 = sd_q30,
      Mean_Q40 = mean_q40,
      SD_Q40 = sd_q40,
      Mean_diff = mean_diff,
      SD_diff = sd_diff,
      Cohens_d = cohens_d_result$cohens_d,
      CI_lower_d = cohens_d_result$ci_lower,
      CI_upper_d = cohens_d_result$ci_upper,
      n_pairs = NA,
      Significance = significance,
      Shapiro_result = shapiro_result,
      Effect_size_magnitude = effect_size_magnitude,
      Mean_Q30_str = mean_q30_str,
      Mean_Q40_str = mean_q40_str,
      Mean_diff_CI = mean_diff_ci_str,
      Cohens_d_CI = cohens_d_ci_str,
      stringsAsFactors = FALSE
    )
  }
  
  # combine all results
  results_df <- do.call(rbind, result_list)
  rownames(results_df) <- NULL
  # return results
  return(results_df)
}