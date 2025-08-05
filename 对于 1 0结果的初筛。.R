
run_evalue_analysis <- function(data, 
                                exposure_vars = NULL,
                                outcome_var = NULL,
                                duration_var = "Duration",
                                output_type = "both", # "OR", "HR", "both"
                                p_threshold = c(0.05, 0.01, 0.001),
                                unique_val_threshold = 10,
                                skewness_threshold = 1) {  # 新增偏度阈值参数
  
  library(moments)
  library(dplyr)
  library(purrr)
  library(broom)
  library(survival)
  library(EValue)
  
  # 创建结果列表
  all_logistic_results <- list()
  all_cox_results <- list()
  
  # 函数：计算变量分组的outcome prevalence并确定rare参数
  calculate_rare_parameter <- function(data, variable_name, outcome_var, is_continuous) {
    tryCatch({
      if (!is_continuous) {
        # 对于分类变量，计算每个分组的prevalence
        prevalence_by_group <- data %>%
          group_by(!!sym(variable_name)) %>%
          summarise(
            total = n(),
            cases = sum(!!sym(outcome_var) == 1, na.rm = TRUE),
            prevalence = cases / total,
            .groups = 'drop'
          )
        
        cat("Prevalence by", variable_name, "groups:\n")
        print(prevalence_by_group)
        
        # 如果任何一个分组的prevalence >= 15%，则rare = FALSE
        max_prevalence <- max(prevalence_by_group$prevalence, na.rm = TRUE)
        rare_param <- max_prevalence < 0.15
        
      } else {
        # 对于连续变量，计算整体prevalence
        overall_prevalence <- sum(data[[outcome_var]] == 1, na.rm = TRUE) / nrow(data)
        cat("Overall prevalence for", variable_name, ":", round(overall_prevalence * 100, 2), "%\n")
        rare_param <- overall_prevalence < 0.15
      }
      
      cat("Setting rare =", rare_param, "for variable", variable_name, "\n\n")
      return(rare_param)
      
    }, error = function(e) {
      cat("Error calculating prevalence for", variable_name, ":", e$message, "\n")
      return(FALSE)  # 默认设置为FALSE（非rare）
    })
  }
  
  # 函数：计算E值（OR）
  calculate_evalue_or <- function(or, ci_lower, ci_upper, rare_param) {
    tryCatch({
      result <- evalues.OR(est = or, lower = ci_lower, upper = ci_upper, rare = rare_param)
      return(list(
        e_value_point = result$point.estimate[1],
        e_value_lower = result$lower.confidence.limit[1]
      ))
    }, error = function(e) {
      return(list(e_value_point = NA, e_value_lower = NA))
    })
  }
  
  # 函数：计算E值（HR）
  calculate_evalue_hr <- function(hr, ci_lower, ci_upper, rare_param) {
    tryCatch({
      result <- evalues.HR(est = hr, lower = ci_lower, upper = ci_upper, rare = rare_param)
      return(list(
        e_value_point = result$point.estimate[1],
        e_value_lower = result$lower.confidence.limit[1]
      ))
    }, error = function(e) {
      return(list(e_value_point = NA, e_value_lower = NA))
    })
  }
  
  # 函数：判断正态分布并创建标准化变量和分组变量
  create_continuous_vars <- function(data, var_name) {
    vals <- data[[var_name]]
    vals_clean <- na.omit(vals)
    
    # 计算偏度
    skewness_val <- moments::skewness(vals_clean)
    cat("Skewness for", var_name, ":", round(skewness_val, 3), "\n")
    
    # 判断是否为正态分布
    is_normal <- abs(skewness_val) <= skewness_threshold
    cat("Distribution type:", ifelse(is_normal, "Normal", "Skewed"), 
        "(threshold:", skewness_threshold, ")\n")
    
    new_vars <- list()
    
    # 1. 创建标准化变量（SD或IQR）
    if (is_normal) {
      # 正态分布：创建每SD增加的变量
      mean_val <- mean(vals, na.rm = TRUE)
      sd_val <- sd(vals, na.rm = TRUE)
      
      standardized_var_name <- paste0(var_name, "_SD")
      new_vars[[standardized_var_name]] <- (vals - mean_val) / sd_val
      
      cat("Created standardized variable:", standardized_var_name, 
          "(Mean:", round(mean_val, 3), ", SD:", round(sd_val, 3), ")\n")
      
    } else {
      # 偏态分布：创建每IQR增加的变量
      q25 <- quantile(vals, 0.25, na.rm = TRUE)
      q75 <- quantile(vals, 0.75, na.rm = TRUE)
      iqr_val <- q75 - q25
      median_val <- median(vals, na.rm = TRUE)
      
      iqr_var_name <- paste0(var_name, "_IQR")
      new_vars[[iqr_var_name]] <- (vals - median_val) / iqr_val
      
      cat("Created IQR-standardized variable:", iqr_var_name, 
          "(Median:", round(median_val, 3), ", IQR:", round(iqr_val, 3), ")\n")
    }
    
    # 2. 创建三分位数分组
    tertile_var_name <- paste0(var_name, "_Tertile")
    tertile_cuts <- quantile(vals, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
    new_vars[[tertile_var_name]] <- cut(vals, 
                                        breaks = tertile_cuts, 
                                        labels = c("T1", "T2", "T3"),
                                        include.lowest = TRUE)
    
    cat("Created tertile variable:", tertile_var_name, "\n")
    cat("Tertile cutpoints:", paste(round(tertile_cuts, 2), collapse = ", "), "\n")
    
    # 3. 创建四分位数分组
    quartile_var_name <- paste0(var_name, "_Quartile")
    quartile_cuts <- quantile(vals, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    new_vars[[quartile_var_name]] <- cut(vals, 
                                         breaks = quartile_cuts, 
                                         labels = c("Q1", "Q2", "Q3", "Q4"),
                                         include.lowest = TRUE)
    
    cat("Created quartile variable:", quartile_var_name, "\n")
    cat("Quartile cutpoints:", paste(round(quartile_cuts, 2), collapse = ", "), "\n")
    
    return(new_vars)
  }
  
  # 函数：运行单个变量的回归分析
  run_single_analysis <- function(data, var_name, var_data, outcome_var, duration_var, 
                                  analysis_type, main_var_name, var_type) {
    
    # 计算rare参数
    is_continuous_var <- var_type %in% c("Standardized", "IQR")
    rare_param <- calculate_rare_parameter(data, var_name, outcome_var, is_continuous_var)
    
    results <- list()
    
    # Logistic回归
    if(analysis_type %in% c("OR", "both")) {
      # 创建临时数据框
      temp_data <- data
      temp_data[[var_name]] <- var_data
      
      formula_str <- paste(outcome_var, "~", var_name)
      logistic_model <- tryCatch({
        glm(as.formula(formula_str), data = temp_data, family = binomial(link = "logit"))
      }, error = function(e) {
        cat("Error fitting logistic model for", var_name, ":", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(logistic_model)) {
        logistic_results <- tryCatch({
          model_summary <- summary(logistic_model)
          coef_df <- as.data.frame(model_summary$coefficients)
          coef_df$Variable <- rownames(coef_df)
          
          coef_df <- coef_df %>%
            filter(Variable != "(Intercept)") %>%
            mutate(
              OR = exp(Estimate),
              CI_lower = exp(Estimate - 1.96 * `Std. Error`),
              CI_upper = exp(Estimate + 1.96 * `Std. Error`),
              p_value_raw = `Pr(>|z|)`
            )
          
          # 计算E值
          if(nrow(coef_df) > 0) {
            evalue_results <- purrr::map_dfr(1:nrow(coef_df), function(i) {
              evalue_result <- calculate_evalue_or(
                coef_df$OR[i], 
                coef_df$CI_lower[i], 
                coef_df$CI_upper[i], 
                rare_param
              )
              
              data.frame(
                Variable = coef_df$Variable[i],
                OR = coef_df$OR[i],
                CI_lower = coef_df$CI_lower[i],
                CI_upper = coef_df$CI_upper[i],
                p_value_raw = coef_df$p_value_raw[i],
                E_value_point = evalue_result$e_value_point,
                E_value_lower = evalue_result$e_value_lower,
                stringsAsFactors = FALSE
              )
            })
            
            # 格式化输出
            evalue_results <- evalue_results %>%
              mutate(
                `OR (95% CI)` = sprintf("%.3f (%.3f-%.3f)", OR, CI_lower, CI_upper),
                `p-value` = ifelse(p_value_raw < 0.001, "<0.001", 
                                   format.pval(p_value_raw, digits = 3, eps = 0.001)),
                `E-value (95% CI)` = ifelse(is.na(E_value_point), "NA", 
                                            sprintf("%.2f (%.2f)", E_value_point, E_value_lower)),
                Rare_Parameter = rare_param,
                Model = "Unadjusted",
                Main_Variable = main_var_name,
                Analysis_Type = "Logistic",
                Variable_Type = var_type
              )
            
            results$logistic <- evalue_results
          }
        }, error = function(e) {
          cat("Error extracting logistic results for", var_name, ":", e$message, "\n")
        })
      }
    }
    
    # Cox回归
    if(analysis_type %in% c("HR", "both")) {
      # 创建临时数据框
      temp_data <- data
      temp_data[[var_name]] <- var_data
      
      cox_formula_str <- paste("Surv(", duration_var, ",", outcome_var, ") ~", var_name)
      cox_model <- tryCatch({
        coxph(as.formula(cox_formula_str), data = temp_data)
      }, error = function(e) {
        cat("Error fitting Cox model for", var_name, ":", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(cox_model)) {
        cox_results <- tryCatch({
          model_summary <- summary(cox_model)
          coef_df <- as.data.frame(model_summary$coefficients)
          coef_df$Variable <- rownames(coef_df)
          
          coef_df <- coef_df %>%
            mutate(
              HR = exp(coef),
              CI_lower = exp(coef - 1.96 * `se(coef)`),
              CI_upper = exp(coef + 1.96 * `se(coef)`),
              p_value_raw = `Pr(>|z|)`
            )
          
          # 计算E值
          if(nrow(coef_df) > 0) {
            evalue_results <- purrr::map_dfr(1:nrow(coef_df), function(i) {
              evalue_result <- calculate_evalue_hr(
                coef_df$HR[i], 
                coef_df$CI_lower[i], 
                coef_df$CI_upper[i], 
                rare_param
              )
              
              data.frame(
                Variable = coef_df$Variable[i],
                HR = coef_df$HR[i],
                CI_lower = coef_df$CI_lower[i],
                CI_upper = coef_df$CI_upper[i],
                p_value_raw = coef_df$p_value_raw[i],
                E_value_point = evalue_result$e_value_point,
                E_value_lower = evalue_result$e_value_lower,
                stringsAsFactors = FALSE
              )
            })
            
            # 格式化输出
            evalue_results <- evalue_results %>%
              mutate(
                `HR (95% CI)` = sprintf("%.3f (%.3f-%.3f)", HR, CI_lower, CI_upper),
                `p-value` = ifelse(p_value_raw < 0.001, "<0.001", 
                                   format.pval(p_value_raw, digits = 3, eps = 0.001)),
                `E-value (95% CI)` = ifelse(is.na(E_value_point), "NA", 
                                            sprintf("%.2f (%.2f)", E_value_point, E_value_lower)),
                Rare_Parameter = rare_param,
                Model = "Unadjusted",
                Main_Variable = main_var_name,
                Analysis_Type = "Cox",
                Variable_Type = var_type
              )
            
            results$cox <- evalue_results
          }
        }, error = function(e) {
          cat("Error extracting Cox results for", var_name, ":", e$message, "\n")
        })
      }
    }
    
    return(results)
  }
  
  # 主循环：处理每个暴露变量
  for(var in exposure_vars) {
    
    cat("\n=== Processing variable:", var, "===\n")
    
    # 判断变量类型 ----------------------------------------------------------------
    is_cont <- FALSE
    if (is.numeric(data[[var]])) {
      unique_vals <- length(unique(na.omit(data[[var]])))
      is_cont <- unique_vals > unique_val_threshold
      cat("Variable type: numeric with", unique_vals, "unique values\n")
      cat("Treated as continuous variable:", is_cont, "(threshold:", unique_val_threshold, ")\n")
    } else if (is.factor(data[[var]]) || is.character(data[[var]])) {
      cat("Variable type: categorical (factor or character)\n")
      is_cont <- FALSE
    } else {
      cat("Variable type:", class(data[[var]]), "treated as categorical by default\n")
      is_cont <- FALSE
    }
    
    # 处理连续变量的极端值 ---------------------------------------------------------
    temp_data <- data
    if (is_cont) {
      vals <- temp_data[[var]]
      mean_val <- mean(vals, na.rm = TRUE)
      sd_val <- sd(vals, na.rm = TRUE)
      lower_bound <- mean_val - 3 * sd_val
      upper_bound <- mean_val + 3 * sd_val
      
      # 标识极端值
      is_outlier <- vals < lower_bound | vals > upper_bound
      is_outlier[is.na(is_outlier)] <- FALSE  # 处理NA值
      n_outliers <- sum(is_outlier, na.rm = TRUE)
      
      if (n_outliers > 0) {
        temp_data <- temp_data[!is_outlier, ]
        cat("Removed", n_outliers, "outliers (3SD rule: ",
            round(lower_bound, 2), "-", round(upper_bound, 2), ")\n")
        cat("Remaining sample size:", nrow(temp_data), "\n")
      } else {
        cat("No outliers detected (3SD rule)\n")
      }
      
      # 为连续变量创建标准化变量和分组变量
      cat("\n--- Creating standardized and grouped variables ---\n")
      new_vars <- create_continuous_vars(temp_data, var)
      
      # 将新变量添加到数据中
      for(new_var_name in names(new_vars)) {
        temp_data[[new_var_name]] <- new_vars[[new_var_name]]
      }
      
      # 对每个新创建的变量进行分析
      for(new_var_name in names(new_vars)) {
        cat("\n--- Analyzing", new_var_name, "---\n")
        
        # 确定变量类型
        if(grepl("_SD$", new_var_name)) {
          var_type <- "Standardized"
        } else if(grepl("_IQR$", new_var_name)) {
          var_type <- "IQR"
        } else if(grepl("_Tertile$", new_var_name)) {
          var_type <- "Tertile"
        } else if(grepl("_Quartile$", new_var_name)) {
          var_type <- "Quartile"
        } else {
          var_type <- "Other"
        }
        
        # 运行分析
        analysis_results <- run_single_analysis(
          temp_data, new_var_name, temp_data[[new_var_name]], 
          outcome_var, duration_var, output_type, var, var_type
        )
        
        # 存储结果
        if(!is.null(analysis_results$logistic)) {
          all_logistic_results[[paste(var, new_var_name, "logistic", sep = "_")]] <- analysis_results$logistic
        }
        if(!is.null(analysis_results$cox)) {
          all_cox_results[[paste(var, new_var_name, "cox", sep = "_")]] <- analysis_results$cox
        }
      }
      
    } else {
      # 分类变量的原有处理逻辑
      # 计算当前变量的rare参数
      rare_param <- calculate_rare_parameter(temp_data, var, outcome_var, is_cont)
      
      # 运行分析
      analysis_results <- run_single_analysis(
        temp_data, var, temp_data[[var]], 
        outcome_var, duration_var, output_type, var, "Categorical"
      )
      
      # 存储结果
      if(!is.null(analysis_results$logistic)) {
        all_logistic_results[[paste(var, "logistic", sep = "_")]] <- analysis_results$logistic
      }
      if(!is.null(analysis_results$cox)) {
        all_cox_results[[paste(var, "cox", sep = "_")]] <- analysis_results$cox
      }
    }
    
    cat("=== Completed variable:", var, "===\n\n")
  }
  
  # 合并结果
  final_results <- list()
  
  if(output_type %in% c("OR", "both") && length(all_logistic_results) > 0) {
    final_logistic_results <- bind_rows(all_logistic_results) %>%
      rename(Effect_Size = OR, `Effect_Size (95% CI)` = `OR (95% CI)`) %>%
      mutate(Effect_Type = "OR")
    final_results$logistic <- final_logistic_results
  }
  
  if(output_type %in% c("HR", "both") && length(all_cox_results) > 0) {
    final_cox_results <- bind_rows(all_cox_results) %>%
      rename(Effect_Size = HR, `Effect_Size (95% CI)` = `HR (95% CI)`) %>%
      mutate(Effect_Type = "HR")
    final_results$cox <- final_cox_results
  }
  
  # 合并所有结果
  if(length(final_results) > 0) {
    all_results_combined <- bind_rows(final_results) %>%
      select(Main_Variable, Analysis_Type, Effect_Type, Model, Variable, Variable_Type,
             `Effect_Size (95% CI)`, `E-value (95% CI)`, `p-value`, 
             Effect_Size, CI_lower, CI_upper, E_value_point, E_value_lower, Rare_Parameter, p_value_raw)
  } else {
    all_results_combined <- data.frame()
  }
  
  # 筛选p值显著的结果
  significant_results <- list()
  
  for(threshold in p_threshold) {
    if(nrow(all_results_combined) > 0) {
      # 使用 p_value_raw 进行筛选
      sig_subset <- all_results_combined %>%
        filter(p_value_raw < threshold) %>%
        mutate(Significance_Level = paste0("p < ", threshold)) %>%
        select(-p_value_raw)  # 移除原始p值列
      
      if(nrow(sig_subset) > 0) {
        significant_results[[paste0("p_", gsub("\\.", "", threshold))]] <- sig_subset
      }
    }
  }
  
  # 从最终结果中移除原始p值列
  if(nrow(all_results_combined) > 0) {
    all_results_combined <- all_results_combined %>%
      select(-p_value_raw)
  }
  
  # 返回结果
  return(list(
    all_results = all_results_combined,
    significant_results = significant_results,
    summary = list(
      total_variables = length(exposure_vars),
      output_type = output_type,
      p_thresholds = p_threshold,
      total_results = nrow(all_results_combined),
      significant_counts = sapply(significant_results, nrow)
    )
  ))
}


