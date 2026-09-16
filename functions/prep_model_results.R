#####################################################################
#
# Data prep for visualization
#
# Author: F.Dhellemmes
# Last upd: 05.01.2026
#
#####################################################################

ci_func <- function(x) {
  c(mean = mean(x), ci_lower = quantile(x, 0.025), ci_upper = quantile(x, 0.975))
}

prep_model_results <- function(fit,response,group,random, feature_idx){
  # Convert arguments to strings
  response_name <- deparse(substitute(response))
  group_name <- deparse(substitute(group))
  random_name <- deparse(substitute(random))
  
  
  # Build the spread_draws expression dynamically: e.g., beta[feature]
  var_expr <- rlang::parse_expr(paste0(response_name, "[", group_name, "]"))
  tmp_main <- tidybayes::spread_draws(fit, !!var_expr)
  
  formula_agg <- as.formula(paste(response_name, "~", group_name))
  
  # Aggregate mean per group
  main_effects <- aggregate(
    formula_agg,
    data = tmp_main,
    FUN = mean
  )
  
  names(main_effects)[2] <- "mean_effect"
  
  
  random_expr <- rlang::parse_expr(paste0("v_",random_name,"[",random_name, ",", group_name, "]"))
  
  tmp_comp <- tidybayes::spread_draws(fit, !!random_expr)
  
  tmp_comp <- tmp_comp[tmp_comp$feature %in% feature_idx, ]
  
  tmp_comp$feature <- tmp_comp$feature - min(feature_idx) +1
  
  formula_effects <- as.formula(paste0("v_",random_name, " ~ ", group_name, " + ", random_name))
  
  comp_effects <- aggregate(
    formula_effects,
    data = tmp_comp,
    FUN = ci_func)
  
  comp_effects<-do.call(data.frame, comp_effects)
  names(comp_effects)[3] <- "mean_effect_offset"
  names(comp_effects)[4] <- "lower_effect_offset"
  names(comp_effects)[5] <- "upper_effect_offset"
  
  effects <- merge(
    main_effects,
    comp_effects,
    by = "feature",
    all.x = TRUE)
  
  names(effects)[3] <- "group"
  
  
  return(effects)
  
}


