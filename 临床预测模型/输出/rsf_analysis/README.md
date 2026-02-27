随机生存森林效果不如Cox回归，存在过拟合风险，main中未保留代码

``` R
# 随机生存森林
library("randomForestSRC")

rsf_output_dir <- file.path(output_dir, "rsf_analysis")
dir.create(rsf_output_dir, recursive = TRUE, showWarnings = FALSE)

# 沿用main.r中的多因素变量
rsf_predictors <- multivariate_predictors
rsf_formula <- as.formula(
  paste("Surv(time, status == 1) ~", paste(rsf_predictors, collapse = " + "))
)

set.seed(412)
rsf_model <- rfsrc(
  formula = rsf_formula,
  data = train_ready[, c("time", "status", rsf_predictors)],
  ntree = 1000,
  nodesize = 15,
  nsplit = 10,
  importance = TRUE,
  na.action = "na.impute",
  seed = 412
)

rsf_train_pred <- predict(rsf_model, newdata = train_ready)
rsf_test_pred <- predict(rsf_model, newdata = test_ready)

calc_cindex <- function(data, risk_score) {
  concordance(Surv(time, status == 1) ~ risk_score, data = data)$concordance
}

rsf_eval <- data.frame(
  dataset = c("Training", "Validation"),
  c_index = c(
    calc_cindex(train_ready, as.numeric(rsf_train_pred$predicted)),
    calc_cindex(test_ready, as.numeric(rsf_test_pred$predicted))
  )
)
write.csv(rsf_eval, file = file.path(rsf_output_dir, "rsf_cindex.csv"), row.names = FALSE)

rsf_vimp <- data.frame(
  variable = names(rsf_model$importance),
  importance = as.numeric(rsf_model$importance),
  stringsAsFactors = FALSE
)
rsf_vimp <- rsf_vimp[order(rsf_vimp$importance, decreasing = TRUE), , drop = FALSE]
write.csv(rsf_vimp, file = file.path(rsf_output_dir, "rsf_variable_importance.csv"), row.names = FALSE)
```

``` R 
# DCA
dca_required_funs <- c("compute_dca_curve", "plot_dca_curve", "collect_dca_rows")
if (!all(dca_required_funs %in% ls())) {
  stop("请先运行main.r中DCA相关函数定义，再执行该代码块。")
}

dca_time_points <- c(6, 12, 18)
dca_output_dir <- file.path(rsf_output_dir, "dca_plots")
dir.create(dca_output_dir, recursive = TRUE, showWarnings = FALSE)

predict_rsf_risk_at_time <- function(rsf_pred, time_point) {
  surv_matrix <- as.matrix(rsf_pred$survival)
  time_grid <- as.numeric(rsf_pred$time.interest)
  if (length(time_grid) == 0 || ncol(surv_matrix) == 0) return(rep(NA_real_, nrow(surv_matrix)))

  # 取最接近目标时间点的生存概率并转换为风险概率
  time_index <- which.min(abs(time_grid - time_point))
  pmin(pmax(1 - as.numeric(surv_matrix[, time_index]), 0), 1)
}

set.seed(412)
rsf_validation_refit_model <- rfsrc(
  formula = rsf_formula,
  data = test_ready[, c("time", "status", rsf_predictors)],
  ntree = 1000,
  nodesize = 15,
  nsplit = 10,
  importance = TRUE,
  na.action = "na.impute",
  seed = 412
)
rsf_validation_refit_pred <- predict(rsf_validation_refit_model, newdata = test_ready)

rsf_dca_rows <- list()
for (time_point in dca_time_points) {
  rsf_training_refit_risk <- predict_rsf_risk_at_time(rsf_train_pred, time_point)
  rsf_validation_refit_risk <- predict_rsf_risk_at_time(rsf_validation_refit_pred, time_point)
  rsf_validation_external_risk <- predict_rsf_risk_at_time(rsf_test_pred, time_point)

  rsf_dca_training_refit <- compute_dca_curve(train_ready, rsf_training_refit_risk, time_point)
  rsf_dca_validation_refit <- compute_dca_curve(test_ready, rsf_validation_refit_risk, time_point)
  rsf_dca_validation_external <- compute_dca_curve(test_ready, rsf_validation_external_risk, time_point)

  plot_dca_curve(
    dca_table = rsf_dca_training_refit,
    time_point = time_point,
    output_file = file.path(dca_output_dir, sprintf("dca_%d_training_refit.png", time_point)),
    cohort_label = "Training cohort",
    strategy_label = "RSF Refit"
  )
  plot_dca_curve(
    dca_table = rsf_dca_validation_refit,
    time_point = time_point,
    output_file = file.path(dca_output_dir, sprintf("dca_%d_validation_refit.png", time_point)),
    cohort_label = "Validation cohort",
    strategy_label = "RSF Refit"
  )
  plot_dca_curve(
    dca_table = rsf_dca_validation_external,
    time_point = time_point,
    output_file = file.path(dca_output_dir, sprintf("dca_%d_validation_external.png", time_point)),
    cohort_label = "Validation cohort",
    strategy_label = "RSF External"
  )

  rsf_dca_rows[[length(rsf_dca_rows) + 1L]] <- transform(
    collect_dca_rows(rsf_dca_training_refit, strategy = "refit", cohort = "training", time_point = time_point),
    model = "rsf"
  )
  rsf_dca_rows[[length(rsf_dca_rows) + 1L]] <- transform(
    collect_dca_rows(rsf_dca_validation_refit, strategy = "refit", cohort = "validation", time_point = time_point),
    model = "rsf"
  )
  rsf_dca_rows[[length(rsf_dca_rows) + 1L]] <- transform(
    collect_dca_rows(rsf_dca_validation_external, strategy = "external", cohort = "validation", time_point = time_point),
    model = "rsf"
  )
}

rsf_dca_summary <- do.call(rbind, rsf_dca_rows)
rsf_dca_summary <- rsf_dca_summary[, c("model", "strategy", "cohort", "time_month", "threshold", "net_benefit")]
write.csv(rsf_dca_summary, file = file.path(dca_output_dir, "dca_net_benefit_summary.csv"), row.names = FALSE)

cox_dca_file <- file.path(output_dir, "dca_plots", "dca_net_benefit_summary.csv")
if (file.exists(cox_dca_file)) {
  cox_dca_summary <- read.csv(cox_dca_file, stringsAsFactors = FALSE)
  cox_dca_summary$model <- "cox"
  cox_dca_summary <- cox_dca_summary[, c("model", "strategy", "cohort", "time_month", "threshold", "net_benefit")]
  write.csv(
    rbind(cox_dca_summary, rsf_dca_summary),
    file = file.path(rsf_output_dir, "dca_compare.csv"),
    row.names = FALSE
  )
}

```

``` R
# 时间依赖ROC
roc_required_funs <- c("plot_roc_combined", "collect_roc_auc_rows")
if (!all(roc_required_funs %in% ls())) {
  stop("请先运行main.r中ROC相关函数定义，再执行该代码块。")
}
if (!requireNamespace("survivalROC", quietly = TRUE)) {
  stop("Package 'survivalROC' is required for time-dependent ROC analysis.")
}

if (!exists("predict_rsf_risk_at_time")) {
  predict_rsf_risk_at_time <- function(rsf_pred, time_point) {
    surv_matrix <- as.matrix(rsf_pred$survival)
    time_grid <- as.numeric(rsf_pred$time.interest)
    if (length(time_grid) == 0 || ncol(surv_matrix) == 0) return(rep(NA_real_, nrow(surv_matrix)))
    time_index <- which.min(abs(time_grid - time_point))
    pmin(pmax(1 - as.numeric(surv_matrix[, time_index]), 0), 1)
  }
}
if (!exists("rsf_validation_refit_pred")) {
  set.seed(412)
  rsf_validation_refit_model <- rfsrc(
    formula = rsf_formula,
    data = test_ready[, c("time", "status", rsf_predictors)],
    ntree = 1000,
    nodesize = 15,
    nsplit = 10,
    importance = TRUE,
    na.action = "na.impute",
    seed = 412
  )
  rsf_validation_refit_pred <- predict(rsf_validation_refit_model, newdata = test_ready)
}

roc_time_points <- c(6, 12, 18)
roc_colors <- c("6" = "#66ccff", "12" = "#EE0000", "18" = "#FFCC66")
roc_output_dir <- file.path(rsf_output_dir, "roc_plots")
dir.create(roc_output_dir, recursive = TRUE, showWarnings = FALSE)

compute_time_roc_list_rsf <- function(data, rsf_pred, time_points = roc_time_points) {
  roc_list <- lapply(time_points, function(time_point) {
    risk_marker <- predict_rsf_risk_at_time(rsf_pred, time_point)

    roc_result <- tryCatch(
      survivalROC::survivalROC(
        Stime = data$time,
        status = data$status,
        marker = risk_marker,
        predict.time = time_point,
        method = "KM"
      ),
      error = function(e) NULL
    )

    if (is.null(roc_result)) return(list(FP = numeric(), TP = numeric(), AUC = NA_real_))
    roc_result
  })

  names(roc_list) <- as.character(time_points)
  roc_list
}

rsf_roc_training_refit <- compute_time_roc_list_rsf(train_ready, rsf_train_pred)
rsf_roc_validation_refit <- compute_time_roc_list_rsf(test_ready, rsf_validation_refit_pred)
rsf_roc_validation_external <- compute_time_roc_list_rsf(test_ready, rsf_test_pred)

plot_roc_combined(
  roc_results = rsf_roc_training_refit,
  output_file = file.path(roc_output_dir, "roc_training_refit.png"),
  time_points = roc_time_points,
  colors = roc_colors
)
plot_roc_combined(
  roc_results = rsf_roc_validation_refit,
  output_file = file.path(roc_output_dir, "roc_validation_refit.png"),
  time_points = roc_time_points,
  colors = roc_colors
)
plot_roc_combined(
  roc_results = rsf_roc_validation_external,
  output_file = file.path(roc_output_dir, "roc_validation_external.png"),
  time_points = roc_time_points,
  colors = roc_colors
)

rsf_roc_summary <- rbind(
  transform(collect_roc_auc_rows(rsf_roc_training_refit, strategy = "refit", cohort = "training"), model = "rsf"),
  transform(collect_roc_auc_rows(rsf_roc_validation_refit, strategy = "refit", cohort = "validation"), model = "rsf"),
  transform(collect_roc_auc_rows(rsf_roc_validation_external, strategy = "external", cohort = "validation"), model = "rsf")
)
rsf_roc_summary <- rsf_roc_summary[, c("model", "strategy", "cohort", "time_month", "auc")]
write.csv(rsf_roc_summary, file = file.path(roc_output_dir, "roc_auc_summary.csv"), row.names = FALSE)

cox_roc_file <- file.path(output_dir, "roc_plots", "roc_auc_summary.csv")
if (file.exists(cox_roc_file)) {
  cox_roc_summary <- read.csv(cox_roc_file, stringsAsFactors = FALSE)
  cox_roc_summary$model <- "cox"
  cox_roc_summary <- cox_roc_summary[, c("model", "strategy", "cohort", "time_month", "auc")]
  write.csv(
    rbind(cox_roc_summary, rsf_roc_summary),
    file = file.path(rsf_output_dir, "roc_compare.csv"),
    row.names = FALSE
  )
}
```
