# -*- coding: utf-8 -*-
### ============================================================
### 第一部分：环境配置和工作目录设置
### ============================================================

library("compareGroups")
library("survival")
library("survminer")
library("MASS")
library("regplot")
library("forestplot")
library("rms")
library("here")
library("quarto")

# 使用 here::here() 获取项目根目录
base_dir <- here::here()

base_dir <- file.path(base_dir, "复现", "临床预测模型")
inputed_dir <- file.path(base_dir, "参考文献及原始数据")
output_dir <- file.path(base_dir, "输出")


### ============================================================
### 第二部分：数据规格定义 - 因子变量配置
### ============================================================

factor_specs <- list(
  age = list(levels = c(1, 2, 3), labels = c("≤ 54", "55-77", "> 78"), label = "Age at diagnosis "),
  sex = list(levels = c(1, 2), labels = c("Female", "Male"), label = "Gender "),
  marital = list(levels = c(1, 2), labels = c("Married", "Unmarried"), label = "Marital status "),
  race = list(levels = c(1, 2, 3), labels = c("White", "Black", "Other"), label = "Race "),
  site = list(levels = c(1, 2, 3, 4, 5), labels = c("Oral cavity", "Oropharynx", "Nasopharynx", "Hypopharynx", "Larynx"), label = "Primary_site "),
  grade = list(levels = c(1, 2, 3, 4), labels = c("I", "II", "III", "IV"), label = "Grade "),
  t = list(levels = c(1, 2, 3, 4), labels = c("T1", "T2", "T3", "T4"), label = "T_stage "),
  n = list(levels = c(1, 2, 3, 4), labels = c("N0", "N1", "N2", "N3"), label = "N_stage "),
  bone = list(levels = c(0, 1), labels = c("No", "Yes"), label = "Bone_met "),
  brain = list(levels = c(0, 1), labels = c("No", "Yes"), label = "Brain_met "),
  liver = list(levels = c(0, 1), labels = c("No", "Yes"), label = "Liver_met "),
  lung = list(levels = c(0, 1), labels = c("No", "Yes"), label = "Lung_met "),
  multi = list(levels = c(0, 1), labels = c("Absence", "Presence"), label = "Multiple primary carcinomas "),
  sur = list(levels = c(0, 1), labels = c("No", "Yes"), label = "Surgery "),
  rad = list(levels = c(0, 1), labels = c("No/unknown", "Yes"), label = "Radiotherapy "),
  che = list(levels = c(0, 1), labels = c("No/unknown", "Yes"), label = "Chemotherapy ")
)

### ============================================================
### 第三部分：辅助函数定义
### ============================================================

# 将数值向量转换为带标签的因子变量
# 将数值转为因子类型: factor()
# 设置对象属性: attr()
# 函数最后一行的表达式自动作为返回值，无需显式使用return()
to_factor_with_label <- function(values, levels, labels, label_text) {

  factor_values <- factor(values, levels = levels, labels = labels)

  attr(factor_values, "label") <- label_text

  factor_values

}

# 批量处理数据框中的所有变量，将数值编码转换为带标签的因子
# $符号访问数据框的列，names()返回列表的所有键名
# [[]]从列表提取元素本身，[]返回包含该元素的子列表
prepare_variables <- function(data) {

  data$time <- as.numeric(data$time)

  vars <- intersect(names(factor_specs), names(data))
  
  data[vars] <- lapply(vars, function(var_name) {

    spec <- factor_specs[[var_name]]
    to_factor_with_label(
      data[[var_name]],
      levels = spec$levels, labels = spec$labels, label_text = spec$label
    )

  })

  data

}

# 生成Cox回归结果表格，将统计结果整理成易读的表格格式
# data = 数据框,
# variables = 变量名向量,
# cox_coef = 系数矩阵(含P值),
# cox_conf = 置信区间矩阵(含HR)
generate_cox_table <- function(data, variables, cox_coef, cox_conf) {

  # 创建空数据框模板，定义结果表的列结构
  # stringsAsFactors=FALSE表示不自动将字符串转为因子
  empty_result <- data.frame(
    Variable = character(), Patients = character(), CI = character(), P = character(), HR = numeric(), LowerCI = numeric(),  UpperCI = numeric(),
    stringsAsFactors = FALSE
  )

  result_rows <- list() # 用列表存储结果行（列表追加元素比数据框更高效）

  # 添加一行数据
  # 通过返回更新后的rows避免使用超级赋值，L后缀表示整数类型
  append_row <- function(rows, variable, patients, ci, p, hr, lower_ci, upper_ci) {

    rows[[length(rows) + 1L]] <- data.frame(
      Variable = variable, Patients = patients, CI = ci, P = p, HR = hr, LowerCI = lower_ci, UpperCI = upper_ci,
      stringsAsFactors = FALSE
    )

    rows

  }

  # 遍历每个变量，为每个变量生成多行结果
  for (var in variables) {
    
    # 获取变量的label属性
    # 如果没有label则使用变量名
    var_label <- attr(data[[var]], "label")
    if (is.null(var_label)) var_label <- var

    # 添加变量标题行（NA_real_是实数类型缺失值）
    # 获取因子变量的所有水平（类别）
    result_rows <- append_row(result_rows, var_label, "", "", "", NA_real_, NA_real_, NA_real_)
    levels_var <- levels(data[[var]]) 

    # seq_along()生成从1到向量长度的整数序列
    # 如: seq_along(c("a","b","c"))，返回1 2 3
    # for循环遍历每个水平
    for (i in seq_along(levels_var)) {

      # 获取当前水平的名称
      level <- levels_var[i]

      # 计算该水平的患者数和百分比
      # 求和: sum()
      # 返回数据框行数: nrow()
      # 四舍五入: round()
      # 格式化字符串: sprintf()，%f表示浮点数，%d表示整数，%%表示百分号
      # na.rm=TRUE表示在计算时忽略NA值
      n_patients <- sum(data[[var]] == level, na.rm = TRUE)
      pct_patients <- round(n_patients / nrow(data) * 100, 2)
      patients_str <- sprintf("%d (%.2f%%)", n_patients, pct_patients) # %d整数，%.2f%%保留2位小数+%符号

      if (i == 1) {

        # 第一个水平作为参考组
        # 连接字符串: paste0()，在水平前添加缩进以区分标题行和水平行
        result_rows <- append_row(result_rows, paste0("    ", level), patients_str, "Reference", "", NA_real_, NA_real_, NA_real_)

      }else {

        # 非参考组：从Cox模型结果中提取统计量
        coef_name <- paste0(var, level)

        # %in%是包含于运算符，rownames()获取行名
        if (coef_name %in% rownames(cox_coef)) {

          # 格式化P值
          # cox_coef包含了P值
          # 如果P值小于0.001则显示为<0.001，否则保留三位小数
          p_value <- cox_coef[coef_name, "Pr(>|z|)"]          
          p_str <- if (p_value < 0.001) "<0.001" else sprintf("%.3f", p_value)

          # 格式化CI
          # cox_conf包含了HR和置信区间
          # exp(coef)是指数化后的系数，即风险比HR
          # lower .95和upper .95是95%置信区间的下限和上限
          # sprintf()格式化字符串，%.2f保留两位小数
          # CI字符串格式为 "HR(lower_ci-upper_ci)"
          # round()函数用于四舍五入，确保表格中显示的HR和CI数值一致
          hr <- cox_conf[coef_name, "exp(coef)"] 
          lower_ci <- cox_conf[coef_name, "lower .95"]
          upper_ci <- cox_conf[coef_name, "upper .95"]
          ci_str <- sprintf("%.2f(%.2f-%.2f)", hr, lower_ci, upper_ci)
          result_rows <- append_row(result_rows, paste0("    ", level), patients_str, ci_str, p_str, round(hr, 2), round(lower_ci, 2), round(upper_ci, 2))
          
        }
      }
    }
  }
  
  # 无结果则返回空数据框
  if (length(result_rows) == 0) return(empty_result) 
  
  # do.call()调用函数并传入列表作为参数
  # rbind按行合并多个数据框
  result_table <- do.call(rbind, result_rows)

  # 移除行名，使用默认数字索引
  # rownames()获取或设置数据框的行名，NULL表示移除行名
  rownames(result_table) <- NULL 

  # 最后生成的结果表格包含了每个变量的标题行和各水平的统计结果，便于后续分析和报告撰写
  result_table

}

# 构建Cox回归公式
build_cox_formula <- function(predictors) {

  # paste() 将多个字符串连接成一个字符串，collapse参数指定连接符
  # 如: c("age", "sex") 变成 "age + sex"
  predictor_string <- paste(predictors, collapse = " + ")

  # Surv() 是survival包的函数，创建生存对象
  # time是生存时间，status==1表示事件发生（如: 死亡）
  # as.formula() 将字符串转换为R的公式对象，~ 是公式符号，左边是因变量，右边是自变量
  as.formula(paste("Surv(time, status == 1) ~", predictor_string)) 

}

# 执行Cox回归并提取结果
# deparse(substitute(data))用于获取传入参数的变量名字符串
# data = 数据框
# predictors = 预测变量向量
# data_name = 数据名称（默认使用传入数据的变量名）
extract_cox_outputs <- function(data, predictors, data_name = deparse(substitute(data))) { 

  # coxph() 是survival包的函数，执行Cox比例风险回归
  # data参数指定数据框
  model <- coxph(build_cox_formula(predictors), data = data)

  # 修改模型对象中保存的数据名称（用于后续引用）
  # as.name() 将字符串转换为R的名称对象
  # summary() 对模型对象进行汇总，提取详细的统计信息
  model$call$data <- as.name(data_name)
  model_summary <- summary(model)

  # 返回包含原始模型对象、系数表（含P值）和置信区间表（含HR）的列表，便于后续分析和表格生成
  # list() 创建命名列表
  # 可以通过$符号访问元素
  list( 
    model = model, # 原始模型对象
    coef = as.data.frame(model_summary$coefficients), # 系数表（包含P值）
    conf = as.data.frame(model_summary$conf.int) # 置信区间表（包含HR）
  )

}

# 分割数据集，制作训练集和测试集，并保存为CSV文件
# data = 原始数据框
# indices = 分割索引
# dataset_label = 数据集标签（训练/测试）
# output_dir = 输出目录
# output_filename = 输出文件名
build_dataset_split <- function(data, indices, dataset_label, output_dir, output_filename) {

  # 提取指定索引的数据行
  # drop=FALSE确保即使只有一行也返回数据框
  raw <- data[indices, , drop = FALSE]

  # 添加dataset列用于标识数据集来源（训练集或测试集），便于后续合并和分析  
  raw$dataset <- dataset_label
  
  # 将原始分割数据保存为CSV文件，便于检查和后续使用
  write.csv(raw, file = file.path(output_dir, output_filename), row.names = FALSE)

  # 返回原始数据和带有标签的数据
  list(
    raw = raw,
    ready = prepare_variables(raw)
  )

}

# 3+1绘图函数
# 统一处理“1张合并图 + 3张单时间点图”的输出流程，减少重复代码
# curves_getter用于从容器中提取指定时间点曲线（内部校准与外部验证结构不同）
# cohort_name和suffix用于构建输出文件名，确保文件命名规范且易于识别
# suffix参数为不同类型的校准图添加区分（如"_external"）
# save_combined_fun和save_single_fun分别是保存合并图和单时间点图的函数，提供灵活性以适应内部校准和外部验证的不同绘图需求
save_calibration_plots <- function( curves_container, cohort_name, suffix, save_combined_fun, save_single_fun, curves_getter) {

  combined_output <- file.path(
    calibration_output_dir,
    sprintf("calibration_%s_combined%s.png", cohort_name, suffix)
  )
  save_combined_fun(curves_container, combined_output)

  for (time_point in calibration_times) {
    time_key <- as.character(time_point)
    single_output <- file.path(
      calibration_output_dir,
      sprintf("calibration_%s_%dmonth%s.png", cohort_name, time_point, suffix)
    )
    save_single_fun(
      curves_getter(curves_container, time_key),
      time_point = time_point,
      output_file = single_output,
      color = calibration_colors[time_key]
    )
  }

}

### ============================================================
### 第四部分：数据读取与分割
### ============================================================

# file.path()根据操作系统自动构建路径，避免手动拼接错误
# read.csv()读取CSV文件，默认第一行为列名，返回数据框
data_file <- file.path(inputed_dir, "data_raw.csv")
data <- read.csv(data_file) 

# set.seed()设置随机种子确保结果可重复
# nrow()返回数据框的行数，此处为患者数量
# sample.int(n)生成1到n的随机排列，如sample.int(5)可能返回3 1 5 2 4
set.seed(412)
n <- nrow(data) 
indices <- sample.int(n) 

# 计算训练集大小（80%），floor()向下取整确保整数索引
# 根据随机索引分割数据集，前80%作为训练集，剩余20%作为测试集
# []表示数据框的行列索引，[a:b]表示从第a行到第b行
train_size <- floor(0.8 * n)
train_indices <- indices[1:train_size]
test_indices <- indices[(train_size + 1):n]

# 提取训练集与测试集并保存
train_split <- build_dataset_split(
  data = data,
  indices = train_indices, dataset_label = "Training",
  output_dir = output_dir, output_filename = "data_train.csv"
)
test_split <- build_dataset_split(
  data = data,
  indices = test_indices, dataset_label = "Validation",
  output_dir = output_dir, output_filename = "data_test.csv"
)

# build_dataset_split()函数返回一个列表，包含原始数据（raw）和准备好的数据（ready）
# raw是分割后的原始数据
# ready是经过prepare_variables()处理后的数据，数值编码转换为带标签的因子
train_raw <- train_split$raw
test_raw <- test_split$raw
train_ready <- train_split$ready
test_ready <- test_split$ready

# 合并训练集和测试集
# rbind()按行合并数据框
combined_raw <- rbind(train_raw, test_raw)
combined_ready <- rbind(train_ready, test_ready)

# 将合并后的原始数据保存为CSV文件，便于检查和后续使用
write.csv(combined_raw, file = file.path(output_dir, "data_combine.csv"), row.names = FALSE)

### ============================================================
### 第五部分：生成描述性统计表
### ============================================================

# descrTable()生成描述性统计表
# 公式格式：分组变量 ~ 要描述的变量
# .表示所有其他变量，-表示排除某些变量
# 此处按dataset分组，描述除ID、status、time外的所有变量
baseline_table <- descrTable(dataset ~ . - ID - status - time, data = combined_ready, show.all = TRUE)

# export2csv()将描述性统计表导出为CSV文件，便于检查和后续使用
# 区别于write.csv()，export2csv()专门处理descrTable对象，格式化输出更适合描述性统计表格
export2csv(baseline_table, file = file.path(output_dir, "baseline_table.csv"))

### ============================================================
### 第六部分：单因素Cox回归分析
### ============================================================

# 单因素分析：每次只分析一个变量对生存的影响
univariate_predictors <- c(
  "age", "sex", "race", "marital", "grade", "site",
  "t", "n", "bone", "brain", "liver", "lung",
  "multi", "sur", "rad", "che"
)

# 单个预测变量的单因素Cox回归 + 结果表生成
build_univariate_table <- function(predictor, data = train_ready) {

  # 执行单因素Cox回归并提取结果
  # extract_cox_outputs()函数执行Cox回归并返回模型对象、P值和置信区间表
  cox_outputs <- extract_cox_outputs(data = data, predictors = predictor)

  # 生成单因素Cox回归结果表
  # generate_cox_table()将这些统计结果整理成易读的表格格式
  # data参数传入原始数据框，variables参数传入当前预测变量，cox_coef和cox_conf分别传入P值和置信区间表
  generate_cox_table(
    data = data, variables = predictor,
    cox_coef = cox_outputs$coef, cox_conf = cox_outputs$conf
  )

}

# lapply(predictors, FUN): 对predictors向量中的每个元素调用FUN函数，返回一个列表，每个元素是一个变量的结果表
# 此处: 对每个单因素预测变量执行build_univariate_table()函数，生成单因素Cox回归结果表
univariate_tables <- lapply(univariate_predictors, build_univariate_table)

# do.call()调用rbind函数将列表中的所有数据框按行合并成一个完整的结果表
# rbind()按行合并数据框
univariate_table <- do.call(rbind, univariate_tables)

# 将单因素Cox回归结果表保存为CSV文件，便于检查和后续使用
# row.names=FALSE表示不保存行名
# na=""表示将NA值保存为空字符串
write.csv(univariate_table, file = file.path(output_dir, "univariate_results.csv"), row.names = FALSE, na = "")

### ============================================================
### 第七部分：多因素Cox回归分析
### ============================================================

# 多因素分析：同时分析多个变量对生存的影响（此处排除了multi变量）
multivariate_predictors <- c(
  "age", "sex", "race", "marital", "grade", "site",
  "t", "n", "bone", "brain", "liver", "lung",
  "sur", "rad", "che"
)

# 执行多因素Cox回归并提取结果
# extract_cox_outputs()函数执行Cox回归并返回模型对象、P值和置信区间表
multivariate_outputs <- extract_cox_outputs(data = train_ready, predictors = multivariate_predictors)

# 生成多因素Cox回归结果表
# generate_cox_table()将这些统计结果整理成易读的表格格式
# data参数传入原始数据框，variables参数传入当前预测变量，cox_coef和cox_conf分别传入P值和置信区间表
multivariate_table <- generate_cox_table(
  data = train_ready, variables = multivariate_predictors,
  cox_coef = multivariate_outputs$coef, cox_conf = multivariate_outputs$conf
)
write.csv(multivariate_table, file = file.path(output_dir, "multivariate_results.csv"), row.names = FALSE, na = "")

### ============================================================
### 第八部分：绘制森林图（Forest Plot）
### ============================================================

# 创建标签矩阵用于森林图文本显示
# rbind()按行合并
# as.matrix()将数据框转为矩阵
# [,1:4]选择前4列
labeltext <- rbind(
  c("Variable", "Patients", "HR", "P"), # 表头
  as.matrix(multivariate_table[, 1:4])
)

# 从多因素分析结果表中提取数据用于森林图绘制
lower_ci <- multivariate_table$LowerCI
upper_ci <- multivariate_table$UpperCI
hr <- multivariate_table$HR

# 识别标题行和参考组行，以便在森林图中进行特殊标记和分组
# nzchar(): 检查字符串是否非空
# trimws(): 去除首尾空白
# !: 逻辑非
# 此处观察到，Patients列为空的是标题行，因此通过检查Patients列是否为空来识别标题行和分组标题行
group_header_rows <- !nzchar(trimws(multivariate_table$Patients))

# c()函数创建一个逻辑向量
# 第一行（表头）标记为TRUE，其他行根据group_header_rows标记
is_summary <- c(TRUE, group_header_rows)

# build_hrzl_line_style()函数根据行号i构建水平线的样式配置
build_hrzl_line_style <- function(i, columns = 1:5, col = "black") {

  line_type <- if (i == 1) 1 else 2
  line_width <- if (i == 1) 2 else 1

  # gpar()是grid包的函数，用于定义图形参数
  # 线条类型（lty）: 1表示实线，2表示虚线
  # 线条宽度（lwd）: 1表示正常，2表示加粗
  # columns参数指定线条应用于哪些列（此处为森林图中的HR列）
  # 颜色（col）: 默认黑色
  gpar(
    lty = line_type, mlwd = line_width, col = col,
    columns = columns
  )

}

# 配置水平分隔线，在标题行和参考组行之间绘制水平线，以区分不同变量的结果
# which()返回满足条件的索引，+1L在标题行下一行画线
# hrzl_lines是一个列表，键名为行号，值为对应行的线条样式配置
# lapply()对group_start_lines中的每个行号调用build_hrzl_line_style()函数，生成对应的线条样式配置
# seq_along(group_start_lines)生成从1到group_start_lines长度的整数序列，作为lapply的索引
group_start_lines <- which(group_header_rows) + 1L 
hrzl_lines <- lapply(seq_along(group_start_lines), build_hrzl_line_style)

# 为列表元素命名，指定线条出现的行号，这样在绘制森林图时就知道在哪些行画水平线
names(hrzl_lines) <- as.character(group_start_lines)

draw_forestplot <- function() {

  # 绘制森林图
  # 除了下面的参数，还有以下参数: legend指定图例位置，title指定图形标题
  forestplot(
    # labeltext提供文本标签，使用之前创建的labeltext矩阵
    labeltext = labeltext,

    # mean、lower、upper提供HR和置信区间数据，使用多因素分析结果中的HR、LowerCI、UpperCI列
    # c(NA,hr)在HR向量前加NA对应表头行
    mean = c(NA, hr), lower = c(NA, lower_ci), upper = c(NA, upper_ci),

    # is.summary标记哪些行是总结行，使用之前定义的is_summary逻辑向量
    is.summary = is_summary,

    # hrzl_lines指定水平线样式，使用之前定义的hrzl_lines列表
    hrzl_lines = hrzl_lines,

    # boxsize调整箱线图的大小
    # ci.vertices.height调整置信区间端点的高度
    # lineheight调整行高
    # colgap调整列间距
    # unit()创建单位对象，lines表示行数，char表示字符宽度
    boxsize = 0.5, ci.vertices.height = 0.2, lineheight = unit(1, "lines"), colgap = unit(3, "char"),  

    # zero指定无效线位置（HR=1）
    # lwd.zero调整无效线宽度
    # lwd.ci调整置信区间线宽
    # lty.ci调整置信区间线型，"solid"表示实线
    zero = 1, lwd.zero = 2, lwd.ci = 1.5, lty.ci = "solid",

    # col指定颜色方案，fpColors()函数定义森林图中颜色方案
    # box表示箱线图颜色
    # summary表示总结行颜色
    # lines表示线条颜色
    # zero表示无效线颜色
    col = fpColors(box = "#66ccff", summary = "black", lines = "black", zero = "lightgray"),

    # txt_gp指定文本样式，fpTxtGp()函数定义森林图中文本样式
    # label调整标签文本大小
    # ticks调整刻度文本大小
    # xlab调整X轴标签文本大小
    # gpar()定义具体的图形参数，cex调整文本大小，0.60表示60%的默认大小，0.8表示80%的默认大小
    txt_gp = fpTxtGp(label = gpar(cex = 0.60), ticks = gpar(cex = 0.60), xlab = gpar(cex = 0.8)),

    # xlab指定X轴标签
    # lwd.xaxis调整X轴线宽
    # xticks指定X轴刻度
    # clip指定剪切范围
    xlab = "Hazard ratio", lwd.xaxis = 2, xticks = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3), clip = c(0.0, 5),

    # align指定文本对齐方式
    # graphwidth调整图形区域宽度，npc表示相对于整个图形设备的比例，0.35表示占35%
    # graph.pos指定图形位置
    # mar调整边距
    align = "l", graphwidth = unit(0.35, "npc"), graph.pos = 3, mar = unit(c(8, 8, 8, 8), "mm") 
  )

}

# 根据标签矩阵的行数动态调整森林图的尺寸，确保图形清晰且内容不拥挤
# nrow()返回矩阵的行数，此处为标签行数
# 根据经验公式计算图形高度，基础高度14英寸，每增加一行增加0.22英寸，外加2英寸的额外空间
# 图形宽度固定为10英寸，足够容纳文本和图形元素
row_count <- nrow(labeltext)
plot_height_in <- max(14, 0.22 * row_count + 2)
plot_width_in <- 10

# 保存森林图为PNG文件
# png()函数打开一个PNG图形设备，指定文件路径、图像尺寸和分辨率
# as.integer()将计算结果转换为整数，300 DPI是出版质量标准分辨率，确保图像清晰
# file.path()构建输出文件路径，"multivariate_forestplot.png"是输出文件名
# draw_forestplot()函数在打开的图形设备上绘制森林图
# dev.off()关闭图形设备完成文件保存（R图形系统标准流程：打开设备->绘图->关闭设备）
png(file.path(output_dir, "multivariate_forestplot.png"),
    width = as.integer(plot_width_in * 300), height = as.integer(plot_height_in * 300),
    res = 300)
draw_forestplot()
dev.off() 

### ============================================================
### 第九部分：绘制列线图（Nomogram）
### ============================================================

# 提取多因素Cox模型的线性预测器（linear predictors），用于后续列线图中示例患者的定位
# 线性预测器是Cox模型中每个患者的风险评分，基于模型中的变量和系数计算得出
multivariate_model <- multivariate_outputs$model
linear_predictors <- multivariate_model$linear.predictors

# 通过计算线性预测器的中位数，找到最接近中位风险水平的患者，作为列线图中的示例观测点
# 这样可以在列线图中直观地展示一个典型患者的风险评分和对应的变量值，增强图形的解释性和实用性
# abs()计算绝对值
# which.min()返回最小值索引
# stats::median()计算中位数，(::表示从stats包调用)
# na.rm=TRUE表示在计算中位数时忽略NA值
nomogram_observation_index <- which.min(abs(linear_predictors - stats::median(linear_predictors, na.rm = TRUE)))

# 提取该中位患者的变量值，作为列线图中的示例观测点
nomogram_observation <- train_ready[nomogram_observation_index, , drop = FALSE]

# 定义列线图绘制函数
nomogram_plot <- function(output_file, width_in = 13, height_in = 16, res = 300) {

  # 根据输入的图像尺寸和分辨率计算像素尺寸，确保输出图像清晰且符合预期大小
  width_px <- as.integer(width_in * res)
  height_px <- as.integer(height_in * res)

  # 保存当前默认图形设备设置，以便在函数退出时恢复，确保不会影响后续的图形输出
  # getOption("device")获取当前默认图形设备函数
  # on.exit()确保在函数退出时执行恢复操作，add=TRUE表示如果函数中已经有其他on.exit()调用，则将这个恢复操作添加到现有的退出操作列表中，而不是覆盖它
  # options(device = old_device)恢复默认图形设备设置
  old_device <- getOption("device")
  on.exit(options(device = old_device), add = TRUE)

  # regplot包会关闭当前设备并打开默认设备，因此临时修改默认设备为PNG输出
  # function(...)创建匿名函数，...表示接受任意参数
  # png()函数打开一个PNG图形设备，指定文件路径、图像尺寸、分辨率和背景色
  # bg="white"确保背景为白色，适合列线图的视觉效果
  options(device = function(...) {
    png(
      filename = output_file,
      width = width_px,
      height = height_px,
      res = res,
      bg = "white",
      ...
    )
  })

  # 绘制列线图
  regplot(

    # reg参数传入Cox模型对象，regplot会基于该模型的变量和系数绘制列线图
    reg = multivariate_model,
    # points参数控制是否显示点数刻度
    points = TRUE,
    # plots参数指定要绘制的图形类型
    # 此处为密度图和箱线图，帮助展示变量分布和异常值
    plots = c("density", "boxes"),
    # rank参数指定变量排序方式
    # 此处按标准差排序，使得变化范围较大的变量排在前面，增强图形的可读性
    rank = "sd",
    # observation参数指定示例患者的变量值，用于在图中标记该患者
    # 此处使用之前提取的nomogram_observation，展示一个典型患者在列线图中的位置
    observation = nomogram_observation,
    # failtime参数指定预测时间点，此处为18、12、6个月，列线图将展示这些时间点的生存概率或失败概率
    failtime = c(18, 12, 6),
    # prfail参数控制是否显示失败概率，此处设置为FALSE，表示显示生存概率而非失败概率
    prfail = FALSE,
    # showP参数控制是否显示P值
    showP = TRUE,
    # subticks参数控制是否显示子刻度，增强刻度的细分和可读性
    subticks = TRUE,
    # clickable参数控制是否生成交互式图形，此处设置为FALSE，生成静态图像，适合保存和打印
    clickable = FALSE,
    # boxcol是箱线图颜色
    # dencol是密度图颜色
    # obscol是示例患者点的颜色
    # spkcol是刻度线颜色
    boxcol = "#66CCFF", dencol = "#FFCC66", obscol = "#EE0000", spkcol = "#00AAFF",
    # droplines参数控制是否显示下垂线
    droplines = TRUE,
    # cexvars是变量标签文本大小
    # cexcats是类别标签文本大小
    # cexscales是刻度标签文本大小
    cexvars = 1.0, cexcats = 1.0, cexscales = 0.9
  )

  # 关闭图形设备保存文件
  dev.off() 

}

# 调用函数生成列线图
nomogram_plot(file.path(output_dir, "multivariate_nomogram.png"))

### ============================================================
### 第十部分：绘制校准曲线
### ============================================================

# calibration_times定义评估时间点（单位：月，numeric向量）
# calibration_colors定义每个时间点曲线颜色（命名字符向量，名称与时间点一一对应）
# calibration_corrected_pch设置bootstrap校正点形状（integer，4表示叉号）
# calibration_bootstrap_B设置bootstrap重采样次数（integer，L后缀表示整型字面量）
calibration_times <- c(6, 12, 18)
calibration_colors <- c("6" = "#66CCFF", "12" = "#EE0000", "18" = "#FFCC66")
calibration_corrected_pch <- 4
calibration_bootstrap_B <- 1000L

# 复用多因素Cox建模公式（formula对象），保证校准分析与主分析变量口径一致
calibration_formula <- build_cox_formula(multivariate_predictors)

# 创建校准曲线输出目录
calibration_output_dir <- file.path(output_dir, "calibration_plots")
dir.create(calibration_output_dir, recursive = TRUE, showWarnings = FALSE)

# 在执行表达式期间临时设置rms::datadist环境
# data是data.frame，必须包含后续cph()/calibrate()会使用的所有变量
# expr是表达式对象（language），会在函数内部被force()强制求值
# dd_name是character，写入.GlobalEnv的数据分布对象名称
# on.exit(add = TRUE)用于“退出即清理”，避免函数中途报错污染环境
with_datadist <- function(data, expr, dd_name = "dd_calibration") {

  # exists()用于检查全局环境中是否已经存在同名对象，避免覆盖用户可能已经定义的datadist对象
  # get()用于获取现有对象的值，以便在函数退出时恢复原状，保持环境整洁
  # inherits = FALSE表示只查当前环境，不向父环境查找
  had_existing_dd <- exists(dd_name, envir = .GlobalEnv, inherits = FALSE)
  if (had_existing_dd) old_dd <- get(dd_name, envir = .GlobalEnv)
  
  # options("datadist")返回命名列表，这里先备份旧值，退出时恢复
  old_options <- options("datadist")

  # assign()向全局环境写入datadist对象，供rms函数读取
  assign(dd_name, datadist(data), envir = .GlobalEnv)
  options(datadist = dd_name)

  # on.exit()中的表达式会在函数退出时执行（无论正常结束还是报错中断）
  on.exit({
    if (had_existing_dd) {
      assign(dd_name, old_dd, envir = .GlobalEnv)
    } else if (exists(dd_name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = dd_name, envir = .GlobalEnv)
    }
    options(datadist = old_options$datadist)
  }, add = TRUE)

  # 强制求值表达式，确保在设置datadist环境后执行校准分析
  force(expr)

}

# 对输入队列执行“重拟合 + bootstrap校准”
fit_internal_calibration <- function(data) {

  # 显式设置生存时间单位，确保rms绘图标签统一显示为month
  # units<-是“替换函数”，语法看起来像赋值，本质是调用`units<-`方法
  units(data$time) <- "month"

  # calibrate()分组大小（integer）
  # 此处按样本量1/3取整，属于常见经验设置
  m_value <- floor(nrow(data) / 3) 

  # 在执行表达式期间临时设置rms::datadist环境，确保cph()和calibrate()函数能够正确访问数据分布信息进行建模和校准
  with_datadist(data, {
    calibration_results <- lapply(calibration_times, function(time_point) {

      # cph()是rms包中的Cox拟合函数；x=TRUE/y=TRUE/surv=TRUE便于后续校准与预测
      fit <- cph(
        calibration_formula,
        data = data,
        x = TRUE,
        y = TRUE,
        surv = TRUE,
        time.inc = time_point
      )

      # calibrate()执行bootstrap校准：
      # cmethod = "KM"使用Kaplan-Meier估计实际结局
      # method = "boot"表示bootstrap重采样校正
      curve <- calibrate(
        fit,
        cmethod = "KM",
        method = "boot",
        u = time_point,
        m = m_value,
        B = calibration_bootstrap_B
      )

      list(fit = fit, curve = curve)

    })

    # names()把时间点转为字符键，方便按名称索引校准结果
    names(calibration_results) <- as.character(calibration_times)

    # 返回一个列表，键名为时间点字符串，值为包含拟合对象和校准曲线的列表
    calibration_results

  })

}

# 保存单个时间点的内部校准图
# calibration_curve：calibrate对象（含原始/校正曲线）
# time_point：时间点（月）
save_internal_single_plot <- function(calibration_curve, time_point, output_file, color) {

  # 创建PNG图形设备
  png(output_file, width = 2200, height = 2200, res = 300)

  # 设备关闭确保，无论函数如何退出（正常或异常），都能正确关闭图形设备，避免资源泄漏和后续图形输出问题
  on.exit(dev.off(), add = TRUE)

  # 绘制校准曲线
  # add=FALSE表示新图
  # conf.int=TRUE显示置信区间
  # subtitles=FALSE不显示子标题
  # riskdist=FALSE不显示风险分布图
  # lwd/lty/col设置线条宽度、类型和颜色
  # errbar.col设置误差条颜色
  # par.corrected设置校正曲线的参数
  # xlim/ylim设置坐标轴范围
  # xlab/ylab设置坐标轴标签，使用sprintf格式化字符串包含时间点信息
  plot(
    calibration_curve,
    add = FALSE,
    conf.int = TRUE,
    subtitles = FALSE,
    riskdist = FALSE,
    lwd = 2,
    lty = 1,
    col = color,
    errbar.col = color,
    par.corrected = list(col = color, lty = 1, lwd = 2, pch = calibration_corrected_pch),
    xlim = c(0.0, 1.0),
    ylim = c(0.0, 1.0),
    xlab = sprintf("Nomogram-Predicted %d-month OS probability", time_point),
    ylab = sprintf("Actual-Predicted %d-month OS probability", time_point)
  )

  # 绘制45度理想校准线
  abline(0, 1, lty = 3, lwd = 2, col = "black") 

}

# 将6/12/18个月内部校准曲线叠加到同一张图
# calibration_results：来自fit_internal_calibration()的返回值，包含各时间点的拟合对象和校准曲线
save_internal_combined_plot <- function(calibration_results, output_file) {

  # 创建PNG图形设备
  first_time <- as.character(calibration_times[1])
  png(output_file, width = 2200, height = 2200, res = 300)
  on.exit(dev.off(), add = TRUE)

  plot(
    calibration_results[[first_time]]$curve,
    add = FALSE,
    conf.int = TRUE,
    subtitles = FALSE,
    riskdist = FALSE,
    lwd = 2,
    lty = 1,
    col = calibration_colors[first_time],
    errbar.col = calibration_colors[first_time],
    par.corrected = list(
      col = calibration_colors[first_time], lty = 1, lwd = 2, pch = calibration_corrected_pch
    ),
    xlim = c(0.0, 1.0),
    ylim = c(0.0, 1.0),
    xlab = "Nomogram-Predicted 6-, 12-, and 18-month OS probability",
    ylab = "Actual-Predicted 6-, 12-, and 18-month OS probability",
    main = ""
  )

  # 循环绘制剩余时间点的校准曲线，add=TRUE表示在同一图上叠加
  for (time_point in calibration_times[-1]) {
    time_key <- as.character(time_point)
    plot(
      calibration_results[[time_key]]$curve,
      add = TRUE,
      conf.int = TRUE,
      subtitles = FALSE,
      riskdist = FALSE,
      lwd = 2,
      lty = 1,
      col = calibration_colors[time_key],
      errbar.col = calibration_colors[time_key],
      par.corrected = list(
        col = calibration_colors[time_key], lty = 1, lwd = 2, pch = calibration_corrected_pch
      )
    )
  }

  # 绘制45度理想校准线
  abline(0, 1, lty = 3, lwd = 2, col = "black")

  # 添加图例，说明不同颜色对应的时间点
  legend(
    x = 0.58,
    y = 0.20,
    legend = c("6-month survival", "12-month survival", "18-month survival"),
    col = unname(calibration_colors[as.character(calibration_times)]),
    lty = 1,
    lwd = 2,
    bty = "n"
  )

}

# 外部校准验证（训练集模型 -> 验证集评估）
# training_results：训练集各时间点的拟合结果（含fit）
# validation_data：验证集数据
# 返回list，名字为"6"/"12"/"18"，元素为val.surv对象
# Surv(time, status == 1)把“1=事件发生”显式转换为逻辑事件指示
compute_external_validation_curves <- function(training_results, validation_data) {

  external_curves <- lapply(calibration_times, function(time_point) {
    fit <- training_results[[as.character(time_point)]]$fit
    val.surv(
      fit = fit,
      newdata = validation_data,
      S = Surv(validation_data$time, validation_data$status == 1),
      u = time_point,
      method = "smoothkm"
    )
  })

  names(external_curves) <- as.character(calibration_times)
  external_curves

}

# 保存外部验证的单时间点校准图
# external_curve：val.surv对象
save_external_single_plot <- function(external_curve, time_point, output_file, color) {

  png(output_file, width = 2200, height = 2200, res = 300)
  on.exit(dev.off(), add = TRUE)

  plot(
    external_curve,
    add = FALSE,
    riskdist = FALSE,
    lwd = 2,
    lty = 1,
    col = color,
    xlim = c(0.0, 1.0),
    ylim = c(0.0, 1.0),
    xlab = sprintf("Nomogram-Predicted %d-month OS probability", time_point),
    ylab = sprintf("Actual-Predicted %d-month OS probability", time_point)
  )

  abline(0, 1, lty = 3, lwd = 2, col = "black")

}

# 将外部验证的6/12/18个月曲线绘制在同一张图
# external_curves：list，compute_external_validation_curves()返回值
save_external_combined_plot <- function(external_curves, output_file) {

  first_time <- as.character(calibration_times[1])
  png(output_file, width = 2200, height = 2200, res = 300)
  on.exit(dev.off(), add = TRUE)

  plot(
    external_curves[[first_time]],
    add = FALSE,
    riskdist = FALSE,
    lwd = 2,
    lty = 1,
    col = calibration_colors[first_time],
    xlim = c(0.0, 1.0),
    ylim = c(0.0, 1.0),
    xlab = "Nomogram-Predicted 6-, 12-, and 18-month OS probability",
    ylab = "Actual-Predicted 6-, 12-, and 18-month OS probability"
  )

  for (time_point in calibration_times[-1]) {
    time_key <- as.character(time_point)
    plot(
      external_curves[[time_key]],
      add = TRUE,
      riskdist = FALSE,
      lwd = 2,
      lty = 1,
      col = calibration_colors[time_key]
    )
  }

  abline(0, 1, lty = 3, lwd = 2, col = "black")
  legend(
    x = 0.58,
    y = 0.20,
    legend = c("6-month survival", "12-month survival", "18-month survival"),
    col = unname(calibration_colors[as.character(calibration_times)]),
    lty = 1,
    lwd = 2,
    bty = "n"
  )

}

# 训练集内部校准：
# fit_internal_calibration()先拟合模型并计算校准曲线
# save_calibration_plots()再统一输出“单时间点 + 合并图”
# curves_getter是函数参数（function类型），用于定义“如何从容器中取曲线”
training_internal_calibration <- fit_internal_calibration(train_ready)
save_calibration_plots(
  curves_container = training_internal_calibration,
  cohort_name = "training",
  suffix = "",
  save_combined_fun = save_internal_combined_plot,
  save_single_fun = save_internal_single_plot,
  curves_getter = function(x, time_key) x[[time_key]]$curve
)

# 验证集内部校准（重拟合）：
# 该步骤用于与外部验证结果做对照，观察“同队列拟合”与“外部外推”的差异
validation_internal_calibration <- fit_internal_calibration(test_ready)
save_calibration_plots(
  curves_container = validation_internal_calibration,
  cohort_name = "validation",
  suffix = "",
  save_combined_fun = save_internal_combined_plot,
  save_single_fun = save_internal_single_plot,
  curves_getter = function(x, time_key) x[[time_key]]$curve
)

# 外部验证校准（训练模型 -> 验证集）：
# compute_external_validation_curves()返回每个时间点的val.surv对象
# save_calibration_plots()负责批量导出图像
# curves_getter在这里直接返回x[[time_key]]，因为外部结果不再嵌套fit/curve二级结构
validation_external_curves <- compute_external_validation_curves(training_internal_calibration, test_ready)
save_calibration_plots(
  curves_container = validation_external_curves,
  cohort_name = "validation",
  suffix = "_external",
  save_combined_fun = save_external_combined_plot,
  save_single_fun = save_external_single_plot,
  curves_getter = function(x, time_key) x[[time_key]]
)

### ============================================================
### 第十一部分：时间依赖ROC曲线
### ============================================================

if (!requireNamespace("survivalROC", quietly = TRUE)) {
  stop("Package 'survivalROC' is required for time-dependent ROC analysis.")
}

# ROC分析固定使用6/12/18个月三个时间点，颜色与图例保持一致，结果输出到独立子目录
roc_time_points <- c(6, 12, 18)
roc_colors <- c("6" = "#66ccff", "12" = "#EE0000", "18" = "#FFCC66")
roc_output_dir <- file.path(output_dir, "roc_plots")
dir.create(roc_output_dir, recursive = TRUE, showWarnings = FALSE)

# 为ROC/DCA复用同一套Cox建模逻辑，避免和前文主模型对象相互覆盖
# data = 输入数据框（训练集或验证集）
# predictors = 建模变量，默认使用多因素分析筛选后的变量
# data_name = 写入model$call$data，便于后续函数在模型调用链中识别数据来源
fit_cox_model_for_roc_dca <- function(data, predictors = multivariate_predictors, data_name = deparse(substitute(data))) {
  model <- coxph(build_cox_formula(predictors), data = data)
  model$call$data <- as.name(data_name)
  model
}

# 逐个时间点计算时间依赖ROC，返回命名列表（名称为月份字符）
# data = 必须包含time和status列
# marker = 风险评分（线性预测值lp）
# time_points = ROC评估时间点向量（单位：月）
compute_time_roc_list <- function(data, marker, time_points = roc_time_points) {
  roc_list <- lapply(time_points, function(time_point) {
    roc_result <- tryCatch(
      survivalROC::survivalROC(
        Stime = data$time,
        status = data$status,
        marker = marker,
        predict.time = time_point,
        method = "KM"
      ),
      error = function(e) NULL
    )

    # 某一时间点计算失败时返回空曲线和NA AUC，保证后续绘图/汇总流程不中断
    if (is.null(roc_result)) {
      return(list(FP = numeric(), TP = numeric(), AUC = NA_real_))
    }

    roc_result
  })

  names(roc_list) <- as.character(time_points)
  roc_list
}

# 将多个时间点ROC曲线叠加在同一张图，图例统一显示AUC
# roc_results = compute_time_roc_list输出结果
# output_file = PNG输出路径
plot_roc_combined <- function(roc_results, output_file, time_points = roc_time_points, colors = roc_colors) {
  png(output_file, width = 2200, height = 2200, res = 300)
  on.exit(dev.off(), add = TRUE)

  plot(
    0, 0,
    type = "n",
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "1-Specificity",
    ylab = "Sensitivity",
    main = ""
  )

  # 仅绘制有效曲线（FP/TP非空）；无效时间点自动跳过
  for (time_point in time_points) {
    time_key <- as.character(time_point)
    roc_curve <- roc_results[[time_key]]
    if (is.null(roc_curve) || length(roc_curve$FP) == 0 || length(roc_curve$TP) == 0) next
    lines(roc_curve$FP, roc_curve$TP, col = colors[time_key], lty = 1, lwd = 2)
  }

  abline(0, 1, col = "black", lty = 2, lwd = 1)

  # 图例文本与曲线顺序严格按time_points一致，便于横向比较
  auc_labels <- sapply(time_points, function(time_point) {
    time_key <- as.character(time_point)
    auc_value <- roc_results[[time_key]]$AUC
    if (is.null(auc_value) || is.na(auc_value)) {
      sprintf("AUC at %d month = NA", time_point)
    } else {
      sprintf("AUC at %d month = %.3f", time_point, auc_value)
    }
  })

  legend(
    x = 0.60,
    y = 0.20,
    legend = unname(auc_labels),
    col = unname(colors[as.character(time_points)]),
    lty = 1,
    lwd = 2,
    cex = 0.8,
    bty = "n"
  )
}

# 将每个时间点AUC整理为长表，便于后续汇总、作图或与DCA结果联动分析
# strategy = 评估策略标签（refit/external）
# cohort = 数据集标签（training/validation）
collect_roc_auc_rows <- function(roc_results, strategy, cohort, time_points = roc_time_points) {
  data.frame(
    strategy = strategy,
    cohort = cohort,
    time_month = time_points,
    auc = sapply(time_points, function(time_point) {
      auc_value <- roc_results[[as.character(time_point)]]$AUC
      if (is.null(auc_value)) NA_real_ else as.numeric(auc_value)
    }),
    stringsAsFactors = FALSE
  )
}

# 针对训练集和验证集分别进行Cox模型重拟合，确保ROC分析的公平性和内部一致性
# 训练集重拟合模型用于评估训练集内部ROC和外部验证ROC；验证集重拟合模型用于评估验证集内部ROC
training_refit_model <- fit_cox_model_for_roc_dca(train_ready)
validation_refit_model <- fit_cox_model_for_roc_dca(test_ready)

# 计算线性预测值lp，作为时间依赖ROC的marker输入
training_refit_lp <- predict(training_refit_model, newdata = train_ready, type = "lp")
validation_refit_lp <- predict(validation_refit_model, newdata = test_ready, type = "lp")
validation_external_lp <- predict(training_refit_model, newdata = test_ready, type = "lp")

# 分别计算三套策略在6/12/18个月下的ROC曲线与AUC
roc_training_refit <- compute_time_roc_list(train_ready, training_refit_lp)
roc_validation_refit <- compute_time_roc_list(test_ready, validation_refit_lp)
roc_validation_external <- compute_time_roc_list(test_ready, validation_external_lp)

# 导出三张ROC图（训练重拟合/验证重拟合/外部验证）
plot_roc_combined(
  roc_results = roc_training_refit,
  output_file = file.path(roc_output_dir, "roc_training_refit.png")
)
plot_roc_combined(
  roc_results = roc_validation_refit,
  output_file = file.path(roc_output_dir, "roc_validation_refit.png")
)
plot_roc_combined(
  roc_results = roc_validation_external,
  output_file = file.path(roc_output_dir, "roc_validation_external.png")
)

# 合并并导出AUC汇总表，每行对应“策略-队列-时间点”
roc_auc_summary <- rbind(
  collect_roc_auc_rows(roc_training_refit, strategy = "refit", cohort = "training"),
  collect_roc_auc_rows(roc_validation_refit, strategy = "refit", cohort = "validation"),
  collect_roc_auc_rows(roc_validation_external, strategy = "external", cohort = "validation")
)
write.csv(roc_auc_summary, file = file.path(roc_output_dir, "roc_auc_summary.csv"), row.names = FALSE)

### ============================================================
### 第十二部分：决策曲线分析（DCA）
### ============================================================

# DCA评估时间点与阈值（与时间依赖ROC保持一致）
dca_time_points <- roc_time_points
dca_thresholds <- seq(0.01, 0.99, by = 0.01)
dca_loess_span <- 0.10
dca_output_dir <- file.path(output_dir, "dca_plots")
dir.create(dca_output_dir, recursive = TRUE, showWarnings = FALSE)

# 基于Cox模型在指定时间点预测事件风险（1-生存概率）
predict_event_risk_at_time <- function(model, data, time_point) {
  surv_values <- tryCatch(
    summary(survfit(model, newdata = data), times = time_point, extend = TRUE)$surv,
    error = function(e) rep(NA_real_, nrow(data))
  )

  if (is.matrix(surv_values)) surv_values <- surv_values[1, ]
  risk_values <- 1 - as.numeric(surv_values)

  if (length(risk_values) != nrow(data)) {
    warning(sprintf("Risk prediction length mismatch at %d months.", time_point))
    return(rep(NA_real_, nrow(data)))
  }

  risk_values
}

# 计算给定样本在目标时间点的KM事件发生率（1-KM生存率）
compute_km_event_rate <- function(time, status, time_point) {
  if (length(time) == 0) return(NA_real_)

  km_fit <- survfit(Surv(time, status == 1) ~ 1)
  surv_value <- tryCatch(
    summary(km_fit, times = time_point, extend = TRUE)$surv,
    error = function(e) NA_real_
  )

  if (length(surv_value) == 0 || is.na(surv_value[1])) return(NA_real_)
  1 - as.numeric(surv_value[1])
}

# 计算DCA曲线（Model / Treat all / Treat none）
compute_dca_curve <- function(data, pred_risk, time_point, thresholds = dca_thresholds, loess_span = dca_loess_span) {
  valid_rows <- complete.cases(data$time, data$status, pred_risk)
  dca_data <- data[valid_rows, , drop = FALSE]
  risk_values <- as.numeric(pred_risk[valid_rows])

  dca_table <- data.frame(
    threshold = thresholds,
    nb_model = NA_real_,
    nb_all = NA_real_,
    nb_none = 0,
    stringsAsFactors = FALSE
  )

  # 总体事件率用于计算Treat all策略净获益
  pd <- compute_km_event_rate(dca_data$time, dca_data$status, time_point)

  for (i in seq_along(thresholds)) {
    pt <- thresholds[i]

    if (!is.na(pd)) {
      dca_table$nb_all[i] <- pd - (1 - pd) * pt / (1 - pt)
    }

    high_risk <- risk_values > pt
    px <- mean(high_risk)

    if (is.na(px) || px == 0) next

    pd_given_x <- compute_km_event_rate(
      time = dca_data$time[high_risk],
      status = dca_data$status[high_risk],
      time_point = time_point
    )
    if (is.na(pd_given_x)) next

    dca_table$nb_model[i] <- pd_given_x * px - (1 - pd_given_x) * px * pt / (1 - pt)
  }

  # 对模型曲线做loess平滑；all/none保持原始折线
  dca_table$nb_model_smooth <- NA_real_
  valid_model_rows <- which(!is.na(dca_table$nb_model))
  if (length(valid_model_rows) >= 3) {
    loess_fit <- loess(nb_model ~ threshold, data = dca_table[valid_model_rows, , drop = FALSE], span = loess_span)
    dca_table$nb_model_smooth <- predict(loess_fit, newdata = data.frame(threshold = thresholds))
  }

  dca_table
}

# 绘制单张DCA图：Model / Treat all / Treat none
plot_dca_curve <- function(dca_table, time_point, output_file, cohort_label, strategy_label) {
  png(output_file, width = 2200, height = 2200, res = 300)
  on.exit(dev.off(), add = TRUE)

  # 固定DCA纵轴范围到[-1, 1]，让9张图的比例尺一致，方便横向比较
  # 同时手动设置0.2步长刻度，使坐标轴更规整
  y_limits <- c(-1, 1)
  y_ticks <- seq(-1, 1, by = 0.2)

  plot(
    dca_table$threshold, dca_table$nb_none,
    type = "n",
    xlim = c(min(dca_table$threshold), max(dca_table$threshold)),
    ylim = y_limits,
    xlab = "Threshold probability",
    ylab = "Net benefit",
    yaxt = "n",
    main = sprintf("%s (%s, %d-month)", cohort_label, strategy_label, time_point)
  )
  axis(2, at = y_ticks, labels = formatC(y_ticks, format = "f", digits = 1), las = 1)

  lines(dca_table$threshold, dca_table$nb_none, lwd = 2, col = "black")
  if (any(is.finite(dca_table$nb_all))) {
    lines(dca_table$threshold, dca_table$nb_all, lwd = 2, col = "#EE0000")
  }

  model_line <- if (any(is.finite(dca_table$nb_model_smooth))) dca_table$nb_model_smooth else dca_table$nb_model
  if (any(is.finite(model_line))) {
    lines(dca_table$threshold, model_line, lwd = 2, col = "#66CCFF")
  }

  legend(
    "topright",
    legend = c("None", "All", "Model"),
    col = c("black", "#EE0000", "#66CCFF"),
    lty = 1,
    lwd = 2,
    cex = 0.8,
    bty = "n"
  )
}

# 整理DCA净获益结果表（用于汇总导出）
collect_dca_rows <- function(dca_table, strategy, cohort, time_point) {
  data.frame(
    strategy = strategy,
    cohort = cohort,
    time_month = time_point,
    threshold = dca_table$threshold,
    net_benefit = dca_table$nb_model,
    stringsAsFactors = FALSE
  )
}

dca_rows <- list()

# 整理DCA净获益结果表（用于汇总导出）
for (time_point in dca_time_points) {
  training_refit_risk <- predict_event_risk_at_time(training_refit_model, train_ready, time_point)
  validation_refit_risk <- predict_event_risk_at_time(validation_refit_model, test_ready, time_point)
  validation_external_risk <- predict_event_risk_at_time(training_refit_model, test_ready, time_point)

  dca_training_refit <- compute_dca_curve(train_ready, training_refit_risk, time_point)
  dca_validation_refit <- compute_dca_curve(test_ready, validation_refit_risk, time_point)
  dca_validation_external <- compute_dca_curve(test_ready, validation_external_risk, time_point)

  plot_dca_curve(
    dca_table = dca_training_refit,
    time_point = time_point,
    output_file = file.path(dca_output_dir, sprintf("dca_%d_training_refit.png", time_point)),
    cohort_label = "Training cohort",
    strategy_label = "Refit"
  )
  plot_dca_curve(
    dca_table = dca_validation_refit,
    time_point = time_point,
    output_file = file.path(dca_output_dir, sprintf("dca_%d_validation_refit.png", time_point)),
    cohort_label = "Validation cohort",
    strategy_label = "Refit"
  )
  plot_dca_curve(
    dca_table = dca_validation_external,
    time_point = time_point,
    output_file = file.path(dca_output_dir, sprintf("dca_%d_validation_external.png", time_point)),
    cohort_label = "Validation cohort",
    strategy_label = "External"
  )

  dca_rows[[length(dca_rows) + 1L]] <- collect_dca_rows(dca_training_refit, strategy = "refit", cohort = "training", time_point = time_point)
  dca_rows[[length(dca_rows) + 1L]] <- collect_dca_rows(dca_validation_refit, strategy = "refit", cohort = "validation", time_point = time_point)
  dca_rows[[length(dca_rows) + 1L]] <- collect_dca_rows(dca_validation_external, strategy = "external", cohort = "validation", time_point = time_point)
}

dca_net_benefit_summary <- do.call(rbind, dca_rows)
write.csv(dca_net_benefit_summary, file = file.path(dca_output_dir, "dca_net_benefit_summary.csv"), row.names = FALSE)

### ============================================================
### 第十三部分：风险分层与生存曲线分析
### ============================================================

# predictors：用于构建 Cox 风险评分模型的 12 个变量
risk_strat_predictors <- c(
  "age", "marital", "site", "t", "n",
  "bone", "brain", "liver", "lung", "sur", "rad", "che"
)

# 根据训练集 Cox 模型的线性预测值分布，经验设定两个切点将患者分为三组风险等级
# cutoffs：总分分层阈值（<=674 为低危，675~758 为中危，>758 为高危）
risk_point_cutoffs <- c(674, 758)

# 输出目录分别保存分组结果表与 KM 生存曲线图
risk_output_dir <- file.path(output_dir, "risk_stratification")
survival_output_dir <- file.path(output_dir, "survival_plots")
dir.create(risk_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(survival_output_dir, recursive = TRUE, showWarnings = FALSE)

# 拟合风险分层 Cox 模型
fit_risk_strat_model <- function(data, predictors = risk_strat_predictors, data_name = deparse(substitute(data))) {
  model <- coxph(build_cox_formula(predictors), data = data)
  model$call$data <- as.name(data_name)
  model
}

# 从 regplot 输出中提取“变量取值 -> Points”映射表
# model：已拟合的 Cox 模型
# observation：用于生成 nomogram 标尺的代表性个体（1 行数据框）
# predictors：仅保留这批变量对应的映射，避免混入其他字段
extract_point_mapping_from_regplot <- function(model, observation, predictors = risk_strat_predictors) {
  temp_plot_file <- tempfile(fileext = ".png")
  # regplot 依赖图形设备，使用临时 PNG 设备并在退出时统一释放
  png(temp_plot_file, width = 2200, height = 2200, res = 300)
  on.exit({
    dev.off()
    if (file.exists(temp_plot_file)) unlink(temp_plot_file)
  }, add = TRUE)

  regplot_result <- regplot(
    reg = model,
    points = TRUE,
    plots = c("density", "boxes"),
    rank = "sd",
    observation = observation,
    failtime = c(18, 12, 6),
    prfail = FALSE,
    showP = TRUE,
    subticks = TRUE,
    clickable = FALSE
  )

  point_mapping <- list()
  mapping_tables <- regplot_result[seq_len(length(regplot_result) - 1L)] # 去掉末尾 Total Points 汇总行

  for (mapping_table in mapping_tables) {
    variable_name <- names(mapping_table)[1]
    if (!variable_name %in% predictors) next

    point_mapping[[variable_name]] <- data.frame(
      Value = trimws(as.character(mapping_table[[1]])),
      Points = as.numeric(mapping_table[[2]]),
      stringsAsFactors = FALSE
    )
  }

  point_mapping
}

# 按映射表计算每个样本的总 Point
# data：待评分数据（可为训练+验证合并队列）
# point_mapping：extract_point_mapping_from_regplot() 的输出
# predictors：按该顺序逐个变量累加分值
score_total_points <- function(data, point_mapping, predictors = risk_strat_predictors) {
  total_points <- rep(0, nrow(data))

  for (predictor in predictors) {
    if (!predictor %in% names(point_mapping)) {
      stop(sprintf("Point mapping missing for predictor: %s", predictor))
    }

    mapping_table <- point_mapping[[predictor]]
    mapping_vector <- setNames(mapping_table$Points, mapping_table$Value)
    predictor_values <- trimws(as.character(data[[predictor]]))
    predictor_points <- as.numeric(mapping_vector[predictor_values])

    if (anyNA(predictor_points)) {
      missing_values <- unique(predictor_values[is.na(predictor_points)])
      warning(sprintf(
        "Point mapping missing levels for predictor %s: %s",
        predictor, paste(missing_values, collapse = ", ")
      ))
      predictor_points[is.na(predictor_points)] <- 0 # 未匹配到映射的水平暂按 0 分处理
    }

    total_points <- total_points + predictor_points
  }

  total_points
}

# 按总 Point 划分低/中/高风险组
# point_vector：总分向量
# cutoffs：两个切点，输出三组有序因子（Low/Medium/High）
assign_risk_group_by_point <- function(point_vector, cutoffs = risk_point_cutoffs) {
  cut(
    point_vector,
    breaks = c(-Inf, cutoffs[1], cutoffs[2], Inf),
    labels = c("Low-risk", "Medium-risk", "High-risk"),
    right = TRUE,
    ordered_result = TRUE
  )
}

# 绘制 Kaplan-Meier 曲线并显示组间 P 值
# data：至少包含 time、status、risk_group 三列
# cohort_name：图例标题（队列名称）
# output_file：PNG 输出路径
plot_km_risk_group <- function(data, cohort_name, output_file) {
  data$risk_group <- droplevels(data$risk_group)
  group_counts <- table(data$risk_group)

  if (length(group_counts[group_counts > 0]) < 2) {
    warning(sprintf("Skip KM plot for %s: fewer than two non-empty risk groups.", cohort_name))
    return(invisible(NULL))
  }

  km_fit <- survfit(Surv(time, status == 1) ~ risk_group, data = data)
  km_plot <- ggsurvplot(
    km_fit,
    data = data,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    xlim = c(0, 120),
    xlab = "Time in months",
    ylim = c(0, 1),
    break.time.by = 24,
    ggtheme = theme_light(),
    risk.table.y.text.col = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = TRUE,
    ncensor.plot = FALSE,
    ncensor.plot.height = 0.25,
    conf.int.style = "step",
    surv.median.line = "none",
    legend.title = cohort_name,
    legend.labs = levels(data$risk_group)
  )

  png(output_file, width = 2600, height = 2200, res = 300)
  on.exit(dev.off(), add = TRUE)
  withCallingHandlers(
    print(km_plot),
    warning = function(w) {
      if (grepl(
        "Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.",
        conditionMessage(w),
        fixed = TRUE
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# 先在训练集拟合模型，再选取线性预测值最接近中位数的样本
# 该样本用于 regplot 打分标尺，避免极端个体导致分值刻度偏移
risk_strat_model <- fit_risk_strat_model(train_ready)
risk_lp <- risk_strat_model$linear.predictors
risk_observation_index <- which.min(abs(risk_lp - stats::median(risk_lp, na.rm = TRUE)))
risk_observation <- train_ready[risk_observation_index, , drop = FALSE]

# 将映射应用到合并队列（训练+验证），计算总分并完成风险分组
risk_point_mapping <- extract_point_mapping_from_regplot(risk_strat_model, risk_observation, predictors = risk_strat_predictors)
risk_combined <- combined_ready
risk_combined$Point <- score_total_points(risk_combined, risk_point_mapping, predictors = risk_strat_predictors)
risk_combined$risk_group <- assign_risk_group_by_point(risk_combined$Point, cutoffs = risk_point_cutoffs)

# 导出逐例分组结果：
# ID / dataset：样本标识与来源队列
# Point / risk_group：总分与风险层级
# time / status：生存结局字段（用于后续 KM 分析）
risk_points_grouped <- data.frame(
  ID = combined_raw$ID,
  dataset = combined_raw$dataset,
  Point = risk_combined$Point,
  risk_group = as.character(risk_combined$risk_group),
  time = combined_raw$time,
  status = combined_raw$status,
  stringsAsFactors = FALSE
)
write.csv(risk_points_grouped, file = file.path(risk_output_dir, "risk_points_grouped.csv"), row.names = FALSE)

# 统计各队列（全体/训练/验证）的风险组人数与占比
# 返回列：risk_group, n, cohort, prop
build_risk_group_counts <- function(data, cohort_name) {
  counts <- as.data.frame(table(data$risk_group), stringsAsFactors = FALSE)
  names(counts) <- c("risk_group", "n")
  counts$cohort <- cohort_name
  counts$prop <- if (sum(counts$n) > 0) counts$n / sum(counts$n) else NA_real_
  counts
}

risk_points_grouped$risk_group <- factor(
  risk_points_grouped$risk_group,
  levels = c("Low-risk", "Medium-risk", "High-risk"),
  ordered = TRUE
)

risk_counts_all <- build_risk_group_counts(risk_points_grouped, cohort_name = "All")
risk_counts_training <- build_risk_group_counts(
  subset(risk_points_grouped, dataset == "Training"),
  cohort_name = "Training"
)
risk_counts_validation <- build_risk_group_counts(
  subset(risk_points_grouped, dataset == "Validation"),
  cohort_name = "Validation"
)
risk_group_counts <- rbind(risk_counts_all, risk_counts_training, risk_counts_validation)
write.csv(risk_group_counts, file = file.path(risk_output_dir, "risk_group_counts.csv"), row.names = FALSE)

# 分别绘制全部、训练、验证队列的 KM 曲线并导出 PNG
plot_km_risk_group(
  data = risk_points_grouped,
  cohort_name = "All cohort",
  output_file = file.path(survival_output_dir, "km_risk_all.png")
)
plot_km_risk_group(
  data = subset(risk_points_grouped, dataset == "Training"),
  cohort_name = "Training cohort",
  output_file = file.path(survival_output_dir, "km_risk_training.png")
)
plot_km_risk_group(
  data = subset(risk_points_grouped, dataset == "Validation"),
  cohort_name = "Validation cohort",
  output_file = file.path(survival_output_dir, "km_risk_validation.png")
)
