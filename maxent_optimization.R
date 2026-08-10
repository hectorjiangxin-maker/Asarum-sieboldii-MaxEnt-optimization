# ============================================================
# ENMeval 2.x + maxent.jar 参数优化（改进版）
# ============================================================

# 清空R当前环境中的全部对象
rm(list = ls())
# 固定随机数种子，有利于论文复现
set.seed(20260731)

# -------------------- 0. 用户参数区 --------------------
WORK_DIR <- "/home/data/t050615/JX/VOC/"
OCC_FILE <- "Occ.csv"

# 分布点列号
SPECIES_COL <- 1
LON_COL <- 2
LAT_COL <- 3

# 环境文件读取，asc或者tif
ENV_PATTERN <- "\\.asc$"

# 分布点坐标读取
OCC_CRS <- "EPSG:4326"

# 创建输出目录
OUT_DIR <- file.path(WORK_DIR, "ENMeval_output_default_vs_optimized")

# 设定FC和RM的范围
FC_SET <- c("L", "LQ", "H", "LQH", "LQHP", "LQHPT")
RM_SET <- seq(0.5, 4.0, by = 0.5)

# 设置抽取背景点
N_BACKGROUND <- 10000

# 验证方法，block将分布点按空间位置分成4组，每轮使用3组训练、1组验证，空间block比随机75%/25%划分更严格。
# "block"、"jackknife"或"auto"
PARTITION_METHOD <- "block"
# block分成4组后，要求每一组至少有3个分布点
MIN_OCC_PER_FOLD <- 3

# 调参阶段不保存48个候选模型的完整预测地图。
RASTER_PREDS <- FALSE

# 并行计算设置
USE_PARALLEL <- FALSE
N_CORES <- max(1, parallel::detectCores(logical = TRUE) - 1)

# 若知道原稿的真实FC和RM，可填入以便直接比较
ORIGINAL_FC <- "LQHP"
ORIGINAL_RM <- 1


# -------------------- 1. 加载软件包 --------------------
required_packages <- c(
  "terra", "ENMeval", "predicts", "rJava",
  "dplyr", "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(terra)
  library(ENMeval)
  library(predicts)
  library(rJava)
  library(dplyr)
  library(ggplot2)
})

# -------------------- 2. 建立目录 --------------------
if (!dir.exists(WORK_DIR)) {
  stop("工作目录不存在：", WORK_DIR)
}

setwd(WORK_DIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("工作目录：", normalizePath(WORK_DIR), "\n")
cat("输出目录：", normalizePath(OUT_DIR), "\n")


# -------------------- 3. 读取环境变量 --------------------
# 找到工作目录下全部环境栅格
env_files <- list.files(
  path = WORK_DIR,
  pattern = ENV_PATTERN,
  full.names = TRUE,
  ignore.case = TRUE
)

# 至少需要两个环境变量
if (length(env_files) < 2) {
  stop(
    "找到的环境栅格少于2个，请检查ENV_PATTERN和工作目录。\n",
    "当前工作目录：", WORK_DIR, "\n",
    "当前匹配模式：", ENV_PATTERN
  )
}

cat("找到以下环境变量文件：\n")
print(basename(env_files))

# 每个环境文件先分别读取为一个SpatRaster
env_list <- lapply(
  env_files,
  terra::rast
)

# 比较第一个环境图层与其余所有环境图层
geom_ok <- terra::compareGeom(
  env_list[[1]],
  env_list[-1],
  crs = TRUE,
  ext = TRUE,
  rowcol = TRUE,
  res = TRUE,
  stopOnError = FALSE,
  messages = TRUE
)

if (!isTRUE(geom_ok)) {
  stop(
    "环境栅格的空间几何信息不一致。\n",
    "请检查各图层的坐标系、范围、分辨率、行列数和栅格原点。"
  )
}

# 几何信息一致后，组合成一个多层SpatRaster
envs <- do.call(
  c,
  env_list
)

# 使用文件名作为环境变量名
names(envs) <- make.names(
  tools::file_path_sans_ext(
    basename(env_files)
  ),
  unique = TRUE
)

# 检查是否读取成功
cat("\n环境变量读取完成：\n")
print(envs)

cat(
  "\n环境变量名称：",
  paste(names(envs), collapse = ", "),
  "\n"
)

cat(
  "环境变量数量：",
  terra::nlyr(envs),
  "\n"
)

cat(
  "坐标系：\n",
  terra::crs(envs, proj = TRUE),
  "\n"
)

cat(
  "分辨率：",
  paste(terra::res(envs), collapse = " × "),
  "\n"
)

cat(
  "空间范围：\n"
)
print(terra::ext(envs))

# -------------------- 4. 读取和清理分布点 --------------------
if (!file.exists(OCC_FILE)) {
  stop("找不到分布点文件：", file.path(WORK_DIR, OCC_FILE))
}

occ_raw <- read.csv(
  OCC_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_index <- max(c(SPECIES_COL, LON_COL, LAT_COL))

if (ncol(occ_raw) < required_index) {
  stop("分布点文件列数不足。")
}

occ_clean <- data.frame(
  species = occ_raw[[SPECIES_COL]],
  lon = suppressWarnings(as.numeric(occ_raw[[LON_COL]])),
  lat = suppressWarnings(as.numeric(occ_raw[[LAT_COL]]))
)

n_original <- nrow(occ_clean)

occ_clean <- occ_clean |>
  filter(
    is.finite(lon),
    is.finite(lat),
    lon >= -180,
    lon <= 180,
    lat >= -90,
    lat <= 90
  ) |>
  distinct(lon, lat, .keep_all = TRUE)

if (nrow(occ_clean) == 0) {
  stop("分布点清理后没有有效坐标。")
}

occ_vect <- terra::vect(
  occ_clean,
  geom = c("lon", "lat"),
  crs = OCC_CRS
)

if (!terra::same.crs(occ_vect, envs)) {
  occ_vect <- terra::project(occ_vect, terra::crs(envs))
}

occs <- as.data.frame(terra::crds(occ_vect))
names(occs) <- c("x", "y")

occ_cells <- terra::cellFromXY(envs[[1]], occs)
keep_inside <- !is.na(occ_cells)

occs <- occs[keep_inside, , drop = FALSE]
occ_cells <- occ_cells[keep_inside]

keep_unique_cell <- !duplicated(occ_cells)
occs <- occs[keep_unique_cell, , drop = FALSE]

occ_env <- terra::extract(envs, occs, ID = FALSE)
keep_complete <- stats::complete.cases(occ_env)

occs <- occs[keep_complete, , drop = FALSE]

if (nrow(occs) == 0) {
  stop("与环境图层匹配后没有有效分布点。")
}

cat("原始记录数：", n_original, "\n")
cat("最终有效分布点数：", nrow(occs), "\n")

write.csv(
  occs,
  file.path(OUT_DIR, "occurrences_cleaned_model_CRS.csv"),
  row.names = FALSE
)


# -------------------- 5. 抽取背景点 --------------------
complete_mask <- terra::app(
  envs,
  fun = function(v) {
    if (all(!is.na(v))) 1 else NA
  }
)

n_valid_cells <- terra::global(
  !is.na(complete_mask),
  fun = "sum",
  na.rm = TRUE
)[1, 1]

n_bg_actual <- min(N_BACKGROUND, as.integer(n_valid_cells))

bg <- terra::spatSample(
  complete_mask,
  size = n_bg_actual,
  method = "random",
  na.rm = TRUE,
  xy = TRUE,
  values = FALSE,
  replace = FALSE,
  as.df = TRUE
)

bg <- bg[, 1:2, drop = FALSE]
names(bg) <- c("x", "y")

if (nrow(bg) == 0) {
  stop("未能抽取有效背景点。")
}

cat("背景点数量：", nrow(bg), "\n")

write.csv(
  bg,
  file.path(OUT_DIR, "background_points.csv"),
  row.names = FALSE
)


# -------------------- 6. 交叉验证方式 --------------------
partition_method <- PARTITION_METHOD

if (PARTITION_METHOD == "auto") {
  partition_method <- if (nrow(occs) < 25) "jackknife" else "block"
}

partition_settings <- NULL
chosen_orientation <- NA_character_
block_groups <- NULL

if (partition_method == "block") {
  orientations <- c("lat_lon", "lon_lat", "lat_lat", "lon_lon")
  
  block_summary_list <- lapply(orientations, function(orientation_i) {
    group_i <- ENMeval::get.block(
      occs = occs,
      bg = bg,
      orientation = orientation_i
    )
    
    fold_n <- tabulate(group_i$occs.grp, nbins = 4)
    
    data.frame(
      orientation = orientation_i,
      fold1 = fold_n[1],
      fold2 = fold_n[2],
      fold3 = fold_n[3],
      fold4 = fold_n[4],
      min_occ_per_fold = min(fold_n),
      sd_occ_per_fold = stats::sd(fold_n)
    )
  })
  
  block_summary <- dplyr::bind_rows(block_summary_list) |>
    arrange(desc(min_occ_per_fold), sd_occ_per_fold)
  
  write.csv(
    block_summary,
    file.path(OUT_DIR, "block_orientation_comparison.csv"),
    row.names = FALSE
  )
  
  chosen_orientation <- block_summary$orientation[1]
  
  block_groups <- ENMeval::get.block(
    occs = occs,
    bg = bg,
    orientation = chosen_orientation
  )
  
  fold_counts <- tabulate(block_groups$occs.grp, nbins = 4)
  
  if (min(fold_counts) < MIN_OCC_PER_FOLD) {
    stop(
      "空间block分区后至少有一个验证折少于",
      MIN_OCC_PER_FOLD, "个分布点。各折数量：",
      paste(fold_counts, collapse = ", ")
    )
  }
  
  partition_settings <- list(
    orientation = chosen_orientation
  )
  
  cat("采用空间block分区：", chosen_orientation, "\n")
  cat("四个空间折点数：", paste(fold_counts, collapse = ", "), "\n")
  
  png(
    file.path(OUT_DIR, "spatial_block_partitions.png"),
    width = 1800,
    height = 1400,
    res = 220
  )
  
  plot(
    envs[[1]],
    main = paste0("Spatial block partitions: ", chosen_orientation)
  )
  
  points(
    occs$x,
    occs$y,
    pch = 21,
    bg = block_groups$occs.grp,
    col = "black",
    cex = 1.2
  )
  
  legend(
    "topright",
    legend = paste("Fold", 1:4),
    pt.bg = 1:4,
    pch = 21,
    bty = "n"
  )
  
  dev.off()
}

if (partition_method == "jackknife") {
  cat("采用jackknife留一法。最终分布点数：", nrow(occs), "\n")
}

if (!partition_method %in% c("block", "jackknife")) {
  stop("当前脚本只支持block、jackknife或auto。")
}


# -------------------- 7. 检查maxent.jar --------------------
maxent_available <- tryCatch(
  predicts::MaxEnt(),
  error = function(e) FALSE
)

if (!isTRUE(maxent_available)) {
  stop(
    "maxent.jar不可用。请检查Java、rJava、predicts及maxent.jar配置。"
  )
}


# -------------------- 8. 批量筛选FC和RM --------------------
validation_bg_setting <- if (
  partition_method == "block"
) "partition" else "full"

cat(
  "开始调参：",
  length(FC_SET), "种FC × ",
  length(RM_SET), "种RM = ",
  length(FC_SET) * length(RM_SET),
  "个候选参数组合。\n"
)

result <- ENMeval::ENMevaluate(
  occs = occs,
  envs = envs,
  bg = bg,
  
  tune.args = list(
    fc = FC_SET,
    rm = RM_SET
  ),
  
  partitions = partition_method,
  partition.settings = partition_settings,
  
  algorithm = "maxent.jar",
  
  other.settings = list(
    pred.type = "cloglog",
    validation.bg = validation_bg_setting,
    abs.auc.diff = TRUE
  ),
  
  doClamp = TRUE,
  raster.preds = RASTER_PREDS,
  
  parallel = USE_PARALLEL,
  numCores = N_CORES,
  
  taxon.name = as.character(occ_clean$species[1]),
  quiet = FALSE
)


# -------------------- 9. 保存候选模型结果 --------------------
results_table <- ENMeval::eval.results(result) |>
  arrange(delta.AICc, AICc)

partition_results <- ENMeval::eval.results.partitions(result)

write.csv(
  results_table,
  file.path(OUT_DIR, "all_candidate_models.csv"),
  row.names = FALSE
)

write.csv(
  partition_results,
  file.path(OUT_DIR, "candidate_models_by_partition.csv"),
  row.names = FALSE
)

ENMeval::saveENMevaluation(
  result,
  file.path(OUT_DIR, "ENMeval_tuning_object.rds")
)


# -------------------- 10. 选择最优模型 --------------------
valid_results <- results_table |>
  filter(is.finite(AICc), is.finite(delta.AICc))

if (nrow(valid_results) == 0) {
  stop(
    "所有候选模型的AICc均无法计算。常见原因是分布点过少，",
    "或候选模型的非零参数数目过多。"
  )
}

equivalent_models <- valid_results |>
  filter(delta.AICc <= 2) |>
  mutate(
    or_rank = ifelse(is.na(or.10p.avg), Inf, or.10p.avg),
    cbi_rank = ifelse(is.na(cbi.val.avg), -Inf, cbi.val.avg),
    aucdiff_rank = ifelse(is.na(auc.diff.avg), Inf, auc.diff.avg)
  ) |>
  arrange(
    delta.AICc,
    or_rank,
    desc(cbi_rank),
    aucdiff_rank,
    ncoef
  )

best_model <- equivalent_models |>
  slice(1) |>
  select(-or_rank, -cbi_rank, -aucdiff_rank)

write.csv(
  equivalent_models |>
    select(-or_rank, -cbi_rank, -aucdiff_rank),
  file.path(OUT_DIR, "models_deltaAICc_le_2.csv"),
  row.names = FALSE
)

write.csv(
  best_model,
  file.path(OUT_DIR, "selected_optimal_model.csv"),
  row.names = FALSE
)

cat("\n最优模型参数：\n")
print(best_model)


# -------------------- 11. 默认模型与优化模型比较 --------------------
# 说明：用户原稿中的默认参数为 LQPH + RM = 1。
# ENMeval 对特征组合使用规范顺序，因此在脚本中写作 LQHP；
# LQPH 与 LQHP 表示相同的四类特征：L、Q、H、P。

# 从48个候选模型中提取原稿默认模型
original_model <- results_table |>
  filter(
    toupper(fc) == toupper(ORIGINAL_FC),
    abs(rm - ORIGINAL_RM) < 1e-12
  ) |>
  slice(1)

if (nrow(original_model) == 0) {
  stop(
    "未在候选模型结果中找到原稿默认模型：FC = ", ORIGINAL_FC,
    ", RM = ", ORIGINAL_RM, "。\n",
    "请检查FC_SET是否包含LQHP，以及RM_SET是否包含1.0。"
  )
}

original_model <- original_model |>
  mutate(model_type = "Default")

optimized_model <- best_model |>
  mutate(model_type = "Optimized")

# 保存包含ENMeval全部评价字段的完整比较表
comparison_full <- bind_rows(
  original_model,
  optimized_model
) |>
  select(model_type, everything())

write.csv(
  comparison_full,
  file.path(OUT_DIR, "default_vs_optimized_full.csv"),
  row.names = FALSE
)

# 为兼容不同ENMeval 2.x小版本，自动识别常用评价指标列名
extract_metric <- function(df, candidate_names) {
  hit <- candidate_names[candidate_names %in% names(df)]
  if (length(hit) == 0) {
    return(rep(NA_real_, nrow(df)))
  }
  suppressWarnings(as.numeric(df[[hit[1]]]))
}

comparison_key <- data.frame(
  Model = comparison_full$model_type,
  FC = comparison_full$fc,
  RM = comparison_full$rm,
  Train_AUC = extract_metric(
    comparison_full,
    c("auc.train", "auc.train.avg")
  ),
  Test_AUC = extract_metric(
    comparison_full,
    c("auc.val.avg", "auc.val")
  ),
  AUCdiff = extract_metric(
    comparison_full,
    c("auc.diff.avg", "auc.diff")
  ),
  OR10 = extract_metric(
    comparison_full,
    c("or.10p.avg", "or.10p")
  ),
  CBI = extract_metric(
    comparison_full,
    c("cbi.val.avg", "cbi.val")
  ),
  AICc = extract_metric(
    comparison_full,
    c("AICc")
  ),
  delta_AICc = extract_metric(
    comparison_full,
    c("delta.AICc")
  ),
  AICc_weight = extract_metric(
    comparison_full,
    c("w.AIC", "AICc.wt")
  ),
  N_nonzero_coefficients = extract_metric(
    comparison_full,
    c("ncoef")
  ),
  stringsAsFactors = FALSE
)

write.csv(
  comparison_key,
  file.path(OUT_DIR, "default_vs_optimized_key_metrics.csv"),
  row.names = FALSE
)

capture.output(
  print(comparison_key, row.names = FALSE),
  file = file.path(OUT_DIR, "default_vs_optimized_key_metrics.txt")
)

cat("\n默认模型与优化模型关键指标比较：\n")
print(comparison_key, row.names = FALSE)

# 自动给出一个仅用于核查的变化方向摘要；论文中的文字仍需结合实际数值谨慎表述
if (nrow(comparison_key) == 2) {
  default_row <- comparison_key[comparison_key$Model == "Default", , drop = FALSE]
  optimized_row <- comparison_key[comparison_key$Model == "Optimized", , drop = FALSE]

  direction_text <- function(new_value, old_value, lower_is_better = FALSE) {
    if (!is.finite(new_value) || !is.finite(old_value)) return("NA")
    if (abs(new_value - old_value) < .Machine$double.eps^0.5) return("unchanged")
    if (lower_is_better) {
      if (new_value < old_value) "improved" else "increased"
    } else {
      if (new_value > old_value) "increased" else "decreased"
    }
  }

  comparison_summary <- data.frame(
    metric = c("Test_AUC", "AUCdiff", "OR10", "AICc", "N_nonzero_coefficients"),
    default_value = c(
      default_row$Test_AUC,
      default_row$AUCdiff,
      default_row$OR10,
      default_row$AICc,
      default_row$N_nonzero_coefficients
    ),
    optimized_value = c(
      optimized_row$Test_AUC,
      optimized_row$AUCdiff,
      optimized_row$OR10,
      optimized_row$AICc,
      optimized_row$N_nonzero_coefficients
    ),
    direction = c(
      direction_text(optimized_row$Test_AUC, default_row$Test_AUC, FALSE),
      direction_text(optimized_row$AUCdiff, default_row$AUCdiff, TRUE),
      direction_text(optimized_row$OR10, default_row$OR10, TRUE),
      direction_text(optimized_row$AICc, default_row$AICc, TRUE),
      direction_text(
        optimized_row$N_nonzero_coefficients,
        default_row$N_nonzero_coefficients,
        TRUE
      )
    )
  )

  write.csv(
    comparison_summary,
    file.path(OUT_DIR, "default_vs_optimized_change_summary.csv"),
    row.names = FALSE
  )
}


# -------------------- 12. 绘图 --------------------
p_delta <- ENMeval::evalplot.stats(
  e = result,
  stats = "delta.AICc",
  x.var = "rm",
  color.var = "fc",
  error.bars = FALSE
)

ggplot2::ggsave(
  filename = file.path(OUT_DIR, "delta_AICc.pdf"),
  plot = p_delta,
  width = 7,
  height = 5
)

p_metrics <- ENMeval::evalplot.stats(
  e = result,
  stats = c("or.10p", "cbi.val", "auc.diff"),
  x.var = "rm",
  color.var = "fc",
  error.bars = TRUE
)

ggplot2::ggsave(
  filename = file.path(
    OUT_DIR,
    "omission_CBI_AUCdiff.pdf"
  ),
  plot = p_metrics,
  width = 9,
  height = 8
)


# -------------------- 13. 保存运行信息 --------------------
run_summary <- data.frame(
  item = c(
    "Occurrence_original",
    "Occurrence_final",
    "Background_n",
    "Partition_method",
    "Block_orientation",
    "FC_candidates",
    "RM_candidates",
    "Raster_predictions_during_tuning",
    "Default_FC",
    "Default_RM",
    "Best_FC",
    "Best_RM"
  ),
  value = c(
    n_original,
    nrow(occs),
    nrow(bg),
    partition_method,
    chosen_orientation,
    paste(FC_SET, collapse = ","),
    paste(RM_SET, collapse = ","),
    RASTER_PREDS,
    ORIGINAL_FC,
    ORIGINAL_RM,
    best_model$fc[1],
    best_model$rm[1]
  )
)

write.csv(
  run_summary,
  file.path(OUT_DIR, "run_summary.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(OUT_DIR, "sessionInfo.txt")
)

cat(
  "\n模型调参完成。\n",
  "全部候选模型：all_candidate_models.csv\n",
  "最优参数：selected_optimal_model.csv\n",
  "默认与优化模型完整比较：default_vs_optimized_full.csv\n",
  "默认与优化模型关键指标：default_vs_optimized_key_metrics.csv\n",
  "输出目录：", OUT_DIR, "\n"
)
