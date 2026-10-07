# ==============================================================================
# 04_maxent_formal.R - MaxEnt正式训练
# ==============================================================================
# 目的：使用最优参数和全部数据进行30次Bootstrap正式训练
# 输入：
#   - data/swd/presence_raw_{species}_{resolution}.csv
#   - data/swd/background_raw_{species}_{resolution}.csv
#   - output/01_mainline/sdm_hsi/maxent/optimization/{species}_{resolution}/{species}_best_model.csv
# 输出：
#   - output/01_mainline/sdm_hsi/maxent/formal/{species}_{resolution}/
#   - output/01_mainline/sdm_hsi/maxent/formal/{species}_{resolution}/maxentResults.csv
# 依赖：config/config.R
# ==============================================================================

# 清空环境
rm(list = ls())

# 加载必要的包
suppressPackageStartupMessages({
  library(dismo)
  library(rJava)
  library(dplyr)
  library(parallel)
})

# 配置项目路径
if (!exists("PROJECT_ROOT")) {
  detect_root <- function() {
    path <- getwd()
    while (length(path) > 0 && path != dirname(path)) {
      if (dir.exists(file.path(path, "config")) && file.exists(file.path(path, "config", "config.R"))) {
        return(path)
      }
      path <- dirname(path)
    }
    return(getwd())
  }
  PROJECT_ROOT <- detect_root()
}
setwd(PROJECT_ROOT)

source(file.path(PROJECT_ROOT, "config/config.R"))

env_species_list <- Sys.getenv("MAXENT_SPECIES_LIST", unset = "")
if (nzchar(env_species_list)) {
  SPECIES_LIST <- trimws(strsplit(env_species_list, ",")[[1]])
  SPECIES_LIST <- SPECIES_LIST[nzchar(SPECIES_LIST)]
}

env_swd_dir <- Sys.getenv("MAXENT_SWD_DATA_DIR", unset = "")
if (nzchar(env_swd_dir)) {
  SWD_DATA_DIR <- if (grepl("^[A-Za-z]:|^/", env_swd_dir)) env_swd_dir else file.path(PROJECT_ROOT, env_swd_dir)
}

env_output_root <- Sys.getenv("MAXENT_OUTPUT_ROOT", unset = "")
if (nzchar(env_output_root)) {
  OUTPUT_DIR <- if (grepl("^[A-Za-z]:|^/", env_output_root)) env_output_root else file.path(PROJECT_ROOT, env_output_root)
  MAXENT_OUTPUT_DIR <- file.path(OUTPUT_DIR, "maxent")
  TEMPORAL_OUTPUT_DIR <- file.path(MAXENT_OUTPUT_DIR, "temporal")
}

env_regions <- Sys.getenv("MAXENT_REGIONS", unset = "")
if (nzchar(env_regions)) {
  keep_regions <- trimws(strsplit(env_regions, ",")[[1]])
  REGIONS <- REGIONS[intersect(names(REGIONS), keep_regions)]
}

# ---- Isolated Global 0.083-degree branch override ----
BRANCH_ROOT <- "D:/ecoMarine/output/01_mainline/inputs/sdm_global_0083_corrected_20260722"
SPECIES_LIST <- c("Acropora", "Lobophora", "Scarus", "Cephalopholis")
env_branch_species <- Sys.getenv("MAXENT_SPECIES_LIST", unset = "")
if (nzchar(env_branch_species)) {
  SPECIES_LIST <- trimws(strsplit(env_branch_species, ",")[[1]])
  SPECIES_LIST <- SPECIES_LIST[nzchar(SPECIES_LIST)]
}
REGIONS <- list()
SWD_DATA_DIR <- "D:/ecoMarine/output/01_mainline/inputs/sdm_global_0083_corrected_20260722/swd"
OUTPUT_DIR <- BRANCH_ROOT
MAXENT_OUTPUT_DIR <- file.path(OUTPUT_DIR, "maxent")
TEMPORAL_OUTPUT_DIR <- file.path(MAXENT_OUTPUT_DIR, "temporal")
GLOBAL_RESOLUTION <- 0.083
env_temporal_years <- Sys.getenv("MAXENT_TEMPORAL_YEARS", unset = "")
if (nzchar(env_temporal_years)) {
  TEMPORAL_YEARS <- as.integer(trimws(strsplit(env_temporal_years, ",")[[1]]))
  TEMPORAL_YEARS <- TEMPORAL_YEARS[!is.na(TEMPORAL_YEARS)]
}
get_final_vars <- function(species, region, mhw_window = NULL) {
  final_file <- file.path("D:/ecoMarine/output/01_mainline/inputs/sdm_global_0083_corrected_20260722", "vif", paste0(species, "_", region), "final_variables.txt")
  if (file.exists(final_file)) {
    vars <- readLines(final_file, warn = FALSE)
    return(vars[nzchar(vars)])
  }
  vars <- setdiff(ALL_ENV_VARS, EXCLUDE_GLOBAL)
  vars <- setdiff(vars, MHW_VARS)
  vars
}
# ---- End branch override ----



# MaxEnt jar 文件路径
MAXENT_JAR <- "D:/maxent/maxent/maxent.jar"
if (!file.exists(MAXENT_JAR)) {
  stop("MaxEnt jar not found: ", MAXENT_JAR)
}
# 初始化 Java 环境
# 硬件配置：i9 10核/20线程，128GB内存
if (is.null(getOption("java.parameters"))) {
  options(java.parameters = "-Xmx16g") # 主进程内存限制
}
.jinit()
.jaddClassPath(MAXENT_JAR)

# ==============================================================================
# 核心函数：运行正式训练
# ==============================================================================

run_formal_training <- function(species_name, resolution_type = "global",
                                region_name = NULL, n_threads = 4) {
  # 1. 确定路径和标识 ----
  if (resolution_type == "regional" && !is.null(region_name)) {
    run_id <- paste0(species_name, "_", region_name)
    pres_file <- file.path(SWD_DATA_DIR, paste0("presence_raw_", species_name, "_", region_name, ".csv"))
    bg_file <- file.path(SWD_DATA_DIR, paste0("background_raw_", species_name, "_", region_name, ".csv"))
  } else {
    run_id <- paste0(species_name, "_", resolution_type)
    pres_file <- file.path(SWD_DATA_DIR, paste0("presence_raw_", species_name, "_global.csv"))
    bg_file <- file.path(SWD_DATA_DIR, paste0("background_raw_", species_name, "_global.csv"))
  }

  opt_dir <- file.path(OUTPUT_DIR, "maxent/optimization", run_id)
  formal_dir <- file.path(OUTPUT_DIR, "maxent/formal", run_id)

  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("开始正式训练 (Bootstrap 10x):", run_id, "\n")
  cat("  时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

  # 2. 检查数据文件 ----
  if (!file.exists(pres_file) || !file.exists(bg_file)) {
    cat("  x 未找到数据文件，跳过\n")
    return(NULL)
  }

  # 2.5 断点续传：检查是否已完成训练 ----
  # 完整训练应包含: maxentResults.csv + 10个 species_*.lambdas 文件
  results_file <- file.path(formal_dir, "maxentResults.csv")
  n_bootstraps <- 10
  lambda_files <- file.path(formal_dir, paste0("species_", 0:(n_bootstraps - 1), ".lambdas"))

  training_complete <- FALSE
  if (file.exists(results_file)) {
    # 检查所有 lambda 文件是否存在
    existing_lambdas <- file.exists(lambda_files)
    n_existing <- sum(existing_lambdas)

    if (all(existing_lambdas)) {
      # 验证结果文件有效性
      tryCatch(
        {
          res <- read.csv(results_file)
          if ("Test.AUC" %in% names(res)) {
            avg_row <- res[res$Species == "species (average)", ]
            if (nrow(avg_row) > 0 && !is.na(avg_row$Test.AUC[1])) {
              training_complete <- TRUE
              cat("  ✓ 训练已完成，跳过\n")
              cat("    Lambda 文件:", n_existing, "/", n_bootstraps, "\n")
              cat("    平均测试 AUC:", round(avg_row$Test.AUC[1], 3), "\n")
              return(list(
                status = "skip",
                run_id = run_id,
                auc = avg_row$Test.AUC[1]
              ))
            }
          }
        },
        error = function(e) {
          cat("  ⚠ 结果文件损坏，将重新训练\n")
        }
      )
    } else {
      cat("  ⚠ 训练不完整 (", n_existing, "/", n_bootstraps, " lambdas)，将重新训练\n")
    }
  }

  # 3. 确定模型参数 (FC, RM) ----
  best_model_file <- file.path(opt_dir, paste0(species_name, "_best_model.csv"))

  # 确定区域键名
  region_key_for_params <- if (resolution_type == "global") "global" else region_name

  # 尝试获取手动参数
  manual_params <- get_manual_params(species_name, region_key_for_params)

  fc <- NULL
  rm <- NULL

  # 逻辑: 优先使用 optimize 的最优结果; 若无最优结果文件，则尝试手动参数
  if (file.exists(best_model_file)) {
    best_params <- read.csv(best_model_file)
    fc <- as.character(best_params$Best_FC[1])
    rm <- as.numeric(best_params$Best_RM[1])
    cat("  最优参数(自动优化): FC =", fc, ", RM =", rm, "\n")
  } else if (!is.null(manual_params)) {
    fc <- manual_params$fc
    rm <- manual_params$rm
    cat("  ! 使用手动参数(无优化结果): FC =", fc, ", RM =", rm, "\n")
  } else {
    cat("  x 未找到最优参数文件，且未配置手动参数，跳过\n")
    return(NULL)
  }

  # 4. 读取并准备数据 ----
  cat("  [1/3] 读取数据...\n")
  pres_data <- read.csv(pres_file)
  bg_data <- read.csv(bg_file)

  # ★ 关键：统一物种名称
  if ("species" %in% names(pres_data)) {
    pres_data$species <- species_name
  }
  if ("species" %in% names(bg_data)) {
    bg_data$species <- "background"
  }

  cat("    原始存在点:", nrow(pres_data), "\n")
  cat("    原始背景点:", nrow(bg_data), "\n")

  # 使用 get_final_vars 获取 VIF+重要性筛选后的最终变量集
  # 确定区域键名
  region_key <- if (resolution_type == "global") "global" else region_name

  cat("  [2/3] 变量准备...\n")

  if (exists("get_final_vars")) {
    selected_vars <- get_final_vars(species_name, region_key)
    if (is.null(selected_vars) || length(selected_vars) == 0) {
      cat("  ⚠ 未找到变量集文件，使用默认配置\n")
      # 回退逻辑：基础变量集，排除 MHW（不参与 MaxEnt 训练）
      selected_vars <- setdiff(ALL_ENV_VARS, EXCLUDE_GLOBAL)
      selected_vars <- setdiff(selected_vars, MHW_VARS)
    } else {
      cat(sprintf("    使用最终变量集: %d 个变量\n", length(selected_vars)))
    }
  } else {
    cat("  ⚠ get_final_vars 未定义，使用默认配置\n")
    # 回退逻辑：基础变量集，排除 MHW（不参与 MaxEnt 训练）
    selected_vars <- setdiff(ALL_ENV_VARS, EXCLUDE_GLOBAL)
    selected_vars <- setdiff(selected_vars, MHW_VARS)
  }

  # 强制排除 MHW 变量（Double Check）
  selected_vars <- setdiff(selected_vars, MHW_VARS)

  cat("    配置变量数:", length(selected_vars), "\n")

  # 检查变量可用性
  available_vars <- intersect(selected_vars, names(pres_data))

  cat("    可用变量数:", length(available_vars), "\n")

  if (length(available_vars) < 3) {
    cat("  x 可用变量过少，跳过\n")
    return(NULL)
  }

  # 提取环境变量（仅保留数值列，不包含元数据）
  pres_df <- pres_data[, available_vars, drop = FALSE]
  bg_df <- bg_data[, available_vars, drop = FALSE]

  # 移除含 NA 的记录
  pres_df <- pres_df[complete.cases(pres_df), ]
  bg_df <- bg_df[complete.cases(bg_df), ]

  cat("    清洗后存在点:", nrow(pres_df), "\n")
  cat("    清洗后背景点:", nrow(bg_df), "\n")
  cat("    使用变量:\n")
  for (v in available_vars) {
    cat("      -", v, "\n")
  }

  # 5. 创建输出目录 ----
  if (!dir.exists(formal_dir)) dir.create(formal_dir, recursive = TRUE)

  # 6. 运行 MaxEnt (Bootstrap) ----
  cat("  [3/3] 训练模型 (10次 Bootstrap)...\n")
  cat("    开始时间:", format(Sys.time(), "%H:%M:%S"), "\n")

  # 构造参数
  args <- c(
    paste0("betamultiplier=", rm),
    "replicates=10",
    "replicatetype=bootstrap",
    "responsecurves=true",
    "jackknife=false",
    "writeplotdata=true",
    "askoverwrite=false",
    "warnings=false",
    "tooltips=false",
    "randomseed=true",
    paste0("threads=", n_threads),
    # 外推与钳制设置（防止预测到训练范围外的极端值）
    "extrapolate=false", # 关闭外推：超出训练范围的环境值不进行外推预测
    "doclamp=true", # 开启钳制：将超出范围的值钳制到训练范围边界
    "fadebyclamping=false", # 策略调整：Fade=OFF (求稳)，禁用渐变
    "maximumbackground=-1", # 解除背景点数量限制（默认10000）
    # 特征类初始化（后续根据 FC 参数覆盖）
    "linear=false", "quadratic=false", "product=false",
    "threshold=false", # 强制禁用 Threshold：避免锯齿状响应曲线
    "hinge=false"
  )

  if (grepl("L", fc)) args[grep("linear=", args)] <- "linear=true"
  if (grepl("Q", fc)) args[grep("quadratic=", args)] <- "quadratic=true"
  if (grepl("P", fc)) args[grep("product=", args)] <- "product=true"
  if (grepl("T", fc)) args[grep("threshold=", args)] <- "threshold=true"
  if (grepl("H", fc)) args[grep("hinge=", args)] <- "hinge=true"

  tryCatch(
    {
      # 合并存在点和背景点（仅环境变量）
      train_data <- rbind(pres_df, bg_df)
      train_p <- c(rep(1, nrow(pres_df)), rep(0, nrow(bg_df)))
      cat("    训练数据总行数:", nrow(train_data), "\n")
      cat("    特征类: FC =", fc, ", 正则化: RM =", rm, "\n")

      gc() # 强制垃圾回收
      model <- maxent(x = train_data, p = train_p, args = args, path = formal_dir)

      cat("    完成时间:", format(Sys.time(), "%H:%M:%S"), "\n")
      cat("  ✓ 训练完成\n")
      cat("  输出目录:", formal_dir, "\n")

      # 简单检查结果文件
      if (file.exists(results_file)) {
        res <- read.csv(results_file)
        if ("Test.AUC" %in% names(res)) {
          avg_auc <- res$Test.AUC[res$Species == "species (average)"]
          if (length(avg_auc) > 0) {
            cat("  平均测试 AUC:", round(avg_auc, 3), "\n")
          }
        }
      }
    },
    error = function(e) {
      cat("  x 训练失败:", e$message, "\n")
      cat("    错误时间:", format(Sys.time(), "%H:%M:%S"), "\n")
    }
  )
}

# ==============================================================================
# 主程序
# ==============================================================================

main <- function() {
  # 1. 收集所有任务
  tasks <- list()

  # 遍历所有物种
  for (species in SPECIES_LIST) {
    # Global Models
    tasks[[length(tasks) + 1]] <- list(
      species = species,
      type = "global",
      region = NULL
    )

    # Regional Models
    for (region_key in names(REGIONS)) {
      check_file <- file.path(SWD_DATA_DIR, paste0("presence_raw_", species, "_", region_key, ".csv"))
      if (file.exists(check_file)) {
        tasks[[length(tasks) + 1]] <- list(
          species = species,
          type = "regional",
          region = region_key
        )
      }
    }
  }

  cat("总任务数:", length(tasks), "\n")
  cat("估计训练时间: ~", round(length(tasks) * 5 / 60, 1), "小时 (假设每个任务5分钟)\n")

  # 2. 并行设置
  # 硬件配置：i9 10核/20线程，128GB内存
  # 设计原则：R并行任务数(5) × MaxEnt线程数(4) = 20线程
  # 内存分配：5个进程 × 16GB = 80GB，留余量给系统
  n_cores <- suppressWarnings(as.integer(Sys.getenv("MAXENT_FORMAL_CORES", unset = "10")))
  if (is.na(n_cores) || n_cores < 1) n_cores <- 1
  cat("使用核心数:", n_cores, "\n")
  cat("\n", paste(rep("-", 70), collapse = ""), "\n\n")

  cl <- makeCluster(n_cores)

  # 3. 导出变量和函数
  clusterExport(cl, c(
    "PROJECT_ROOT", "SWD_DATA_DIR", "OUTPUT_DIR", "MAXENT_JAR",
    "run_formal_training", "get_manual_params", "get_final_vars",
    "SPECIES_LIST", "REGIONS", "ALL_ENV_VARS",
    "MAXENT_MANUAL_PARAMS", "BASE_ENV_VARS", "MHW_VARS",
    "EXCLUDE_GLOBAL"
  ),
  envir = environment()
  )

  # 4. 初始化 Worker 环境
  clusterEvalQ(cl, {
    library(dismo)
    library(rJava)
    library(dplyr)

    # 设置 Java 内存 (并行模式下，每个进程分配 16GB)
    # 硬件配置：i9 10核/20线程，128GB内存，4个并行进程各16GB
    options(java.parameters = "-Xmx16g")
    .jinit()
    .jaddClassPath(MAXENT_JAR)

    # 加载配置 (确保所有辅助函数可用)
    source(file.path(PROJECT_ROOT, "config/config.R"))

    # Rebind after config.R: config.R defines a mainline get_final_vars() and
    # would otherwise overwrite the isolated branch function in PSOCK workers.
    get_final_vars <- function(species, region, mhw_window = NULL) {
      final_file <- file.path("D:/ecoMarine/output/01_mainline/inputs/sdm_global_0083_corrected_20260722", "vif", paste0(species, "_", region), "final_variables.txt")
      if (file.exists(final_file)) {
        vars <- readLines(final_file, warn = FALSE)
        return(vars[nzchar(vars)])
      }
      stop("Branch final variable file not found: ", final_file)
    }

    env_species_list <- Sys.getenv("MAXENT_SPECIES_LIST", unset = "")
    if (nzchar(env_species_list)) {
      SPECIES_LIST <- trimws(strsplit(env_species_list, ",")[[1]])
      SPECIES_LIST <- SPECIES_LIST[nzchar(SPECIES_LIST)]
    }

    env_swd_dir <- Sys.getenv("MAXENT_SWD_DATA_DIR", unset = "")
    if (nzchar(env_swd_dir)) {
      SWD_DATA_DIR <- if (grepl("^[A-Za-z]:|^/", env_swd_dir)) env_swd_dir else file.path(PROJECT_ROOT, env_swd_dir)
    }

    env_output_root <- Sys.getenv("MAXENT_OUTPUT_ROOT", unset = "")
    if (nzchar(env_output_root)) {
      OUTPUT_DIR <- if (grepl("^[A-Za-z]:|^/", env_output_root)) env_output_root else file.path(PROJECT_ROOT, env_output_root)
      MAXENT_OUTPUT_DIR <- file.path(OUTPUT_DIR, "maxent")
      TEMPORAL_OUTPUT_DIR <- file.path(MAXENT_OUTPUT_DIR, "temporal")
    }

    env_regions <- Sys.getenv("MAXENT_REGIONS", unset = "")
    if (nzchar(env_regions)) {
      keep_regions <- trimws(strsplit(env_regions, ",")[[1]])
      REGIONS <- REGIONS[intersect(names(REGIONS), keep_regions)]
    }
  })

  # 5. 并行执行任务
  results <- parLapply(cl, tasks, function(task) {
    tryCatch(
      {
        res <- run_formal_training(
          species_name = task$species,
          resolution_type = task$type,
          region_name = task$region,
          n_threads = 1 # 并行模式下，MaxEnt 内部使用单线程，避免资源争抢
        )
        # 区分跳过和新训练
        if (!is.null(res) && is.list(res) && !is.null(res$status) && res$status == "skip") {
          return(paste("Skip:", task$species, task$type))
        }
        return(paste("Success:", task$species, task$type))
      },
      error = function(e) {
        return(paste("Error:", task$species, task$type, "-", e$message))
      }
    )
  })

  stopCluster(cl)

  # 6. 打印结果摘要
  cat("\n", paste(rep("#", 70), collapse = ""), "\n")
  cat("训练完成统计\n")
  cat("结束时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat(paste(rep("#", 70), collapse = ""), "\n\n")

  results_vec <- unlist(results)

  success_count <- sum(grepl("^Success:", results_vec))
  skip_count <- sum(grepl("^Skip:", results_vec))
  error_count <- sum(grepl("^Error:", results_vec))

  cat("新训练:", success_count, "\n")
  cat("已跳过:", skip_count, "(断点续传)\n")
  cat("失败:", error_count, "\n")
  cat("总计:", length(results_vec), "\n\n")

  if (error_count > 0) {
    cat("失败任务:\n")
    error_tasks <- results_vec[grepl("^Error:", results_vec)]
    for (err in error_tasks) {
      cat("  ", err, "\n")
    }
  }
}

main()
