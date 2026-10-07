# ==============================================================================
# 05_maxent_projection.R - 统一投影脚本（时序 + 未来情景 + MESS）
# ==============================================================================
# 目的：
#   1) 基于正式 MaxEnt 模型执行统一投影流程
#   2) 以逐年时序投影为核心（1993-2022），生成年度 HSI 序列
#   3) 从时序结果导出 temporal_mean（替代原 baseline）
#   4) 生成未来情景投影（ssp245_2050 / ssp585_2050）
#   5) 生成 MESS（全变量 + 非热变量）与二值掩膜
#   6) 统一输出 bootstrap 聚合结果（均值 + SD）
#
# 输入：
#   - output/01_mainline/sdm_hsi/maxent/formal/{species}_{scope}/species_*.lambdas
#   - data/processed/temporal_scenarios_shared/temporal_{year}_{resolution}/
#   - data/processed/scenarios_shared/{scenario}_{resolution}/
#   - data/swd/background_raw_{species}_{scope}.csv
#
# 输出：
#   - output/01_mainline/sdm_hsi/maxent/temporal/{species}_{scope}/hsi_{year}.tif
#   - output/01_mainline/sdm_hsi/maxent/temporal/{species}_{scope}/hsi_{year}_unconstrained.tif
#   - output/01_mainline/sdm_hsi/maxent/temporal/{species}_{scope}/hsi_{year}_sd.tif
#   - output/01_mainline/sdm_hsi/maxent/temporal/{species}_{scope}/temporal_mean.tif (时序均值)
#   - output/01_mainline/sdm_hsi/maxent/temporal/{species}_{scope}/temporal_mean_unconstrained.tif
#   - output/01_mainline/sdm_hsi/maxent/temporal/{species}_{scope}/temporal_sd.tif (年际SD)
#   - output/01_mainline/sdm_hsi/maxent/projection/{species}_{scope}/{scenario}_hsi_mean.tif
#   - output/01_mainline/sdm_hsi/maxent/projection/{species}_{scope}/{scenario}_hsi_mean_unconstrained.tif
#   - output/01_mainline/sdm_hsi/maxent/projection/{species}_{scope}/{scenario}_hsi_sd.tif
#   - output/01_mainline/sdm_hsi/maxent/mess/{run_id}/mess_*.tif
#   - output/01_mainline/sdm_hsi/maxent/mess/{run_id}/mess_*_binary.tif
#   - output/01_mainline/sdm_hsi/maxent/mess/all_mess_summary.csv
#   - output/01_mainline/sdm_hsi/maxent/projection/all_projection_summary.csv
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(terra)
  library(raster)
  library(dplyr)
  library(parallel)
})

# ------------------------------------------------------------------------------
# 0) 项目路径与环境初始化
# ------------------------------------------------------------------------------
if (!exists("PROJECT_ROOT")) {
  detect_root <- function() {
    path <- getwd()
    while (length(path) > 0 && path != dirname(path)) {
      if (dir.exists(file.path(path, "config")) &&
          file.exists(file.path(path, "config", "config.R"))) {
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
TEMPORAL_OVERRIDE_ROOT <- Sys.getenv("MAXENT_TEMPORAL_OVERRIDE_ROOT", unset = "")
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


if (!exists("SCENARIOS_SHARED_DIR")) {
  SCENARIOS_SHARED_DIR <- file.path(PROCESSED_DATA_DIR, "scenarios_shared")
}

MAXENT_JAR <- "D:/maxent/maxent/maxent.jar"
if (!file.exists(MAXENT_JAR)) {
  stop("MaxEnt jar not found: ", MAXENT_JAR)
}

get_default_workers <- function() {
  # Windows 上 detectCores(logical=FALSE) 常返回 NA，需回退到逻辑核心数
  cpu_physical <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.na(cpu_physical) || cpu_physical < 1) {
    cpu_logical <- suppressWarnings(parallel::detectCores(logical = TRUE))
    cpu_physical <- if (!is.na(cpu_logical) && cpu_logical >= 2) as.integer(cpu_logical / 2) else 10L
  }
  cat(sprintf("  检测到物理核心数: %d\n", cpu_physical))
  max(1L, min(8L, as.integer(cpu_physical - 2L)))
}

MAXENT_PARALLEL <- tolower(Sys.getenv("MAXENT_PARALLEL", unset = "true")) %in% c("1", "true", "yes")
MAXENT_WORKERS <- suppressWarnings(as.integer(Sys.getenv("MAXENT_WORKERS", unset = NA)))
if (is.na(MAXENT_WORKERS) || MAXENT_WORKERS < 1) {
  MAXENT_WORKERS <- get_default_workers()
}
MAXENT_JAVA_XMX_GB <- suppressWarnings(as.integer(Sys.getenv("MAXENT_JAVA_XMX_GB", unset = "7")))
if (is.na(MAXENT_JAVA_XMX_GB) || MAXENT_JAVA_XMX_GB < 2) {
  MAXENT_JAVA_XMX_GB <- 7L
}

if (!MAXENT_PARALLEL) {
  MAXENT_WORKERS <- 1L
}

master_java_heap_gb <- if (MAXENT_WORKERS > 1) MAXENT_JAVA_XMX_GB else 16L
options(java.parameters = c(
  sprintf("-Xmx%dg", master_java_heap_gb),
  "-XX:+UseG1GC",
  "-XX:ParallelGCThreads=4"
))

suppressPackageStartupMessages({
  library(rJava)
  library(dismo)
})

try(.jinit(), silent = TRUE)
.jaddClassPath(MAXENT_JAR)
terra::terraOptions(memfrac = 0.7, threads = if (MAXENT_WORKERS > 1) 1 else 4)

PROJECTION_DIR <- file.path(MAXENT_OUTPUT_DIR, "projection")
MESS_DIR <- file.path(MAXENT_OUTPUT_DIR, "mess")
if (!dir.exists(PROJECTION_DIR)) dir.create(PROJECTION_DIR, recursive = TRUE)
if (!dir.exists(MESS_DIR)) dir.create(MESS_DIR, recursive = TRUE)
if (!dir.exists(TEMPORAL_OUTPUT_DIR)) dir.create(TEMPORAL_OUTPUT_DIR, recursive = TRUE)

env_scenarios <- Sys.getenv("MAXENT_SCENARIOS", unset = "ssp126_2050,ssp245_2050,ssp585_2050")
if (tolower(trimws(env_scenarios)) %in% c("none", "temporal_only", "temporal-only")) {
  SCENARIOS <- character(0)
} else {
  SCENARIOS <- trimws(strsplit(env_scenarios, ",")[[1]])
  SCENARIOS <- SCENARIOS[nzchar(SCENARIOS)]
}
DEFAULT_MHW_WINDOW <- MAXENT_CANONICAL_WINDOW

# ------------------------------------------------------------------------------
# 跳过本地缓存设置
# SKIP_EXISTING_SCENARIOS = TRUE: 若某情景的 hsi_mean.tif 已存在则跳过该情景投影
# SKIP_EXISTING_TEMPORAL  = TRUE: 若 temporal_mean.tif 已存在则跳过时序投影
# 可通过环境变量覆盖: MAXENT_SKIP_EXISTING=false Rscript 05_maxent_projection.R
# ------------------------------------------------------------------------------
SKIP_EXISTING_SCENARIOS <- tolower(Sys.getenv("MAXENT_SKIP_EXISTING", unset = "true")) %in%
  c("1", "true", "yes")
SKIP_EXISTING_TEMPORAL <- tolower(Sys.getenv("MAXENT_SKIP_TEMPORAL", unset = "true")) %in%
  c("1", "true", "yes")

get_run_id <- function(species_name, scope) {
  paste0(species_name, "_", scope)
}

get_mess_run_id <- function(species_name, scope, mhw_window = DEFAULT_MHW_WINDOW) {
  paste0(species_name, "_", scope, "_mhw", mhw_window, "d")
}

get_resolution_type <- function(scope) {
  if (scope == "global") "regional" else "regional"
}

get_region_extent <- function(scope) {
  if (scope == "global") return(NULL)
  if (!scope %in% names(REGIONS)) return(NULL)
  rb <- REGIONS[[scope]]
  ext(rb$xmin, rb$xmax, rb$ymin, rb$ymax)
}

safe_stat_mean <- function(r) {
  out <- tryCatch(global(r, "mean", na.rm = TRUE)[[1]], error = function(e) NA_real_)
  as.numeric(out)
}

safe_stat_sd <- function(r) {
  out <- tryCatch(global(r, "sd", na.rm = TRUE)[[1]], error = function(e) NA_real_)
  as.numeric(out)
}

# ------------------------------------------------------------------------------
# 1) 生理约束（inline，用于时序投影）
# ------------------------------------------------------------------------------
compute_physio_coral <- function(thetao_rast, omega_rast = NULL, mhw_rast = NULL) {
  pp <- PHYSIO_PARAMS$coral

  # f(T): 不对称高斯 TPC — 慢性温度适宜性
  sigma <- ifel(thetao_rast <= pp$T_opt, pp$sigma_L, pp$sigma_R)
  f_t <- exp(-((thetao_rast - pp$T_opt)^2) / (2 * sigma^2))

  phi <- f_t

  # g(Ω): Sigmoid 碳酸钙约束 — 酸化胁迫
  if (!is.null(omega_rast)) {
    g_omega <- 1 / (1 + exp(-pp$k_omega * (omega_rast - pp$omega_mid)))
    g_omega <- clamp(g_omega, lower = 0, upper = 1)
    phi <- phi * g_omega
  }

  clamp(phi, lower = 0, upper = 1)
}

compute_physio_algae <- function(thetao_rast) {
  pp <- PHYSIO_PARAMS$algae
  t_k <- thetao_rast + 273.15

  arrhenius <- exp((pp$E_a / pp$k_B) * (1 / pp$T_ref - 1 / t_k))
  deactivation <- 1 + exp((pp$E_h / pp$k_B) * (1 / pp$T_h - 1 / t_k))
  r_t <- arrhenius / deactivation

  t_opt_k <- pp$T_h * (pp$E_h - pp$E_a) /
    (pp$E_h - pp$E_a + pp$k_B * pp$T_h * log(pp$E_h / pp$E_a))
  arr_opt <- exp((pp$E_a / pp$k_B) * (1 / pp$T_ref - 1 / t_opt_k))
  deact_opt <- 1 + exp((pp$E_h / pp$k_B) * (1 / pp$T_h - 1 / t_opt_k))
  r_opt <- arr_opt / deact_opt

  phi <- r_t / r_opt
  clamp(phi, lower = 0, upper = 1)
}

# --- Thermal Niche Correction (TNC) for algae ---
# 用文献 TPC 替换 MaxEnt 因采样偏差导致的过窄热维度响应
# 方法: α(T) = TPC_lit(T) / TPC_MaxEnt_norm(T)
# 依据: Kearney & Porter (2009) Ecol Lett; Bush et al. (2018) Glob Ecol Biogeogr
load_maxent_thetao_response <- function(species_name, scope) {
  dat_path <- file.path(OUTPUT_DIR, "maxent/formal",
                        paste0(species_name, "_", scope),
                        "plots/species_thetao.dat")
  if (!file.exists(dat_path)) {
    cat(sprintf("    TNC: MaxEnt thetao response not found: %s\n", dat_path))
    return(NULL)
  }
  dat <- read.csv(dat_path)
  list(
    interp = approxfun(dat$x, dat$y, rule = 2),
    peak = max(dat$y)
  )
}

compute_tnc_algae <- function(thetao_rast, species_name, scope) {
  pp <- PHYSIO_PARAMS$algae
  if (!isTRUE(pp$tnc_enabled)) return(NULL)

  # 1) 文献 TPC (Sharpe-Schoolfield with T_h=35°C)
  phi_lit <- compute_physio_algae(thetao_rast)

  # 2) MaxEnt 边际 thetao 响应 (归一化到 [0,1])
  maxent_resp <- load_maxent_thetao_response(species_name, scope)
  if (is.null(maxent_resp)) return(NULL)

  thetao_vals <- values(thetao_rast)
  maxent_norm <- maxent_resp$interp(thetao_vals) / maxent_resp$peak
  maxent_norm <- pmax(maxent_norm, pp$tnc_floor)

  # 3) 修正因子 α = TPC_lit / TPC_MaxEnt_norm
  phi_lit_vals <- values(phi_lit)
  alpha <- phi_lit_vals / maxent_norm
  alpha <- pmax(alpha, 0)
  alpha <- pmin(alpha, pp$tnc_cap)

  alpha_rast <- thetao_rast
  values(alpha_rast) <- alpha
  cat(sprintf("    TNC: α range [%.3f, %.3f], mean=%.3f\n",
              min(alpha, na.rm = TRUE), max(alpha, na.rm = TRUE),
              mean(alpha, na.rm = TRUE)))
  alpha_rast
}

compute_physio_fish <- function(thetao_rast, o2_rast, species_name) {
  pp <- PHYSIO_PARAMS$fish
  t_k <- thetao_rast + 273.15

  t_opt_mmr <- switch(species_name,
    Scarus = pp$T_opt_MMR_Scarus,
    Cephalopholis = pp$T_opt_MMR_Cephalopholis,
    pp$T_opt_MMR_default
  )
  fas <- switch(species_name,
    Scarus = pp$FAS_Scarus,
    Cephalopholis = pp$FAS_Cephalopholis,
    pp$FAS_default
  )
  gamma <- switch(species_name,
    Scarus = pp$gamma_Scarus,
    Cephalopholis = pp$gamma_Cephalopholis,
    pp$gamma_default
  )

  smr <- exp((pp$E_smr / pp$k_B) * (1 / pp$T_ref - 1 / t_k))
  t_opt_k <- t_opt_mmr + 273.15
  smr_at_opt <- exp((pp$E_smr / pp$k_B) * (1 / pp$T_ref - 1 / t_opt_k))

  mmr_peak <- fas * smr_at_opt
  gauss_t <- exp(-((thetao_rast - t_opt_mmr)^2) / (2 * pp$sigma_MMR^2))
  do_limit <- o2_rast / (pp$K_crit + o2_rast)
  mmr <- mmr_peak * gauss_t * do_limit

  aas <- mmr - smr
  aas <- ifel(aas < 0, 0, aas)

  aas_max <- (fas - 1) * smr_at_opt
  aas_ratio <- aas / aas_max
  aas_ratio <- clamp(aas_ratio, lower = 0, upper = 1)

  phi <- aas_ratio^gamma
  clamp(phi, lower = 0, upper = 1)
}

get_physio_type <- function(species_name) {
  if (species_name %in% c("Acropora")) return("coral")
  if (species_name %in% c("Lobophora")) return("algae")
  if (species_name %in% c("Scarus", "Cephalopholis")) return("fish")
  return(NULL)
}

compute_physio_constraint <- function(species_name, env_rasters, scope = NULL) {
  physio_type <- get_physio_type(species_name)
  if (is.null(physio_type)) return(NULL)

  if (physio_type == "coral") {
    if (is.null(env_rasters[["thetao"]])) return(NULL)
    return(compute_physio_coral(env_rasters[["thetao"]], env_rasters[["omega_arag"]],
                                env_rasters[["mhw_annual_cumulative"]]))
  }

  if (physio_type == "algae") {
    if (is.null(env_rasters[["thetao"]])) return(NULL)
    # TNC: 用文献 TPC 修正 MaxEnt 热 niche 偏差
    if (!is.null(scope)) {
      tnc <- compute_tnc_algae(env_rasters[["thetao"]], species_name, scope)
      if (!is.null(tnc)) return(tnc)
    }
    return(compute_physio_algae(env_rasters[["thetao"]]))
  }

  if (physio_type == "fish") {
    if (is.null(env_rasters[["thetao"]]) || is.null(env_rasters[["o2"]])) return(NULL)
    return(compute_physio_fish(env_rasters[["thetao"]], env_rasters[["o2"]], species_name))
  }

  NULL
}

# ------------------------------------------------------------------------------
# 2) 数据加载
# ------------------------------------------------------------------------------
load_env_stack <- function(dynamic_dir, static_dir, model_vars, region_extent = NULL,
                           dynamic_override_dir = NULL) {
  static_vars <- if (exists("STATIC_VARS")) {
    STATIC_VARS
  } else {
    c("bathymetry", "slope", "rugosity", "distance_to_coast", "log_shift_slope")
  }

  raster_list <- list()
  missing_vars <- c()

  for (var_name in model_vars) {
    fpath <- if (var_name %in% static_vars) {
      file.path(static_dir, paste0(var_name, ".tif"))
    } else {
      file.path(dynamic_dir, paste0(var_name, ".tif"))
    }
    if (!(var_name %in% static_vars) && !is.null(dynamic_override_dir)) {
      override_path <- file.path(dynamic_override_dir, paste0(var_name, ".tif"))
      if (file.exists(override_path)) fpath <- override_path
    }

    if (!file.exists(fpath)) {
      missing_vars <- c(missing_vars, var_name)
      next
    }
    raster_list[[var_name]] <- rast(fpath)
  }

  if (length(missing_vars) > 0) {
    stop("缺失变量: ", paste(missing_vars, collapse = ", "))
  }

  # 对齐 extent：以第一个栅格为 template，resample 不匹配的栅格
  template_r <- raster_list[[1]]
  for (vn in names(raster_list)[-1]) {
    if (!compareGeom(raster_list[[vn]], template_r, stopOnError = FALSE)) {
      raster_list[[vn]] <- resample(raster_list[[vn]], template_r, method = "bilinear")
    }
  }

  env_stack <- rast(raster_list)
  if (!is.null(region_extent)) {
    env_stack <- crop(env_stack, region_extent)
  }
  env_stack
}

load_scenario_rasters <- function(scenario, resolution_type, model_vars, region_extent = NULL) {
  dynamic_dir <- file.path(SCENARIOS_SHARED_DIR, paste0(scenario, "_", resolution_type))
  static_dir <- file.path(SCENARIOS_SHARED_DIR, paste0("static_", resolution_type))
  # Future dynamic regional rasters have no duplicated static directory in
  # scenarios_shared; reuse the same-resolution immutable temporal static layers.
  if (!dir.exists(static_dir)) {
    temporal_static_dir <- file.path(TEMPORAL_SHARED_DIR, paste0("static_", resolution_type))
    if (dir.exists(temporal_static_dir)) static_dir <- temporal_static_dir
  }

  if (!dir.exists(dynamic_dir)) stop("情景目录不存在: ", dynamic_dir)
  if (!dir.exists(static_dir)) stop("静态目录不存在: ", static_dir)

  load_env_stack(dynamic_dir, static_dir, model_vars, region_extent)
}

load_temporal_rasters <- function(year, resolution_type, model_vars, region_extent = NULL) {
  dynamic_dir <- file.path(TEMPORAL_SHARED_DIR, paste0("temporal_", year, "_", resolution_type))
  static_dir <- file.path(TEMPORAL_SHARED_DIR, paste0("static_", resolution_type))
  override_dir <- NULL
  if (nzchar(TEMPORAL_OVERRIDE_ROOT)) {
    candidate <- file.path(TEMPORAL_OVERRIDE_ROOT, paste0("temporal_", year, "_", resolution_type))
    if (dir.exists(candidate)) override_dir <- candidate
  }
  if (!dir.exists(static_dir)) {
    static_dir <- file.path(SCENARIOS_SHARED_DIR, paste0("static_", resolution_type))
  }

  if (!dir.exists(dynamic_dir)) stop("时序目录不存在: ", dynamic_dir)
  if (!dir.exists(static_dir)) stop("静态目录不存在: ", static_dir)

  load_env_stack(dynamic_dir, static_dir, model_vars, region_extent, override_dir)
}

load_temporal_physio_env <- function(year, resolution_type, species_name, template_rast) {
  temporal_dir <- file.path(TEMPORAL_SHARED_DIR, paste0("temporal_", year, "_", resolution_type))
  baseline_dir <- file.path(SCENARIOS_SHARED_DIR, paste0("baseline_", resolution_type))

  out <- list(thetao = NULL, o2 = NULL, omega_arag = NULL, mhw_annual_cumulative = NULL)

  thetao_path <- file.path(temporal_dir, "thetao.tif")
  if (file.exists(thetao_path)) {
    r <- rast(thetao_path)
    if (!compareGeom(r, template_rast, stopOnError = FALSE)) r <- resample(r, template_rast, method = "bilinear")
    out[["thetao"]] <- r
  }

  if (identical(get_physio_type(species_name), "fish")) {
    o2_path <- file.path(temporal_dir, "o2.tif")
    if (file.exists(o2_path)) {
      r <- rast(o2_path)
      if (!compareGeom(r, template_rast, stopOnError = FALSE)) r <- resample(r, template_rast, method = "bilinear")
      out[["o2"]] <- r
    }
  }

  if (identical(get_physio_type(species_name), "coral")) {
    omega_temporal <- file.path(temporal_dir, "omega_arag.tif")
    omega_baseline <- file.path(baseline_dir, "omega_arag.tif")

    omega_path <- if (file.exists(omega_temporal)) omega_temporal else omega_baseline
    if (file.exists(omega_path)) {
      r <- rast(omega_path)
      if (!compareGeom(r, template_rast, stopOnError = FALSE)) r <- resample(r, template_rast, method = "bilinear")
      out[["omega_arag"]] <- r
    }

  }

  out
}

load_scenario_physio_env <- function(scenario, resolution_type, species_name, template_rast) {
  # 从情景环境目录加载生理约束所需的原始环境变量 (thetao, o2, omega_arag, mhw)
  # 与 load_temporal_physio_env() 完全对称，确保 inline 计算一致性
  scenario_dir <- file.path(SCENARIOS_SHARED_DIR, paste0(scenario, "_", resolution_type))
  baseline_dir <- file.path(SCENARIOS_SHARED_DIR, paste0("baseline_", resolution_type))

  out <- list(thetao = NULL, o2 = NULL, omega_arag = NULL, mhw_annual_cumulative = NULL)

  # thetao: 所有物种的生理约束都依赖温度
  thetao_path <- file.path(scenario_dir, "thetao.tif")
  if (file.exists(thetao_path)) {
    r <- rast(thetao_path)
    if (!compareGeom(r, template_rast, stopOnError = FALSE)) r <- resample(r, template_rast, method = "bilinear")
    out[["thetao"]] <- r
  }

  # o2: 鱼类 OCLTT AAS 模型需要溶解氧
  if (identical(get_physio_type(species_name), "fish")) {
    o2_path <- file.path(scenario_dir, "o2.tif")
    if (file.exists(o2_path)) {
      r <- rast(o2_path)
      if (!compareGeom(r, template_rast, stopOnError = FALSE)) r <- resample(r, template_rast, method = "bilinear")
      out[["o2"]] <- r
    }
  }

  # omega_arag + MHW: 珊瑚 TPC × Ω × 白化惩罚
  if (identical(get_physio_type(species_name), "coral")) {
    omega_scenario <- file.path(scenario_dir, "omega_arag.tif")
    omega_baseline <- file.path(baseline_dir, "omega_arag.tif")

    omega_path <- if (file.exists(omega_scenario)) omega_scenario else omega_baseline
    if (file.exists(omega_path)) {
      r <- rast(omega_path)
      if (!compareGeom(r, template_rast, stopOnError = FALSE)) r <- resample(r, template_rast, method = "bilinear")
      out[["omega_arag"]] <- r
    }

  }

  out
}

load_lambda_contents <- function(run_id) {
  formal_dir <- file.path(OUTPUT_DIR, "maxent/formal", run_id)
  if (!dir.exists(formal_dir)) stop("未找到正式模型目录: ", formal_dir)

  lambda_files <- Sys.glob(file.path(formal_dir, "species_*.lambdas"))
  if (length(lambda_files) == 0) stop("未找到 lambda 文件: ", formal_dir)

  lapply(lambda_files, readLines)
}

# ------------------------------------------------------------------------------
# 3) MaxEnt 预测与 bootstrap 聚合
# ------------------------------------------------------------------------------
predict_with_lambdas <- function(env_stack, lambda_contents) {
  rstack <- raster::stack(env_stack)
  n_layers <- raster::nlayers(rstack)

  pred_list <- list()
  for (i in seq_along(lambda_contents)) {
    pred_i <- tryCatch({
      me <- new("MaxEnt")
      me@lambdas <- lambda_contents[[i]]
      me@hasabsence <- TRUE
      dummy <- as.data.frame(matrix(0, nrow = 1, ncol = n_layers))
      names(dummy) <- names(rstack)
      me@presence <- dummy
      me@absence <- dummy

      pred <- dismo::predict(
        me,
        rstack,
        args = c("doclamp=true", "extrapolate=false", "fadebyclamping=true")
      )
      rast(pred)
    }, error = function(e) {
      cat(sprintf("      bootstrap %d 失败: %s\n", i, e$message))
      NULL
    })

    if (!is.null(pred_i)) pred_list[[length(pred_list) + 1]] <- pred_i
  }

  if (length(pred_list) == 0) {
    stop("全部 bootstrap 投影失败")
  }

  pred_stack <- rast(pred_list)
  list(
    per_bootstrap = pred_list,
    mean = mean(pred_stack, na.rm = TRUE),
    sd = app(pred_stack, sd, na.rm = TRUE)
  )
}

aggregate_pred_list <- function(pred_list) {
  pred_stack <- rast(pred_list)
  list(
    mean = mean(pred_stack, na.rm = TRUE),
    sd = app(pred_stack, sd, na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# 4) MESS
# ------------------------------------------------------------------------------
get_training_background_df <- function(species_name, scope) {
  bg_file <- file.path(SWD_DATA_DIR, paste0("background_raw_", species_name, "_", scope, ".csv"))
  if (!file.exists(bg_file)) {
    stop("未找到 MESS 背景训练数据: ", bg_file)
  }
  read.csv(bg_file, stringsAsFactors = FALSE)
}

calc_mess_raster <- function(env_stack, training_df, vars_use) {
  if (length(vars_use) == 0) return(NULL)

  vars_use <- intersect(vars_use, names(env_stack))
  vars_use <- intersect(vars_use, names(training_df))
  if (length(vars_use) == 0) return(NULL)

  train_x <- training_df[, vars_use, drop = FALSE]
  train_x <- train_x[stats::complete.cases(train_x), , drop = FALSE]
  if (nrow(train_x) < 10) return(NULL)

  x_raster <- raster::stack(env_stack[[vars_use]])
  mess_r <- dismo::mess(x_raster, train_x, full = FALSE)
  rast(mess_r)
}

mess_binary <- function(mess_r) {
  ifel(mess_r >= 0, 1, 0)
}

write_mess_with_alias <- function(mess_dir_primary, mess_dir_alias, filename, raster_obj) {
  if (!dir.exists(mess_dir_primary)) dir.create(mess_dir_primary, recursive = TRUE)
  if (!dir.exists(mess_dir_alias)) dir.create(mess_dir_alias, recursive = TRUE)

  out_primary <- file.path(mess_dir_primary, filename)
  out_alias <- file.path(mess_dir_alias, filename)

  writeRaster(raster_obj, out_primary, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  invisible(file.copy(out_primary, out_alias, overwrite = TRUE))
}

run_mess_for_scenario <- function(species_name, scope, scenario, env_stack, model_vars) {
  run_id <- get_run_id(species_name, scope)
  mess_run_id <- get_mess_run_id(species_name, scope)

  mess_dir_primary <- file.path(MESS_DIR, run_id)
  mess_dir_alias <- file.path(MESS_DIR, mess_run_id)

  train_df <- get_training_background_df(species_name, scope)

  full_vars <- intersect(model_vars, names(env_stack))
  nonthermal_exclude <- c("thetao", "bottomT")
  nonthermal_vars <- setdiff(full_vars, nonthermal_exclude)

  mess_full <- calc_mess_raster(env_stack, train_df, full_vars)
  mess_nonthermal <- calc_mess_raster(env_stack, train_df, nonthermal_vars)

  stats_row <- data.frame(
    run_id = run_id,
    mess_run_id = mess_run_id,
    species = species_name,
    region = scope,
    scenario = scenario,
    mess_mean = NA_real_,
    mess_sd = NA_real_,
    novelty_fraction = NA_real_,
    novelty_nonthermal_fraction = NA_real_,
    stringsAsFactors = FALSE
  )

  if (!is.null(mess_full)) {
    write_mess_with_alias(mess_dir_primary, mess_dir_alias, paste0("mess_", scenario, ".tif"), mess_full)
    mess_full_bin <- mess_binary(mess_full)
    write_mess_with_alias(mess_dir_primary, mess_dir_alias, paste0("mess_", scenario, "_binary.tif"), mess_full_bin)

    stats_row$mess_mean <- safe_stat_mean(mess_full)
    stats_row$mess_sd <- safe_stat_sd(mess_full)
    stats_row$novelty_fraction <- 1 - safe_stat_mean(mess_full_bin)
  }

  if (!is.null(mess_nonthermal)) {
    write_mess_with_alias(mess_dir_primary, mess_dir_alias, paste0("mess_nonthermal_", scenario, ".tif"), mess_nonthermal)
    mess_nonthermal_bin <- mess_binary(mess_nonthermal)
    write_mess_with_alias(mess_dir_primary, mess_dir_alias, paste0("mess_nonthermal_", scenario, "_binary.tif"), mess_nonthermal_bin)

    stats_row$novelty_nonthermal_fraction <- 1 - safe_stat_mean(mess_nonthermal_bin)
  }

  stats_row
}

# ------------------------------------------------------------------------------
# 5) 单情景投影（ssp245/ssp585 未来情景）
# ------------------------------------------------------------------------------
run_scenario_projection <- function(species_name, scope, scenario, lambda_contents, model_vars) {
  run_id <- get_run_id(species_name, scope)
  out_dir <- file.path(PROJECTION_DIR, run_id)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  resolution_type <- get_resolution_type(scope)
  region_extent <- get_region_extent(scope)

  cat(sprintf("    场景投影: %s | %s\n", run_id, scenario))

  env_stack <- load_scenario_rasters(scenario, resolution_type, model_vars, region_extent)
  pred_obj <- predict_with_lambdas(env_stack, lambda_contents)

  # 从情景环境栅格 inline 计算生理约束（与时序路径统一）
  physio_env <- load_scenario_physio_env(scenario, resolution_type, species_name, pred_obj$mean)
  phi <- compute_physio_constraint(species_name, physio_env, scope = scope)

  constrained_bs <- lapply(pred_obj$per_bootstrap, function(p) {
    if (is.null(phi)) return(p)
    if (!compareGeom(phi, p, stopOnError = FALSE)) {
      p <- resample(p, phi, method = "bilinear")
    }
    # 幂次衰减约束: HSI × Φ^γ (γ < 1 软化约束)
    # TNC 修正时 phi 可能 >1 (提升 MaxEnt 采样偏差导致的低估)
    gamma <- if (exists("PHYSIO_GAMMA")) PHYSIO_GAMMA else 1.0
    upper <- if (get_physio_type(species_name) == "algae" &&
                 isTRUE(PHYSIO_PARAMS$algae$tnc_enabled)) {
      PHYSIO_PARAMS$algae$tnc_cap
    } else {
      1.0
    }
    phi_eff <- clamp(phi, lower = 0, upper = upper) ^ gamma
    clamp(p * phi_eff, lower = 0, upper = 1)
  })

  constrained_agg <- aggregate_pred_list(constrained_bs)

  out_unc <- file.path(out_dir, paste0(scenario, "_hsi_mean_unconstrained.tif"))
  out_mean <- file.path(out_dir, paste0(scenario, "_hsi_mean.tif"))
  out_sd <- file.path(out_dir, paste0(scenario, "_hsi_sd.tif"))
  out_clamp <- file.path(out_dir, paste0("clamp_", scenario, "_binary.tif"))

  writeRaster(pred_obj$mean, out_unc, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  writeRaster(constrained_agg$mean, out_mean, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  writeRaster(constrained_agg$sd, out_sd, overwrite = TRUE, gdal = c("COMPRESS=LZW"))

  # 当前脚本未显式提取 clamp 元信息，先输出占位（全1，表示默认可用域）
  clamp_bin <- ifel(!is.na(constrained_agg$mean), 1, NA)
  writeRaster(clamp_bin, out_clamp, overwrite = TRUE, gdal = c("COMPRESS=LZW"))

  mess_stats <- run_mess_for_scenario(species_name, scope, scenario, env_stack, model_vars)

  list(
    scenario = scenario,
    hsi_mean = safe_stat_mean(constrained_agg$mean),
    hsi_sd_global = safe_stat_mean(constrained_agg$sd),
    unconstrained_mean = safe_stat_mean(pred_obj$mean),
    mess_stats = mess_stats,
    mean_raster = constrained_agg$mean
  )
}

# ------------------------------------------------------------------------------
# 6) 时序投影（逐年）
# ------------------------------------------------------------------------------
run_temporal_projection <- function(species_name, scope, lambda_contents, model_vars) {
  run_id <- get_run_id(species_name, scope)
  out_dir <- file.path(TEMPORAL_OUTPUT_DIR, run_id)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  resolution_type <- get_resolution_type(scope)
  region_extent <- get_region_extent(scope)

  available_years <- c()
  for (yr in TEMPORAL_YEARS) {
    td <- file.path(TEMPORAL_SHARED_DIR, paste0("temporal_", yr, "_", resolution_type))
    if (dir.exists(td)) available_years <- c(available_years, yr)
  }

  if (length(available_years) == 0) {
    cat(sprintf("    时序数据不可用: %s\n", run_id))
    return(list(
      has_temporal = FALSE,
      available_years = integer(0),
      temporal_mean = NULL,
      temporal_sd = NULL,
      temporal_mean_unconstrained = NULL
    ))
  }

  cat(sprintf("    时序投影: %s | 可用年份 %d\n", run_id, length(available_years)))

  annual_constrained <- list()
  annual_unconstrained <- list()
  names_constrained <- c()

  for (yr in available_years) {
    cat(sprintf("      年份 %d ...\n", yr))

    out_hsi <- file.path(out_dir, paste0("hsi_", yr, ".tif"))
    out_hsi_unc <- file.path(out_dir, paste0("hsi_", yr, "_unconstrained.tif"))
    out_hsi_sd <- file.path(out_dir, paste0("hsi_", yr, "_sd.tif"))

    if (SKIP_EXISTING_TEMPORAL &&
        file.exists(out_hsi) && file.exists(out_hsi_unc) && file.exists(out_hsi_sd)) {
      cached_hsi <- tryCatch(rast(out_hsi), error = function(e) NULL)
      cached_unc <- tryCatch(rast(out_hsi_unc), error = function(e) NULL)
      if (!is.null(cached_hsi) && !is.null(cached_unc)) {
        annual_constrained[[length(annual_constrained) + 1]] <- cached_hsi
        annual_unconstrained[[length(annual_unconstrained) + 1]] <- cached_unc
        names_constrained <- c(names_constrained, as.character(yr))
        next
      }
    }

    env_stack <- tryCatch(
      load_temporal_rasters(yr, resolution_type, model_vars, region_extent),
      error = function(e) {
        cat(sprintf("        跳过 %d: %s\n", yr, e$message))
        NULL
      }
    )
    if (is.null(env_stack)) next

    pred_obj <- tryCatch(
      predict_with_lambdas(env_stack, lambda_contents),
      error = function(e) {
        cat(sprintf("        跳过 %d: %s\n", yr, e$message))
        NULL
      }
    )
    if (is.null(pred_obj)) next

    physio_env <- load_temporal_physio_env(yr, resolution_type, species_name, pred_obj$mean)
    phi <- compute_physio_constraint(species_name, physio_env, scope = scope)

    constrained_bs <- lapply(pred_obj$per_bootstrap, function(p) {
      if (is.null(phi)) return(p)
      if (!compareGeom(phi, p, stopOnError = FALSE)) {
        p <- resample(p, phi, method = "bilinear")
      }
      # 幂次衰减约束: HSI × Φ^γ (γ < 1 软化约束)
      # TNC 修正时 phi 可能 >1 (提升 MaxEnt 采样偏差导致的低估)
      gamma <- if (exists("PHYSIO_GAMMA")) PHYSIO_GAMMA else 1.0
      upper <- if (get_physio_type(species_name) == "algae" &&
                   isTRUE(PHYSIO_PARAMS$algae$tnc_enabled)) {
        PHYSIO_PARAMS$algae$tnc_cap
      } else {
        1.0
      }
      phi_eff <- clamp(phi, lower = 0, upper = upper) ^ gamma
      clamp(p * phi_eff, lower = 0, upper = 1)
    })

    constrained_agg <- aggregate_pred_list(constrained_bs)

    writeRaster(constrained_agg$mean, out_hsi, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    writeRaster(pred_obj$mean, out_hsi_unc, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    writeRaster(constrained_agg$sd, out_hsi_sd, overwrite = TRUE, gdal = c("COMPRESS=LZW"))

    annual_constrained[[length(annual_constrained) + 1]] <- constrained_agg$mean
    annual_unconstrained[[length(annual_unconstrained) + 1]] <- pred_obj$mean
    names_constrained <- c(names_constrained, as.character(yr))

    gc(verbose = FALSE)
  }

  if (length(annual_constrained) == 0) {
    return(list(
      has_temporal = FALSE,
      available_years = integer(0),
      temporal_mean = NULL,
      temporal_sd = NULL,
      temporal_mean_unconstrained = NULL
    ))
  }

  if (isTRUE(ANNUAL_ONLY)) {
    cat(sprintf("    Annual-only mode: wrote %d annual raster set(s); aggregation deferred.\n", length(annual_constrained)))
    return(list(
      has_temporal = TRUE,
      available_years = as.integer(names_constrained),
      temporal_mean = NULL,
      temporal_sd = NULL,
      temporal_mean_unconstrained = NULL
    ))
  }

  names(annual_constrained) <- names_constrained
  names(annual_unconstrained) <- names_constrained

  # 导出时序均值（替代旧 baseline）
  temporal_stack <- rast(annual_constrained)
  temporal_mean <- mean(temporal_stack, na.rm = TRUE)
  temporal_sd <- app(temporal_stack, sd, na.rm = TRUE)

  temporal_unc_stack <- rast(annual_unconstrained)
  temporal_mean_unc <- mean(temporal_unc_stack, na.rm = TRUE)

  writeRaster(temporal_mean,
    file.path(out_dir, "temporal_mean.tif"),
    overwrite = TRUE, gdal = c("COMPRESS=LZW")
  )
  writeRaster(temporal_sd,
    file.path(out_dir, "temporal_sd.tif"),
    overwrite = TRUE, gdal = c("COMPRESS=LZW")
  )
  writeRaster(temporal_mean_unc,
    file.path(out_dir, "temporal_mean_unconstrained.tif"),
    overwrite = TRUE, gdal = c("COMPRESS=LZW")
  )

  cat(sprintf("    已导出 temporal_mean.tif (%d 年聚合)\n", length(annual_constrained)))

  list(
    has_temporal = TRUE,
    available_years = as.integer(names_constrained),
    temporal_mean = temporal_mean,
    temporal_sd = temporal_sd,
    temporal_mean_unconstrained = temporal_mean_unc
  )
}

# ------------------------------------------------------------------------------
# 7) 单 run 主流程
# ------------------------------------------------------------------------------
run_unified_projection <- function(species_name, scope) {
  run_id <- get_run_id(species_name, scope)
  out_dir <- file.path(PROJECTION_DIR, run_id)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("统一投影:", run_id, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")

  model_vars <- get_final_vars(species_name, scope)
  if (length(model_vars) == 0) {
    stop("模型变量为空: ", run_id)
  }
  cat("  模型变量数:", length(model_vars), "\n")

  lambda_contents <- load_lambda_contents(run_id)
  cat("  bootstrap 模型数:", length(lambda_contents), "\n")

  # 1) 时序核心：30 年逐年 HSI + temporal_mean 导出
  temporal_out_path <- file.path(TEMPORAL_OUTPUT_DIR, run_id, "temporal_mean.tif")
  if (!ANNUAL_ONLY && SKIP_EXISTING_TEMPORAL && file.exists(temporal_out_path)) {
    cat(sprintf("  [跳过] 时序投影已有缓存: %s\n", temporal_out_path))
    # 构造兼容返回对象，供后续逻辑读取已有结果
    temporal_mean_cached <- tryCatch(rast(temporal_out_path), error = function(e) NULL)
    temporal_sd_path <- file.path(TEMPORAL_OUTPUT_DIR, run_id, "temporal_sd.tif")
    temporal_sd_cached <- tryCatch(
      if (file.exists(temporal_sd_path)) rast(temporal_sd_path) else NULL,
      error = function(e) NULL
    )
    temporal_unc_path <- file.path(TEMPORAL_OUTPUT_DIR, run_id, "temporal_mean_unconstrained.tif")
    temporal_unc_cached <- tryCatch(
      if (file.exists(temporal_unc_path)) rast(temporal_unc_path) else NULL,
      error = function(e) NULL
    )
    temporal_res <- list(
      has_temporal = !is.null(temporal_mean_cached),
      available_years = integer(0),
      temporal_mean = temporal_mean_cached,
      temporal_sd = temporal_sd_cached,
      temporal_mean_unconstrained = temporal_unc_cached
    )
  } else {
    temporal_res <- run_temporal_projection(species_name, scope, lambda_contents, model_vars)
  }

  if (isTRUE(ANNUAL_ONLY)) {
    return(list(
      run_id = run_id,
      species = species_name,
      region = scope,
      temporal_mean_hsi = NA_real_,
      projection_summary = data.frame(),
      mess_summary = data.frame()
    ))
  }

  scenario_rows <- list()
  mess_rows <- list()

  # 2) 将 temporal_mean 复制到 projection 目录（供下游 SEM/coupling 使用）
  if (!is.null(temporal_res$temporal_mean)) {
    writeRaster(temporal_res$temporal_mean,
      file.path(out_dir, "temporal_mean_hsi.tif"),
      overwrite = TRUE, gdal = c("COMPRESS=LZW")
    )
    # 向后兼容: 下游脚本 (fig_03, SEM, validation 等) 仍引用 baseline_hsi_mean.tif
    file.copy(file.path(out_dir, "temporal_mean_hsi.tif"),
              file.path(out_dir, "baseline_hsi_mean.tif"),
              overwrite = TRUE)

    if (!is.null(temporal_res$temporal_mean_unconstrained)) {
      writeRaster(temporal_res$temporal_mean_unconstrained,
        file.path(out_dir, "temporal_mean_hsi_unconstrained.tif"),
        overwrite = TRUE, gdal = c("COMPRESS=LZW")
      )
      file.copy(file.path(out_dir, "temporal_mean_hsi_unconstrained.tif"),
                file.path(out_dir, "baseline_hsi_mean_unconstrained.tif"),
                overwrite = TRUE)
    }
    if (!is.null(temporal_res$temporal_sd)) {
      writeRaster(temporal_res$temporal_sd,
        file.path(out_dir, "temporal_sd_hsi.tif"),
        overwrite = TRUE, gdal = c("COMPRESS=LZW")
      )
      file.copy(file.path(out_dir, "temporal_sd_hsi.tif"),
                file.path(out_dir, "baseline_hsi_sd.tif"),
                overwrite = TRUE)
    }
  }

  # 3) SSP 未来情景投影
  for (sc in SCENARIOS) {
    # 跳过本地缓存：若两个关键输出均已存在则直接读取统计量，跳过重新训练
    cached_mean_path <- file.path(out_dir, paste0(sc, "_hsi_mean.tif"))
    cached_unc_path  <- file.path(out_dir, paste0(sc, "_hsi_mean_unconstrained.tif"))
    if (SKIP_EXISTING_SCENARIOS &&
        file.exists(cached_mean_path) && file.exists(cached_unc_path)) {
      cat(sprintf("  [跳过] %s 情景已有缓存\n", sc))
      cached_mean_r <- tryCatch(rast(cached_mean_path), error = function(e) NULL)
      cached_unc_r  <- tryCatch(rast(cached_unc_path),  error = function(e) NULL)
      sc_res <- list(
        scenario     = sc,
        hsi_mean     = if (!is.null(cached_mean_r)) safe_stat_mean(cached_mean_r) else NA_real_,
        hsi_sd_global = NA_real_,
        unconstrained_mean = if (!is.null(cached_unc_r)) safe_stat_mean(cached_unc_r) else NA_real_,
        mess_stats   = data.frame(
                          run_id = run_id,
                          mess_run_id = run_id,
                          species = species_name,
                          region = scope,
                          scenario = sc,
                          mess_mean = NA_real_,
                          mess_sd = NA_real_,
                          novelty_fraction = NA_real_,
                          novelty_nonthermal_fraction = NA_real_,
                          stringsAsFactors = FALSE),
        mean_raster  = cached_mean_r
      )
    } else {
      sc_res <- run_scenario_projection(species_name, scope, sc, lambda_contents, model_vars)
    }
    scenario_rows[[length(scenario_rows) + 1]] <- data.frame(
      run_id = run_id,
      species = species_name,
      region = scope,
      scenario = sc,
      hsi_mean = sc_res$hsi_mean,
      hsi_sd_global = sc_res$hsi_sd_global,
      unconstrained_mean = sc_res$unconstrained_mean,
      stringsAsFactors = FALSE
    )
    mess_rows[[length(mess_rows) + 1]] <- sc_res$mess_stats
  }

  # 4) 汇总
  projection_summary <- bind_rows(scenario_rows)

  # 添加时序均值统计到摘要
  if (!is.null(temporal_res$temporal_mean)) {
    scenario_rows[[length(scenario_rows) + 1]] <- data.frame(
      run_id = run_id,
      species = species_name,
      region = scope,
      scenario = "temporal_mean",
      hsi_mean = safe_stat_mean(temporal_res$temporal_mean),
      hsi_sd_global = safe_stat_mean(temporal_res$temporal_sd),
      unconstrained_mean = safe_stat_mean(temporal_res$temporal_mean_unconstrained),
      stringsAsFactors = FALSE
    )
    projection_summary <- bind_rows(scenario_rows)
  }

  write.csv(projection_summary,
    file.path(out_dir, "projection_summary.csv"),
    row.names = FALSE
  )

  list(
    run_id = run_id,
    species = species_name,
    region = scope,
    temporal_mean_hsi = if (!is.null(temporal_res$temporal_mean)) safe_stat_mean(temporal_res$temporal_mean) else NA_real_,
    projection_summary = projection_summary,
    mess_summary = bind_rows(mess_rows)
  )
}

# ------------------------------------------------------------------------------
# 8) 全部运行
# ------------------------------------------------------------------------------
cat("\n开始统一投影流程（时序 + 情景 + MESS）...\n")
cat("项目根目录:", PROJECT_ROOT, "\n")
cat("时序输入:", TEMPORAL_SHARED_DIR, "\n")
cat("情景输入:", SCENARIOS_SHARED_DIR, "\n")
cat("输出目录:", PROJECTION_DIR, "\n\n")

all_rows <- list()
all_proj_summaries <- list()
all_mess_summaries <- list()

scopes_all <- c("global", names(REGIONS))
task_grid <- expand.grid(
  species = SPECIES_LIST,
  scope = scopes_all,
  stringsAsFactors = FALSE
)

FORCE_RERUN <- tolower(Sys.getenv("FORCE_RERUN", unset = "false")) %in% c("1", "true", "yes")
ANNUAL_ONLY <- tolower(Sys.getenv("MAXENT_ANNUAL_ONLY", unset = "false")) %in% c("1", "true", "yes")

cat(sprintf("并行配置: workers=%d, per-worker JVM=-Xmx%dg, 强制重跑=%s\n",
  MAXENT_WORKERS, MAXENT_JAVA_XMX_GB, ifelse(FORCE_RERUN, "是", "否")))

# 检查单个 run 是否已完成：projection_summary.csv 存在 + temporal 年份数 >= 可用年份数
check_task_complete <- function(sp, scope) {
  run_id <- paste0(sp, "_", scope)

  # 1) projection_summary.csv 必须存在
  proj_csv <- file.path(PROJECTION_DIR, run_id, "projection_summary.csv")
  if (!file.exists(proj_csv)) return(FALSE)

  # 2) 时序 HSI 文件数 >= 可用年份数
  resolution_type <- get_resolution_type(scope)
  available_years <- 0L
  for (yr in TEMPORAL_YEARS) {
    td <- file.path(TEMPORAL_SHARED_DIR, paste0("temporal_", yr, "_", resolution_type))
    if (dir.exists(td)) available_years <- available_years + 1L
  }

  if (available_years > 0) {
    temporal_dir <- file.path(TEMPORAL_OUTPUT_DIR, run_id)
    if (!dir.exists(temporal_dir)) return(FALSE)
    n_hsi <- length(Sys.glob(file.path(temporal_dir, "hsi_[0-9]*.tif")))
    if (n_hsi < available_years) return(FALSE)
  }

  # 3) 所有情景的 hsi_mean.tif 文件必须存在（新增情景时可强制重新运行）
  for (sc in SCENARIOS) {
    mean_path <- file.path(PROJECTION_DIR, run_id, paste0(sc, "_hsi_mean.tif"))
    if (!file.exists(mean_path)) return(FALSE)
  }

  TRUE
}

run_one_task <- function(sp, scope) {
  run_id <- paste0(sp, "_", scope)

  # 跳过已完成的任务
  if (!ANNUAL_ONLY && !FORCE_RERUN && check_task_complete(sp, scope)) {
    cat(sprintf("  ✓ 跳过（已完成）: %s\n", run_id))

    # 读取已有结果摘要返回，使汇总 CSV 完整
    proj_csv <- file.path(PROJECTION_DIR, run_id, "projection_summary.csv")
    existing_summary <- tryCatch(read.csv(proj_csv, stringsAsFactors = FALSE), error = function(e) NULL)

    return(list(
      ok = TRUE,
      run_id = run_id,
      species = sp,
      scope = scope,
      error = NA_character_,
      elapsed_sec = 0,
      skipped = TRUE,
      result = list(
        run_id = run_id,
        species = sp,
        region = scope,
        temporal_mean_hsi = if (!is.null(existing_summary) && "hsi_mean" %in% names(existing_summary))
          existing_summary$hsi_mean[existing_summary$scenario == "temporal_mean"][1] else NA_real_,
        projection_summary = existing_summary,
        mess_summary = data.frame()
      )
    ))
  }

  t0 <- Sys.time()

  out <- tryCatch({
    res <- run_unified_projection(sp, scope)
    list(
      ok = TRUE,
      run_id = run_id,
      species = sp,
      scope = scope,
      error = NA_character_,
      result = res
    )
  }, error = function(e) {
    list(
      ok = FALSE,
      run_id = run_id,
      species = sp,
      scope = scope,
      error = conditionMessage(e),
      result = NULL
    )
  })

  out$elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out
}

results <- list()

# 显式列出需要导出到 worker 的函数和标量变量
# 不导出 rJava/terra 外部指针对象（不可跨进程序列化）
EXPORT_FUNCTIONS <- c(
  # 本脚本定义的函数
  "run_one_task", "run_unified_projection", "run_temporal_projection",
  "run_scenario_projection", "run_mess_for_scenario",
  "predict_with_lambdas", "aggregate_pred_list",
  "load_env_stack", "load_scenario_rasters", "load_temporal_rasters",
  "load_temporal_physio_env", "load_scenario_physio_env", "load_lambda_contents",
  "calc_mess_raster", "mess_binary", "write_mess_with_alias",
  "get_training_background_df",
  "compute_physio_coral", "compute_physio_algae", "compute_physio_fish",
  "load_maxent_thetao_response", "compute_tnc_algae",
  "get_physio_type", "compute_physio_constraint",
  "get_run_id", "get_mess_run_id", "get_resolution_type", "get_region_extent",
  "safe_stat_mean", "safe_stat_sd",
  "check_task_complete"
)

EXPORT_SCALARS <- c(
  # 脚本级路径常量
  "PROJECT_ROOT", "MAXENT_JAR", "MAXENT_JAVA_XMX_GB",
  "PROJECTION_DIR", "MESS_DIR", "TEMPORAL_OUTPUT_DIR",
  "SCENARIOS",
  "DEFAULT_MHW_WINDOW",
  "FORCE_RERUN", "ANNUAL_ONLY",
  "SKIP_EXISTING_SCENARIOS", "SKIP_EXISTING_TEMPORAL",
  # 并行任务表
  "task_grid"
)

if (MAXENT_WORKERS > 1) {
  cat("启用并行投影（Windows PSOCK）...\n")
  cl <- parallel::makeCluster(MAXENT_WORKERS, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # 每个 worker 独立初始化：加载包 + source config.R + 初始化 JVM
  parallel::clusterExport(cl,
    varlist = c(EXPORT_FUNCTIONS, EXPORT_SCALARS),
    envir = .GlobalEnv
  )
    parallel::clusterEvalQ(cl, {
    options(java.parameters = c(
      sprintf("-Xmx%dg", MAXENT_JAVA_XMX_GB),
      "-XX:+UseG1GC",
      "-XX:ParallelGCThreads=1"
    ))
    suppressPackageStartupMessages({
      library(rJava)
      library(dismo)
      library(terra)
      library(raster)
      library(dplyr)
    })
    try(.jinit(), silent = TRUE)
    .jaddClassPath(MAXENT_JAR)
    terra::terraOptions(memfrac = 0.6, threads = 1)

    # Worker 自行 source config.R 获取 REGIONS / PHYSIO_PARAMS / get_final_vars 等
      setwd(PROJECT_ROOT)
      source(file.path(PROJECT_ROOT, "config/config.R"))

      # config.R restores mainline output roots. In this isolated branch the
      # workers must instead load the branch formal lambdas and write only
      # branch outputs; PROJECTION_DIR was explicitly exported from master.
      MAXENT_OUTPUT_DIR <- dirname(PROJECTION_DIR)
      OUTPUT_DIR <- dirname(MAXENT_OUTPUT_DIR)
      TEMPORAL_OUTPUT_DIR <- file.path(MAXENT_OUTPUT_DIR, "temporal")
      branch_root_worker <- OUTPUT_DIR
      # Do not retain config.R's broad default predictor list in workers.
      # Formal branch models were fitted with species-specific VIF-selected
      # variables; using any additional layer would invalidate projection.
      get_final_vars <- function(species, region, mhw_window = NULL) {
        selected <- file.path(branch_root_worker, "vif", paste0(species, "_", region), "final_variables.txt")
        if (!file.exists(selected)) stop("Missing branch final-variable file: ", selected)
        values <- readLines(selected, warn = FALSE)
        values[nzchar(values)]
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

      if (!exists("SCENARIOS_SHARED_DIR")) {
        SCENARIOS_SHARED_DIR <<- file.path(PROCESSED_DATA_DIR, "scenarios_shared")
      }
    NULL
  })

  cat(sprintf("  任务数: %d, 分配到 %d workers\n", nrow(task_grid), MAXENT_WORKERS))

  results <- parallel::parLapply(cl, seq_len(nrow(task_grid)), function(i) {
    sp <- task_grid$species[i]
    scope <- task_grid$scope[i]
    run_one_task(sp, scope)
  })
} else {
  cat("并行已关闭，使用串行投影。\n")
  results <- lapply(seq_len(nrow(task_grid)), function(i) {
    sp <- task_grid$species[i]
    scope <- task_grid$scope[i]
    run_one_task(sp, scope)
  })
}

if (isTRUE(ANNUAL_ONLY)) {
  failures <- Filter(function(x) is.null(x) || !isTRUE(x$ok), results)
  if (length(failures) > 0) stop("Annual-only projection task failed")
  cat("\nAnnual-only projection completed; shared summaries intentionally deferred.\n")
  quit(save = "no", status = 0)
}

for (one in results) {
  if (is.null(one)) next

  if (isTRUE(one$ok) && !is.null(one$result)) {
    res <- one$result
    all_rows[[length(all_rows) + 1]] <- data.frame(
      run_id = res$run_id,
      species = res$species,
      region = res$region,
      temporal_mean_hsi = res$temporal_mean_hsi,
      stringsAsFactors = FALSE
    )
    all_proj_summaries[[length(all_proj_summaries) + 1]] <- res$projection_summary
    all_mess_summaries[[length(all_mess_summaries) + 1]] <- res$mess_summary
  } else {
    cat(sprintf("\n✗ 失败: %s | %s (%.1fs)\n", one$run_id, one$error, one$elapsed_sec))
  }
  gc(verbose = FALSE)
}

fail_rows <- lapply(results, function(one) {
  if (is.null(one) || isTRUE(one$ok)) return(NULL)
  data.frame(
    run_id = one$run_id,
    species = one$species,
    scope = one$scope,
    error = one$error,
    elapsed_sec = one$elapsed_sec,
    stringsAsFactors = FALSE
  )
})
fail_rows <- Filter(Negate(is.null), fail_rows)

if (length(fail_rows) > 0) {
  fail_df <- bind_rows(fail_rows)
  fail_path <- file.path(PROJECTION_DIR, "parallel_failures.csv")
  write.csv(fail_df, fail_path, row.names = FALSE)
  cat(sprintf("\n⚠ 已写出失败任务清单: %s\n", fail_path))
}

if (length(all_rows) == 0) {
  stop("全部任务失败，未生成任何投影结果。请检查 output/01_mainline/sdm_hsi/maxent/projection/parallel_failures.csv")
}

if (length(all_rows) > 0) {
  all_projection_summary <- bind_rows(all_rows)
  write.csv(all_projection_summary,
    file.path(PROJECTION_DIR, "all_projection_summary.csv"),
    row.names = FALSE
  )
  cat("\n✓ 已写出 all_projection_summary.csv\n")
}

if (length(all_proj_summaries) > 0) {
  all_proj_detail <- bind_rows(all_proj_summaries)
  write.csv(all_proj_detail,
    file.path(PROJECTION_DIR, "all_projection_summary_detailed.csv"),
    row.names = FALSE
  )
}

if (length(all_mess_summaries) > 0) {
  combined_mess <- bind_rows(Filter(function(x) nrow(x) > 0, all_mess_summaries))
  if (nrow(combined_mess) > 0 && all(c("run_id", "scenario") %in% names(combined_mess))) {
    all_mess_summary <- distinct(combined_mess, run_id, scenario, .keep_all = TRUE)
  } else {
    all_mess_summary <- combined_mess
  }
  write.csv(all_mess_summary,
    file.path(MESS_DIR, "all_mess_summary.csv"),
    row.names = FALSE
  )
  cat("✓ 已写出 all_mess_summary.csv\n")
}

cat("\n统一投影流程完成。\n")
cat("后续可运行：08_future_mhw_risk.R\n")
