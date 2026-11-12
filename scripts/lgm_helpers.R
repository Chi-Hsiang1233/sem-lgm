# ============================================================
# lgm_helpers.R —— 所有分析腳本共用的函式
# 其他腳本以 source(file.path(script_dir(), "lgm_helpers.R")) 載入
# ============================================================

script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}

# 解析 --key value 與 --flag 形式的命令列參數
parse_args <- function(args = commandArgs(TRUE)) {
  out <- list(); i <- 1
  while (i <= length(args)) {
    k <- args[i]
    if (!startsWith(k, "--")) stop("無法解析參數：", k)
    key <- gsub("-", "_", sub("^--", "", k))
    if (i < length(args) && !startsWith(args[i + 1], "--")) {
      out[[key]] <- args[i + 1]; i <- i + 2
    } else {
      out[[key]] <- TRUE; i <- i + 1
    }
  }
  out
}
arg_or <- function(a, key, default = NULL) if (is.null(a[[key]])) default else a[[key]]
need <- function(a, keys) {
  miss <- keys[vapply(keys, function(k) is.null(a[[k]]), TRUE)]
  if (length(miss)) stop("缺少必要參數：", paste0("--", gsub("_", "-", miss), collapse = ", "), call. = FALSE)
}
split_csv <- function(x) if (is.null(x) || isTRUE(x)) NULL else trimws(strsplit(x, ",")[[1]])

require_pkgs <- function(...) {
  p <- c(...)
  miss <- p[!vapply(p, requireNamespace, TRUE, quietly = TRUE)]
  if (length(miss)) {
    stop("缺少 R 套件：", paste(miss, collapse = ", "),
         "\n請執行：Rscript scripts/check_env.R --install", call. = FALSE)
  }
}

# 讀取 csv / tsv / sav / dta / xlsx / rds
read_any <- function(path, na = NULL) {
  ext <- tolower(tools::file_ext(path))
  na_str <- c("NA", "", na)
  d <- switch(ext,
    csv = utils::read.csv(path, na.strings = na_str, fileEncoding = "UTF-8-BOM", check.names = FALSE),
    tsv = utils::read.delim(path, na.strings = na_str, fileEncoding = "UTF-8-BOM", check.names = FALSE),
    txt = utils::read.delim(path, na.strings = na_str, fileEncoding = "UTF-8-BOM", check.names = FALSE),
    sav = { require_pkgs("haven"); as.data.frame(haven::zap_labels(haven::read_sav(path))) },
    dta = { require_pkgs("haven"); as.data.frame(haven::zap_labels(haven::read_dta(path))) },
    xlsx = { require_pkgs("readxl"); as.data.frame(readxl::read_excel(path)) },
    xls = { require_pkgs("readxl"); as.data.frame(readxl::read_excel(path)) },
    rds = readRDS(path),
    stop("不支援的檔案格式：", ext, call. = FALSE))
  if (!is.null(na)) {
    nn <- suppressWarnings(as.numeric(na))
    if (!is.na(nn)) d[] <- lapply(d, function(x) { if (is.numeric(x)) x[x == nn] <- NA; x })
  }
  d
}

read_model <- function(path) paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

ensure_dir <- function(p) { dir.create(p, showWarnings = FALSE, recursive = TRUE); p }
wtab <- function(x, dir, name, digits = 6) {
  num <- vapply(x, is.numeric, TRUE)
  x[num] <- lapply(x[num], round, digits)
  utils::write.csv(x, file.path(dir, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(x)
}

# 判斷一組模型可共同使用的指標類型（比較模型時不可混用 robust 與 scaled）
common_kind <- function(fits) {
  fm <- lapply(fits, lavaan::fitMeasures)
  ok <- function(n) all(vapply(fm, function(f) n %in% names(f) && !is.na(f[[n]]), TRUE))
  if (!ok("chisq.scaled")) return("standard")
  if (ok("cfi.robust")) "robust" else "scaled"
}

# 從 lavaan 物件取出配適指標。
# 有校正檢定統計量時 robust > scaled；一般 ML（含自助法 SE）只用未校正指標。
fit_table <- function(fit, kind = NULL) {
  fm <- lavaan::fitMeasures(fit)
  get1 <- function(stems) {
    for (s in stems) if (s %in% names(fm) && !is.na(fm[[s]])) return(unname(fm[[s]]))
    NA_real_
  }
  if (is.null(kind)) kind <- common_kind(list(fit))
  suffix <- switch(kind, robust = c(".robust", ".scaled", ""), scaled = c(".scaled", ""), standard = "")
  chain <- function(st) paste0(st, suffix)
  test <- if (kind == "standard") "" else ".scaled"
  data.frame(
    n_obs = lavaan::lavInspect(fit, "ntotal"),
    npar = get1("npar"),
    chisq = get1(paste0("chisq", test)),
    df = get1(paste0("df", test)),
    pvalue = get1(paste0("pvalue", test)),
    cfi = get1(chain("cfi")),
    tli = get1(chain("tli")),
    rmsea = get1(chain("rmsea")),
    rmsea_lo = get1(chain("rmsea.ci.lower")),
    rmsea_hi = get1(chain("rmsea.ci.upper")),
    srmr = get1(c("srmr", "srmr_within")),
    srmr_between = get1("srmr_between"),
    aic = get1("aic"),
    bic = get1("bic"),
    index_type = kind
  )
}

# 估計問題檢查：Heywood、非正定、標準化負荷 > 1、收斂
# 只評估自由估計的負荷（成長模型與 LCS 的固定負荷不是測量品質問題）；
# 低負荷以「每題最大負荷」判斷，避免把 ESEM 的小跨負荷或雙因子的一般因子負荷誤判。
diagnose <- function(fit) {
  msg <- character()
  if (!lavaan::lavInspect(fit, "converged")) return("[嚴重] 模型未收斂，結果不可解釋。")
  pe <- lavaan::parameterEstimates(fit)
  neg <- pe[pe$op == "~~" & pe$lhs == pe$rhs & !is.na(pe$est) & pe$est < 0, ]
  if (nrow(neg)) msg <- c(msg, paste0("[嚴重] 負變異數（Heywood case）：", paste(unique(neg$lhs), collapse = ", ")))
  ss <- lavaan::standardizedSolution(fit, se = FALSE)
  free_load <- pe[pe$op == "=~" & !is.na(pe$se) & pe$se > 0, ]
  keep <- paste(free_load$lhs, free_load$rhs)
  L <- ss[ss$op == "=~" & !is.na(ss$est.std) & paste(ss$lhs, ss$rhs) %in% keep, ]
  big <- L[abs(L$est.std) > 1, ]
  if (nrow(big)) msg <- c(msg, paste0("[嚴重] 標準化負荷絕對值 > 1：", paste(unique(paste(big$lhs, big$rhs, sep = "=~")), collapse = ", ")))
  if (nrow(L)) {
    allL <- ss[ss$op == "=~" & !is.na(ss$est.std) & ss$rhs %in% L$rhs, ]
    mx <- tapply(abs(allL$est.std), allL$rhs, max)
    low <- names(mx)[mx < .40]
    if (length(low)) msg <- c(msg, paste0("[注意] 主要標準化負荷 < .40 的指標：", paste(low, collapse = ", ")))
  }
  lv <- lavaan::lavNames(fit, "lv")
  hi <- ss[ss$op == "~~" & ss$lhs != ss$rhs & !is.na(ss$est.std) & abs(ss$est.std) > .90 & ss$lhs %in% lv & ss$rhs %in% lv, ]
  if (nrow(hi)) msg <- c(msg, paste0("[注意] 潛在相關 > .90：", paste(unique(paste(hi$lhs, hi$rhs, sep = "~~")), collapse = ", ")))
  tryCatch(withCallingHandlers(lavaan::lavInspect(fit, "post.check"),
    warning = function(w) { msg <<- c(msg, paste0("[注意] post.check：", gsub("\\s+", " ", conditionMessage(w)))); invokeRestart("muffleWarning") }),
    error = function(e) NULL)
  if (!length(msg)) msg <- "未偵測到估計問題（收斂、無負變異數、標準化負荷皆在合理範圍）。"
  msg
}

fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < .001, "< .001", paste0("= ", sub("^0", "", formatC(p, 3, format = "f")))))
fmt2 <- function(x, d = 3) sub("^(-?)0\\.", "\\1.", formatC(x, d, format = "f"))

# 一行的配適摘要
fit_sentence <- function(ft, estimator = "") {
  if (isTRUE(ft$df == 0)) return(sprintf("配適：恰好識別（df = 0），配適指標沒有檢驗意義（%s）", estimator))
  between <- if (!is.null(ft$srmr_between) && !is.na(ft$srmr_between)) sprintf(", SRMR(between) = %s", fmt2(ft$srmr_between)) else ""
  sprintf("配適：χ²(%s) = %.2f, p %s, CFI = %s, TLI = %s, RMSEA = %s [%s, %s], SRMR = %s%s（%s，%s 指標）",
          format(round(ft$df, 2)), ft$chisq, fmt_p(ft$pvalue), fmt2(ft$cfi), fmt2(ft$tli), fmt2(ft$rmsea),
          fmt2(ft$rmsea_lo), fmt2(ft$rmsea_hi), fmt2(ft$srmr), between, estimator, ft$index_type)
}

# ---- 選模指標與門檻：全部步驟只引用這裡的數值 ------------------------------------------
# 內容與 references/09-fit-criteria.md、使用手冊「估計與模型比較」的兩張門檻表相同；改門檻時三處一起改。
# 整體配適門檻（單一模型是否可接受）：
#   RMSEA < .05 良好、.05–.08 可接受、.08–.10 勉強、≥ .10 不佳（Grimm et al., 2017；Browne & Cudeck, 1993；MacCallum et al., 1996）；
#   CFI、TLI > .95 良好、.90–.95 可接受、< .90 不佳（Grimm et al., 2017；Hu & Bentler, 1999）；SRMR ≤ .08（Hu & Bentler, 1999）。
#   自由度 ≤ 3 時 RMSEA 不納入總判定、也不單獨觸發提醒（Kenny et al., 2015）。
# 模型比較門檻：LRT p < .05（Grimm et al., 2017）；BIC 差距 < 2 微弱、2–6 有證據、6–10 強、> 10 非常強（Raftery, 1995）；
#   AIC 差距 ≤ 2 無法區分、> 10 幾乎不支持較差的模型（Burnham & Anderson, 2002）；CFI 差距 ≤ .01 沒有實質差異（Cheung & Rensvold, 2002）。
# 都是經驗法則，不是硬性標準。
FIT_RULES <- list(
  rmsea_good = .05, rmsea_ok = .08, rmsea_poor = .10,   # < .05 良好；< .08 可接受；< .10 勉強；≥ .10 不佳
  cfi_good = .95, cfi_poor = .90,            # CFI 與 TLI 共用：> .95 良好；.90–.95 可接受；< .90 不佳
  srmr_ok = .08,
  rmsea_min_df = 4,                          # 自由度 ≤ 3（< 4）時 RMSEA 僅供參考
  alpha = .05,                               # LRT
  bic_cut = c(2, 6, 10),                     # Raftery (1995) 的分界
  bic_lab = c("微弱", "有證據", "強", "非常強"),
  aic_indist = 2, aic_strong = 10,
  dcfi = .01,
  resid_mean_sd = .10                        # 單一模型的提醒：各波殘差平均最大絕對值 ≥ 各波 SD 平均的 .10 倍，視為有系統性偏離
)
# 各指標的等級
rmsea_level <- function(x) { r <- FIT_RULES; ifelse(is.na(x), NA, ifelse(x < r$rmsea_good, "良好", ifelse(x < r$rmsea_ok, "可接受", ifelse(x < r$rmsea_poor, "勉強", "不佳")))) }
cfi_level <- function(x) { r <- FIT_RULES; ifelse(is.na(x), NA, ifelse(x > r$cfi_good, "良好", ifelse(x >= r$cfi_poor, "可接受", "不佳"))) }
srmr_level <- function(x) ifelse(is.na(x), NA, ifelse(x <= FIT_RULES$srmr_ok, "良好", "未達"))
rmsea_counts <- function(df) is.na(df) | df >= FIT_RULES$rmsea_min_df
# 整體配適的總判定：所有指標都良好 → 良好；任一指標不佳 → 不佳；其餘 → 可接受（括號列出未達良好的指標）。
# 自由度 ≤ 3 時 RMSEA 不納入，並加註。
fit_verdict <- function(cfi, tli, rmsea, srmr, df = NA) {
  n <- max(length(cfi), length(df)); rep_n <- function(x) if (length(x) == n) x else rep(x, length.out = n)
  cfi <- rep_n(cfi); tli <- rep_n(tli); rmsea <- rep_n(rmsea); srmr <- rep_n(srmr); df <- rep_n(df)
  vapply(seq_len(n), function(k) {
    if (!is.na(df[k]) && df[k] == 0) return("恰好識別")
    if (is.na(cfi[k]) || is.na(rmsea[k])) return("無法判定")
    use_r <- rmsea_counts(df[k])
    lv <- c(CFI = cfi_level(cfi[k]), TLI = cfi_level(tli[k]), RMSEA = if (use_r) rmsea_level(rmsea[k]) else NA, SRMR = srmr_level(srmr[k]))
    lv <- lv[!is.na(lv)]
    v <- if (any(lv == "不佳")) "不佳" else if (all(lv == "良好")) "良好" else "可接受"
    notes <- c(if (v == "可接受") paste0("未達良好：", paste(names(lv)[lv != "良好"], collapse = "、")),
               if (!use_r) "自由度小，RMSEA 僅供參考")
    if (length(notes)) sprintf("%s（%s）", v, paste(notes, collapse = "；")) else v
  }, "")
}
# 確認點 1 直接指定單一模型時的配適提醒（CFI < .90 或 RMSEA ≥ .10；自由度 ≤ 3 時不因 RMSEA 提醒）
fit_poor_flags <- function(ft) {
  r <- FIT_RULES; f <- character()
  if (isTRUE(ft$df == 0)) return(f)
  if (isTRUE(ft$cfi < r$cfi_poor)) f <- c(f, sprintf("CFI = %s（< %s）", fmt2(ft$cfi), fmt2(r$cfi_poor, 2)))
  if (isTRUE(ft$rmsea >= r$rmsea_poor) && rmsea_counts(ft$df)) f <- c(f, sprintf("RMSEA = %s（≥ %s）", fmt2(ft$rmsea), fmt2(r$rmsea_poor, 2)))
  f
}
# BIC 差距的證據強度（Raftery, 1995）
bic_evidence <- function(d) { d <- abs(d); FIT_RULES$bic_lab[findInterval(d, FIT_RULES$bic_cut) + 1] }
aic_evidence <- function(d) { d <- abs(d); ifelse(d <= FIT_RULES$aic_indist, "無法區分", ifelse(d > FIT_RULES$aic_strong, "幾乎不支持較差的模型", "有差異")) }
# 兩個模型的比較結論（模型比較門檻）：simple 是較簡單（限制較多）的模型
compare_verdict <- function(f_simple, f_complex, lab_simple, lab_complex, nested = TRUE) {
  ts <- fit_table(f_simple); tc <- fit_table(f_complex)
  if (anyNA(c(ts$bic, tc$bic))) return("無法比較（模型未收斂）")
  dB <- ts$bic - tc$bic; dC <- round(tc$cfi - ts$cfi, 3) + 0
  out <- character()
  if (nested) {
    p <- tryCatch({ l <- as.data.frame(lavaan::lavTestLRT(f_simple, f_complex)); l[2, "Pr(>Chisq)"] }, error = function(e) NA)
    out <- c(out, if (is.na(p)) "LRT 無法計算" else sprintf("LRT p %s（%s）", fmt_p(p), if (p < FIT_RULES$alpha) paste0(lab_complex, "顯著較好") else "差異不顯著"))
  }
  out <- c(out, sprintf("ΔBIC = %.1f，支持%s（%s）", abs(dB), if (dB > 0) lab_complex else lab_simple, bic_evidence(dB)))
  if (!is.na(dC)) out <- c(out, sprintf("ΔCFI = %s%s", fmt2(dC), if (abs(dC) <= FIT_RULES$dcfi) "（≤ .01，沒有實質差異）" else ""))
  paste(out, collapse = "；")
}
# 門檻的一行說明（報告用）
fit_rules_text <- function() {
  r <- FIT_RULES
  sprintf("整體配適門檻：RMSEA < %s 良好、%s–%s 可接受、%s–%s 勉強、≥ %s 不佳；CFI、TLI > %s 良好、%s–%s 可接受、< %s 不佳；SRMR ≤ %s。自由度 ≤ %d 時 RMSEA 只供參考，不納入判定（Grimm et al., 2017；Browne & Cudeck, 1993；MacCallum et al., 1996；Hu & Bentler, 1999；Kenny et al., 2015）。這些是經驗法則，不是硬性標準。",
          fmt2(r$rmsea_good, 2), fmt2(r$rmsea_good, 2), fmt2(r$rmsea_ok, 2), fmt2(r$rmsea_ok, 2), fmt2(r$rmsea_poor, 2), fmt2(r$rmsea_poor, 2),
          fmt2(r$cfi_good, 2), fmt2(r$cfi_poor, 2), fmt2(r$cfi_good, 2), fmt2(r$cfi_poor, 2), fmt2(r$srmr_ok, 2), r$rmsea_min_df - 1)
}
# 將配適指標與參考值對照（只作參考，不是判定規則）
fit_flags <- function(ft) {
  if (isTRUE(ft$df == 0)) return("df = 0，不做配適判讀。")
  r <- FIT_RULES
  c(sprintf("CFI %s：%s", fmt2(ft$cfi), cfi_level(ft$cfi)),
    sprintf("TLI %s：%s", fmt2(ft$tli), cfi_level(ft$tli)),
    sprintf("RMSEA %s：%s%s", fmt2(ft$rmsea), rmsea_level(ft$rmsea), if (rmsea_counts(ft$df)) "" else "（自由度小，僅供參考）"),
    sprintf("SRMR %s：%s", fmt2(ft$srmr), if (srmr_level(ft$srmr) == "良好") "良好" else sprintf("> %s", fmt2(r$srmr_ok, 2))),
    sprintf("整體配適：%s（經驗法則，不是硬性標準）", fit_verdict(ft$cfi, ft$tli, ft$rmsea, ft$srmr, ft$df)))
}

# 自動路徑圖僅供檢查模型設定；正式論文的圖請另外繪製
save_diagram <- function(fit, file, layout = "tree2", structural = FALSE) {
  if (!requireNamespace("semPlot", quietly = TRUE)) return(invisible(FALSE))
  ok <- tryCatch({
    grDevices::png(file, width = 2400, height = 1600, res = 200)
    on.exit(grDevices::dev.off(), add = TRUE)
    semPlot::semPaths(fit, what = "std", layout = layout, edge.label.cex = .8, sizeMan = 4, sizeLat = 7,
                      residuals = FALSE, intercepts = FALSE, nCharNodes = 0, rotation = 2, fade = FALSE, structural = structural)
    TRUE
  }, error = function(e) { message("路徑圖略過：", conditionMessage(e)); FALSE })
  invisible(ok)
}

# ---- 縱貫 skill 專用 ---------------------------------------------------------
# 依欄名中的波次編號排序；可給逗號清單或正規表示式（以 ^ 開頭）
wave_vars <- function(x, data = NULL) {
  if (is.null(x)) return(NULL)
  if (startsWith(x, "^")) {
    v <- grep(x, names(data), value = TRUE)
    if (!length(v)) stop("找不到符合的欄位：", x, call. = FALSE)
    n <- as.numeric(sub(".*?([0-9]+)$", "\\1", v))
    return(v[order(n)])
  }
  split_csv(x)
}
check_cols <- function(v, data) {
  miss <- setdiff(v, names(data))
  if (length(miss)) stop("資料中找不到欄位：", paste(miss, collapse = ", "), call. = FALSE)
  invisible(v)
}
# 捕捉警告後估計
quiet_fit <- function(expr) {
  warn <- character()
  f <- withCallingHandlers(expr, warning = function(w) { warn <<- c(warn, gsub("\\s+", " ", conditionMessage(w))); invokeRestart("muffleWarning") })
  attr(f, "lsem_warn") <- unique(warn)
  f
}
# 模型是否可採用：有收斂且無負變異數／標準化解異常
admissible <- function(fit) {
  if (!lavaan::lavInspect(fit, "converged")) return(FALSE)
  !any(grepl("^\\[嚴重\\]", diagnose(fit)))
}
.lgm_state <- new.env()
.lgm_state$data_dirs <- character()
# 寫出重現模型所需的檔案：model.lav（完整語法）、call.txt（估計設定）；
# 每個步驟資料夾另存 analysis_data.rds（實際使用的資料）與 cases.txt（模型的順序），供 export_qmd.R 產生 qmd
save_replication <- function(fit, dir, sub, name) {
  if (!is.null(attr(fit, "syntax"))) writeLines(attr(fit, "syntax"), file.path(sub, "model.lav"), useBytes = TRUE)
  fa <- attr(fit, "fit_args")
  if (!is.null(fa)) writeLines(sprintf("%s=%s", names(fa), vapply(fa, function(v) paste(v, collapse = ","), "")), file.path(sub, "call.txt"), useBytes = TRUE)
  if (is.null(name) || is.null(attr(fit, "fit_data"))) return(invisible())
  key <- normalizePath(dir, mustWork = FALSE)
  if (!key %in% .lgm_state$data_dirs) {
    saveRDS(attr(fit, "fit_data"), file.path(dir, "analysis_data.rds"))
    writeLines(character(), file.path(dir, "cases.txt"))
    .lgm_state$data_dirs <- c(.lgm_state$data_dirs, key)
  }
  cat(name, file = file.path(dir, "cases.txt"), sep = "\n", append = TRUE)
}
save_case <- function(fit, dir, name = NULL) {
  sub <- if (is.null(name)) dir else ensure_dir(file.path(dir, name))
  save_replication(fit, dir, sub, name)
  writeLines(c(diagnose(fit), if (length(attr(fit, "lsem_warn"))) paste0("[lavaan 警告] ", attr(fit, "lsem_warn"))), file.path(sub, "diagnostics.txt"))
  if (!lavaan::lavInspect(fit, "converged")) return(invisible(NULL))
  wtab(fit_table(fit), sub, "fit_indices.csv")
  wtab(lavaan::parameterEstimates(fit, standardized = TRUE), sub, "parameters.csv")
  sink(file.path(sub, "summary.txt")); on.exit(sink(), add = TRUE)
  print(summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE))
  invisible(NULL)
}
# ---- 分析決定紀錄（確認點）：<results>/decisions.csv，一個確認點一列 ----------------------
# 欄位：point（1 篩檢、2 成長模型、3 預測變數）、checkpoint、recommendation、choice、differs（TRUE/FALSE/NA）、
#       reason（兩者不同時使用者的理由）、mode（使用者確認／自動）、updated
DECISION_POINTS <- c(`1` = "確認點 1：資料篩檢之後", `2` = "確認點 2：選擇成長模型", `3` = "確認點 3：是否加入預測變數", `4` = "遺失值的處理")
read_decisions <- function(dir) {
  f <- file.path(dir, "decisions.csv")
  if (!file.exists(f)) return(data.frame(point = integer(), checkpoint = character(), recommendation = character(), choice = character(),
                                         differs = logical(), reason = character(), mode = character(), updated = character(), stringsAsFactors = FALSE))
  x <- utils::read.csv(f, fileEncoding = "UTF-8", stringsAsFactors = FALSE, na.strings = "")
  for (k in c("recommendation", "choice", "reason", "mode")) x[[k]] <- as.character(x[[k]])
  x
}
# 只更新有給的欄位；同一個確認點只保留一列
write_decision <- function(dir, point, recommendation = NULL, choice = NULL, differs = NULL, reason = NULL, mode = NULL) {
  x <- read_decisions(dir); point <- as.integer(point)
  i <- match(point, x$point)
  if (is.na(i)) { x[nrow(x) + 1, "point"] <- point; i <- nrow(x); x$checkpoint[i] <- DECISION_POINTS[[as.character(point)]] }
  if (!is.null(recommendation)) x$recommendation[i] <- recommendation
  if (!is.null(choice)) x$choice[i] <- choice
  if (!is.null(differs)) x$differs[i] <- differs
  if (!is.null(reason)) x$reason[i] <- reason
  if (!is.null(mode)) x$mode[i] <- mode
  x$updated[i] <- format(Sys.time(), "%Y-%m-%d %H:%M")
  x <- x[order(x$point), ]
  utils::write.csv(x, file.path(dir, "decisions.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(x)
}

# prev：具名字元向量，模型名 → 要比較的前一個（較簡單的）模型名；有給時加上 vs_prev 欄（依模型比較門檻的比較結論）
# labels：顯示用的模型名稱（比較結論中使用）
fit_compare_table <- function(fits, kind = NULL, prev = NULL, labels = NULL) {
  ok <- vapply(fits, lavaan::lavInspect, TRUE, "converged")
  if (is.null(kind)) kind <- common_kind(fits[ok])
  tab <- do.call(rbind, lapply(names(fits), function(n) {
    if (!ok[[n]]) return(data.frame(model = n, chisq = NA, df = NA, pvalue = NA, cfi = NA, tli = NA, rmsea = NA, srmr = NA, aic = NA, bic = NA, admissible = FALSE))
    f <- fit_table(fits[[n]], kind)
    data.frame(model = n, chisq = f$chisq, df = f$df, pvalue = f$pvalue, cfi = f$cfi, tli = f$tli, rmsea = f$rmsea, srmr = f$srmr, aic = f$aic, bic = f$bic, admissible = admissible(fits[[n]]))
  }))
  tab$fit_verdict <- fit_verdict(tab$cfi, tab$tli, tab$rmsea, tab$srmr, tab$df)
  if (!is.null(prev)) {
    lab <- function(n) if (!is.null(labels) && n %in% names(labels)) labels[[n]] else n
    tab$vs_prev <- vapply(tab$model, function(n) {
      p <- if (n %in% names(prev)) prev[[n]] else NA
      if (is.na(p) || !p %in% names(fits) || !ok[[n]] || !ok[[p]]) return("—")
      paste0("vs ", lab(p), "：", compare_verdict(fits[[p]], fits[[n]], lab(p), lab(n)))
    }, "")
  }
  tab
}
# 表格中的模型名稱用中文（CSV 的 model 欄維持英文代碼，供程式讀取）
MODEL_LAB <- c(none = "無成長", linear = "線性", quad = "二次", basis = "潛在基底", free = "殘差變異自由估計", equal = "殘差變異相等")
model_label <- function(x) ifelse(x %in% names(MODEL_LAB), MODEL_LAB[x], x)
compare_md <- function(tab) {
  vp <- "vs_prev" %in% names(tab)
  tab$model <- model_label(tab$model)
  c(paste0("| 模型 | χ² | df | CFI | TLI | RMSEA | SRMR | AIC | BIC | 可採用 | 整體配適 |", if (vp) " 比較結論 |" else ""),
    paste0("|---|---|---|---|---|---|---|---|---|---|---|", if (vp) "---|" else ""),
    paste0(sprintf("| %s | %.2f | %.0f | %.3f | %.3f | %.3f | %.3f | %.1f | %.1f | %s | %s |", tab$model, tab$chisq, tab$df, tab$cfi, tab$tli, tab$rmsea, tab$srmr, tab$aic, tab$bic,
                   ifelse(tab$admissible, "是", "否（見 diagnostics.txt）"), tab$fit_verdict), if (vp) paste0(" ", tab$vs_prev, " |") else ""))
}
par_line <- function(pe, label = NULL) {
  lab <- if (is.null(label)) paste(pe$lhs, pe$op, pe$rhs) else label
  sprintf("- %s：b = %.3f, SE = %.3f, p %s, 95%% CI [%.3f, %.3f], β = %.3f", lab, pe$est, pe$se, fmt_p(pe$pvalue), pe$ci.lower, pe$ci.upper, pe$std.all)
}
# 類別色依固定順序使用（取自已驗證色盲區辨的預設色盤）
PAL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
line_chart <- function(file, x, series, xlab, ylab, main = NULL, band = NULL) {
  # series：具名 list，每個元素是與 x 等長的數值向量；圖例一律顯示，並於線尾直接標示
  grDevices::png(file, 1800, 1200, res = 220, bg = "#fcfcfb")
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- graphics::par(mar = c(4.2, 4.2, if (is.null(main)) 3 else 3.6, 8), las = 1, bty = "n", col.axis = "#555555", col.lab = "#333333", family = "sans")
  on.exit(graphics::par(op), add = TRUE)
  n <- names(series); yr <- range(unlist(series), na.rm = TRUE)
  graphics::plot(NA, xlim = range(x), ylim = yr, xlab = xlab, ylab = ylab, main = main, axes = FALSE, cex.main = 1, font.main = 1)
  graphics::abline(h = pretty(yr), col = "#e6e6e6", lwd = 1)
  graphics::axis(1, at = x, col = "#aaaaaa", tcl = -.25); graphics::axis(2, col = NA, tcl = 0)
  # 名為 Observed 的系列最後才畫，並用空心圓標示，才不會被重疊的模型線蓋住
  obs_k <- match("Observed", names(series))
  for (k in c(setdiff(seq_along(series), obs_k), obs_k[!is.na(obs_k)])) {
    if (identical(n[k], "Observed")) {
      graphics::points(x, series[[k]], pch = 21, bg = NA, col = PAL[k], cex = 1.9, lwd = 2.2)
    } else {
      graphics::lines(x, series[[k]], col = PAL[k], lwd = 2.4)
      graphics::points(x, series[[k]], pch = 21, bg = PAL[k], col = "white", cex = 1.25)
    }
  }
  last <- vapply(series, function(v) v[length(v)], 0)
  # 線尾標示：互相太近時往上下推開，避免重疊
  ord <- order(last); pos <- last[ord]; gap <- diff(graphics::par("usr")[3:4]) * .045
  for (j in seq_along(pos)[-1]) if (pos[j] - pos[j - 1] < gap) pos[j] <- pos[j - 1] + gap
  graphics::par(xpd = NA)
  if (length(series) > 1) {
    # 超過 4 條線時線尾標示容易擠在一起，只留圖例
    if (length(series) <= 4) graphics::text(max(x) + diff(range(x)) * .025, pos, n[ord], adj = 0, cex = .78, col = "#333333")
    is_obs <- n == "Observed"
    graphics::legend("top", inset = -.09, legend = n, col = PAL[seq_along(n)], lwd = ifelse(is_obs, NA, 2.4), pch = ifelse(is_obs, 1, 19), horiz = TRUE, bty = "n", cex = .78, text.col = "#333333")
  }
  invisible(NULL)
}
