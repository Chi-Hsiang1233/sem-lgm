# ============================================================
# lgm_core.R —— 潛在成長模型的語法產生、估計與診斷（各步驟腳本共用）
# 依 Grimm, Ram, & Estabrook (2017) 的模型設定
# ============================================================

# 載入同資料夾的 lgm_helpers.R
# 搜尋順序：LGM_SCRIPTS_DIR（在 qmd / 互動式 R 中先設定）→ Rscript 的 --file → 工作目錄 → 工作目錄/scripts
.lgm_dir <- (function() {
  cand <- character()
  if (exists("LGM_SCRIPTS_DIR")) cand <- c(cand, get("LGM_SCRIPTS_DIR"))
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) cand <- c(cand, dirname(normalizePath(f[1], mustWork = FALSE)))
  cand <- c(cand, getwd(), file.path(getwd(), "scripts"))
  hit <- cand[file.exists(file.path(cand, "lgm_helpers.R"))]
  if (!length(hit)) stop("找不到 lgm_helpers.R；在互動式 R 中請先設定 LGM_SCRIPTS_DIR <- \"<skill>/scripts\"", call. = FALSE)
  normalizePath(hit[1])
})()
source(file.path(.lgm_dir, "lgm_helpers.R"), encoding = "UTF-8")

GF_NAMES <- c("i", "s", "q")

# 截距位置（anchor）：時間編碼中哪一個時間點當作 0，決定截距的意義
#   first（預設）= 第一波的水準；last = 最後一波的水準；mid = 追蹤期中點；也可給數值（例如 2）
#   移動截距位置不會改變模型配適，只改變截距平均、截距變異與截距–斜率共變的意義
resolve_anchor <- function(x, times) {
  if (is.null(x) || (length(x) == 1 && is.na(x)) || identical(x, "first")) return(times[1])
  if (identical(x, "last")) return(times[length(times)])
  if (x %in% c("mid", "center")) return(mean(range(times)))
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) stop("無法解析 --anchor：", x, "（可用 first、last、mid 或數值）", call. = FALSE)
  v
}
anchor_label <- function(anchor, times) {
  if (isTRUE(all.equal(anchor, times[1]))) "第一波" else if (isTRUE(all.equal(anchor, times[length(times)]))) "最後一波" else sprintf("時間 = %s", format(anchor))
}

# 讀取步驟二的 chosen.txt（--from），回傳具名字元向量；沒有時回傳 character()
# （必須是 character() 而不是 c()/NULL：character()["form"] 是 NA，NULL["form"] 會讓 if() 出錯）
read_prev <- function(a) {
  if (is.null(a$from) || !file.exists(file.path(a$from, "chosen.txt"))) return(character())
  kv <- strsplit(readLines(file.path(a$from, "chosen.txt"), encoding = "UTF-8"), "=")
  setNames(vapply(kv, `[`, "", 2), vapply(kv, `[`, "", 1))
}

# 取得設定值：命令列 > chosen.txt > 預設
prev_or <- function(a, prev, key, default = NULL) {
  if (!is.null(a[[key]])) return(a[[key]])
  if (!is.na(prev[key])) return(prev[[key]])
  default
}

# 時間負荷：時間減去 anchor（預設第一個時間點），因此 anchor 那一波的負荷為 0
time_loadings <- function(times, form, anchor = times[1]) {
  t0 <- times - anchor
  switch(form,
    none = list(),
    linear = list(s = t0),
    quad = list(s = t0, q = t0^2),
    basis = list(s = c(t0[1], rep(NA, length(t0) - 2), t0[length(t0)])),
    stop("未知的軌跡形式：", form, "（可用 none、linear、quad、basis）", call. = FALSE))
}

# 產生成長模型語法
# resid: free（殘差變異自由估計）| equal（殘差變異相等，對應多層次模型的同質變異）
# cov_paths: free | zero（預測變數路徑固定為 0，作為階層比較的基準模型）
# tvc_mode: equal（跨波相等）| free（各波自由）；tvc_cor = TRUE 讓時變共變數與成長因子相關
lgm_syntax <- function(y, times, form = "linear", resid = "free", covs = NULL, cov_paths = "free",
                       tvc = NULL, tvc_mode = "equal", tvc_cor = FALSE, distal = NULL, fix_var = NULL,
                       anchor = times[1]) {
  Tn <- length(y)
  stopifnot(length(times) == Tn, is.null(tvc) || length(tvc) == Tn)
  L <- time_loadings(times, form, anchor)
  s <- paste("i =~", paste0("1*", y, collapse = " + "))
  for (g in names(L)) {
    ld <- ifelse(is.na(L[[g]]), paste0("NA*", y), paste0(format(L[[g]], trim = TRUE, scientific = FALSE, drop0trailing = TRUE), "*", y))
    s <- c(s, paste(g, "=~", paste(ld, collapse = " + ")))
  }
  gf <- c("i", names(L))
  if (resid == "equal") s <- c(s, sprintf("%s ~~ theta*%s", y, y))
  for (v in fix_var) s <- c(s, sprintf("%s ~~ 0*%s", v, v), sprintf("%s ~~ 0*%s", v, setdiff(gf, v)))
  if (length(covs)) {
    pre <- if (cov_paths == "zero") "0*" else ""
    s <- c(s, paste0(gf, " ~ ", vapply(gf, function(g) paste0(pre, covs, collapse = " + "), "")))
  }
  if (length(tvc)) {
    lab <- if (tvc_mode == "free") paste0("btvc", seq_len(Tn)) else rep("btvc", Tn)
    s <- c(s, sprintf("%s ~ %s*%s", y, lab, tvc))
    if (tvc_cor) s <- c(s, as.vector(outer(gf, tvc, function(g, c) sprintf("%s ~~ %s", g, c))))
  }
  if (length(distal)) s <- c(s, paste0(distal, " ~ ", paste(gf, collapse = " + ")), paste0(distal, " ~ 1"))
  paste(s, collapse = "\n")
}

exo_free <- function(exo) {
  s <- sprintf("%s ~ 1", exo)
  if (length(exo) > 1) { p <- utils::combn(exo, 2); s <- c(s, sprintf("%s ~~ %s", p[1, ], p[2, ])) }
  s
}

# 估計：外生變數有遺失或時變共變數要與成長因子相關時，改用 fixed.x = FALSE，
# 並讓這些變數的平均自由估計（growth() 預設把觀察變數截距固定為 0）
fit_lgm <- function(syntax, data, exo = NULL, estimator = "MLR", force_random_x = FALSE, ...) {
  random_x <- force_random_x || (length(exo) && anyNA(data[, exo, drop = FALSE]))
  # fixed.x = FALSE 時，growth() 不會自動估計外生變數的平均與彼此的共變（會固定為 0），要明確列出
  if (random_x && length(exo)) syntax <- paste(c(syntax, exo_free(exo)), collapse = "\n")
  f <- quiet_fit(lavaan::growth(syntax, data, estimator = estimator, missing = "fiml", fixed.x = !random_x, ...))
  attr(f, "random_x") <- random_x
  attr(f, "syntax") <- syntax
  # 匯出 qmd（export_qmd.R）用：記下估計設定與實際使用的資料，save_case() 會寫到輸出資料夾
  attr(f, "fit_args") <- c(list(estimator = estimator, missing = "fiml", fixed.x = !random_x), list(...))
  attr(f, "fit_data") <- data
  f
}

# 成長因子的平均、變異、相關（可多群組）
growth_params <- function(fit) {
  pe <- lavaan::parameterEstimates(fit, standardized = TRUE)
  if (is.null(pe$group)) pe$group <- 1
  gf <- intersect(GF_NAMES, lavaan::lavNames(fit, "lv"))
  list(mean = pe[pe$op == "~1" & pe$lhs %in% gf, ],
       var = pe[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% gf, ],
       cor = pe[pe$op == "~~" & pe$lhs != pe$rhs & pe$lhs %in% gf & pe$rhs %in% gf, ],
       reg = pe[pe$op == "~", ],
       load = pe[pe$op == "=~" & pe$lhs != "i", ],
       all = pe)
}

lrt_line <- function(f0, f1, label) {
  l <- tryCatch(as.data.frame(lavaan::lavTestLRT(f0, f1)), error = function(e) NULL)
  if (is.null(l) || nrow(l) < 2) return(sprintf("- %s：無法計算 LRT", label))
  sprintf("- %s：Δχ²(%s) = %.2f, p %s", label, l[2, "Df diff"], l[2, "Chisq diff"], fmt_p(l[2, "Pr(>Chisq)"]))
}
lrt_p <- function(f0, f1) { l <- tryCatch(as.data.frame(lavaan::lavTestLRT(f0, f1)), error = function(e) NULL); if (is.null(l)) NA else l[2, "Pr(>Chisq)"] }

# 殘差診斷：以因子分數得到個人預測軌跡與殘差
# - 各波殘差平均偏離 0 且有系統趨勢 → 未建模的非線性
# - 各波殘差 SD 差異大 → 殘差變異跨波不等
# - 成長因子分數的偏態峰度 → 常態假設；極端分數 → 可能的離群個案
# 因子分數有縮減（shrinkage），只作為檢查工具，不要拿來做後續推論
lgm_diagnostics <- function(fit, data, y, times, out, tag = "") {
  yhat <- tryCatch(lavaan::lavPredict(fit, type = "ov"), error = function(e) NULL)
  fs <- tryCatch(lavaan::lavPredict(fit, type = "lv"), error = function(e) NULL)
  if (is.null(yhat) || is.null(fs)) return("（無法計算因子分數，略過殘差診斷）")
  idx <- lavaan::lavInspect(fit, "case.idx")
  Y <- as.matrix(data[idx, y]); R <- Y - yhat[, y, drop = FALSE]
  rs <- data.frame(time = times, mean_resid = colMeans(R, na.rm = TRUE), sd_resid = apply(R, 2, sd, na.rm = TRUE))
  wtab(rs, out, paste0("residuals_by_wave", tag, ".csv"))
  gf <- intersect(GF_NAMES, colnames(fs))
  sk <- function(v) mean((v - mean(v))^3) / sd(v)^3; ku <- function(v) mean((v - mean(v))^4) / sd(v)^4 - 3
  fsd <- data.frame(factor = gf, mean = colMeans(fs[, gf, drop = FALSE]), sd = apply(fs[, gf, drop = FALSE], 2, sd),
                    skew = apply(fs[, gf, drop = FALSE], 2, sk), kurtosis = apply(fs[, gf, drop = FALSE], 2, ku))
  wtab(fsd, out, paste0("factor_score_distribution", tag, ".csv"))
  md <- stats::mahalanobis(fs[, gf, drop = FALSE], colMeans(fs[, gf, drop = FALSE]), stats::cov(fs[, gf, drop = FALSE]))
  flag <- which(stats::pchisq(md, length(gf), lower.tail = FALSE) < .001)
  if (length(flag)) wtab(data.frame(row = idx[flag], mahalanobis = md[flag], fs[flag, gf, drop = FALSE]), out, paste0("factor_score_outliers", tag, ".csv"))

  # 圖 1：個人預測軌跡（抽樣）＋ 平均預測軌跡
  set.seed(1); pick <- sample(nrow(yhat), min(60, nrow(yhat)))
  grDevices::png(file.path(out, paste0("predicted_trajectories", tag, ".png")), 1800, 1200, res = 220, bg = "#fcfcfb")
  graphics::par(mar = c(4.2, 4.2, 2.4, 1), las = 1, bty = "n", col.axis = "#555555", col.lab = "#333333")
  graphics::plot(NA, xlim = range(times), ylim = range(yhat[pick, y], na.rm = TRUE), xlab = "Time", ylab = "Model-predicted score", axes = FALSE,
                 main = "Individual predicted trajectories (sample of 60)", cex.main = .95, font.main = 1)
  graphics::abline(h = pretty(range(yhat[pick, y], na.rm = TRUE)), col = "#e6e6e6")
  graphics::axis(1, at = times, col = "#aaaaaa", tcl = -.25); graphics::axis(2, col = NA, tcl = 0)
  for (k in pick) graphics::lines(times, yhat[k, y], col = "#2a78d633")
  graphics::lines(times, colMeans(yhat[, y]), col = "#2a78d6", lwd = 3)
  grDevices::dev.off()

  # 圖 2：殘差散布與各波平均（應接近 0、且沒有趨勢）
  grDevices::png(file.path(out, paste0("residual_plot", tag, ".png")), 1800, 1200, res = 220, bg = "#fcfcfb")
  graphics::par(mar = c(4.2, 4.2, 2.4, 1), las = 1, bty = "n", col.axis = "#555555", col.lab = "#333333")
  jit <- stats::runif(length(R), -.08, .08) * diff(range(times))
  graphics::plot(rep(times, each = nrow(R)) + jit, as.vector(R), pch = 16, cex = .45, col = "#2a78d640", xlab = "Time", ylab = "Residual (observed − predicted)", axes = FALSE,
                 main = "Residuals by time (line = mean residual)", cex.main = .95, font.main = 1)
  graphics::abline(h = 0, col = "#888888", lty = 2)
  graphics::axis(1, at = times, col = "#aaaaaa", tcl = -.25); graphics::axis(2, col = NA, tcl = 0)
  graphics::lines(times, rs$mean_resid, col = "#eb6834", lwd = 2.5); graphics::points(times, rs$mean_resid, pch = 21, bg = "#eb6834", col = "white", cex = 1.3)
  grDevices::dev.off()

  trend <- stats::cor(times, rs$mean_resid)
  c(sprintf("- 各波殘差平均：%s（應接近 0 且沒有規律）", paste(sprintf("%.2f", rs$mean_resid), collapse = "、")),
    sprintf("- 各波殘差 SD：%s（差異大代表殘差變異跨波不等）", paste(sprintf("%.2f", rs$sd_resid), collapse = "、")),
    sprintf("- 成長因子分數：%s", paste(sprintf("%s 偏態 %.2f、峰度 %.2f", fsd$factor, fsd$skew, fsd$kurtosis), collapse = "；")),
    sprintf("- 成長因子分數的多變量離群值（p < .001）：%d 人%s", length(flag), if (length(flag)) "，見 factor_score_outliers.csv" else ""),
    "- 因子分數有縮減，只作為檢查工具，不要拿來做後續推論。")
}
