# step1_screen.R —— 步驟一：資料篩檢：各波描述、各波信度、軌跡圖、遺失值檢查（遺失比例、涵蓋率、MCAR 檢定、流失者比較）、穩定性、ICC
# Rscript step1_screen.R --data wide.csv --y "^dep_t" [--times 0,6,12,24] [--baseline age_c,female,treat] [--aux a1,a2] [--items "^dep_t[0-9]+_i"] [--na -99] --out results/screen
#   報告順序：各波信度 → 遺失機制檢驗 → 各波描述與軌跡（含「從篩檢結果看成長形式」）→ 跨波穩定性與 ICC → 簡易判讀
#   --times：實際測量時間，只用來描述「每單位時間」的平均變化（不給時以波為單位）
#   --items：各波的題目欄位（正規表示式或逗號清單）。比對到的那一段裡的第一個數字就是波次，例如 dep_t2_i3 → 第 2 波；
#            依波次排序後與 --y 的各波對應。沒有提供時略過信度，報告寫「只有總分，無法計算內部一致性」。
#   --aux：輔助變數（不進分析模型，但與遺失或結果有關），一併納入流失者比較與 MCAR 檢定；多重插補時也會用到。
# 多個重複變數時，分別執行。
source(file.path((function() { a <- commandArgs(FALSE); dirname(normalizePath(sub("^--file=", "", a[grep("^--file=", a)][1]))) })(), "lgm_helpers.R"))
a <- parse_args(); need(a, c("data", "y"))
out <- ensure_dir(arg_or(a, "out", "results/step1"))
d <- read_any(a$data, na = arg_or(a, "na"))
y <- check_cols(wave_vars(a$y, d), d); base <- split_csv(a$baseline); Tn <- length(y); Y <- d[, y]
aux <- check_cols(split_csv(a$aux), d)
has_lavaan <- requireNamespace("lavaan", quietly = TRUE)
rep <- c("# 步驟一：資料篩檢", "", sprintf("N = %d；%d 波：%s", nrow(d), Tn, paste(y, collapse = ", ")))

# 各波描述
sk <- function(v) { v <- v[!is.na(v)]; mean((v - mean(v))^3) / sd(v)^3 }
ku <- function(v) { v <- v[!is.na(v)]; mean((v - mean(v))^4) / sd(v)^4 - 3 }
w <- data.frame(wave = y, n = colSums(!is.na(Y)), retention_pct = 100 * colSums(!is.na(Y)) / sum(!is.na(Y[[1]])),
                mean = colMeans(Y, na.rm = TRUE), sd = apply(Y, 2, sd, na.rm = TRUE),
                skew = apply(Y, 2, sk), kurtosis = apply(Y, 2, ku))
wtab(w, out, "wave_descriptives.csv", 3)
seg_desc <- c("", "**各波描述**", "", "| 波次 | n | 保留率% | M | SD | 偏態 | 峰度 |", "|---|---|---|---|---|---|---|",
         sprintf("| %s | %d | %.1f | %.2f | %.2f | %.2f | %.2f |", w$wave, w$n, w$retention_pct, w$mean, w$sd, w$skew, w$kurtosis))
# 從篩檢結果看成長形式：只描述各波平均的走勢與 SD 的變化，不下結論（形式在步驟二決定）
tm <- if (!is.null(a$times)) as.numeric(split_csv(a$times)) else seq_len(Tn) - 1
if (length(tm) != Tn) stop("--times 的個數和 --y 的波數不同。", call. = FALSE)
rate <- diff(w$mean) / diff(tm); unit <- if (!is.null(a$times)) "每單位時間" else "每波"
shape <- if (Tn < 3) "只有兩波，看不出走勢的形狀" else if (max(abs(rate)) < .1 * mean(w$sd)) "各段的變化都很小，走勢接近水平" else
  if (!(all(rate > 0) || all(rate < 0))) "變化方向不一致，走勢有升有降" else if (diff(range(rate)) <= .25 * mean(abs(rate))) "各段的變化速度相近，走勢接近直線" else
  if (all(diff(abs(rate)) < 0)) "變化速度逐段變小，走勢逐漸趨緩（彎曲）" else if (all(diff(abs(rate)) > 0)) "變化速度逐段變大，走勢逐漸加速（彎曲）" else "變化速度時快時慢，走勢不規則"
sd_txt <- if (w$sd[Tn] > 1.2 * w$sd[1] && all(diff(w$sd) > 0)) "逐波擴大" else if (w$sd[Tn] > 1.2 * w$sd[1]) "整體變大，但不是每一波都變大" else if (w$sd[Tn] < w$sd[1] / 1.2) "逐漸變小" else "大致不變"
seg_shape <- c("", "**從篩檢結果看成長形式**", "",
               sprintf("- 各波平均依序為 %s；%s的變化為 %s：%s。", paste(sprintf("%.2f", w$mean), collapse = "、"), unit, paste(sprintf("%+.2f", rate), collapse = "、"), shape),
               sprintf("- 標準差從 %.2f 到 %.2f，%s。", w$sd[1], w$sd[Tn], sd_txt),
               "- 這裡只描述觀察到的平均與標準差；成長形式要在步驟二用模型比較決定。")

# 各波信度（選用：需要題目層級資料）
# α 與 Feldt (1965) 95% CI；ω 由單因子模型（lavaan，因子變異固定為 1）計算 (Σλ)² / ((Σλ)² + Σθ)；
# 校正後題目–總分相關 = 題目與「其他題目總和」的相關。各波以該波題目都有作答的人計算。
rel <- NULL; rel_low_items <- character()
if (!is.null(a$items)) {
  icol <- if (startsWith(a$items, "^")) grep(a$items, names(d), value = TRUE) else check_cols(split_csv(a$items), d)
  if (!length(icol)) stop("找不到符合 --items 的欄位：", a$items, call. = FALSE)
  hit_part <- if (startsWith(a$items, "^")) regmatches(icol, regexpr(a$items, icol)) else icol
  wnum <- suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", ifelse(grepl("[0-9]", hit_part), hit_part, icol))))
  if (anyNA(wnum)) stop("無法從 --items 的欄位名稱判斷波次：", paste(icol[is.na(wnum)], collapse = ", "), call. = FALSE)
  sets <- split(icol, wnum); sets <- sets[order(as.numeric(names(sets)))]
  lab <- if (length(sets) == Tn) y else paste0("第 ", names(sets), " 組")
  rel <- do.call(rbind, lapply(seq_along(sets), function(k) {
    X <- d[, sets[[k]], drop = FALSE]; X <- X[stats::complete.cases(X), , drop = FALSE]; ki <- ncol(X); n <- nrow(X)
    if (ki < 2 || n < 10) return(data.frame(wave = lab[k], n_items = ki, n = n, alpha = NA, alpha_lo = NA, alpha_hi = NA, omega = NA, min_citc = NA, min_citc_item = NA))
    tot <- rowSums(X); al <- ki / (ki - 1) * (1 - sum(apply(X, 2, stats::var)) / stats::var(tot))
    df1 <- n - 1; df2 <- (n - 1) * (ki - 1)
    ci <- c(1 - (1 - al) * stats::qf(.975, df1, df2), 1 - (1 - al) * stats::qf(.025, df1, df2))
    om <- NA
    if (has_lavaan && ki >= 3) {
      cf <- tryCatch(suppressWarnings(lavaan::cfa(paste("f =~", paste(names(X), collapse = " + ")), X, std.lv = TRUE)), error = function(e) NULL)
      if (!is.null(cf) && lavaan::lavInspect(cf, "converged")) {
        pe <- lavaan::parameterEstimates(cf); ld <- pe$est[pe$op == "=~"]; th <- pe$est[pe$op == "~~" & pe$lhs %in% names(X) & pe$lhs == pe$rhs]
        if (all(th > 0)) om <- sum(ld)^2 / (sum(ld)^2 + sum(th))
      }
    }
    citc <- vapply(names(X), function(v) stats::cor(X[[v]], tot - X[[v]]), 0)
    low <- names(citc)[citc < .30]; if (length(low)) rel_low_items <<- c(rel_low_items, sprintf("%s（%s）", paste(low, collapse = "、"), lab[k]))
    data.frame(wave = lab[k], n_items = ki, n = n, alpha = al, alpha_lo = ci[1], alpha_hi = ci[2], omega = om,
               min_citc = min(citc), min_citc_item = names(citc)[which.min(citc)])
  }))
  wtab(rel, out, "reliability.csv", 3)
  seg_rel <- c("", "**各波信度**", "", "| 波次 | 題數 | n | α | α 95% CI | ω | 最低題目–總分相關 |", "|---|---|---|---|---|---|---|",
           sprintf("| %s | %d | %d | %s | [%s, %s] | %s | %s（%s） |", rel$wave, rel$n_items, rel$n, fmt2(rel$alpha), fmt2(rel$alpha_lo), fmt2(rel$alpha_hi),
                   ifelse(is.na(rel$omega), "—", fmt2(rel$omega)), fmt2(rel$min_citc), rel$min_citc_item),
           "", "α：Cronbach's α（95% CI 依 Feldt 法）；ω：McDonald's ω（單因子模型）；題目–總分相關已校正（不含該題）。各波以該波題目都有作答的人計算。",
           if (length(sets) != Tn) sprintf("⚠ 題目分成 %d 組，但 --y 有 %d 波，無法逐波對應，請確認 --items。", length(sets), Tn))
} else {
  seg_rel <- c("", "**各波信度**：只有總分，無法計算內部一致性（有題目資料時用 `--items` 指定各波的題目欄位）。")
}

# 遺失型態與單調性
pat <- apply(!is.na(Y), 1, function(r) paste(ifelse(r, "X", "."), collapse = ""))
pt <- sort(table(pat), decreasing = TRUE)
wtab(data.frame(pattern = names(pt), n = as.integer(pt), pct = round(100 * as.numeric(pt) / nrow(d), 1)), out, "missing_patterns.csv")
mono <- vapply(names(pt), function(p) !grepl("\\.X", p), TRUE)
rep_body <- rep; rep <- character()   # 遺失機制檢驗的段落先收在 rep，最後再依順序組合
rep <- c(rep, "", "**遺失型態**", "",
         if (all(complete.cases(Y))) sprintf("- 沒有遺失值（%d 人全部完整作答）", nrow(d)) else
         c(sprintf("- 共 %d 種型態；完整作答 %.1f%%；單調遺失（中斷後不再回來）型態占 %.1f%%", length(pt), 100 * mean(complete.cases(Y)), 100 * sum(pt[mono]) / nrow(d)),
           sprintf("- 最常見：%s", paste(sprintf("%s (%d)", names(pt)[seq_len(min(4, length(pt)))], pt[seq_len(min(4, length(pt)))]), collapse = "、"))))

# 各變數的遺失比例、共變數涵蓋率、Little (1988) 的 MCAR 檢定
mvars <- unique(c(y, base, aux)); mvars <- mvars[vapply(d[mvars], is.numeric, TRUE)]
M <- d[, mvars, drop = FALSE]
mp <- data.frame(variable = mvars, role = ifelse(mvars %in% y, "重複測量", ifelse(mvars %in% aux, "輔助變數", "基線變數")),
                 n = colSums(!is.na(M)), missing_pct = 100 * colMeans(is.na(M)))
wtab(mp, out, "missing_by_variable.csv", 2)
cov_m <- crossprod(!is.na(as.matrix(M))) / nrow(M)
wtab(cbind(variable = mvars, as.data.frame(round(cov_m, 3))), out, "covariance_coverage.csv", 3)
mcar <- NULL
if (!anyNA(M)) {
  rep <- c(rep, "", sprintf("**遺失比例與 MCAR 檢定**：%s 都沒有遺失值；涵蓋率都是 100%%，不需要做 MCAR 檢定。", paste(mvars, collapse = "、")))
} else {
  lowtri <- which(lower.tri(cov_m), arr.ind = TRUE); mn <- lowtri[which.min(cov_m[lowtri]), ]
  rep <- c(rep, "", "**各變數的遺失比例**", "", "| 變數 | 角色 | n | 遺失 % |", "|---|---|---|---|",
           sprintf("| %s | %s | %d | %.1f |", mp$variable, mp$role, mp$n, mp$missing_pct), "",
           sprintf("- 共變數涵蓋率（每一對變數同時有資料的比例）：最低 %.1f%%（%s 與 %s），完整矩陣見 covariance_coverage.csv。涵蓋率太低（例如 < 10%%）時，該共變數幾乎只靠模型推估。",
                   100 * cov_m[mn[1], mn[2]], mvars[mn[1]], mvars[mn[2]]))
  # Little's MCAR test：以 ML（FIML）估計的平均與共變數為基準，比較各遺失型態的觀察平均
  if (has_lavaan) {
    mc <- tryCatch({
      Mk <- M[rowSums(!is.na(M)) > 0, , drop = FALSE]
      sat <- lavaan::lavCor(Mk, missing = "ml", output = "fit")
      im <- lavaan::lavInspect(sat, "implied"); S <- im$cov[mvars, mvars]; mu <- im$mean[mvars]
      R <- !is.na(Mk); key <- apply(R, 1, function(r) paste(as.integer(r), collapse = ""))
      d2 <- 0; dfm <- 0
      for (g in unique(key)) {
        o <- R[match(g, key), ]; rows <- Mk[key == g, o, drop = FALSE]
        dv <- colMeans(rows) - mu[o]
        d2 <- d2 + nrow(rows) * as.numeric(t(dv) %*% solve(S[o, o, drop = FALSE]) %*% dv); dfm <- dfm + sum(o)
      }
      dfm <- dfm - length(mvars)
      data.frame(chisq = d2, df = dfm, p = stats::pchisq(d2, dfm, lower.tail = FALSE), n_patterns = length(unique(key)), variables = paste(mvars, collapse = ","))
    }, error = function(e) NULL)
    if (!is.null(mc)) {
      mcar <- mc; wtab(mc, out, "mcar_test.csv", 4)
      rep <- c(rep, sprintf("- Little 的 MCAR 檢定（變數：%s）：χ²(%d) = %.2f, p %s。", paste(mvars, collapse = "、"), mc$df, mc$chisq, fmt_p(mc$p)),
               "  顯著只代表遺失不是完全隨機（MCAR）；樣本大時很小的差異也會顯著；這個檢定無法分辨 MAR 與 MNAR。")
    } else rep <- c(rep, "- Little 的 MCAR 檢定：無法計算（共變數矩陣無法估計）。")
  } else rep <- c(rep, "- Little 的 MCAR 檢定：需要 lavaan 套件，已略過。")
}

# 遺失值分析：完成者 vs 有遺失者
dropout <- is.na(Y[[Tn]])
if (any(dropout) && !all(dropout)) {
  cand <- unique(c(y[1], base, aux))
  at <- do.call(rbind, lapply(cand, function(v) {
    x <- suppressWarnings(as.numeric(d[[v]])); a1 <- x[!dropout]; a0 <- x[dropout]
    if (sum(!is.na(a0)) < 3) return(NULL)
    tt <- stats::t.test(a0, a1); sp <- sqrt((var(a0, na.rm = TRUE) + var(a1, na.rm = TRUE)) / 2)
    data.frame(variable = v, mean_dropout = mean(a0, na.rm = TRUE), mean_completer = mean(a1, na.rm = TRUE), d = (mean(a0, na.rm = TRUE) - mean(a1, na.rm = TRUE)) / sp, p = tt$p.value)
  }))
  wtab(at, out, "attrition_comparison.csv", 3)
  rep <- c(rep, "", sprintf("**遺失值分析**（最後一波遺失 %d 人，%.1f%%）", sum(dropout), 100 * mean(dropout)), "",
           "| 變數 | 有遺失者 M | 完成者 M | d | p |", "|---|---|---|---|---|",
           sprintf("| %s | %.2f | %.2f | %.2f | %s |", at$variable, at$mean_dropout, at$mean_completer, at$d, sub("= ", "", fmt_p(at$p))))
  # 前一波分數是否預測下一波遺失（MAR 的直接證據）
  if (Tn >= 3) {
    rows <- list()
    for (t in 2:Tn) {
      inprev <- !is.na(Y[[t - 1]])
      if (sum(is.na(Y[[t]][inprev])) < 5) next
      m <- stats::glm(is.na(Y[[t]][inprev]) ~ scale(Y[[t - 1]][inprev]), family = stats::binomial)
      co <- summary(m)$coefficients[2, ]
      rows[[length(rows) + 1]] <- data.frame(transition = paste(y[t - 1], "→", y[t]), OR_per_SD = exp(co[1]), p = co[4])
    }
    if (length(rows)) {
      tr <- do.call(rbind, rows); wtab(tr, out, "dropout_by_prior_score.csv", 3)
      rep <- c(rep, "", "前一波分數（每 1 SD）對下一波遺失的勝算比：", "", sprintf("- %s：OR = %.2f, p %s", tr$transition, tr$OR_per_SD, fmt_p(tr$p)),
               "", "若前一波分數預測遺失，代表遺失不是 MCAR，列刪除會偏誤；FIML 在 MAR 下仍可用。資料無法排除 MNAR，應做敏感度分析。")
    }
  }
}

# 穩定性
seg_miss <- rep; rep <- rep_body
R <- cor(Y, use = "pairwise")
wtab(cbind(wave = y, as.data.frame(round(R, 3))), out, "stability_correlations.csv", 3)
lag1 <- R[cbind(2:Tn, 1:(Tn - 1))]
far <- R[Tn, 1]
seg_stab <- c("", "**跨波穩定性**", "", sprintf("- 相鄰波相關：%s", paste(sprintf("%.2f", lag1), collapse = "、")), sprintf("- 第 1 波與最後一波相關：%.2f", far),
         sprintf("- 相關隨間隔%s（單純 AR(1) 會隨間隔快速衰減；不太衰減暗示有穩定特質，應考慮 RI-CLPM 或成長模型）", if (far < min(lag1) - .1) "明顯下降" else "幾乎不下降"))

# 個體間比例（ICC）：以單因子變異數分析估計
long <- data.frame(id = rep(seq_len(nrow(Y)), Tn), v = unlist(Y, use.names = FALSE)); long <- long[!is.na(long$v), ]
ni <- table(long$id); k <- (length(long$v) - sum(ni^2) / length(long$v)) / (length(ni) - 1)
msb <- sum(ni * (tapply(long$v, long$id, mean) - mean(long$v))^2) / (length(ni) - 1)
msw <- sum((long$v - ave(long$v, long$id))^2) / (length(long$v) - length(ni))
icc <- (msb - msw) / (msb + (k - 1) * msw)
seg_stab <- c(seg_stab, "", sprintf("- ICC（個體間變異比例）= %.3f；個體內變異比例 = %.3f", icc, 1 - icc))
# 報告順序：信度 → 遺失機制檢驗 → 各波描述與軌跡（含成長形式的描述）→ 跨波穩定性與 ICC → 簡易判讀
rep <- c(rep, seg_rel, seg_miss, seg_desc, "", "軌跡圖：trajectories.png（抽樣 60 人的個體軌跡與整體平均）。", seg_shape, seg_stab)

# 簡易判讀：給使用者確認資料後，再決定要不要進入步驟二
ok <- character(); warn <- character()
if (nrow(d) < 100) {
  warn <- c(warn, sprintf("樣本只有 %d 人，模型可能不穩定（容易不收斂或出現負變異數），建議只估較簡單的形式。", nrow(d)))
} else ok <- c(ok, sprintf("樣本 %d 人、%d 波，可以估計成長模型%s。", nrow(d), Tn, if (Tn == 3) "（只有 3 波，只能估線性）" else ""))
codes <- c(-999, -99, -98, -9, -1, 99, 999, 9999)
hit <- codes[vapply(codes, function(cd) any(unlist(Y) == cd, na.rm = TRUE), TRUE)]
if (length(hit)) warn <- c(warn, sprintf("分數中出現 %s，可能是遺失值代碼；若是，請用 --na 指定後重跑。", paste(hit, collapse = "、")))
miss_cell <- 100 * mean(is.na(as.matrix(Y)))
if (miss_cell == 0) {
  ok <- c(ok, "沒有遺失值，不需要做 MCAR 檢定或多重插補。")
} else {
  msg <- sprintf("重複測量有 %.1f%% 的格子遺失，最後一波保留 %.1f%%；分析會用 FIML 保留所有人，前提是 MAR。", miss_cell, w$retention_pct[Tn])
  if (w$retention_pct[Tn] < 70) warn <- c(warn, paste0(msg, "最後一波遺失超過三成，結果的代表性要特別說明。")) else ok <- c(ok, msg)
  sig_at <- if (exists("at") && !is.null(at)) at$variable[at$p < .05] else character()
  sig_tr <- if (exists("tr")) tr$transition[tr$p < .05] else character()
  if (length(sig_at) || length(sig_tr)) warn <- c(warn, sprintf("遺失和已觀察到的資料有關（%s），遺失不是完全隨機（MCAR）。FIML 在 MAR 下仍可用，報告時要說明。", paste(c(if (length(sig_at)) paste("有遺失者與完成者在", paste(sig_at, collapse = "、"), "不同"), if (length(sig_tr)) paste("前一波分數預測遺失：", paste(sig_tr, collapse = "、"))), collapse = "；")))
  else if (any(dropout)) ok <- c(ok, "有遺失者與完成者沒有顯著差異，前一波分數也不預測遺失，沒有違反 MAR 的跡象。")
  if (!is.null(mcar)) {
    s <- sprintf("Little 的 MCAR 檢定 χ²(%d) = %.2f, p %s", mcar$df, mcar$chisq, fmt_p(mcar$p))
    if (mcar$p < .05) warn <- c(warn, paste0(s, "。顯著只代表遺失不是完全隨機（MCAR）；樣本大時很小的差異也會顯著；這個檢定無法分辨 MAR 與 MNAR。",
                                              "FIML 在 MAR 下仍可用；有與遺失或結果相關的變數時用 `--aux` 納入，並可用多重插補做敏感度分析。"))
    else ok <- c(ok, paste0(s, "，沒有證據顯示遺失不是完全隨機；但檢定力有限，也無法分辨 MAR 與 MNAR。"))
  }
}
if (!is.null(rel)) {
  lowr <- rel$wave[(!is.na(rel$alpha) & rel$alpha < .70) | (!is.na(rel$omega) & rel$omega < .70)]
  spread <- suppressWarnings(max(c(diff(range(rel$alpha, na.rm = TRUE)), if (any(!is.na(rel$omega))) diff(range(rel$omega, na.rm = TRUE))), na.rm = TRUE))
  if (length(lowr)) warn <- c(warn, sprintf("%s 的 α 或 ω 低於 .70，分數的測量誤差較大；先檢查題目，或考慮二階成長模型。", paste(lowr, collapse = "、")))
  if (is.finite(spread) && spread > .10) warn <- c(warn, sprintf("各波信度最多相差 %.2f（> .10），變異的改變可能部分來自誤差變異的改變，也可能是測量恆等性的警訊。", spread))
  if (length(rel_low_items)) warn <- c(warn, sprintf("校正後題目–總分相關 < .30 的題目：%s，考慮刪修。", paste(rel_low_items, collapse = "；")))
  if (!length(lowr) && !(is.finite(spread) && spread > .10) && !length(rel_low_items))
    ok <- c(ok, sprintf("各波 α 介於 %s–%s%s，信度足夠且跨波相近。", fmt2(min(rel$alpha, na.rm = TRUE), 2), fmt2(max(rel$alpha, na.rm = TRUE), 2),
                        if (any(!is.na(rel$omega))) sprintf("、ω 介於 %s–%s", fmt2(min(rel$omega, na.rm = TRUE), 2), fmt2(max(rel$omega, na.rm = TRUE), 2)) else ""))
}
bad_sk <- w$wave[abs(w$skew) > 2 | abs(w$kurtosis) > 7]
if (length(bad_sk)) {
  warn <- c(warn, sprintf("%s 偏態或峰度過大（|偏態| > 2 或 |峰度| > 7），預設的 MLR 可部分校正，必要時考慮轉換。", paste(bad_sk, collapse = "、")))
} else ok <- c(ok, "各波分配大致對稱（|偏態| ≤ 2、|峰度| ≤ 7），預設的 MLR 足以處理。")
edge <- vapply(Y, function(v) { v <- v[!is.na(v)]; max(mean(v == min(v)), mean(v == max(v))) }, 0)
if (any(edge > .15)) warn <- c(warn, sprintf("%s 有超過 15%% 的人落在最低或最高分，可能有地板或天花板效應，變化會被低估。", paste(y[edge > .15], collapse = "、")))
dm <- diff(w$mean)
ok <- c(ok, sprintf("平均從 %.2f 變到 %.2f，%s。", w$mean[1], w$mean[Tn], if (all(dm > 0)) "逐波上升" else if (all(dm < 0)) "逐波下降" else if (max(abs(dm)) < .1 * mean(w$sd)) "變化很小" else "方向不一致，可能不是線性"))
if (w$sd[Tn] > 1.2 * w$sd[1]) ok <- c(ok, "SD 逐漸變大，個人的變化速度可能不同（斜率有個別差異）。")
ok <- c(ok, sprintf("相鄰波相關 %.2f–%.2f，%s；ICC = %.2f，表示 %.0f%% 的變異來自個人之間的差異。", min(lag1), max(lag1),
                    if (far < min(lag1) - .1) "相隔越遠相關越低（個人變化速度不同時常見）" else "相隔較遠時相關仍高，適合用成長模型", icc, 100 * icc))
if (length(hit)) {   # 遺失值代碼沒處理時，其他統計量都不可信
  ok <- character()
  warn <- c(sprintf("分數中出現 %s，可能是遺失值代碼；若是，請用 --na 指定後重跑。", paste(hit, collapse = "、")),
            "代碼被當成分數時，平均、SD、相關與偏態都不可信，其他判讀等處理後再看。")
}
rep <- c(rep, "", "## 簡易判讀", "", if (length(ok)) paste0("- ✔ ", ok), if (length(warn)) paste0("- ⚠ ", warn), "",
         if (length(warn)) "**結論**：上面有 ⚠ 的項目，建議先確認或處理，再進入步驟二。" else "**結論**：資料沒有明顯問題，可以進入步驟二。",
         "", if (is.null(rel)) "只有單一分數時無法檢驗測量恆等性（同一量尺跨時間是否相同），報告時要說明這是假設。"
             else "有題目資料：建議進一步檢驗縱貫測量恆等性（同一組題目在各波的負荷與截距是否相同），並考慮二階成長模型（因子曲線，見 references/06-extensions.md）。")

# 圖：抽樣個體軌跡與整體平均（單一系列標題命名，不需圖例）
grDevices::png(file.path(out, "trajectories.png"), 1800, 1200, res = 220, bg = "#fcfcfb")
graphics::par(mar = c(4.2, 4.2, 2.4, 1), las = 1, bty = "n", col.axis = "#555555", col.lab = "#333333")
set.seed(1); idx <- sample(nrow(Y), min(60, nrow(Y)))
graphics::plot(NA, xlim = c(1, Tn), ylim = range(Y, na.rm = TRUE), xlab = "Wave", ylab = "Score", axes = FALSE, main = "Individual trajectories (sample of 60) and mean", cex.main = .95, font.main = 1)
graphics::abline(h = pretty(range(Y, na.rm = TRUE)), col = "#e6e6e6")
graphics::axis(1, at = 1:Tn, col = "#aaaaaa", tcl = -.25); graphics::axis(2, col = NA, tcl = 0)
for (i in idx) graphics::lines(1:Tn, as.numeric(Y[i, ]), col = "#2a78d633", lwd = 1)
graphics::lines(1:Tn, w$mean, col = "#2a78d6", lwd = 3); graphics::points(1:Tn, w$mean, pch = 21, bg = "#2a78d6", col = "white", cex = 1.4)
grDevices::dev.off()
writeLines(rep, file.path(out, "screen_report.md"))
cat(paste(rep, collapse = "\n"), "\n")
