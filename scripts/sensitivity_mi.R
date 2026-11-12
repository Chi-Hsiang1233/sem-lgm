# sensitivity_mi.R —— 預測模型的多重插補（主結果）與 FIML 對照
#   流程中的位置：確認點 3 加入的預測變數有遺失值、使用者選多重插補時，由 run_pipeline.R --missing mi 在步驟 3–5 之後執行。
#   使用者主動要求時，也可以直接對無條件模型（step2/chosen_fit.rds）執行，作為插補對照（不在流程中）。
# Rscript sensitivity_mi.R --final lgm_results/step2/chosen_fit.rds [--aux a1,a2] [--m 20] [--seed 20261112] [--maxit 10] --out lgm_results/sensitivity_mi
#   --final：已選定的最終模型（步驟腳本存的 *_fit.rds；裡面有模型語法、估計設定與實際使用的資料）
#   --aux：輔助變數，只放進插補模型，不進分析模型
#   --m：插補組數（預設 max(20, 有遺失者的百分比)）
# 做法：
#   1. 用 mice 以寬格式插補：插補模型包含分析模型的所有觀察變數（含分組變數）與 --aux；
#      連續變數用 pmm，二元變數用 logreg；固定亂數種子。
#   2. 每組資料用和主分析相同的語法與估計設定（估計法、fixed.x 等）估計。
#   3. 以 Rubin (1987) 規則合併：估計值取平均；總變異 T = Ū + (1 + 1/m)B；
#      缺失資訊比例 FMI = (r + 2/(ν + 3)) / (r + 1)，r = (1 + 1/m)B/Ū，ν = (m − 1)(1 + 1/r)²。
#      （lavaan 本身沒有多重插補的合併函式（另見 lavaan.mi 套件）；本腳本直接計算 Rubin 規則，不需要額外套件）
#   4. 以插補結果為主，並列同一模型的 FIML 估計值、SE 與差異；差異超過 FIML 的 1 個 SE 時標為注意。
#   5. 插補下做不到的模型比較（配適指標、LRT、解釋變異量、時變共變數模型比較、多群組比較）註明「插補下未提供」。
# 沒有遺失值時不執行，只寫出說明。
source(file.path((function() { a <- commandArgs(FALSE); dirname(normalizePath(sub("^--file=", "", a[grep("^--file=", a)][1]))) })(), "lgm_core.R"), encoding = "UTF-8")
a <- parse_args(); need(a, "final")
out <- ensure_dir(arg_or(a, "out", "results/sensitivity_mi"))
done <- function(msg) {
  rep <- c("# 多重插補（主結果）與 FIML 對照", "", msg)
  writeLines(rep, file.path(out, "mi_report.md")); cat(paste(rep, collapse = "\n"), "\n"); quit(status = 0)
}
fit <- readRDS(a$final)
syn <- attr(fit, "syntax"); fa <- attr(fit, "fit_args"); dat <- attr(fit, "fit_data")
if (is.null(syn) || is.null(fa) || is.null(dat)) stop("--final 的模型檔沒有記錄語法、估計設定或資料；請用目前版本的步驟腳本重新產生。", call. = FALSE)
aux <- split_csv(a$aux); check_cols(aux, dat)
grp <- fa$group
vars <- unique(c(lavaan::lavNames(fit, "ov"), grp, aux))
X <- dat[, vars, drop = FALSE]
pct <- 100 * mean(!stats::complete.cases(X))
if (!anyNA(X)) done(sprintf("分析模型的變數（%s）都沒有遺失值，多重插補與 FIML 的結果會相同，所以不執行。", paste(vars, collapse = "、")))
if (!requireNamespace("mice", quietly = TRUE)) done("需要 R 套件 mice 才能做多重插補，已略過。安裝方式：Rscript scripts/check_env.R（列出缺少的選用套件）。其他步驟不受影響。")

m <- if (!is.null(a$m)) as.integer(a$m) else max(20L, as.integer(ceiling(pct)))
seed <- as.integer(arg_or(a, "seed", 20261112)); maxit <- as.integer(arg_or(a, "maxit", 10))

# ---- 1. 插補 ----------------------------------------------------------------
binary <- vapply(X, function(v) is.numeric(v) && length(unique(stats::na.omit(v))) == 2, TRUE)
Xi <- X; for (v in names(X)[binary]) Xi[[v]] <- factor(X[[v]])
meth <- ifelse(colSums(is.na(Xi)) == 0, "", ifelse(binary, "logreg", "pmm"))
imp <- suppressWarnings(mice::mice(Xi, m = m, method = meth, seed = seed, maxit = maxit, printFlag = FALSE))
comp <- lapply(seq_len(m), function(k) {
  dk <- mice::complete(imp, k)
  for (v in names(X)[binary]) dk[[v]] <- as.numeric(as.character(dk[[v]]))
  dk
})

# ---- 2. 每組資料估計（與主分析相同的語法與設定） --------------------------------
fits <- lapply(comp, function(dk) tryCatch(quiet_fit(do.call(lavaan::growth, c(list(model = syn, data = dk), fa))), error = function(e) NULL))
okk <- vapply(fits, function(f) !is.null(f) && lavaan::lavInspect(f, "converged"), TRUE)
if (sum(okk) < 2) done(sprintf("插補後的資料有 %d 組無法收斂，無法合併。", sum(!okk)))

# ---- 3. 主要參數與 Rubin 規則 ---------------------------------------------------
lv <- lavaan::lavNames(fit, "lv")
key_of <- function(pe) paste(pe$lhs, pe$op, pe$rhs, if (is.null(pe$group)) 1 else pe$group)
main_rows <- function(pe) (pe$op == "~1" & pe$lhs %in% lv) | (pe$op == "~~" & pe$lhs %in% lv & pe$rhs %in% lv) | pe$op == "~"
pe0 <- lavaan::parameterEstimates(fit); pe0 <- pe0[main_rows(pe0) & !is.na(pe0$se) & pe0$se > 0, ]
pe0 <- pe0[order(if (is.null(pe0$group)) rep(1, nrow(pe0)) else pe0$group, match(pe0$op, c("~1", "~~", "~")), pe0$lhs != pe0$rhs), ]   # 平均 → 變異 → 共變 → 迴歸
keys <- key_of(pe0)
Q <- U <- matrix(NA_real_, sum(okk), length(keys), dimnames = list(NULL, keys))
for (j in seq_len(sum(okk))) {
  pk <- lavaan::parameterEstimates(fits[okk][[j]]); idx <- match(keys, key_of(pk))
  Q[j, ] <- pk$est[idx]; U[j, ] <- pk$se[idx]^2
}
mm <- nrow(Q); qbar <- colMeans(Q); ubar <- colMeans(U); B <- apply(Q, 2, stats::var)
Tv <- ubar + (1 + 1 / mm) * B; r <- (1 + 1 / mm) * B / ubar
nu <- ifelse(r > 0, (mm - 1) * (1 + 1 / r)^2, Inf); fmi <- ifelse(r > 0, (r + 2 / (nu + 3)) / (r + 1), 0)

lab <- function(pe) {
  nm <- c(i = "截距", s = "斜率", q = "二次項")
  f <- function(x) ifelse(x %in% names(nm), nm[x], x)
  ifelse(pe$op == "~1", paste0(f(pe$lhs), "平均"),
  ifelse(pe$op == "~~" & pe$lhs == pe$rhs, paste0(f(pe$lhs), "變異"),
  ifelse(pe$op == "~~", paste0(f(pe$lhs), "–", f(pe$rhs), "共變"), paste0(f(pe$lhs), " ← ", pe$rhs))))
}
tab <- data.frame(parameter = lab(pe0), lavaan = paste(pe0$lhs, pe0$op, pe0$rhs),
                  fiml_est = pe0$est, fiml_se = pe0$se, mi_est = qbar, mi_se = sqrt(Tv), fmi = fmi,
                  diff = qbar - pe0$est, diff_in_fiml_se = (qbar - pe0$est) / pe0$se, row.names = NULL)
if (!is.null(pe0$group)) tab$group <- pe0$group
tab$flag <- abs(tab$diff_in_fiml_se) > 1
wtab(tab, out, "mi_comparison.csv", 4)
writeLines(c(sprintf("final=%s", normalizePath(a$final)), sprintf("variables=%s", paste(vars, collapse = ",")), sprintf("aux=%s", paste(aux, collapse = ",")),
             sprintf("m=%d", m), sprintf("seed=%d", seed), sprintf("maxit=%d", maxit),
             sprintf("methods=%s", paste(sprintf("%s:%s", names(meth), ifelse(meth == "", "（無遺失）", meth)), collapse = ",")),
             sprintf("estimator=%s", fa$estimator), sprintf("fixed.x=%s", fa$fixed.x)), file.path(out, "mi_settings.txt"))

# ---- 4. 報告 ----------------------------------------------------------------
# 模型種類與插補下做不到的比較
fb <- basename(a$final)
what <- if (grepl("tvc", fb)) "時變共變數模型（步驟 4）" else if (grepl("conditional", fb)) "時不變預測變數的條件模型（步驟 3）" else "無條件模型（插補對照）"
not_avail <- paste(c("整體配適指標與模型間的 LRT、BIC 比較", if (grepl("conditional", fb)) "階層式 LRT（預測變數的整體效果）與解釋變異量（R²）", if (grepl("tvc", fb)) "三個時變共變數模型的比較", if (!is.null(grp)) "多群組 M1–M4 的比較"), collapse = "、")
rep <- c("# 多重插補（主結果）與 FIML 對照", "",
         sprintf("- 模型：%s（%s）。", basename(a$final), what),
         sprintf("- 插補：mice，%d 組（有遺失者 %.1f%%）、迭代 %d 次、亂數種子 %d；連續變數用 pmm，二元變數用 logreg。", m, pct, maxit, seed),
         sprintf("- 插補模型的變數：%s%s。", paste(setdiff(vars, aux), collapse = "、"), if (length(aux)) sprintf("；輔助變數 %s", paste(aux, collapse = "、")) else "（沒有指定輔助變數）"),
         sprintf("- 每組資料用和主分析相同的語法與設定估計（%s），%d 組都收斂%s，以 Rubin 規則合併。", fa$estimator, sum(okk), if (any(!okk)) sprintf("；%d 組未收斂，已排除", sum(!okk)) else ""),
         "", "| 參數 | 多重插補估計值 | 多重插補 SE | FMI | FIML 估計值 | FIML SE | 差異（插補 − FIML） | 差異 / FIML SE |", "|---|---|---|---|---|---|---|---|",
         sprintf("| %s%s | %.3f | %.3f | %.2f | %.3f | %.3f | %+.3f | %+.2f%s |", tab$parameter, if (!is.null(tab$group)) sprintf("（組 %s）", tab$group) else "",
                 tab$mi_est, tab$mi_se, tab$fmi, tab$fiml_est, tab$fiml_se, tab$diff, tab$diff_in_fiml_se, ifelse(tab$flag, " ⚠", "")),
         "", "FMI（缺失資訊比例）：因為遺失而損失的資訊占多少；越大代表這個參數越依賴遺失資料的處理方式。",
         "", sprintf("**插補下未提供**：%s；這些以 FIML 的結果為準。", not_avail),
         "", "## 簡易判讀", "")
if (any(tab$flag)) {
  rep <- c(rep, sprintf("- ⚠ %s 的差異超過 FIML 的 1 個 SE：結論對遺失資料的處理方式敏感，報告時兩種結果都列出並說明。", paste(tab$parameter[tab$flag], collapse = "、")))
} else rep <- c(rep, sprintf("- ✔ 所有參數的差異都在 FIML 的 1 個 SE 以內（最大 %.2f 個 SE），結論對遺失資料的處理方式不敏感。", max(abs(tab$diff_in_fiml_se))))
rep <- c(rep, sprintf("- FMI 介於 %.2f–%.2f。", min(tab$fmi), max(tab$fmi)),
         if (length(aux)) "- 多重插補多用了輔助變數；兩種方法的差異可能來自輔助變數提供的資訊，而不只是方法本身。要讓 FIML 也用到輔助變數，可用飽和相關法（auxiliary variables as saturated correlates）。",
         "- FIML 與多重插補都假設 MAR；兩者一致不代表遺失不是 MNAR。MNAR 無法由資料檢驗，只能在限制中說明，必要時做 MNAR 模型的敏感度分析。")
writeLines(rep, file.path(out, "mi_report.md"))
cat(paste(rep, collapse = "\n"), "\n")
