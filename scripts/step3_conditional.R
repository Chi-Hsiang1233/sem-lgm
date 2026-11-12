# step3_conditional.R —— 加入時不變預測變數
# Rscript step3_conditional.R --data d.csv --y "^t[0-9]$" --covs x1,x2 [--form linear] [--residuals free] [--times ...]
#     [--center] [--interaction x1:x2] [--focal x2] [--distal z] [--anchor first|last|mid] --out results/step3
#   --center：預測變數先減去平均，讓截距平均代表「預測變數在平均值時」的起始水準
#   --interaction：加入中心化後的乘積項（調節）
#   --focal：畫預測軌跡時用哪個預測變數（預設第一個）
#   形式與殘差預設讀取步驟二的 chosen.txt（--from results/step2），否則用 linear / free
source(file.path((function() { a <- commandArgs(FALSE); dirname(normalizePath(sub("^--file=", "", a[grep("^--file=", a)][1]))) })(), "lgm_core.R"), encoding = "UTF-8")
a <- parse_args(); need(a, c("data", "y", "covs"))
require_pkgs("lavaan"); suppressPackageStartupMessages(library(lavaan))
out <- ensure_dir(arg_or(a, "out", "results/step3"))
d <- read_any(a$data, na = arg_or(a, "na"))
y <- check_cols(wave_vars(a$y, d), d); Tn <- length(y)
prev <- read_prev(a)
form <- arg_or(a, "form", if (!is.na(prev["form"])) prev[["form"]] else "linear")
resid <- arg_or(a, "residuals", if (!is.na(prev["residuals"])) prev[["residuals"]] else "free")
times <- if (!is.null(a$times)) as.numeric(split_csv(a$times)) else if (!is.na(prev["times"])) as.numeric(split_csv(prev[["times"]])) else 0:(Tn - 1)
anchor <- resolve_anchor(prev_or(a, prev, "anchor"), times)
est <- arg_or(a, "estimator", "MLR")
covs <- check_cols(split_csv(a$covs), d); distal <- split_csv(a$distal)
cmeans <- vapply(covs, function(v) mean(d[[v]], na.rm = TRUE), 0)
if (isTRUE(a$center)) for (v in covs) d[[v]] <- d[[v]] - cmeans[[v]]
for (it in split_csv(a$interaction)) {
  xy <- strsplit(it, ":")[[1]]; nm <- paste0(xy[1], "_x_", xy[2])
  d[[nm]] <- (d[[xy[1]]] - mean(d[[xy[1]]], na.rm = TRUE)) * (d[[xy[2]]] - mean(d[[xy[2]]], na.rm = TRUE)); covs <- c(covs, nm)
}
rep <- c("# 步驟三：時不變預測變數", "", sprintf("形式 %s；殘差 %s；時間 %s；預測變數 %s%s。", form, resid, paste(times, collapse = ", "), paste(covs, collapse = ", "), if (isTRUE(a$center)) "（已中心化）" else ""),
         "", sprintf("預測變數平均：%s。%s", paste(sprintf("%s = %.3f", names(cmeans), cmeans), collapse = "；"),
                     if (!isTRUE(a$center) && any(abs(cmeans) > .25 * vapply(names(cmeans), function(v) sd(read_any(a$data, na = arg_or(a, "na"))[[v]], na.rm = TRUE), 0))) "預測變數平均離 0 較遠，截距平均是「預測變數 = 0」時的值；若 0 沒有意義，請加 --center。" else "截距與斜率的平均是預測變數都等於 0 時的值。"))

# 無條件模型（比較變異）、基準模型（路徑固定 0）與條件模型
f0 <- fit_lgm(lgm_syntax(y, times, form, resid, anchor = anchor), d, estimator = est)
syn_zero <- lgm_syntax(y, times, form, resid, covs = covs, cov_paths = "zero", anchor = anchor)
syn_full <- lgm_syntax(y, times, form, resid, covs = covs, distal = distal, anchor = anchor)
writeLines(syn_full, file.path(out, "model_conditional.lav"))
fz <- fit_lgm(syn_zero, d, exo = covs, estimator = est)
f1 <- fit_lgm(syn_full, d, exo = covs, estimator = est)
if (!length(distal)) save_case(fz, out, "baseline_zero")
save_case(f1, out, "conditional")
if (isTRUE(attr(f1, "random_x"))) rep <- c(rep, "", "預測變數有遺失：已改用 fixed.x = FALSE，讓 FIML 保留這些個案。")
rep <- c(rep, "", "## 配適", "", fit_sentence(fit_table(f1), est), "",
         if (length(distal)) "（含遠端結果，與基準模型的變數不同，不做以下 LRT）" else lrt_line(fz, f1, "加入預測變數的整體效果（基準模型：路徑固定為 0）"),
         "", "加入時不變預測變數通常不會明顯改變 CFI、RMSEA 等整體指標；預測變數是否有用，要看上面的 LRT 與個別路徑（Grimm, 2017）。")

gp <- growth_params(f1); gp0 <- growth_params(f0)
reg <- gp$reg[gp$reg$lhs %in% GF_NAMES, ]
rep <- c(rep, "", "## 預測變數的效果", "", "| 路徑 | b | SE | 95% CI | p | β |", "|---|---|---|---|---|---|",
         sprintf("| %s ← %s | %.3f | %.3f | [%.3f, %.3f] | %s | %.3f |", reg$lhs, reg$rhs, reg$est, reg$se, reg$ci.lower, reg$ci.upper, sub("= ", "", fmt_p(reg$pvalue)), reg$std.all))
if (length(distal)) {
  dr <- gp$reg[gp$reg$lhs %in% distal, ]
  rep <- c(rep, "", "**遠端結果**", "", sprintf("- %s ← %s：b = %.3f（95%% CI [%.3f, %.3f]），β = %.3f", dr$lhs, dr$rhs, dr$est, dr$ci.lower, dr$ci.upper, dr$std.all))
}
# 解釋變異量：無條件變異 vs 條件殘差變異
gf <- intersect(gp$var$lhs, gp0$var$lhs)
v0 <- setNames(gp0$var$est, gp0$var$lhs)[gf]; v1 <- setNames(gp$var$est, gp$var$lhs)[gf]
r2 <- tryCatch(lavInspect(f1, "r2")[gf], error = function(e) rep(NA, length(gf)))
pr2 <- (v0 - v1) / v0
wtab(data.frame(factor = gf, var_unconditional = v0, resid_var_conditional = v1, pseudo_R2 = pr2, model_R2 = r2), out, "explained_variance.csv")
rep <- c(rep, "", "## 解釋變異量", "", "| 成長因子 | 無條件變異 | 條件殘差變異 | 解釋比例（變異縮減） | 模型 R² |", "|---|---|---|---|---|",
         sprintf("| %s | %.3f | %.3f | %.1f%% | %.3f |", gf, v0, v1, 100 * pr2, r2))
gc <- gp$cor; gc0 <- gp0$cor
if (nrow(gc) && nrow(gc0)) rep <- c(rep, "", sprintf("成長因子相關：無條件 %s；加入預測變數後的殘差相關 %s。", paste(sprintf("%s-%s r = %.2f", gc0$lhs, gc0$rhs, gc0$std.all), collapse = "；"), paste(sprintf("%s-%s r = %.2f", gc$lhs, gc$rhs, gc$std.all), collapse = "；")))

# 預測軌跡：焦點預測變數在平均 ±1 SD（其他預測變數固定在平均）
focal <- arg_or(a, "focal", covs[1])
if ("s" %in% gf && focal %in% covs) {
  mu <- vapply(covs, function(v) mean(d[[v]], na.rm = TRUE), 0); sdv <- sd(d[[focal]], na.rm = TRUE)
  cf <- function(g, v) { r <- reg[reg$lhs == g & reg$rhs == v, "est"]; if (length(r)) r else 0 }
  m <- setNames(gp$mean$est, gp$mean$lhs)
  L <- c(list(i = rep(1, Tn)), time_loadings(times, form, anchor))
  if (form == "basis") L$s <- gp$load[gp$load$lhs == "s", "est"]
  pred <- function(val) { xv <- mu; xv[focal] <- val; Reduce(`+`, lapply(names(L), function(g) (m[[g]] + sum(vapply(covs, function(v) cf(g, v) * xv[[v]], 0))) * L[[g]])) }
  series <- list(pred(mu[[focal]] + sdv), pred(mu[[focal]]), pred(mu[[focal]] - sdv))
  names(series) <- sprintf("%s %s", focal, c("+1 SD", "mean", "-1 SD"))
  line_chart(file.path(out, "predicted_by_covariate.png"), times, series, "Time", "Predicted score", sprintf("Predicted trajectories by %s (others at mean)", focal))
  rep <- c(rep, "", sprintf("預測軌跡圖：%s 在平均 ±1 SD，其他預測變數固定在平均（predicted_by_covariate.png）。", focal))
}
rep <- c(rep, "", "## 本步驟依賴的假設", "", "- 所有人遵循同一種軌跡形式，不論預測變數的值。", "- 成長因子的殘差變異與共變，以及觀察分數的殘差變異，不隨預測變數改變。", "- 這些假設可用多群組模型（步驟五）檢驗類別型預測變數（Grimm, 2017）。")
saveRDS(f1, file.path(out, "conditional_fit.rds"))
writeLines(rep, file.path(out, "conditional_report.md"))
cat(paste(rep, collapse = "\n"), "\n")
