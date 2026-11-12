# step4_tvc.R —— 時變共變數
# Rscript step4_tvc.R --data d.csv --y "^t[0-9]$" --tvc "^c[0-9]$" [--covs x1,x2] [--form linear] [--residuals free] [--from results/step2] [--anchor first|last|mid] --out results/step4
# 依序比較三個模型：
#   (1) 效果跨波相等、時變共變數與成長因子不相關（多層次模型的預設）
#   (2) 效果各波自由
#   (3) 效果自由，且時變共變數可與成長因子相關（SEM 才能放寬；兩個變數都是內生時較合理）
source(file.path((function() { a <- commandArgs(FALSE); dirname(normalizePath(sub("^--file=", "", a[grep("^--file=", a)][1]))) })(), "lgm_core.R"), encoding = "UTF-8")
a <- parse_args(); need(a, c("data", "y", "tvc"))
require_pkgs("lavaan"); suppressPackageStartupMessages(library(lavaan))
out <- ensure_dir(arg_or(a, "out", "results/step4"))
d <- read_any(a$data, na = arg_or(a, "na"))
y <- check_cols(wave_vars(a$y, d), d); tvc <- check_cols(wave_vars(a$tvc, d), d); Tn <- length(y)
if (length(tvc) != Tn) stop("時變共變數的波數要與結果變數相同。", call. = FALSE)
prev <- read_prev(a)
form <- arg_or(a, "form", if (!is.na(prev["form"])) prev[["form"]] else "linear")
resid <- arg_or(a, "residuals", if (!is.na(prev["residuals"])) prev[["residuals"]] else "free")
times <- if (!is.null(a$times)) as.numeric(split_csv(a$times)) else if (!is.na(prev["times"])) as.numeric(split_csv(prev[["times"]])) else 0:(Tn - 1)
anchor <- resolve_anchor(prev_or(a, prev, "anchor"), times)
est <- arg_or(a, "estimator", "MLR"); covs <- split_csv(a$covs)
exo <- c(covs, tvc)
# 三個模型都用 fixed.x = FALSE，才能做巢套比較（模型 3 需要時變共變數是隨機變數）
mk <- function(mode, cor) fit_lgm(lgm_syntax(y, times, form, resid, covs = covs, tvc = tvc, tvc_mode = mode, tvc_cor = cor, anchor = anchor), d, exo = exo, estimator = est, force_random_x = TRUE)
fits <- list(`相等-不相關` = mk("equal", FALSE), `自由-不相關` = mk("free", FALSE), `自由-相關` = mk("free", TRUE))
for (n in names(fits)) { save_case(fits[[n]], out, make.names(n)); writeLines(attr(fits[[n]], "syntax"), file.path(out, paste0("model_", make.names(n), ".lav"))) }
tab <- fit_compare_table(fits, prev = setNames(head(names(fits), -1), tail(names(fits), -1))); wtab(tab, out, "comparison.csv")
p1 <- lrt_p(fits[[1]], fits[[2]]); p2 <- lrt_p(fits[[2]], fits[[3]])
rep <- c("# 步驟四：時變共變數", "", sprintf("形式 %s；殘差 %s；時變共變數 %s%s。", form, resid, paste(tvc, collapse = ", "), if (length(covs)) sprintf("；時不變預測 %s", paste(covs, collapse = ", ")) else ""),
         "", compare_md(tab), "", lrt_line(fits[[1]], fits[[2]], "效果跨波相等 → 各波自由"), lrt_line(fits[[2]], fits[[3]], "不相關 → 與成長因子相關"))
bic <- setNames(tab$bic, tab$model)
pick <- if (!is.na(p2) && p2 < .05) 3 else if (!is.na(p1) && p1 < .05) 2 else 1
rep <- c(rep, "", sprintf("**依 LRT 採用：%s**。BIC 最低的是 %s%s", names(fits)[pick], names(bic)[which.min(bic)],
                          if (names(bic)[which.min(bic)] != names(fits)[pick]) "；兩者不一致時，請把兩個結果都報告，並說明採用的理由（例如少做一個假設、效果差異的實質大小）。" else "，兩者一致。"))
for (k in seq_along(fits)) {
  gp <- growth_params(fits[[k]]); tv <- gp$reg[gp$reg$rhs %in% tvc, ]; tv <- tv[!duplicated(tv$label) | !nzchar(tv$label), ]
  if (names(fits)[k] == "相等-不相關") tv <- tv[1, ]
  rep <- c(rep, "", sprintf("**%s**", names(fits)[k]), "", sprintf("- %s：b = %.3f（95%% CI [%.3f, %.3f]）", if (nrow(tv) == 1) "時變共變數（各波相同）" else paste(tv$lhs, "←", tv$rhs), tv$est, tv$ci.lower, tv$ci.upper),
           sprintf("- 平均 %s：%.3f；變異 %s：%.3f", gp$mean$lhs, gp$mean$est, gp$var$lhs, gp$var$est))
  if (names(fits)[k] == "自由-相關") {
    cc <- gp$all[gp$all$op == "~~" & gp$all$lhs %in% GF_NAMES & gp$all$rhs %in% tvc, ]
    rep <- c(rep, sprintf("- 成長因子與時變共變數的相關：%s", paste(sprintf("%s-%s %.2f", cc$lhs, cc$rhs, cc$std.all), collapse = "、")))
  }
}
rep <- c(rep, "", "## 解讀", "", "- 加入時變共變數後，截距與斜率是「控制時變共變數之後」的成長，解讀要隨之改變。",
         "- 效果是同一波的關聯，不代表因果方向；要問時間順序，可改用前一波的值（滯後效果）。",
         "- 若時變共變數本身也在改變、且與結果互相影響，平行歷程模型（sem-longitudinal 的 parallel_growth.R）通常更合適（Grimm, 2017）。")
saveRDS(fits[[pick]], file.path(out, "tvc_fit.rds"))
writeLines(rep, file.path(out, "tvc_report.md"))
cat(paste(rep, collapse = "\n"), "\n")
