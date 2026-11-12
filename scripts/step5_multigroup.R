# step5_multigroup.R —— 多群組成長模型（M1–M4 依序比較）
# Rscript step5_multigroup.R --data d.csv --y "^t[0-9]$" --group g [--form linear] [--from results/step2] [--anchor first|last|mid] --out results/step5
#   或用 --split x1 以中位數切分連續變數（只建議作為檢查步驟三假設的工具，會損失資訊）
# 模型：M1 全部相等 → M2 平均不同 → M3 + 變異共變不同 → M4 + 殘差變異不同
source(file.path((function() { a <- commandArgs(FALSE); dirname(normalizePath(sub("^--file=", "", a[grep("^--file=", a)][1]))) })(), "lgm_core.R"), encoding = "UTF-8")
a <- parse_args(); need(a, c("data", "y"))
if (is.null(a$group) && is.null(a$split)) stop("需要 --group 或 --split。", call. = FALSE)
require_pkgs("lavaan"); suppressPackageStartupMessages(library(lavaan))
out <- ensure_dir(arg_or(a, "out", "results/step5"))
d <- read_any(a$data, na = arg_or(a, "na"))
y <- check_cols(wave_vars(a$y, d), d); Tn <- length(y)
prev <- read_prev(a)
form <- arg_or(a, "form", if (!is.na(prev["form"])) prev[["form"]] else "linear")
times <- if (!is.null(a$times)) as.numeric(split_csv(a$times)) else if (!is.na(prev["times"])) as.numeric(split_csv(prev[["times"]])) else 0:(Tn - 1)
anchor <- resolve_anchor(prev_or(a, prev, "anchor"), times)
est <- arg_or(a, "estimator", "MLR")
if (!is.null(a$split)) {
  med <- stats::median(d[[a$split]], na.rm = TRUE)
  d$.grp <- ifelse(d[[a$split]] > med, paste0(a$split, " 高"), paste0(a$split, " 低")); g <- ".grp"
  note <- sprintf("以 %s 的中位數（%.3f）切成兩組。切分連續變數會損失資訊，這裡只用來檢查步驟三的假設。", a$split, med)
} else { g <- a$group; note <- sprintf("分組變數：%s。", g) }
d <- d[!is.na(d[[g]]), ]
syn <- lgm_syntax(y, times, form, "free", anchor = anchor)
eqs <- list(M1 = c("means", "lv.variances", "lv.covariances", "residuals"), M2 = c("lv.variances", "lv.covariances", "residuals"), M3 = "residuals", M4 = character())
fits <- lapply(eqs, function(eq) fit_lgm(syn, d, estimator = est, group = g, group.equal = eq))
for (n in names(fits)) save_case(fits[[n]], out, n)
lab <- c(M1 = "全部相等", M2 = "平均不同", M3 = "＋變異共變不同", M4 = "＋殘差不同")
tab <- fit_compare_table(fits, prev = c(M2 = "M1", M3 = "M2", M4 = "M3")); wtab(tab, out, "comparison.csv")
rep <- c("# 步驟五：多群組成長模型", "", note, sprintf("各組人數：%s。", paste(sprintf("%s %d", names(table(d[[g]])), as.integer(table(d[[g]]))), collapse = "；")), "",
         sub("^\\| (M[1-4]) ", "| \\1 ", compare_md(transform(tab, model = paste(tab$model, lab[tab$model])))), "", "依序比較（Grimm, 2017）：", "",
         lrt_line(fits$M1, fits$M2, "M1 → M2：平均軌跡是否不同"), lrt_line(fits$M2, fits$M3, "M2 → M3：個別差異（變異共變）是否不同"), lrt_line(fits$M3, fits$M4, "M3 → M4：殘差變異是否不同"))
gp <- growth_params(fits$M4); gl <- lavInspect(fits$M4, "group.label")
rep <- c(rep, "", "**各組估計（M4，最少限制）**", "", sprintf("- %s 平均 %s：%.3f（SE = %.3f）", gl[gp$mean$group], gp$mean$lhs, gp$mean$est, gp$mean$se),
         sprintf("- %s 變異 %s：%.3f", gl[gp$var$group], gp$var$lhs, gp$var$est),
         "", "M2 對 M1 是兩個成長因子平均的整體檢定；與步驟三分別檢定截距、斜率的作法不同，結論可能不一致（Grimm, 2017）。")
writeLines(rep, file.path(out, "multigroup_report.md"))
cat(paste(rep, collapse = "\n"), "\n")
