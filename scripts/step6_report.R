# step6_report.R —— 彙整各步驟結果，產生完整分析報告、各模型的評估指標與最終模型摘要
# Rscript step6_report.R --results results [--final results/step3/conditional_fit.rds] --out results/report
# 會依序收集 step1–step5 的 *_report.md 與比較表（輸出 fit_all.csv）；--final 指定作為最終模型的 fit.rds（預設：最後一個存在的步驟）
source(file.path((function() { a <- commandArgs(FALSE); dirname(normalizePath(sub("^--file=", "", a[grep("^--file=", a)][1]))) })(), "lgm_core.R"), encoding = "UTF-8")
a <- parse_args(); need(a, "results")
require_pkgs("lavaan"); suppressPackageStartupMessages(library(lavaan))
root <- a$results; out <- ensure_dir(arg_or(a, "out", file.path(root, "report")))
steps <- c(step1 = "screen_report.md", step2 = "unconditional_report.md", step3 = "conditional_report.md", step4 = "tvc_report.md", step5 = "multigroup_report.md")
body <- character()
for (s in names(steps)) { f <- file.path(root, s, steps[[s]]); if (file.exists(f)) body <- c(body, readLines(f, encoding = "UTF-8"), "") }
cand <- file.path(root, c("step4/tvc_fit.rds", "step3/conditional_fit.rds", "step2/chosen_fit.rds"))
final <- arg_or(a, "final", cand[file.exists(cand)][1])
if (is.na(final) || !file.exists(final)) stop("找不到最終模型，請用 --final 指定。", call. = FALSE)
f <- readRDS(final); ft <- fit_table(f); gp <- growth_params(f)
# 各模型的評估指標：收集步驟二到五的比較表
rd <- function(...) { p <- file.path(root, ...); if (file.exists(p)) utils::read.csv(p, fileEncoding = "UTF-8", stringsAsFactors = FALSE) }
cols <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr", "aic", "bic")
form_lab <- c(none = "無成長", linear = "線性", quad = "二次", basis = "潛在基底")
chosen <- if (file.exists(file.path(root, "step2", "chosen.txt"))) { kv <- strsplit(readLines(file.path(root, "step2", "chosen.txt"), encoding = "UTF-8"), "="); setNames(vapply(kv, `[`, "", 2), vapply(kv, `[`, "", 1)) } else character()
parts <- list()
if (!is.null(x <- rd("step2", "comparison.csv"))) parts$s2 <- data.frame(step = "步驟二", model = ifelse(x$model %in% names(form_lab), form_lab[x$model], x$model), x[cols])
if (!is.null(x <- rd("step2", "residual_structure.csv")) && "equal" %in% x$model)
  parts$s2r <- data.frame(step = "步驟二", model = paste0(form_lab[chosen["form"]], "＋殘差變異相等"), x[x$model == "equal", cols])
if (!is.null(x <- rd("step3", "conditional", "fit_indices.csv"))) parts$s3 <- data.frame(step = "步驟三", model = "條件模型", x[cols])
if (!is.null(x <- rd("step4", "comparison.csv"))) parts$s4 <- data.frame(step = "步驟四", model = x$model, x[cols])
if (!is.null(x <- rd("step5", "comparison.csv"))) parts$s5 <- data.frame(step = "步驟五", model = x$model, x[cols])
fit_all <- do.call(rbind, parts); rownames(fit_all) <- NULL
idx_md <- character()
if (!is.null(fit_all) && nrow(fit_all)) {
  wtab(fit_all, out, "fit_all.csv")
  idx_md <- c("## 各模型的評估指標", "",
              "| 步驟 | 模型 | χ² | df | p | CFI | TLI | RMSEA | SRMR | AIC | BIC |", "|---|---|---|---|---|---|---|---|---|---|---|",
              sprintf("| %s | %s | %.2f | %s | %s | %s | %s | %s | %s | %.1f | %.1f |", fit_all$step, fit_all$model, fit_all$chisq, format(fit_all$df),
                      fmt2(fit_all$pvalue), fmt2(fit_all$cfi), fmt2(fit_all$tli), fmt2(fit_all$rmsea), fmt2(fit_all$srmr), fit_all$aic, fit_all$bic), "",
              "- χ²、df、p：模型與資料的差異；樣本大時容易顯著，只作參考。",
              paste0("- ", fit_rules_text()),
              sprintf("- AIC、BIC 越低越好，只能在同一步驟的模型之間比較（不同步驟加入的變數不同）。BIC 差距 < %d 微弱、%d–%d 有證據、%d–%d 強、> %d 非常強（Raftery, 1995）。",
                      FIT_RULES$bic_cut[1], FIT_RULES$bic_cut[1], FIT_RULES$bic_cut[2], FIT_RULES$bic_cut[2], FIT_RULES$bic_cut[3], FIT_RULES$bic_cut[3]), "")
}

gf <- c(i = "截距", s = "斜率", q = "二次項")
summ <- c("## 最終模型摘要", "",
          sprintf("- 配適：χ²(%s) = %.2f, p %s, CFI = %s, TLI = %s, RMSEA = %s [%s, %s], SRMR = %s",
                  format(ft$df), ft$chisq, fmt_p(ft$pvalue), fmt2(ft$cfi), fmt2(ft$tli), fmt2(ft$rmsea), fmt2(ft$rmsea_lo), fmt2(ft$rmsea_hi), fmt2(ft$srmr)),
          sprintf("- 平均 %s：%.3f [%.3f, %.3f]", gf[gp$mean$lhs], gp$mean$est, gp$mean$ci.lower, gp$mean$ci.upper))
reg <- gp$reg[gp$reg$lhs %in% GF_NAMES, ]
if (nrow(reg)) summ <- c(summ, sprintf("- %s ← %s：b = %.3f [%.3f, %.3f]，β = %.2f", gf[reg$lhs], reg$rhs, reg$est, reg$ci.lower, reg$ci.upper, reg$std.all))
tv <- gp$reg[!(gp$reg$lhs %in% GF_NAMES), ]
if (nrow(tv)) summ <- c(summ, sprintf("- 時變共變數：b = %.3f 至 %.3f", min(tv$est), max(tv$est)))
writeLines(c("# 潛在成長模型分析報告", "", sprintf("最終模型：%s", final), "", body, idx_md, summ), file.path(out, "analysis_report.md"))
cat(paste(c(idx_md, summ), collapse = "\n"), "\n")
