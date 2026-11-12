# run_all_demos.R —— 以 demo_growth.csv 跑完整個分析流程
# 從 skill 根目錄執行：Rscript scripts/run_all_demos.R --out results [--only step2,step3]
a <- commandArgs(TRUE)
get_arg <- function(k, d = NULL) { i <- match(k, a); if (is.na(i)) d else a[i + 1] }
out <- get_arg("--out", "results"); only <- get_arg("--only"); only <- if (!is.null(only)) strsplit(only, ",")[[1]]
rs <- file.path(R.home("bin"), "Rscript"); dat <- "assets/data/demo_growth.csv"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Y <- c("--data", dat, "--y", "^t[0-9]$")
cases <- list(
  step1 = list("step1_screen.R", c(Y, "--baseline", "x1,x2")),
  step2 = list("step2_unconditional.R", c(Y, "--times", "0,1,2,3")),
  step3 = list("step3_conditional.R", c(Y, "--covs", "x1,x2", "--focal", "x2", "--from", file.path(out, "step2"))),
  step3_interaction = list("step3_conditional.R", c(Y, "--covs", "x1,x2", "--interaction", "x1:x2", "--center", "--from", file.path(out, "step2"))),
  step4 = list("step4_tvc.R", c(Y, "--tvc", "^c[0-9]$", "--covs", "x1,x2", "--from", file.path(out, "step2"))),
  step5 = list("step5_multigroup.R", c(Y, "--split", "x2", "--from", file.path(out, "step2"))),
  step6 = list("step6_report.R", c("--results", out, "--final", file.path(out, "step4", "tvc_fit.rds")), noout = TRUE),
  qmd = list("export_qmd.R", c("--results", out, "--title", "潛在成長模型分析結果（Demo.growth）"), noout = TRUE),
  # 以下用兩份「為示範而模擬」的資料（scripts/make_demo_data.R 產生），輸出放在子資料夾，不影響上面 Demo.growth 的結果
  missing_step1 = list("step1_screen.R", c("--data", "assets/data/demo_growth_missing.csv", "--y", "^t[0-9]$", "--baseline", "x1,x2", "--aux", "aux",
                                           "--out", file.path(out, "missing_demo", "step1")), noout = TRUE),
  missing_step2 = list("step2_unconditional.R", c("--data", "assets/data/demo_growth_missing.csv", "--y", "^t[0-9]$", "--times", "0,1,2,3", "--choose", "linear",
                                                  "--out", file.path(out, "missing_demo", "step2")), noout = TRUE),
  missing_mi = list("sensitivity_mi.R", c("--final", file.path(out, "missing_demo", "step2", "chosen_fit.rds"), "--aux", "aux",
                                          "--out", file.path(out, "missing_demo", "sensitivity_mi")), noout = TRUE, pkg = "mice"),
  # 預測變數有遺失（demo_growth_missing_x.csv，為示範而模擬）：條件模型先用 FIML，再用多重插補（主結果）並列 FIML 對照
  missing_x_step2 = list("step2_unconditional.R", c("--data", "assets/data/demo_growth_missing_x.csv", "--y", "^t[0-9]$", "--times", "0,1,2,3", "--choose", "linear",
                                                    "--out", file.path(out, "missing_x_demo", "step2")), noout = TRUE),
  missing_x_step3 = list("step3_conditional.R", c("--data", "assets/data/demo_growth_missing_x.csv", "--y", "^t[0-9]$", "--covs", "x1,x2",
                                                  "--from", file.path(out, "missing_x_demo", "step2"), "--out", file.path(out, "missing_x_demo", "step3")), noout = TRUE),
  missing_x_mi = list("sensitivity_mi.R", c("--final", file.path(out, "missing_x_demo", "step3", "conditional_fit.rds"), "--aux", "aux",
                                            "--out", file.path(out, "missing_x_demo", "sensitivity_mi")), noout = TRUE, pkg = "mice"),
  # 確認點 1 直接指定線性：只估計線性（未做形式比較），估計值應與 step2 採用線性時相同
  direct_linear = list("step2_unconditional.R", c(Y, "--times", "0,1,2,3", "--forms", "linear", "--out", file.path(out, "direct_linear", "step2")), noout = TRUE),
  items_step1 = list("step1_screen.R", c("--data", "assets/data/demo_items.csv", "--y", "^s_t[0-9]$", "--items", "^s_t[0-9]+_i",
                                         "--out", file.path(out, "items_demo", "step1")), noout = TRUE)
)
res <- data.frame(case = character(), status = character(), seconds = numeric())
for (n in names(cases)) {
  if (!is.null(only) && !(n %in% only)) next
  cs <- cases[[n]]
  if (!is.null(cs$pkg) && !requireNamespace(cs$pkg, quietly = TRUE)) { res[nrow(res) + 1, ] <- list(n, "SKIP(pkg)", 0); next }
  ca <- cs[[2]]; if (!isTRUE(cs$noout)) ca <- c(ca, "--out", file.path(out, n))
  t0 <- Sys.time(); log <- file.path(out, paste0(n, ".log"))
  st <- system2(rs, c(shQuote(file.path("scripts", cs[[1]])), vapply(ca, shQuote, "")), stdout = log, stderr = log)
  res[nrow(res) + 1, ] <- list(n, if (st == 0) "OK" else paste0("FAIL(", st, ")"), round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
  cat(sprintf("%-20s %-10s %6.1fs\n", n, res$status[nrow(res)], res$seconds[nrow(res)]))
}
sf <- file.path(out, "demo_run_status.csv"); all_res <- res
if (!is.null(only) && file.exists(sf)) { old <- utils::read.csv(sf); all_res <- rbind(old[!(old$case %in% res$case), ], res) }
utils::write.csv(all_res, sf, row.names = FALSE)
cat(sprintf("\nOK %d / %d\n", sum(res$status == "OK"), nrow(res)))
if (any(grepl("FAIL", res$status))) quit(status = 1)
