# check_env.R —— 檢查潛在成長模型 skill 需要的 R 套件
# Rscript check_env.R                     只檢查
# Rscript check_env.R --install           安裝缺少的必要套件
# Rscript check_env.R --install-optional  一併安裝缺少的選用套件
args <- commandArgs(TRUE)
need <- c(lavaan = "所有模型；也用於 Little 的 MCAR 檢定與 McDonald's ω", haven = "讀取 SPSS/Stata 檔", readxl = "讀取 Excel 檔")
# 選用套件：缺少時只略過對應的功能，其他步驟照常執行
opt <- c(mice = "預測變數有遺失時的多重插補（run_pipeline.R --missing mi、sensitivity_mi.R）",
         semPlot = "自動路徑圖（只供檢查模型設定；缺少時略過路徑圖）",
         psych = "非必要：可用 psych::alpha() 交叉核對步驟一的 α 與題目–總分相關")
# 不需要的套件：MCAR 檢定（naniar）、ω（semTools）、多重插補合併（lavaan.mi）都已在腳本內以 lavaan 與 Rubin 規則直接計算
has <- function(p) requireNamespace(p, quietly = TRUE)
show <- function(v) for (p in names(v)) cat(sprintf("  %-8s %-10s %s\n", p, if (has(p)) as.character(utils::packageVersion(p)) else "未安裝", v[[p]]))
cat(R.version.string, "\n\n必要套件：\n"); show(need)
cat("\n選用套件（缺少時對應功能略過，不影響其他步驟）：\n"); show(opt)
miss_opt <- names(opt)[!vapply(names(opt), has, TRUE)]
if (length(miss_opt) && !("--install-optional" %in% args)) cat("\n要安裝選用套件：Rscript check_env.R --install-optional\n")
if (any(c("--install", "--install-optional") %in% args)) {
  pk <- c(names(need), if ("--install-optional" %in% args) names(opt))
  miss <- pk[!vapply(pk, has, TRUE)]
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org") else cat("\n沒有需要安裝的套件。\n")
}
