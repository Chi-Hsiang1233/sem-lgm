# export_qmd.R —— 把分析結果匯出成 Quarto 文件（.qmd）：一段程式碼、一段輸出
# Rscript export_qmd.R --results lgm_results [--name lgm_analysis] [--title "標題"] [--render]
#   --results：run_pipeline.R（或各步驟腳本）的輸出資料夾；qmd 會寫在這個資料夾裡，路徑都用相對路徑
#   --render：寫完後用 Quarto 渲染成 HTML（需要 Quarto；RStudio 內建的也可以）
# qmd 內容：
#   - 各步驟的結果摘要（*_report.md）與圖
#   - 每個模型的 lavaan 程式碼（語法、估計設定與資料都和腳本實際使用的相同），程式碼預設摺疊
#   - summary() 的輸出包在可展開的區塊裡；模型比較（LRT、配適指標表）直接顯示
source(file.path((function() { a <- commandArgs(FALSE); dirname(normalizePath(sub("^--file=", "", a[grep("^--file=", a)][1]))) })(), "lgm_helpers.R"), encoding = "UTF-8")
a <- parse_args(); need(a, "results")
root <- normalizePath(a$results, mustWork = TRUE)
name <- arg_or(a, "name", "lgm_analysis")
title <- arg_or(a, "title", "潛在成長模型分析結果")
qmd <- file.path(root, paste0(name, ".qmd"))

rl <- function(p) if (file.exists(p)) readLines(p, encoding = "UTF-8", warn = FALSE) else character()
chunk <- function(code, opts = character()) c("```{r}", if (length(opts)) paste0("#| ", opts), code, "```", "")
q <- function(x) sprintf("\"%s\"", gsub("\"", "\\\\\"", x))
cvec <- function(x) if (length(x) == 1) q(x) else sprintf("c(%s)", paste(q(x), collapse = ", "))

# ---- 各步驟資料夾 -----------------------------------------------------------
dirs <- list.dirs(root, recursive = FALSE, full.names = FALSE)
dirs <- dirs[grepl("^step[1-5]", dirs)]
# 只納入實際執行過的步驟：有 pipeline_status.csv（run_pipeline.R）時，只留這次流程中狀態為 OK 的步驟
st_file <- file.path(root, "pipeline_status.csv")
if (file.exists(st_file)) {
  st <- utils::read.csv(st_file, fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  dirs <- dirs[dirs %in% st$step[st$status == "OK"]]
}
dirs <- dirs[order(as.integer(substr(dirs, 5, 5)), nchar(dirs), dirs)]
report_of <- function(dd) { f <- list.files(file.path(root, dd), "_report\\.md$"); if (length(f)) file.path(dd, f[1]) else NA }
step_title <- function(dd) {
  r <- report_of(dd); h <- if (!is.na(r)) grep("^# ", rl(file.path(root, r)), value = TRUE)[1] else NA
  t <- if (is.na(h)) dd else sub("^# ", "", h)
  sfx <- sub("^step[0-9]+_?", "", dd); sfx <- if (sfx == "interaction") "交互作用" else sfx
  if (nzchar(sfx)) sprintf("%s（%s）", t, sfx) else t
}

form_lab <- c(none = "無成長", linear = "線性", quad = "二次", basis = "潛在基底")
case_label <- function(n) {
  if (n %in% names(form_lab)) return(paste0(form_lab[[n]], "模型（殘差變異自由估計）"))
  if (grepl("_equal_resid$", n)) { b <- sub("_equal_resid$", "", n); return(paste0(if (b %in% names(form_lab)) form_lab[[b]] else b, "模型（殘差變異相等）")) }
  fixed <- c(baseline_zero = "基準模型（預測路徑固定為 0）", conditional = "條件模型",
             M1 = "M1 全部相等", M2 = "M2 平均不同", M3 = "M3 ＋變異共變不同", M4 = "M4 ＋殘差不同")
  if (n %in% names(fixed)) return(fixed[[n]])
  paste0("時變共變數模型：", gsub("\\.", "、", n))
}
obj_name <- function(step, n, k) {
  s <- gsub("[^A-Za-z0-9]+", "_", n)
  if (grepl("[^ -~]", n) || !nzchar(gsub("_", "", s))) s <- paste0("m", k)
  paste0("f", step, "_", s)
}
read_call <- function(p) {
  kv <- regmatches(rl(p), regexpr("=", rl(p)), invert = TRUE)
  setNames(lapply(kv, `[`, 2), vapply(kv, `[`, "", 1))
}
growth_call <- function(mobj, dobj, cl) {
  args <- c(sprintf("data = %s", dobj))
  for (k in names(cl)) {
    v <- cl[[k]]; if (is.na(v) || !nzchar(v)) next
    args <- c(args, if (v %in% c("TRUE", "FALSE")) sprintf("%s = %s", k, v) else sprintf("%s = %s", k, cvec(strsplit(v, ",")[[1]])))
  }
  sprintf("growth(%s,\n       %s)", mobj, paste(args, collapse = ", "))
}
fig_cap <- c(trajectories.png = "抽樣個體軌跡與整體平均", implied_vs_observed.png = "各模型的隱含平均與樣本平均",
             residual_plot.png = "殘差散布與各波平均殘差", predicted_trajectories.png = "個人預測軌跡",
             predicted_by_covariate.png = "預測變數不同水準的預測軌跡")

# ---- 開頭 ------------------------------------------------------------------
out <- c("---", sprintf("title: %s", q(title)), "date: today", "lang: zh-TW", "format:", "  html:",
         "    toc: true", "    toc-depth: 3", "    toc-expand: true", "    number-sections: true", "    embed-resources: true",
         "    theme: cosmo", "    code-fold: true", "    code-summary: \"顯示程式碼\"", "    code-tools: true",
         "    include-in-header:", "      text: |", "        <style>",
         "        .report { border-left: 4px solid #2780e3; background: rgba(39,128,227,.06); padding: .8rem 1rem; margin: 1rem 0; border-radius: 4px; }",
         "        .report table { font-size: .88rem; }",
         "        details.output-fold > summary { cursor: pointer; color: #2780e3; margin: .3rem 0; }",
         "        </style>",
         "execute:", "  echo: true", "  warning: false", "  message: false", "  error: true", "---", "")
out <- c(out, chunk(c(
  "# summary() 的輸出包在可展開的 <details> 裡（區塊選項 fold_output: true）",
  "local({",
  "  default_output <- knitr::knit_hooks$get(\"output\")",
  "  knitr::knit_hooks$set(output = function(x, options) {",
  "    o <- default_output(x, options)",
  "    if (isTRUE(options$fold_output)) paste0(\"\\n\\n<details class=\\\"output-fold\\\"><summary>顯示 summary 輸出</summary>\\n\\n\", o, \"\\n\\n</details>\\n\\n\") else o",
  "  })",
  "})",
  "# 讀入腳本產生的摘要（*_report.md）",
  "show_report <- function(path, from = NULL) {",
  "  if (!file.exists(path)) { cat(\"\\n> 找不到 `\", path, \"`：該步驟可能執行失敗，請看對應的 .log。\\n\\n\", sep = \"\"); return(invisible()) }",
  "  x <- readLines(path, encoding = \"UTF-8\", warn = FALSE)",
  "  if (!is.null(from)) { i <- grep(from, x); x <- if (length(i)) x[i[1]:length(x)] else character() }",
  "  x <- x[!grepl(\"^# \", x)]",
  "  x <- sub(\"^## (.*)$\", \"**\\\\1**\", x)",
  "  cat(\"\\n::: {.report}\\n\\n\", paste(x, collapse = \"\\n\"), \"\\n\\n:::\\n\\n\", sep = \"\")",
  "}"), c("label: setup", "include: false")))
# ---- 分析決定紀錄（確認點）：decisions.csv 存在時放在報告開頭 ----------------------
if (file.exists(file.path(root, "decisions.csv"))) {
  dc <- read_decisions(root)
  cell <- function(x) ifelse(is.na(x) | !nzchar(x), "—", gsub("\\|", "／", x))
  why <- ifelse(is.na(dc$differs), cell(dc$reason), ifelse(dc$differs, ifelse(is.na(dc$reason) | !nzchar(dc$reason), "未說明", dc$reason), "—"))
  out <- c(out, "# 分析決定紀錄", "",
           "每個確認點先呈現結果與系統的建議，再由使用者決定。「理由」只在使用者的選擇和建議不同時填寫。", "",
           "| 確認點 | 系統的建議 | 使用者的選擇 | 理由 |", "|---|---|---|---|",
           sprintf("| %s | %s | %s | %s |", dc$checkpoint, cell(dc$recommendation), cell(dc$choice), cell(why)), "",
           if (any(grepl("^自動", dc$mode))) "這次分析設定為「不用問我，直接跑完」：確認點沒有詢問使用者，直接採用系統建議。" else "",
           "")
}
out <- c(out, "# 分析設定", "",
         "- 每個模型都用和分析腳本相同的語法、估計設定與資料重新估計，所以數值與藍色摘要一致。",
         "- 點「顯示程式碼」展開程式碼，點「顯示 summary 輸出」展開 `summary()` 的完整輸出；右上角的「程式碼」選單可以一次展開或摺疊全部程式碼。",
         "- 藍色區塊是分析腳本產生的摘要。", "",
         chunk(c("library(lavaan)",
                 "show_summary <- function(fit) summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)",
                 "# 配適指標表：有穩健版本（MLR）時用穩健值",
                 "fit_compare <- function(fits) {",
                 "  rob <- c(chisq = \"chisq.scaled\", df = \"df.scaled\", pvalue = \"pvalue.scaled\", cfi = \"cfi.robust\", tli = \"tli.robust\", rmsea = \"rmsea.robust\")",
                 "  keys <- c(\"chisq\", \"df\", \"pvalue\", \"cfi\", \"tli\", \"rmsea\", \"srmr\", \"aic\", \"bic\")",
                 "  tab <- t(sapply(fits, function(f) {",
                 "    fm <- tryCatch(fitMeasures(f), error = function(e) NULL)",
                 "    sapply(keys, function(k) { r <- rob[k]",
                 "      if (is.null(fm)) NA else if (!is.na(r) && r %in% names(fm) && !is.na(fm[[r]])) fm[[r]] else if (k %in% names(fm)) fm[[k]] else NA })",
                 "  }))",
                 "  knitr::kable(tab, digits = 3, col.names = c(\"χ²\", \"df\", \"p\", \"CFI\", \"TLI\", \"RMSEA\", \"SRMR\", \"AIC\", \"BIC\"))",
                 "}"), c("label: functions", "code-summary: \"顯示共用函式\"")))

# ---- 各步驟 ----------------------------------------------------------------
y_names <- character()
s2case <- file.path(root, "step2", rl(file.path(root, "step2", "cases.txt"))[1], "model.lav")
if (file.exists(s2case)) {
  l1 <- grep("^i =~", rl(s2case), value = TRUE)[1]
  y_names <- sub("^1\\*", "", trimws(strsplit(sub("^i =~", "", l1), "\\+")[[1]]))
}

for (dd in dirs) {
  step <- substr(dd, 5, 5)
  out <- c(out, paste("#", step_title(dd)), "")
  # 步驟一：描述統計
  if (step == "1" && length(y_names) && file.exists(file.path(root, "step2", "analysis_data.rds"))) {
    out <- c(out, "## 各波描述統計", "", chunk(c(
      "d <- readRDS(\"step2/analysis_data.rds\")",
      sprintf("y <- %s", cvec(y_names)),
      "data.frame(波次 = y, N = colSums(!is.na(d[y])), 平均 = colMeans(d[y], na.rm = TRUE),",
      "           SD = sapply(d[y], sd, na.rm = TRUE), `遺失 %` = 100 * colMeans(is.na(d[y])), check.names = FALSE) |>",
      "  knitr::kable(digits = 2, row.names = FALSE)",
      "round(cor(d[y], use = \"pairwise.complete.obs\"), 2)")))
  }
  r <- report_of(dd)
  if (!is.na(r)) out <- c(out, "## 結果摘要", "", chunk(sprintf("show_report(%s)", q(r)), c("echo: false", "results: asis")))
  pngs <- list.files(file.path(root, dd), "\\.png$")
  if (length(pngs)) {
    out <- c(out, "## 圖", "")
    for (p in pngs) out <- c(out, sprintf("![%s](%s/%s){width=80%%}", if (p %in% names(fig_cap)) fig_cap[[p]] else p, dd, p), "")
  }
  cases <- unique(rl(file.path(root, dd, "cases.txt")))
  cases <- cases[file.exists(file.path(root, dd, cases, "model.lav"))]
  if (!length(cases)) next
  dobj <- paste0("d", step, if (grepl("_", dd)) paste0("_", gsub("[^A-Za-z0-9]", "", sub("^step[0-9]+_", "", dd))))
  out <- c(out, "## lavaan 程式碼與結果", "", "讀入這一步實際使用的資料（已處理遺失值代碼、中心化、交互作用項或分組）：", "",
           chunk(sprintf("%s <- readRDS(%s)", dobj, q(file.path(dd, "analysis_data.rds")))))
  objs <- character()
  for (k in seq_along(cases)) {
    n <- cases[k]; fo <- obj_name(paste0(step, if (grepl("_", dd)) "i"), n, k); mo <- sub("^f", "m", fo)
    syn <- rl(file.path(root, dd, n, "model.lav")); cl <- read_call(file.path(root, dd, n, "call.txt"))
    out <- c(out, paste("###", case_label(n)), "",
             chunk(c(sprintf("%s <- '", mo), paste0("  ", syn), "'", sprintf("%s <- %s", fo, growth_call(mo, dobj, cl)), sprintf("show_summary(%s)", fo)),
                   "fold_output: true"))
    objs[case_label(n)] <- fo
  }
  # 模型比較
  pairs <- list()
  has <- function(x) x %in% cases
  ob <- function(x) objs[[case_label(x)]]
  if (step == "2") {
    for (p in list(c("none", "linear"), c("linear", "quad"), c("linear", "basis"))) if (all(has(p))) pairs[[length(pairs) + 1]] <- p
    eq <- grep("_equal_resid$", cases, value = TRUE)
    for (e in eq) if (has(sub("_equal_resid$", "", e))) pairs[[length(pairs) + 1]] <- c(e, sub("_equal_resid$", "", e))
  } else if (step == "3" && all(has(c("baseline_zero", "conditional")))) {
    pairs[[1]] <- c("baseline_zero", "conditional")
  } else if (step %in% c("4", "5") && length(cases) > 1) {
    for (k in seq_len(length(cases) - 1)) pairs[[k]] <- cases[k:(k + 1)]
  }
  if (length(objs) > 1 || length(pairs)) {
    code <- c("# 配適指標", sprintf("fit_compare(list(%s))", paste(sprintf("`%s` = %s", names(objs), objs), collapse = ", ")))
    if (length(pairs)) code <- c(code, "", "# 巢套比較（MLR 時自動使用尺度校正的差異檢定；看 Chisq diff 與 Pr(>Chisq)）",
                                 vapply(pairs, function(p) sprintf("lavTestLRT(%s, %s)   # %s → %s", ob(p[1]), ob(p[2]), case_label(p[1]), case_label(p[2])), ""))
    if (step == "3" && has("conditional")) code <- c(code, "", "# 解釋變異量（R²）", sprintf("round(lavInspect(%s, \"r2\"), 3)", ob("conditional")))
    out <- c(out, "### 模型比較", "", chunk(code))
  }
}

# ---- 選用：多重插補敏感度分析 ---------------------------------------------------
mi_ran <- !file.exists(st_file) || "step_mi" %in% st$step[st$status == "OK"]
if (mi_ran && file.exists(file.path(root, "sensitivity_mi", "mi_report.md"))) {
  out <- c(out, "# 預測模型：多重插補（主結果）與 FIML 對照", "")
  if (file.exists(file.path(root, "sensitivity_mi", "mi_comparison.csv"))) {
    out <- c(out, "加入的預測變數有遺失值，選擇以多重插補估計預測模型。插補設定見 `sensitivity_mi/mi_settings.txt`；下表以多重插補（Rubin 規則合併）的結果為主，並列同一模型的 FIML 結果與差異。", "",
             chunk(c("mi <- read.csv(\"sensitivity_mi/mi_comparison.csv\", fileEncoding = \"UTF-8\")",
                     "knitr::kable(mi[c(\"parameter\", \"mi_est\", \"mi_se\", \"fmi\", \"fiml_est\", \"fiml_se\", \"diff\", \"diff_in_fiml_se\")], digits = 3,",
                     "             col.names = c(\"參數\", \"多重插補估計值\", \"多重插補 SE\", \"FMI\", \"FIML 估計值\", \"FIML SE\", \"差異（插補 − FIML）\", \"差異 / FIML SE\"))")),
             chunk("show_report(\"sensitivity_mi/mi_report.md\", from = \"^\\\\*\\\\*插補下未提供\")", c("echo: false", "results: asis")))
  } else out <- c(out, chunk("show_report(\"sensitivity_mi/mi_report.md\")", c("echo: false", "results: asis")))
}

# ---- 步驟六：各模型的評估指標 -------------------------------------------------
if (file.exists(file.path(root, "report", "fit_all.csv"))) {
  out <- c(out, "# 各模型的評估指標", "",
           chunk(c("fa <- read.csv(\"report/fit_all.csv\", fileEncoding = \"UTF-8\", check.names = FALSE)", "knitr::kable(fa, digits = 3, row.names = FALSE)")),
           paste0("- ", fit_rules_text()),
           "- AIC、BIC 只能在同一步驟的模型之間比較。", "")
  if (any(grepl("^## 最終模型摘要", rl(file.path(root, "report", "analysis_report.md")))))
    out <- c(out, "# 最終模型摘要", "", chunk("show_report(\"report/analysis_report.md\", from = \"^## 最終模型摘要\")", c("echo: false", "results: asis")))
}

con <- file(qmd, open = "w", encoding = "UTF-8"); writeLines(out, con); close(con)
cat("已寫出：", qmd, "\n", sep = "")

# ---- 渲染 ------------------------------------------------------------------
if (isTRUE(a$render)) {
  cand <- c(Sys.which("quarto"), Sys.getenv("QUARTO_PATH"),
            "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe",
            "C:/Program Files/Posit/RStudio/resources/app/bin/quarto/bin/quarto.exe",
            "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto",
            "/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto")
  qb <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(qb)) {
    cat("找不到 Quarto，未渲染。安裝 Quarto 後執行：quarto render ", qmd, "\n", sep = "")
  } else {
    Sys.setenv(QUARTO_R = R.home("bin"))
    old <- setwd(root); st <- system2(qb, c("render", shQuote(basename(qmd)))); setwd(old)
    if (st == 0) cat("已渲染：", sub("\\.qmd$", ".html", qmd), "\n", sep = "") else { cat("渲染失敗（狀態 ", st, "）。\n", sep = ""); quit(status = 1) }
  }
}
