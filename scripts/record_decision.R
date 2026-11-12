# record_decision.R —— 把某個確認點的建議與使用者的選擇寫進分析決定紀錄（<results>/decisions.csv）
# 一般情況下 run_pipeline.R 會自動記錄；只有不必執行任何分析的選擇（例如在確認點 1 選「停止」）、
# 或要補上使用者的理由時，才需要直接呼叫：
#   Rscript record_decision.R --results lgm_results --point 1 --choice 停止 [--recommendation 繼續] [--reason "樣本太少"]
#   --point：1（資料篩檢之後）、2（選擇成長模型）、3（是否加入預測變數）
#   --differs true|false：建議與選擇是否不同（不給時：兩者文字不同就算不同；確認點 3 不比較）
source(file.path((function() { a <- commandArgs(FALSE); dirname(normalizePath(sub("^--file=", "", a[grep("^--file=", a)][1]))) })(), "lgm_helpers.R"), encoding = "UTF-8")
a <- parse_args(); need(a, c("results", "point"))
dir <- ensure_dir(a$results); p <- as.integer(a$point)
if (!p %in% 1:3) stop("--point 只能是 1、2 或 3", call. = FALSE)
old <- read_decisions(dir); old <- old[old$point == p, ]
rec <- arg_or(a, "recommendation", if (nrow(old)) old$recommendation else NULL)
ch <- arg_or(a, "choice")
dif <- if (!is.null(a$differs)) tolower(a$differs) == "true" else if (p == 3 || is.null(rec) || is.null(ch)) NULL else !identical(rec, ch)
x <- write_decision(dir, p, recommendation = arg_or(a, "recommendation"), choice = ch, differs = dif, reason = arg_or(a, "reason"),
                    mode = arg_or(a, "mode", "使用者確認"))
print(x[x$point == p, c("checkpoint", "recommendation", "choice", "differs", "reason", "mode")], row.names = FALSE)
