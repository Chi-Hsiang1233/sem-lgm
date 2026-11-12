# describe_vars.R —— 沒有資料說明時，依名稱與數值彙整變數的資料型態與可能角色
# 用法（R 內）：source("scripts/describe_vars.R"); describe_vars(read.csv("data.csv"))
# 同一字首、結尾是編號的欄位（如 t1–t4）合併成一列；範圍是整個系列的最小值～最大值。
describe_vars <- function(d) {
  nm <- names(d)
  stem <- ifelse(grepl("[0-9]+$", nm), sub("[0-9]+$", "", nm), nm)   # 去掉結尾的編號
  single <- as.vector(table(stem)[stem]) == 1; stem[single] <- nm[single]   # 只有一欄的維持原名
  groups <- split(nm, factor(stem, levels = unique(stem)))
  type_of <- function(x) {
    ok <- x[!is.na(x)]; k <- length(unique(ok))
    if (!is.numeric(x)) "類別（文字）" else if (k == 2) "二元" else if (k <= 7) "類別／次序"
    else if (all(ok == round(ok))) "整數" else "連續"
  }
  rows <- lapply(groups, function(v) {
    x <- unlist(d[v], use.names = FALSE); ok <- x[!is.na(x)]
    num <- all(vapply(d[v], is.numeric, TRUE))
    n <- vapply(d[v], function(z) sum(!is.na(z)), 0)
    role <- if (length(v) == 1 && num && length(unique(ok)) == nrow(d) && all(ok == round(ok))) "ID（每列唯一）"
            else if (length(v) >= 3) "重複測量（結果或時變共變數）"
            else "不隨時間改變（時不變預測或分組）"
    data.frame(變數 = if (length(v) > 1) paste0(v[1], "–", v[length(v)]) else v,
               欄數 = length(v),
               資料型態 = paste(unique(vapply(d[v], type_of, "")), collapse = "／"),
               有效人數 = if (min(n) == max(n)) format(min(n)) else sprintf("%d～%d", min(n), max(n)),
               `遺失 %` = round(100 * mean(is.na(x)), 1),
               `範圍／類別` = if (num) sprintf(if (all(ok == round(ok))) "%.0f ～ %.0f" else "%.2f ～ %.2f", min(ok), max(ok)) else paste(head(unique(ok), 4), collapse = "、"),
               可能的角色 = role, check.names = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}
