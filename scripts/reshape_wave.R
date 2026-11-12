# reshape_wave.R —— 長格式 ↔ 寬格式（SEM 需要寬格式：每人一列）
# 長→寬：Rscript reshape_wave.R --to wide --data long.csv --id id --time visit --vars score,stress [--keep ses,age] --output wide.csv
#         產生 score_t1、score_t2……。time 應為波次編號；時間不規則時請改用多層次成長模型（見 references/06-extensions.md）。
# 寬→長：Rscript reshape_wave.R --to long --data wide.csv --id id --stems score_t,stress_t --output long.csv
source(file.path((function() { a <- commandArgs(FALSE); dirname(normalizePath(sub("^--file=", "", a[grep("^--file=", a)][1]))) })(), "lgm_helpers.R"))
a <- parse_args(); need(a, c("to", "data", "id", "output"))
d <- read_any(a$data, na = arg_or(a, "na"))
if (a$to == "wide") {
  need(a, c("time", "vars"))
  vars <- split_csv(a$vars); keep <- split_csv(a$keep)
  dup <- anyDuplicated(d[, c(a$id, a$time)])
  if (dup) stop("同一人同一波有重複資料（第 ", dup, " 列），請先處理。", call. = FALSE)
  ids <- unique(d[[a$id]]); waves <- sort(unique(d[[a$time]]))
  out <- setNames(data.frame(ids), a$id)
  if (length(keep)) {
    first <- d[match(ids, d[[a$id]]), keep, drop = FALSE]
    varying <- keep[vapply(keep, function(k) any(tapply(d[[k]], d[[a$id]], function(v) length(unique(v[!is.na(v)]))) > 1), TRUE)]
    if (length(varying)) warning("這些 --keep 變數會隨時間變動，只保留第一筆：", paste(varying, collapse = ", "), call. = FALSE)
    out <- cbind(out, first)
  }
  for (v in vars) for (w in waves) {
    s <- d[d[[a$time]] == w, ]
    out[[paste0(v, "_t", w)]] <- s[[v]][match(ids, s[[a$id]])]
  }
} else {
  need(a, "stems")
  stems <- split_csv(a$stems)
  cols <- lapply(stems, function(s) grep(paste0("^", s, "[0-9]+$"), names(d), value = TRUE))
  waves <- sort(unique(as.integer(sub(".*?([0-9]+)$", "\\1", unlist(cols)))))
  other <- setdiff(names(d), unlist(cols))
  out <- do.call(rbind, lapply(waves, function(w) {
    r <- d[, other, drop = FALSE]; r$wave <- w
    for (s in stems) { cc <- paste0(s, w); r[[sub("_t$|_$", "", s)]] <- if (cc %in% names(d)) d[[cc]] else NA }
    r
  }))
  out <- out[order(out[[a$id]], out$wave), ]
}
utils::write.csv(out, a$output, row.names = FALSE, na = "", fileEncoding = "UTF-8")
cat(sprintf("已寫入 %s：%d 列，%d 欄\n", a$output, nrow(out), ncol(out)))
