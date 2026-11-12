# make_demo_data.R —— 產生兩份「為示範而模擬」的資料（不是真實資料），固定亂數種子，可重現
# 從 skill 根目錄執行：Rscript scripts/make_demo_data.R
#   assets/data/demo_growth_missing.csv：在 Demo.growth 上製造 MAR 流失（前一波分數越低，下一波越可能流失；流失後不再回來），
#                                        並加上一個模擬的輔助變數 aux（與後期分數有關，用來示範 --aux）
#   assets/data/demo_items.csv：每波 4 題的題目資料（1–5 分），只用來示範各波信度；第 4 波刻意放一題鑑別度低的題目
#   assets/data/demo_growth_missing_x.csv：在 demo_growth_missing.csv 上再讓 x1、x2 有遺失，示範預測變數有遺失時的處理
# 原本的 demo_growth.csv 不會被修改。
src <- "assets/data/demo_growth.csv"
if (!file.exists(src)) stop("請在 skill 根目錄執行（找不到 ", src, "）", call. = FALSE)
d <- utils::read.csv(src)
y <- c("t1", "t2", "t3", "t4")

# ---- 1. demo_growth_missing.csv ---------------------------------------------
set.seed(20261112)
z <- scale(d[y])
# 輔助變數：與第 3、4 波分數有關（流失者沒有這兩波，但 aux 有資料，可補充資訊）
d$aux <- round(as.numeric(scale(.6 * z[, "t3"] + .6 * z[, "t4"] + stats::rnorm(nrow(d), 0, .8))), 3)
m <- d
gone <- rep(FALSE, nrow(d))
for (t in 2:4) {
  p <- stats::plogis(-2.75 - 1.0 * z[, t - 1])          # 前一波（已觀察到的）分數越低，越可能在這一波流失
  gone <- gone | (stats::runif(nrow(d)) < p)            # 單調流失：一旦流失就不再回來
  m[gone, y[t]] <- NA
}
utils::write.csv(m, "assets/data/demo_growth_missing.csv", row.names = FALSE, na = "")
cat(sprintf("demo_growth_missing.csv：N = %d，有遺失者 %d 人（%.1f%%）；各波遺失：%s\n", nrow(m), sum(!stats::complete.cases(m[y])),
            100 * mean(!stats::complete.cases(m[y])), paste(sprintf("%s %d", y, colSums(is.na(m[y]))), collapse = "、")))

# ---- 2. demo_items.csv ------------------------------------------------------
set.seed(20261113)
lam <- list(c(.80, .75, .75, .70), c(.80, .75, .75, .70), c(.80, .75, .70, .70), c(.75, .70, .65, .15))   # 第 4 波第 4 題鑑別度低
it <- data.frame(id = d$id)
for (w in 1:4) {
  for (j in 1:4) {
    raw <- lam[[w]][j] * z[, w] + stats::rnorm(nrow(d), 0, sqrt(1 - lam[[w]][j]^2))
    it[[sprintf("s_t%d_i%d", w, j)]] <- pmin(5, pmax(1, round(3 + 1.1 * raw)))
  }
}
for (w in 1:4) it[[sprintf("s_t%d", w)]] <- rowSums(it[sprintf("s_t%d_i%d", w, 1:4)])
it <- it[c("id", sprintf("s_t%d", 1:4), unlist(lapply(1:4, function(w) sprintf("s_t%d_i%d", w, 1:4))))]
utils::write.csv(it, "assets/data/demo_items.csv", row.names = FALSE)
cat(sprintf("demo_items.csv：N = %d，4 波 × 4 題（1–5 分），總分 s_t1–s_t4\n", nrow(it)))

# ---- 3. demo_growth_missing_x.csv：預測變數也有遺失 ----------------------------------
# 以 demo_growth_missing.csv 為基礎（結果變數已有 MAR 流失），再讓時不變預測變數有遺失：
#   x1：完全隨機遺失約 10%；x2：輔助變數 aux 越低越容易遺失（MAR），約 15%。用來示範確認點 3 之後的遺失值追問。
set.seed(20261114)
mx <- m
mx$x1[stats::runif(nrow(mx)) < .10] <- NA
mx$x2[stats::runif(nrow(mx)) < stats::plogis(-1.9 - .8 * mx$aux)] <- NA
utils::write.csv(mx, "assets/data/demo_growth_missing_x.csv", row.names = FALSE, na = "")
cat(sprintf("demo_growth_missing_x.csv：x1 遺失 %.1f%%、x2 遺失 %.1f%%（其餘同 demo_growth_missing.csv）\n", 100 * mean(is.na(mx$x1)), 100 * mean(is.na(mx$x2))))
