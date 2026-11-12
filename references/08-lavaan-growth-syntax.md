# lavaan 成長模型語法速查

```r
library(lavaan)

# 線性
m <- '
  i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
  s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
'
fit <- growth(m, data = d, estimator = "MLR", missing = "fiml")
summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

# 不等距（0、1、3、6 個月）
s =~ 0*t1 + 1*t2 + 3*t3 + 6*t4

# 截距放在最後一波（--anchor last）：時間 − 最後一個時間點
s =~ -3*t1 + -2*t2 + -1*t3 + 0*t4

# 二次
q =~ 0*t1 + 1*t2 + 4*t3 + 9*t4

# 潛在基底（首尾固定，中間自由）
s =~ 0*t1 + NA*t2 + NA*t3 + 3*t4


# 殘差變異相等
t1 ~~ e*t1;  t2 ~~ e*t2;  t3 ~~ e*t3;  t4 ~~ e*t4

# 固定某個變異為 0（例如二次項沒有個別差異）
q ~~ 0*q;  q ~~ 0*i;  q ~~ 0*s

# 時不變預測變數
i + s ~ x1 + x2

# 時變共變數（效果相等）
t1 ~ b*c1;  t2 ~ b*c2;  t3 ~ b*c3;  t4 ~ b*c4

# 時變共變數與成長因子相關（需 fixed.x = FALSE，並自由估計外生變數的平均與共變）
i ~~ c1 + c2 + c3 + c4
s ~~ c1 + c2 + c3 + c4

# 遠端結果
z ~ i + s
z ~ 1

# 多群組
fit <- growth(m, data = d, group = "g", group.equal = c("means", "lv.variances", "lv.covariances", "residuals"))
```

## 容易誤會的預設

- `growth()` 把觀察變數截距固定為 0，並估計成長因子平均。
- `fixed.x = TRUE`（預設）時，外生變數有遺失的人會被刪除。
- `fixed.x = FALSE` 時，`growth()` 會把外生變數的平均固定為 0，也**不會**自動估計外生變數之間的共變；要自行加入 `x ~ 1` 與 `x1 ~~ x2` 等。本 skill 的腳本已自動處理。
- 加上標籤（如 `b*`）會讓同標籤的參數相等；多群組時同一標籤也會跨組相等。

## 取出結果

```r
parameterEstimates(fit, standardized = TRUE)
fitMeasures(fit, c("chisq.scaled", "df.scaled", "cfi.robust", "tli.robust", "rmsea.robust", "srmr"))
lavInspect(fit, "r2")
lavPredict(fit)                 # 因子分數（只用於診斷）
lavPredict(fit, type = "ov")    # 個人預測值
lavTestLRT(fit0, fit1)          # 巢套比較
```
