# 延伸模型與工具分工

本 skill 專注於單一結果變數的潛在成長模型。下列延伸若使用者另外安裝了 `sem-longitudinal` skill，就改用它的腳本；沒有安裝時，依下表與 `08-lavaan-growth-syntax.md` 手寫 lavaan 語法，並先告知使用者這部分不是本 skill 已驗證的流程：

| 需求 | 模型 | 去哪裡 |
|---|---|---|
| 每人測量時間不同 | 多層次成長、定義變數 | `sem-longitudinal` 的 `mlm_growth.R` |
| 有不同軌跡類型 | 成長混合（GMM）、LCGA | `sem-longitudinal` 的 `gmm.R` |
| 兩個歷程一起變化 | 平行歷程（多變量成長） | `sem-longitudinal` 的 `parallel_growth.R` |
| 有飽和或 S 型趨勢 | 指數、Jenss–Bayley、邏輯成長 | 需手寫語法 |
| 有轉折點（介入、升學） | 分段（piecewise）模型 | 需手寫語法 |
| 二元或次序的重複測量 | 類別結果的成長模型 | lavaan `ordered` 並設定閾值不變 |
| 多題量表的成長 | 因子曲線（二階成長） | `sem-longitudinal` 的 `factor_curve.R` |
| 改變取決於前一波水準 | 潛在改變分數 | `sem-longitudinal` 的 `lcs.R` |
| 個體內的交叉影響 | RI-CLPM | `sem-longitudinal` 的 `panel.R` |

## 二元或次序結果的提醒

- 以潛在反應傾向（latent response propensity）建模，閾值通常假設跨時間不變。
- 比例的改變因此才能解讀為人的改變，而不是變數定義的改變。
