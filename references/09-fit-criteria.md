# 選模指標與門檻

和使用手冊「估計與模型比較」一節的兩張門檻表相同。腳本中的數值集中定義在 `scripts/lgm_helpers.R` 的 `FIT_RULES`，所有步驟都引用它（`fit_verdict()`、`fit_poor_flags()`、`compare_verdict()`、`bic_evidence()`）；改門檻時，`FIT_RULES`、本檔與手冊三處一起改。

## 整體配適門檻（單一模型是否可接受）

| 指標 | 門檻 | 出處 |
|---|---|---|
| χ² | 不顯著代表完美配適；樣本大時通常不切實際，只作參考 | Grimm et al. (2017) |
| RMSEA | < .05 良好；.05–.08 可接受；.08–.10 勉強；≥ .10 不佳 | Grimm et al. (2017)；Browne & Cudeck (1993)；MacCallum et al. (1996) |
| CFI、TLI | > .95 良好；.90–.95 可接受；< .90 不佳 | Grimm et al. (2017)；Hu & Bentler (1999) |
| SRMR | ≤ .08 | Hu & Bentler (1999) |

自由度 ≤ 3 時，RMSEA 不納入總判定，只供參考（Kenny et al., 2015）。

腳本的「整體配適」總判定（`fit_verdict()`）：

- **良好**：所有納入的指標都良好（CFI、TLI > .95；RMSEA < .05；SRMR ≤ .08）。
- **不佳**：任一指標不佳（CFI 或 TLI < .90，或 RMSEA ≥ .10）。
- **可接受**：其餘情況，括號中列出未達良好的指標，例如「可接受（未達良好：RMSEA）」；RMSEA 勉強（.08–.10）也屬於這一類。
- 自由度 ≤ 3 時加註「自由度小，RMSEA 僅供參考」；確認點 1 的配適提醒也不單獨因 RMSEA 觸發。
- df = 0 時標為「恰好識別」，不判定。

## 模型比較門檻

| 指標 | 門檻 | 適用 | 出處 |
|---|---|---|---|
| LRT（Δχ²） | p < .05 表示較複雜的模型顯著較好 | 巢套模型 | Grimm et al. (2017) |
| BIC 差距 | < 2 微弱；2–6 有證據；6–10 強；> 10 非常強 | 巢套與非巢套 | Raftery (1995) |
| AIC 差距 | ≤ 2 無法區分；> 10 幾乎不支持較差的模型 | 巢套與非巢套 | Burnham & Anderson (2002) |
| CFI 差距 | ≤ .01 視為沒有實質差異 | 輔助判斷 | Cheung & Rensvold (2002) |

比較表的「比較結論」欄（`compare_verdict()`）寫出：和前一個較簡單的模型比的 LRT 結論、BIC 支持哪一個模型與證據強度、ΔCFI 是否 ≤ .01。

## 使用時要注意

- 這些是經驗法則，不是硬性標準，要和參數是否合理（負變異數、相關超過 ±1）、理論能否解釋一起判斷。
- CFI 與 TLI 取決於基準模型；一般的獨立模型不適合追蹤資料，應使用合適的基準模型（Widaman & Thompson, 2003；Grimm et al., 2017）。步驟二的「B. 以潛在基底模型為基準的比較」就是這個做法（見 `03-unconditional.md` §2）。
- 檢定變異 = 0（例如無成長 vs 線性，或固定斜率變異為 0）時有邊界值問題：變異不能小於 0，LRT 的分配是卡方的混合分配，一般 p 值偏保守；常用近似是把 p 除以 2（Stoel et al., 2006）。

## 在流程中的用法

| 地方 | 用到的門檻 |
|---|---|
| 確認點 2 的建議（`step2_unconditional.R`） | 無成長 → 線性：LRT p < .05（或無成長不可採用）。線性 → 二次：二次可採用、LRT p < .05，且二次的 BIC 低至少 2（Raftery 的「有證據」）。建議理由寫出用到的指標與門檻 |
| 確認點 1 直接指定單一模型的提醒（`form_check.txt`） | CFI < .90 或 RMSEA ≥ .10（自由度 ≤ 3 時不因 RMSEA 提醒）；不收斂或有負變異數；各波殘差平均最大絕對值 ≥ 各波 SD 平均的 .10 倍（系統性偏離） |
| 比較表的判定欄（步驟二、四、五） | 「整體配適」依整體配適門檻；「比較結論」依模型比較門檻 |
| 報告的參考值說明（`step6_report.R`、`export_qmd.R`） | `fit_rules_text()` |

## 文獻

- Browne, M. W., & Cudeck, R. (1993). Alternative ways of assessing model fit. In K. A. Bollen & J. S. Long (Eds.), *Testing structural equation models* (pp. 136–162). Sage.
- Burnham, K. P., & Anderson, D. R. (2002). *Model selection and multimodel inference* (2nd ed.). Springer.
- Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit indexes for testing measurement invariance. *Structural Equation Modeling, 9*(2), 233–255.
- Grimm, K. J., Ram, N., & Estabrook, R. (2017). *Growth modeling: Structural equation and multilevel modeling approaches*. Guilford.
- Hu, L., & Bentler, P. M. (1999). Cutoff criteria for fit indexes in covariance structure analysis. *Structural Equation Modeling, 6*(1), 1–55.
- Kenny, D. A., Kaniskan, B., & McCoach, D. B. (2015). The performance of RMSEA in models with small degrees of freedom. *Sociological Methods & Research, 44*(3), 486–507.
- MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis and determination of sample size for covariance structure modeling. *Psychological Methods, 1*(2), 130–149.
- Raftery, A. E. (1995). Bayesian model selection in social research. *Sociological Methodology, 25*, 111–163.
- Widaman, K. F., & Thompson, J. S. (2003). On specifying the null model for incremental fit indices in structural equation modeling. *Psychological Methods, 8*(1), 16–37.
