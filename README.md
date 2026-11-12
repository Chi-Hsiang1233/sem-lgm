# sem-lgm：用 SEM 分析長期追蹤資料的 Claude Skill

`sem-lgm` 是一個 [Claude Skill](https://docs.claude.com/)，讓 Claude 用結構方程模型（SEM，R 的 lavaan 套件）建立**潛在成長模型**（latent growth model, LGM），分析長期追蹤資料「怎麼變」以及「為什麼每個人變得不一樣」。

你用一般的中文描述資料和研究問題，Claude 會依固定流程執行內建的 R 腳本、在三個確認點停下來讓你決定，最後匯出一份「一段程式碼、一段輸出」的 Quarto 報告（qmd／HTML）。

分析流程依照 Grimm, Ram, & Estabrook (2017) *Growth Modeling: Structural Equation and Multilevel Modeling Approaches*。

## 功能

- **資料篩檢**：各波信度（α、ω）、遺失機制檢驗（遺失比例、涵蓋率、Little 的 MCAR 檢定、流失者比較、輔助變數）、各波描述與軌跡、跨波穩定性與 ICC。
- **無條件成長模型**：比較無成長、線性、二次（潛在基底作參考），也可以直接指定模型；殘差結構（殘差變異相等或自由估計）、假設診斷、時間編碼（截距位置）敏感度。
- **預測模型**：時不變預測變數（含交互作用、解釋變異量、遠端結果）、時變共變數、多群組 M1–M4。
- **遺失值處理**：成長模型一律用 FIML；加入的預測變數有遺失值時，可選 FIML 或多重插補（mice、Rubin 規則，並列 FIML 對照）。
- **選模門檻集中定義**：整體配適（RMSEA、CFI、TLI、SRMR）與模型比較（LRT、BIC、AIC、ΔCFI）的門檻與出處，見 `references/09-fit-criteria.md`。
- **報告**：分析決定紀錄（每個確認點的建議與你的選擇）、各步驟摘要與圖、每個模型的 lavaan 程式碼與 `summary()` 輸出（程式碼預設摺疊、可展開）。

## 需求

- [R](https://www.r-project.org/)（開發時使用 4.6.1）
- 必要套件：`lavaan`、`haven`（SPSS／Stata 檔）、`readxl`（Excel 檔）
- 選用套件：`mice`（多重插補）、`semPlot`（自動路徑圖）、`psych`（交叉核對 α）
- 匯出 HTML 報告需要 [Quarto](https://quarto.org/)（RStudio 內建的也可以）

檢查與安裝套件：

```bash
Rscript scripts/check_env.R                     # 只檢查
Rscript scripts/check_env.R --install           # 安裝缺少的必要套件
Rscript scripts/check_env.R --install-optional  # 一併安裝選用套件
```

## 安裝

**Claude Code**：把整個資料夾放到 `~/.claude/skills/sem-lgm/`，輸入 `/skills` 確認已載入。

```bash
git clone <本 repo 的網址> ~/.claude/skills/sem-lgm
```

**Claude 應用程式（claude.ai）**：把資料夾壓縮成 zip（最上層是 `sem-lgm/`），在設定的 Skills 中上傳。上傳前請先刪掉 `SKILL.md` 開頭的 `argument-hint` 這一行，那是 Claude Code 專用的欄位，claude.ai 不接受。

## 使用

在放資料的資料夾啟動 Claude Code，指令以 `/sem-lgm` 開頭：

```text
/sem-lgm 檢查 R 環境與需要的套件，然後用範例資料示範一次完整流程，每個確認點都停下來讓我選。
```

```text
/sem-lgm 分析 @data.csv
- 重複測量：dep_t1、dep_t2、dep_t3、dep_t4
- 測量時間：0、6、12、24 個月
- 截距放在：第一波
- 遺失值代碼：-99
- 可能的預測變數：age、sex、stress_t1–stress_t4、group
我想知道：憂鬱是否隨時間下降，以及年齡和性別能不能預測起點與下降速度
```

也可以不透過 Claude，直接執行腳本：

```bash
Rscript scripts/run_pipeline.R --data d.csv --y "^score_t" --times 0,1,2,4 --out lgm_results --stage screen
Rscript scripts/run_pipeline.R --data d.csv --y "^score_t" --times 0,1,2,4 --out lgm_results --stage model
Rscript scripts/run_pipeline.R --data d.csv --y "^score_t" --times 0,1,2,4 --out lgm_results --stage predictors --choose linear --covs age,sex --render
```

完整說明見使用手冊 `docs/sem-lgm-guide.html`。

## 檔案結構

```
sem-lgm/
├── SKILL.md                         # Claude 讀取的主要指示（流程、確認點、判讀原則）
├── README.md
├── references/                      # 方法與判讀文件，Claude 在對應步驟讀取
│   ├── 00-flowcharts.md             # 分析流程圖與各步驟的決策樹
│   ├── 01-workflow-assumptions.md   # 流程總覽與假設總表
│   ├── 02-preliminaries.md          # 步驟一：資料準備、信度、遺失值
│   ├── 03-unconditional.md          # 步驟二：軌跡形式、殘差結構、診斷、時間編碼
│   ├── 04-predictors.md             # 步驟三、四：預測變數
│   ├── 05-multigroup.md             # 步驟五：多群組
│   ├── 06-extensions.md             # 延伸模型與分工
│   ├── 07-reporting-troubleshooting.md  # 報告清單與問題排除
│   ├── 08-lavaan-growth-syntax.md   # lavaan 語法速查
│   └── 09-fit-criteria.md           # 選模指標與門檻（含出處）
├── scripts/                         # R 腳本
│   ├── run_pipeline.R               # 完整流程（分三段或一次跑完），最後匯出 qmd
│   ├── step1_screen.R … step6_report.R  # 各步驟分析
│   ├── sensitivity_mi.R             # 預測變數有遺失時的多重插補
│   ├── export_qmd.R                 # 匯出 qmd／HTML 報告
│   ├── lgm_helpers.R, lgm_core.R    # 共用函式（門檻定義在 FIT_RULES）
│   ├── record_decision.R            # 分析決定紀錄
│   ├── reshape_wave.R               # 長格式 ↔ 寬格式
│   ├── describe_vars.R              # 資料型態表
│   ├── check_env.R                  # 檢查／安裝 R 套件
│   ├── run_all_demos.R              # 用範例資料跑完整流程
│   └── make_demo_data.R             # 產生示範用的模擬資料
├── assets/data/                     # 範例資料（說明見 assets/data/README.md）
│   ├── demo_growth.csv              # lavaan 內建的 Demo.growth
│   ├── demo_growth_missing.csv      # 模擬：結果變數 MAR 流失＋輔助變數
│   ├── demo_growth_missing_x.csv    # 模擬：再加上預測變數的遺失
│   └── demo_items.csv               # 模擬：每波 4 題（示範信度）
└── docs/
    ├── sem-lgm-guide.qmd            # 使用手冊原始檔
    └── sem-lgm-guide.html           # 渲染好的使用手冊
```

## 輸出

每次分析的結果放在工作資料夾的 `lgm_results/`（不會寫進 skill 資料夾）：

- `stepN/`：各步驟的摘要（`*_report.md`）、CSV、圖，以及每個模型的語法、參數與 `summary()`。
- `decisions.csv`：分析決定紀錄。
- `report/analysis_report.md`、`report/fit_all.csv`：彙整報告與各模型的評估指標。
- `lgm_analysis.qmd`（加 `--render` 時另有 `.html`）：完整分析報告。

## 限制

- 假設所有人在相同的時間點測量；每個人測量時間不同時，請改用多層次成長模型。
- 只處理單一變數的成長軌跡；成長混合、平行歷程、因子曲線、潛在改變分數、RI-CLPM、分段模型等延伸模型沒有實作（見 `references/06-extensions.md`）。
- 範例資料中有三份是為示範而模擬的資料，不是真實資料。
- 腳本需要在有 R 的環境執行；沒有 R 時，Claude 只會說明要跑的指令，不會編造結果。

## 參考文獻

- Grimm, K. J., Ram, N., & Estabrook, R. (2017). *Growth modeling: Structural equation and multilevel modeling approaches*. Guilford.
- Rosseel, Y. (2012). lavaan: An R package for structural equation modeling. *Journal of Statistical Software, 48*(2), 1–36.

完整的參考文獻列在使用手冊與 `references/09-fit-criteria.md`。

## AI 協作說明

本 skill 的腳本、文件與使用手冊在作者的指導下由 Anthropic 的 Claude 協助撰寫；範例結果都是實際執行 R（lavaan）產生。使用本 skill 所做的分析、解讀與報告，由使用者自行負責。
