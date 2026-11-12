# run_pipeline.R —— 用自己的資料跑完整流程（步驟 1→2→3→4→5→6，依提供的參數自動略過不需要的步驟）
# 從 skill 根目錄或任何地方執行：
#   Rscript <skill>/scripts/run_pipeline.R --data d.csv --y "^score_t" --out lgm_results
#     [--times 0,1,2,4] [--anchor first|last|mid] [--na -99] [--estimator MLR]
#     [--baseline age,sex]            步驟一遺失值分析要比較的基線變數（預設同 --covs）
#     [--items "^dep_t[0-9]+_i"]      步驟一：各波的題目欄位，計算各波信度（α、ω、題目–總分相關）
#     [--aux a1,a2]                   輔助變數：步驟一的遺失值分析與多重插補模型會納入
#     [--missing fiml|mi] [--m 20] [--seed 1]   預測模型的遺失值處理（確認點 3 之後的追問；--mi 等同 --missing mi）：
#                                     步驟 2 一律用 FIML。加入的預測變數（含時變共變數、分組變數、遠端結果）有遺失值、又沒有給 --missing 時，
#                                     列出遺失比例後停下，請使用者選 FIML（預設）或多重插補；選 mi 時在步驟 3–5 之後用 sensitivity_mi.R
#                                     對最終的預測模型做多重插補（主結果），並列 FIML 對照（需要 mice）。沒有加入預測變數時不插補。
#     [--reason4 "理由"]              使用者選多重插補（和預設 FIML 不同）時的理由
#     [--covs age,sex] [--center] [--interaction age:sex] [--focal age] [--distal z]   步驟三（--distal 要和 --covs 一起用）
#     [--tvc "^stress_t"]             步驟四
#     [--group g | --split x]         步驟五
#     [--choose none|linear|quad] [--residuals free|equal]   確認點 2 使用者選的成長模型（basis 只在使用者主動要求時用）
#     [--render] [--title "標題"] [--name lgm_analysis] [--no-qmd]
#                                     最後匯出 <out>/lgm_analysis.qmd（一段程式碼、一段輸出）；--render 渲染成 HTML
#
# 分段執行（三個確認點；每段跑完先讓使用者確認，再跑下一段，參數保持相同）：
#   --stage screen       步驟一（資料篩檢）               → 確認點 1：a 繼續，比較所有模型（預設）／b 直接指定成長模型／
#                                                          c 先修正設定後重跑篩檢／d 停止
#   --stage model        步驟二（估計無成長、線性、二次，並列出潛在基底作參考）→ 確認點 2：選成長模型
#     [--forms linear | --forms linear,quad]   確認點 1 選 b：只估計指定的形式。只指定一個時略過確認點 2（不做形式比較），
#                                              有提醒（step2/form_check.txt）時先問使用者是否改為比較所有模型；--stage predictors 也要帶相同的 --forms
#   --stage predictors   依 --choose 採用的模型（和步驟二不同時自動重跑步驟二），再跑使用者選的預測模型（步驟三、四、五；
#                        預測變數有遺失值時先問 --missing）、選用的多重插補、步驟六，最後匯出 qmd → 確認點 3 的選擇由 --covs、--tvc、--group/--split、--distal 決定
#   --reason1/--reason2/--reason3 "理由"   使用者的選擇和建議不同時的理由（寫進分析決定紀錄）
#   --rec3 "建議"        確認點 3 的系統建議（不給時寫「依研究問題決定」）
# 一次跑完（不分段）：不給 --stage（也不給 --steps），或加 --auto（「不用問我，直接跑完」）；確認點都採用系統建議，紀錄中註明為自動。
# 進階：--steps 1,2,3,4,5,6 只跑指定步驟。
# 每一步的輸出在 <out>/stepN/，執行紀錄在 <out>/stepN.log，狀態表在 <out>/pipeline_status.csv，
# 分析決定紀錄在 <out>/decisions.csv（export_qmd.R 會放在報告開頭）
a <- commandArgs(TRUE)
has <- function(k) k %in% a
get_arg <- function(k, d = NULL) { i <- match(k, a); if (is.na(i) || i == length(a) || startsWith(a[i + 1], "--")) d else a[i + 1] }
if (is.null(get_arg("--data")) || is.null(get_arg("--y"))) stop("需要 --data 與 --y。用法見本檔開頭的註解。", call. = FALSE)

sdir <- (function() { x <- commandArgs(FALSE); dirname(normalizePath(sub("^--file=", "", x[grep("^--file=", x)][1]))) })()
source(file.path(sdir, "lgm_helpers.R"), encoding = "UTF-8")
rs <- file.path(R.home("bin"), "Rscript")
out <- get_arg("--out", "lgm_results"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
stage <- get_arg("--stage")
if (!is.null(stage) && !stage %in% c("screen", "model", "predictors")) stop("--stage 只能是 screen、model 或 predictors", call. = FALSE)
steps <- if (!is.null(stage)) switch(stage, screen = 1L, model = 2L, predictors = 3:6) else as.integer(strsplit(get_arg("--steps", "1,2,3,4,5,6"), ",")[[1]])
auto <- has("--auto") || (is.null(stage) && is.null(get_arg("--steps")))   # 一次跑完：確認點採用系統建議
mode <- if (auto) "自動（未詢問，採用系統建議）" else "使用者確認"

pass <- function(keys) unlist(lapply(keys, function(k) { v <- get_arg(k); if (!is.null(v)) c(k, v) }))
flag <- function(keys) keys[vapply(keys, has, TRUE)]
base <- c("--data", get_arg("--data"), "--y", get_arg("--y"), pass(c("--na", "--estimator")))
from <- c("--from", file.path(out, "step2"))
covs <- get_arg("--covs"); tvc <- get_arg("--tvc"); choose <- get_arg("--choose")
grp <- c(pass("--group"), pass("--split"))
fl <- c(none = "無成長", linear = "線性", quad = "二次", basis = "潛在基底")
kv <- function(f) { if (!file.exists(f)) return(character()); x <- strsplit(readLines(f, encoding = "UTF-8"), "=", fixed = TRUE); setNames(vapply(x, function(z) paste(z[-1], collapse = "="), ""), vapply(x, `[`, "", 1)) }

# 確認點 3 之前：使用者選的成長模型和步驟二目前採用的不同時，先重跑步驟二
if (!(2 %in% steps) && any(steps >= 3) && !is.null(choose)) {
  cur <- kv(file.path(out, "step2", "chosen.txt"))["form"]
  if (is.na(cur) || cur != choose || (!is.null(get_arg("--residuals")) && kv(file.path(out, "step2", "chosen.txt"))["residuals"] != get_arg("--residuals")))
    steps <- sort(unique(c(2L, steps)))
}

plan <- list(
  `1` = list("step1_screen.R", c(base, pass(c("--times", "--baseline", "--items", "--aux")), if (is.null(get_arg("--baseline")) && !is.null(covs)) c("--baseline", covs)), TRUE),
  `2` = list("step2_unconditional.R", c(base, pass(c("--times", "--anchor", "--forms", "--choose", "--residuals"))), TRUE),
  `3` = list("step3_conditional.R", c(base, from, pass(c("--covs", "--interaction", "--focal", "--distal")), flag("--center")), !is.null(covs)),
  `4` = list("step4_tvc.R", c(base, from, pass(c("--tvc", "--covs"))), !is.null(tvc)),
  `5` = list("step5_multigroup.R", c(base, from, grp), length(grp) > 0)
)
res <- data.frame(step = character(), status = character(), seconds = numeric(), stringsAsFactors = FALSE)
run <- function(n, script, args) {
  log <- file.path(out, paste0("step", n, ".log")); t0 <- Sys.time()
  st <- system2(rs, c(shQuote(file.path(sdir, script)), vapply(args, shQuote, "")), stdout = log, stderr = log)
  res[nrow(res) + 1, ] <<- list(paste0("step", n), if (st == 0) "OK" else paste0("FAIL(", st, ")，見 ", log), round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
  cat(sprintf("step%-3s %-8s %6.1fs  %s\n", n, if (st == 0) "OK" else "FAIL", res$seconds[nrow(res)], script))
  st == 0
}

# ---- 分析決定紀錄 -------------------------------------------------------------
dec <- read_decisions(out)
prev_choice <- function(p) { r <- dec$choice[dec$point == p]; if (length(r) && !is.na(r)) r else NA_character_ }
REC1 <- "繼續，比較所有模型"
forms_arg <- get_arg("--forms"); forms_sel <- if (!is.null(forms_arg)) split_csv(forms_arg) else NULL
single_form <- length(forms_sel) == 1
log1_after_screen <- function() {   # 確認點 1 的建議：依簡易判讀
  rp <- file.path(out, "step1", "screen_report.md"); if (!file.exists(rp)) return(invisible())
  w <- sum(grepl("^- ⚠", readLines(rp, encoding = "UTF-8")))
  rec <- if (w) sprintf("先確認 %d 個 ⚠ 項目再%s", w, REC1) else REC1
  ch <- if (auto) REC1 else if (!is.na(prev_choice(1))) "先修正設定後重跑篩檢（待確認）" else "待確認"
  write_decision(out, 1, recommendation = rec, choice = ch, differs = if (auto) w > 0 else NA, reason = if (auto && w) "未詢問（自動執行）" else NA_character_, mode = mode)
}
continue_point1 <- function() {      # 進入步驟二 = 使用者在確認點 1 選擇繼續（a）或直接指定成長模型（b）
  if (!nrow(dec) || !(1 %in% dec$point) || auto) return(invisible())
  pc <- prev_choice(1); rec <- dec$recommendation[dec$point == 1]
  go <- if (is.null(forms_sel)) REC1 else if (single_form) sprintf("直接指定成長模型：%s（未做形式比較）", fl[[forms_sel]])
        else sprintf("直接指定成長模型：只比較%s", paste(fl[forms_sel], collapse = "、"))
  ch <- if (!is.na(pc) && startsWith(pc, "先修正設定後重跑篩檢")) paste("先修正設定後重跑篩檢 →", go) else go
  dif <- rec != REC1   # 建議是否繼續；直接指定模型仍是「繼續」，選擇欄會寫出指定的模型
  write_decision(out, 1, choice = ch, differs = dif, reason = get_arg("--reason1", if (dif) "未說明" else NA_character_), mode = mode)
}
form_flags <- function() {           # 直接指定單一模型時，步驟二寫出的提醒
  f <- file.path(out, "step2", "form_check.txt"); if (!file.exists(f)) return(character())
  x <- sub("^flag=", "", readLines(f, encoding = "UTF-8")); x[nzchar(x)]
}
log2_after_model <- function() {    # 確認點 2 的建議：步驟二的 recommendation.txt
  r <- kv(file.path(out, "step2", "recommendation.txt")); if (!length(r)) return(invisible())
  rec <- sprintf("%s（%s）", fl[[r[["form"]]]], r[["reason"]])
  used <- kv(file.path(out, "step2", "chosen.txt"))
  if (identical(r["compared"], c(compared = "FALSE"))) {   # 確認點 1 直接指定單一模型：略過確認點 2
    fg <- form_flags()
    ch <- sprintf("%s（殘差%s）", fl[[used[["form"]]]], if (used[["residuals"]] == "equal") "變異相等" else "變異自由估計")
    pending <- identical(stage, "model") && length(fg)
    write_decision(out, 2, recommendation = "—（使用者於確認點 1 直接指定，未做形式比較）",
                   choice = if (pending) paste0(ch, "（待確認：已提醒配適問題）") else ch, differs = FALSE,
                   reason = if (length(fg)) paste0("已提醒：", paste(fg, collapse = "；"), if (pending) "" else "；使用者仍採用") else NA_character_, mode = mode)
    return(invisible())
  }
  if (identical(stage, "model") && is.null(choose)) {
    write_decision(out, 2, recommendation = rec, choice = "待確認", differs = NA, reason = NA_character_, mode = mode)
  } else {
    ch <- sprintf("%s（殘差%s）", fl[[used[["form"]]]], if (used[["residuals"]] == "equal") "變異相等" else "變異自由估計")
    dif <- used[["form"]] != r[["form"]]
    write_decision(out, 2, recommendation = rec, choice = ch, differs = dif, reason = if (dif) get_arg("--reason2", "未說明") else NA_character_, mode = mode)
  }
}
log3_predictors <- function() {
  nice <- function(x) { v <- tryCatch(wave_vars(x, utils::read.csv(get_arg("--data"), nrows = 1, check.names = FALSE)), error = function(e) split_csv(x))
                        if (length(v) > 2 && startsWith(x, "^")) paste0(v[1], "–", v[length(v)]) else paste(v, collapse = "、") }
  sel <- c(if (!is.null(covs)) sprintf("時不變預測變數（%s%s%s）", nice(covs), if (has("--center")) "，中心化" else "",
                                      if (!is.null(get_arg("--interaction"))) paste0("，交互作用 ", get_arg("--interaction")) else ""),
           if (!is.null(tvc)) sprintf("時變共變數（%s）", nice(tvc)),
           if (length(grp)) sprintf("多群組（%s）", paste(grp[c(FALSE, TRUE)], collapse = "、")),
           if (!is.null(get_arg("--distal")) && !is.null(covs)) sprintf("遠端結果（%s）", get_arg("--distal")))
  write_decision(out, 3, recommendation = get_arg("--rec3", if (auto) "依提供的參數執行" else "依研究問題決定"),
                 choice = if (length(sel)) paste(sel, collapse = "；") else "不加入", differs = NA,
                 reason = get_arg("--reason3", NA_character_), mode = mode)
}

# ---- 預測變數的遺失值（確認點 3 之後的追問）；分析決定紀錄第 4 列「遺失值的處理」 -----------------
miss_arg <- get_arg("--missing", if (has("--mi")) "mi" else NULL)
if (!is.null(miss_arg) && !miss_arg %in% c("fiml", "mi")) stop("--missing 只能是 fiml 或 mi", call. = FALSE)
pred_on <- !is.null(covs) || !is.null(tvc) || length(grp) > 0
use_mi <- FALSE
if (any(steps >= 3)) {
  mp <- numeric()
  if (pred_on) {
    d0 <- read_any(get_arg("--data"), na = get_arg("--na"))
    pv <- unique(c(split_csv(covs), if (!is.null(tvc)) wave_vars(tvc, d0), grp[c(FALSE, TRUE)], split_csv(get_arg("--distal"))))
    pv <- intersect(pv, names(d0)); mp <- 100 * colMeans(is.na(d0[pv])); mp <- mp[mp > 0]
  }
  mtxt <- if (length(mp)) paste(sprintf("%s 遺失 %.1f%%", names(mp), mp), collapse = "、") else ""
  if (!pred_on) {
    if (identical(miss_arg, "mi")) cat("沒有加入預測變數，流程中不插補；需要時直接執行 sensitivity_mi.R 對無條件模型做插補對照。\n")
    write_decision(out, 4, recommendation = "FIML（步驟 2 一律用 FIML）", choice = "FIML（沒有加入預測變數，不插補）", differs = NA, reason = NA_character_, mode = mode)
  } else if (!length(mp)) {
    write_decision(out, 4, recommendation = "FIML（預設）", choice = "FIML（加入的預測變數沒有遺失值，不詢問）", differs = NA, reason = NA_character_, mode = mode)
  } else if (auto) {
    write_decision(out, 4, recommendation = "FIML（預設）", choice = "FIML", differs = FALSE,
                   reason = sprintf("預測變數有遺失（%s）；未詢問（自動執行），採用 FIML", mtxt), mode = mode)
  } else if (is.null(miss_arg)) {
    write_decision(out, 4, recommendation = "FIML（預設）", choice = "待確認", differs = NA, reason = paste0("有遺失的預測變數：", mtxt), mode = mode)
    cat(sprintf("\n加入的預測變數有遺失值：%s。\n請先問使用者「遺失值怎麼處理？」：FIML（預設）｜多重插補；再用相同的參數加 --missing fiml 或 --missing mi 重跑 --stage predictors（選多重插補時可加 --aux）。\n", mtxt))
    quit(status = 0)
  } else {
    use_mi <- miss_arg == "mi"
    write_decision(out, 4, recommendation = "FIML（預設）", choice = if (use_mi) "多重插補（主結果，並列 FIML 對照）" else "FIML",
                   differs = use_mi, reason = sprintf("有遺失的預測變數：%s%s", mtxt, if (use_mi) paste0("；理由：", get_arg("--reason4", "未說明")) else ""), mode = mode)
  }
  dec <- read_decisions(out)
}

# ---- 執行 --------------------------------------------------------------------
for (n in names(plan)) {
  if (!(as.integer(n) %in% steps)) next
  p <- plan[[n]]
  if (n == "2") continue_point1()
  if (!isTRUE(p[[3]])) { res[nrow(res) + 1, ] <- list(paste0("step", n), "略過（未提供參數）", 0); next }
  ok <- run(n, p[[1]], c(p[[2]], "--out", file.path(out, paste0("step", n))))
  if (n == "1" && ok) { log1_after_screen(); dec <- read_decisions(out) }
  if (n == "2" && ok) { log2_after_model(); dec <- read_decisions(out) }
  if (!ok && n == "2") { cat("步驟二失敗，後續步驟需要它的結果，停止。\n"); break }
}
if (any(steps >= 3)) {
  if (!(2 %in% steps) && file.exists(file.path(out, "step2", "recommendation.txt"))) log2_after_model()
  log3_predictors()
}
# 預測模型的多重插補（確認點 3 之後選 --missing mi 才執行；用最終的預測模型，主結果是插補，並列 FIML 對照）
if (any(steps >= 3)) {
  mi_cand <- file.path(out, c(if (!is.null(tvc)) "step4/tvc_fit.rds", if (!is.null(covs)) "step3/conditional_fit.rds"))
  mi_cand <- mi_cand[file.exists(mi_cand)]
  if (use_mi && length(mi_cand)) {
    invisible(run("_mi", "sensitivity_mi.R", c("--final", mi_cand[1], pass(c("--aux", "--m", "--seed")), "--out", file.path(out, "sensitivity_mi"))))
  } else {
    if (use_mi) cat("只加入多群組時，插補下未提供 M1–M4 的比較；這一步改用 FIML。\n")
    res[nrow(res) + 1, ] <- list("step_mi", "略過（未選多重插補）", 0)   # 讓匯出報告不納入舊的插補結果
  }
}
if (6 %in% steps && file.exists(file.path(out, "step2", "chosen_fit.rds"))) {
  cand <- file.path(out, c(if (!is.null(tvc)) "step4/tvc_fit.rds", if (!is.null(covs)) "step3/conditional_fit.rds", "step2/chosen_fit.rds"))
  invisible(run("6", "step6_report.R", c("--results", out, "--final", cand[file.exists(cand)][1])))
}
# 狀態表：分段執行時保留之前步驟的狀態；這次略過的步驟標為略過（匯出報告時不納入舊的結果）
sf <- file.path(out, "pipeline_status.csv"); all_res <- res
if (file.exists(sf)) { old <- utils::read.csv(sf, fileEncoding = "UTF-8", stringsAsFactors = FALSE); all_res <- rbind(old[!(old$step %in% res$step), ], res) }
write_status <- function(x) utils::write.csv(x[order(x$step %in% c("step_mi", "step_qmd"), x$step == "step_qmd", x$step), ], sf, row.names = FALSE, fileEncoding = "UTF-8")
write_status(all_res)
# 匯出 qmd（一段程式碼、一段輸出）；--render 時一併渲染成 HTML。分段時只在最後一段（含步驟三以後）匯出
if (!has("--no-qmd") && (any(steps >= 3) || (auto && any(steps >= 2))) && file.exists(file.path(out, "step2", "cases.txt"))) {
  invisible(run("_qmd", "export_qmd.R", c("--results", out, pass(c("--name", "--title")), flag("--render"))))
  all_res <- rbind(all_res[all_res$step != "step_qmd", ], res[res$step == "step_qmd", ]); write_status(all_res)
}

if (identical(steps, 1L)) {
  cat(sprintf("\n步驟一完成（確認點 1）。先讀 %s/step1/screen_report.md 的「從篩檢結果看成長形式」與最後的「簡易判讀」，請使用者選擇：a 繼續，比較所有模型（--stage model）／b 直接指定成長模型（--stage model --forms linear 或 linear,quad）／c 先修正設定後重跑篩檢／d 停止。\n", out))
} else if (identical(steps, 2L)) {
  r <- kv(file.path(out, "step2", "recommendation.txt"))
  if (identical(r["compared"], c(compared = "FALSE"))) {
    fg <- form_flags()
    cat(sprintf("\n步驟二完成。使用者於確認點 1 直接指定%s，未做形式比較，略過確認點 2。%s\n", fl[[r[["form"]]]],
                if (length(fg)) sprintf("\n⚠ 請先提醒使用者：%s。詢問是否改為比較所有模型（不帶 --forms 重跑 --stage model）；使用者堅持採用時，用 --stage predictors --forms %s --choose %s 繼續（分析決定紀錄會記下提醒內容）。",
                                        paste(fg, collapse = "；"), r[["form"]], r[["form"]])
                else sprintf("配適與診斷沒有問題；接著是確認點 3，再用 --stage predictors --forms %s --choose %s 繼續。", r[["form"]], r[["form"]])))
  } else
  cat(sprintf("\n步驟二完成（確認點 2）。建議：%s（%s）。比較表在 %s/step2/comparison.csv 與 unconditional_report.md；請使用者選擇%s，再用 --stage predictors --choose <形式> 繼續%s。\n",
              fl[[r[["form"]]]], r[["reason"]], out, paste(fl[intersect(c("none", "linear", "quad"), split_csv(r[["forms"]]))], collapse = "／"),
              if (!is.null(forms_sel)) sprintf("（也要帶 --forms %s）", forms_arg) else ""))
} else cat(sprintf("\n完成：OK %d / 執行 %d。先讀 %s/report/analysis_report.md 與各步驟的 *_report.md。\n",
                   sum(res$status == "OK"), sum(!grepl("^略過", res$status)), out))
if (any(grepl("^FAIL", res$status))) quit(status = 1)
