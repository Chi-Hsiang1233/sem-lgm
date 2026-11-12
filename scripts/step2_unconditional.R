# step2_unconditional.R —— 無條件成長模型：比較軌跡形式與殘差結構，並做假設診斷
# Rscript step2_unconditional.R --data d.csv --y "^t[0-9]$" [--times 0,1,2,3] [--forms none,linear,quad,basis] [--na -99] --out results/step2
#   --forms：要估計的形式（預設 none,linear,quad,basis）。確認點 1 使用者直接指定模型時傳入，例如 --forms linear 或 --forms linear,quad。
#            只給一個形式：只估計該模型，不產生形式比較（表、LRT、潛在基底的參考比較），報告註明「未做形式比較」，
#            另寫 form_check.txt（配適、收斂、負變異數、殘差平均的提醒）；殘差結構、假設診斷、時間編碼敏感度照常。
#   --choose linear：指定採用的形式（none、linear、quad；使用者主動要求時也可以是 basis）。
#                    沒有指定時採用建議（門檻見 lgm_helpers.R 的 FIT_RULES）：線性相對無成長的 LRT 不顯著（p ≥ .05）→ 無成長；
#                    否則線性，除非二次可採用、LRT 顯著且 BIC 低至少 2（Raftery 的「有證據」）→ 二次。潛在基底只作參考，不列為建議。
#   產出 recommendation.txt（建議的形式與理由；compared=FALSE 代表只估計了一個形式），供確認點 2 與分析決定紀錄使用
#   --anchor first|last|mid|數值：截距放在哪個時間點（預設 first = 第一波）
#   --no-sensitivity：不做時間編碼（截距位置）敏感度分析
# 產出：comparison.csv、residual_structure.csv、baseline_comparison.csv（以潛在基底為基準）、各模型子資料夾、診斷圖、unconditional_report.md、chosen.txt
source(file.path((function() { a <- commandArgs(FALSE); dirname(normalizePath(sub("^--file=", "", a[grep("^--file=", a)][1]))) })(), "lgm_core.R"), encoding = "UTF-8")
a <- parse_args(); need(a, c("data", "y"))
require_pkgs("lavaan"); suppressPackageStartupMessages(library(lavaan))
out <- ensure_dir(arg_or(a, "out", "results/step2"))
d <- read_any(a$data, na = arg_or(a, "na"))
y <- check_cols(wave_vars(a$y, d), d); Tn <- length(y)
times <- if (!is.null(a$times)) as.numeric(split_csv(a$times)) else 0:(Tn - 1)
anchor <- resolve_anchor(a$anchor, times)
est <- arg_or(a, "estimator", "MLR")
forms_given <- !is.null(a$forms)
forms <- split_csv(arg_or(a, "forms", "none,linear,quad,basis"))
bad_form <- setdiff(forms, c("none", "linear", "quad", "basis"))
if (length(bad_form)) stop("--forms 只能是 none、linear、quad、basis：", paste(bad_form, collapse = ", "), call. = FALSE)
if (Tn < 4) {                                              # 3 波只能估線性
  if (forms_given && "quad" %in% forms) message("只有 3 波，不估計二次模型（自由參數會超過資料提供的資訊）。")
  forms <- setdiff(forms, c("quad", "basis"))
  if (!length(forms)) stop("只有 3 波，不能估計二次或潛在基底模型；請改指定無成長或線性。", call. = FALSE)
}
single <- length(forms) == 1                               # 確認點 1 直接指定單一模型：不做形式比較
fl <- c(none = "無成長", linear = "線性", quad = "二次", basis = "潛在基底")
for (f in c("baseline_comparison.csv", "form_check.txt")) unlink(file.path(out, f))   # 清掉上一次執行留下的檔案
rep <- c("# 步驟二：無條件成長模型", "", sprintf("N = %d；%d 波；測量時間 %s；截距位置 %s（時間負荷 %s）；估計法 %s（遺失以 FIML 處理）。", nrow(d), Tn, paste(times, collapse = ", "),
                 anchor_label(anchor, times), paste(times - anchor, collapse = ", "), est),
         if (forms_given) sprintf("估計的形式：%s（使用者於確認點 1 指定%s）。", paste(fl[forms], collapse = "、"), if (single) "，未做形式比較" else "，只比較這些模型"))

# A. 軌跡形式（殘差變異自由估計）
fits <- list()
for (fm in forms) {
  syn <- lgm_syntax(y, times, fm, anchor = anchor)
  writeLines(syn, file.path(out, paste0("model_", fm, ".lav")))
  fits[[fm]] <- fit_lgm(syn, d, estimator = est)
  save_case(fits[[fm]], out, fm)
}
# 比較結論：每個形式和前一個較簡單的形式比（潛在基底和線性比）
seq_forms <- intersect(c("none", "linear", "quad"), forms)
prev <- c(setNames(head(seq_forms, -1), tail(seq_forms, -1)), if (all(c("linear", "basis") %in% forms)) c(basis = "linear"))
tab <- fit_compare_table(fits, prev = if (single) NULL else prev, labels = fl); wtab(tab, out, "comparison.csv")
# 建議（確認點 2）：只在無成長、線性、二次之中建議；潛在基底只作參考。門檻見 FIT_RULES（模型比較門檻）
adm_ok <- function(k) k %in% tab$model[tab$admissible]
bic_of <- function(k) tab$bic[tab$model == k]
R <- FIT_RULES
if (single) {
  rec <- forms; why <- "使用者於確認點 1 直接指定，未做形式比較"
} else {
  # 依無成長 → 線性 → 二次的順序，只在指定的形式中往較複雜的形式走
  rec <- seq_forms[1]; why <- character()
  for (nx in seq_forms[-1]) {
    p <- lrt_p(fits[[rec]], fits[[nx]]); dB <- bic_of(rec) - bic_of(nx)
    if (nx == "quad") {   # 往二次：可採用、LRT 顯著且 BIC 至少低 2
      move <- adm_ok("quad") && !is.na(p) && p < R$alpha && isTRUE(dB >= R$bic_cut[1])
      why <- c(why, if (!adm_ok("quad")) "二次模型不可採用（見 diagnostics.txt）" else
        sprintf("二次相對%s：LRT p %s（門檻 p < %s）；BIC 支持%s，差 %.1f（Raftery：%s；採用二次的門檻是二次的 BIC 低至少 %d）→ %s", fl[[rec]], fmt_p(p), fmt2(R$alpha, 2),
                if (isTRUE(dB > 0)) "二次" else fl[[rec]], abs(dB), bic_evidence(dB), R$bic_cut[1], if (move) "採用二次" else "沒有明顯改善"))
    } else {              # 往線性：LRT 顯著，或無成長不可採用
      move <- is.na(p) || p < R$alpha || !adm_ok(rec)
      why <- c(why, sprintf("%s相對%s：LRT p %s（門檻 p < %s）→ %s", fl[[nx]], fl[[rec]], fmt_p(p), fmt2(R$alpha, 2),
                            if (is.na(p)) "LRT 無法計算，採用較複雜的形式" else if (p < R$alpha) "顯著改善" else if (move) sprintf("不顯著，但%s不可採用", fl[[rec]]) else "沒有顯著改善"))
    }
    if (move) rec <- nx else break
  }
  if (Tn < 4) why <- c(why, "只有 3 波，不估計二次模型（自由參數會超過資料提供的資訊）")
  if (!adm_ok(rec)) why <- c(why, sprintf("⚠ %s模型也不可採用，請看 diagnostics.txt", fl[[rec]]))
  why <- paste(why, collapse = "；")
}
writeLines(c(paste0("form=", rec), paste0("reason=", why), paste0("compared=", !single), paste0("forms=", paste(forms, collapse = ","))), file.path(out, "recommendation.txt"))
chosen <- arg_or(a, "choose", rec)
if (!chosen %in% forms) stop(sprintf("--choose %s 不在這次估計的形式（%s）中；請在 --forms 加入它，或不給 --forms。", chosen, paste(forms, collapse = ",")), call. = FALSE)
if (single) {
  ft1 <- fit_table(fits[[chosen]])
  rep <- c(rep, "", sprintf("## A. 軌跡形式：%s（使用者於確認點 1 直接指定，未做形式比較）", fl[[chosen]]), "",
           sprintf("%s（殘差變異自由估計）的%s", fl[[chosen]], fit_sentence(ft1, est)), "", sprintf("%s（殘差變異自由估計）的整體配適：%s。%s", fl[[chosen]], fit_verdict(ft1$cfi, ft1$tli, ft1$rmsea, ft1$srmr, ft1$df), fit_rules_text()), "",
           "沒有和其他形式比較，所以沒有比較表、LRT 與潛在基底的參考比較；需要時不指定模型、重新執行這一步。")
} else {
  rep <- c(rep, "", "## A. 軌跡形式", "", compare_md(tab), "",
           sprintf("**建議：%s**（%s）。%s", fl[[rec]], why, if ("basis" %in% forms) "潛在基底只作參考（判斷線性是否足夠），不列為建議；使用者主動要求時才採用。" else ""),
           "「整體配適」依整體配適門檻、「比較結論」依模型比較門檻（門檻與出處見 references/09-fit-criteria.md）；這些是經驗法則，要和參數是否合理、理論能否解釋一起判斷。")
}
rep <- c(rep, sprintf("**採用：%s**（%s）。", fl[[chosen]], if (!is.null(a$choose)) "使用者指定" else if (single) "使用者於確認點 1 指定" else "採用建議"),
         if (!adm_ok(chosen)) sprintf("⚠ 採用的模型（%s）不可採用：有負變異數、未收斂或標準化負荷 > 1，見 %s/diagnostics.txt。", chosen, chosen))
# 巢套比較：線性 vs 二次；線性 vs 潛在基底
nest <- c(if (all(c("none", "linear") %in% forms)) lrt_line(fits$none, fits$linear, "無成長 → 線性"),
          if (all(c("linear", "quad") %in% forms)) lrt_line(fits$linear, fits$quad, "線性 → 二次"),
          if (all(c("none", "quad") %in% forms) && !"linear" %in% forms) lrt_line(fits$none, fits$quad, "無成長 → 二次"),
          if (all(c("linear", "basis") %in% forms)) lrt_line(fits$linear, fits$basis, "線性 → 潛在基底"))
if (length(nest)) rep <- c(rep, "", "巢套比較（尺度校正 LRT）：", "", nest)
if ("basis" %in% forms && lavInspect(fits$basis, "converged")) {
  l <- growth_params(fits$basis)$load
  rep <- c(rep, "", sprintf("潛在基底的時間負荷：%s（線性時間為 %s；接近代表線性已足夠）。", paste(sprintf("%.2f", l$est), collapse = ", "), paste(times - anchor, collapse = ", ")))
}
for (fm in forms) {
  dg <- readLines(file.path(out, fm, "diagnostics.txt"), encoding = "UTF-8"); bad <- dg[grepl("嚴重", dg)]
  if (length(bad)) rep <- c(rep, "", sprintf("⚠ %s：%s", fl[[fm]], paste(sub("\\[嚴重\\] ", "", bad), collapse = "；")))
}

seg_A <- rep; rep <- character()   # 報告順序：A 形式 → B 以潛在基底為基準 → C 殘差 → D 採用模型 → E 假設診斷 → F 時間編碼（最後組合）
# C. 殘差結構（以採用的形式比較：殘差變異自由估計 vs 相等）
f_free <- fits[[chosen]]
syn_eq <- lgm_syntax(y, times, chosen, resid = "equal", anchor = anchor)
f_eq <- fit_lgm(syn_eq, d, estimator = est); save_case(f_eq, out, paste0(chosen, "_equal_resid"))
rtab <- fit_compare_table(list(free = f_free, equal = f_eq), prev = c(free = "equal"), labels = c(free = "殘差變異自由估計", equal = "殘差變異相等"))
wtab(rtab, out, "residual_structure.csv")
p_res <- lrt_p(f_eq, f_free)
resid_choice <- arg_or(a, "residuals", if (!is.na(p_res) && p_res >= .05) "equal" else "free")
rep <- c(rep, "", "## C. 殘差結構", "", compare_md(rtab), "", lrt_line(f_eq, f_free, "殘差變異相等 → 自由估計"),
         sprintf("**採用：%s**（%s）。等同限制不被拒絕時，可採較簡約的相等殘差；兩者的成長參數通常差異很小，結論不受影響時任一皆可。", if (resid_choice == "equal") "殘差變異相等" else "殘差變異自由估計",
                 if (!is.null(a$residuals)) "使用者指定" else if (resid_choice == "equal") "LRT 不顯著" else "LRT 顯著"))

seg_res <- rep; rep <- character()
# D. 採用模型的參數
fin <- if (resid_choice == "equal") f_eq else f_free
gp <- growth_params(fin)
rep <- c(rep, "", sprintf("## D. 採用模型（%s＋%s）", fl[[chosen]], if (resid_choice == "equal") "殘差變異相等" else "殘差變異自由估計"), "", fit_sentence(fit_table(fin), est), "",
         sprintf("- 平均 %s：%.3f（SE = %.3f, 95%% CI [%.3f, %.3f], p %s）", gp$mean$lhs, gp$mean$est, gp$mean$se, gp$mean$ci.lower, gp$mean$ci.upper, fmt_p(gp$mean$pvalue)),
         sprintf("- 變異 %s：%.3f（SE = %.3f, p %s）", gp$var$lhs, gp$var$est, gp$var$se, fmt_p(gp$var$pvalue)),
         if (nrow(gp$cor)) sprintf("- 相關 %s, %s：r = %.3f（p %s）", gp$cor$lhs, gp$cor$rhs, gp$cor$std.all, fmt_p(gp$cor$pvalue)))
rep <- c(rep, "", "變異的檢定在參數邊界上（變異不能小於 0），一般 z 檢定偏保守；個別差異是否存在，可比較固定變異為 0 的模型。")

seg_adopt <- rep; rep <- character()
# E. 假設診斷
rep <- c(rep, "", "## E. 假設診斷", "", lgm_diagnostics(fin, d, y, times, out))
obs <- as.numeric(lavInspect(fin, "sampstat")$mean[y])
use <- setdiff(forms, "none"); if (!length(use)) use <- forms   # 只估計無成長時也要畫
use <- use[vapply(use, function(k) lavInspect(fits[[k]], "converged"), TRUE)]
imp <- setNames(lapply(use, function(k) as.numeric(lavInspect(fits[[k]], "implied")$mean[y])), use)
line_chart(file.path(out, "implied_vs_observed.png"), times, c(list(Observed = obs), imp), "Time", "Mean", "Observed (FIML) vs model-implied means")
if (length(use)) rep <- c(rep, sprintf("- 隱含平均與樣本平均的最大差距：%s", paste(sprintf("%s %.2f", fl[use], vapply(imp, function(v) max(abs(v - obs)), 0)), collapse = "；")))

# 確認點 1 直接指定單一模型時的提醒（門檻見 FIT_RULES）：配適不佳、不收斂或負變異數、各波殘差平均有系統性偏離
if (single) {
  ftf <- fit_table(fin); chk <- fit_poor_flags(ftf)
  dg <- unique(c(diagnose(fits[[chosen]]), diagnose(fin))); sev <- sub("^\\[嚴重\\] ", "", dg[grepl("^\\[嚴重\\]", dg)])
  if (length(sev)) chk <- c(chk, sev)
  rb <- file.path(out, "residuals_by_wave.csv")
  if (file.exists(rb)) {
    rs <- utils::read.csv(rb, fileEncoding = "UTF-8"); sdy <- mean(apply(d[y], 2, sd, na.rm = TRUE))
    mx <- max(abs(rs$mean_resid))
    if (mx >= R$resid_mean_sd * sdy)
      chk <- c(chk, sprintf("各波殘差平均有系統性偏離（%s；最大 %.2f，≥ 各波 SD 平均的 %s 倍）", paste(sprintf("%.2f", rs$mean_resid), collapse = "、"), mx, fmt2(R$resid_mean_sd, 2)))
  }
  writeLines(if (length(chk)) paste0("flag=", chk) else "flag=", file.path(out, "form_check.txt"))
  rep <- c(rep, "", "## 指定模型的檢查", "",
           if (length(chk)) c("⚠ 直接指定的模型有下列問題，建議改為比較所有模型（不指定模型重新執行這一步）；使用者堅持採用時照做，並記在分析決定紀錄：", "", paste0("- ⚠ ", chk))
           else sprintf("✔ 採用的模型（%s＋%s）收斂、沒有負變異數，整體配適%s（CFI = %s、RMSEA = %s），各波殘差平均沒有系統性偏離。", fl[[chosen]], if (resid_choice == "equal") "殘差變異相等" else "殘差變異自由估計", fit_verdict(ftf$cfi, ftf$tli, ftf$rmsea, ftf$srmr, ftf$df), fmt2(ftf$cfi), fmt2(ftf$rmsea)))
}

seg_diag <- rep; rep <- character()
# F. 時間編碼敏感度：把截距移到另一端重估，配適應完全相同，但截距相關參數的意義改變
if (!isTRUE(a$no_sensitivity) && chosen != "none") {
  alt <- if (isTRUE(all.equal(anchor, times[Tn]))) times[1] else times[Tn]
  f_alt <- fit_lgm(lgm_syntax(y, times, chosen, resid = resid_choice, anchor = alt), d, estimator = est)
  pick_par <- function(f) {
    pe <- parameterEstimates(f, standardized = TRUE)
    g <- function(l, o, r) { k <- pe[pe$lhs == l & pe$op == o & pe$rhs == r, ]; if (nrow(k)) k[1, ] else NULL }
    slope <- intersect("s", lavNames(f, "lv"))[1]
    rows <- list(`截距平均` = g("i", "~1", ""), `截距變異` = g("i", "~~", "i"),
                 `斜率平均` = if (!is.na(slope)) g(slope, "~1", ""), `截距–斜率相關` = if (!is.na(slope)) g("i", "~~", slope))
    rows <- rows[!vapply(rows, is.null, TRUE)]
    data.frame(parameter = names(rows), est = vapply(rows, function(r) r$est, 0), std = vapply(rows, function(r) r$std.all, 0),
               p = vapply(rows, function(r) r$pvalue, 0), row.names = NULL)
  }
  if (lavInspect(f_alt, "converged")) {
    s1 <- pick_par(fin); s2 <- pick_par(f_alt)
    sens <- data.frame(parameter = s1$parameter, anchor_main = s1$est, anchor_alt = s2$est[match(s1$parameter, s2$parameter)],
                       std_main = s1$std, std_alt = s2$std[match(s1$parameter, s2$parameter)])
    wtab(sens, out, "time_coding_sensitivity.csv")
    show <- function(k) if (grepl("相關", sens$parameter[k])) sprintf("r = %.3f → %.3f", sens$std_main[k], sens$std_alt[k]) else sprintf("%.3f → %.3f", sens$anchor_main[k], sens$anchor_alt[k])
    chi <- c(fit_table(fin)$chisq, fit_table(f_alt)$chisq)   # 和 D 節相同：用 fit_table() 的 χ²（MLR 時是校正後的 χ²）
    rep <- c(rep, "", "## F. 時間編碼敏感度", "",
             sprintf("截距由%s移到%s重新估計（χ²%s %.2f vs %.2f，%s）：", anchor_label(anchor, times), anchor_label(alt, times),
                     if (fit_table(fin)$index_type != "standard") "（校正後，和 D 節相同）" else "", chi[1], chi[2],
                     if (abs(diff(chi)) < 1e-4) "配適完全相同" else "配適不同，請檢查模型設定"), "",
             sprintf("- %s：%s", sens$parameter, vapply(seq_len(nrow(sens)), show, "")), "",
             "斜率平均不變；截距平均、截距變異與截距–斜率相關會隨截距位置改變。報告時要寫明時間 0 的位置，並依研究問題（關心起點或終點）選擇。")
  }
}

seg_time <- rep; rep <- character()
# B. 以潛在基底模型為基準：不預設形狀，再檢查較簡單的模型是否明顯變差
#   無成長、線性（含殘差變異相等）是潛在基底的巢套模型 → LRT；二次不是巢套 → 只比 ΔAIC、ΔBIC
if ("basis" %in% forms && isTRUE(lavInspect(fits$basis, "converged"))) {
  base <- fits$basis; t0 <- times - anchor
  lab_form <- function(k) { eq <- grepl("_equal_resid$", k); b <- sub("_equal_resid$", "", k); paste0(ifelse(b %in% names(fl), fl[b], b), ifelse(eq, "＋殘差變異相等", "")) }
  others <- fits[setdiff(names(fits), "basis")]
  others[[paste0(chosen, "_equal_resid")]] <- f_eq
  nested_forms <- c("none", "linear", "linear_equal_resid")
  ic <- function(f) fitMeasures(f, c("df", "aic", "bic"))
  b <- ic(base)
  btab <- do.call(rbind, lapply(names(others), function(k) {
    f <- others[[k]]; m <- ic(f); nested <- k %in% nested_forms
    l <- if (nested) tryCatch(as.data.frame(lavTestLRT(f, base)), error = function(e) NULL)
    ok <- !is.null(l) && nrow(l) == 2
    data.frame(model = k, df = m[["df"]], d_chisq = if (ok) l[2, "Chisq diff"] else NA, d_df = if (ok) l[2, "Df diff"] else NA,
               p = if (ok) l[2, "Pr(>Chisq)"] else NA, d_aic = m[["aic"]] - b[["aic"]], d_bic = m[["bic"]] - b[["bic"]])
  }))
  wtab(btab, out, "baseline_comparison.csv")
  ld <- growth_params(base)$load$est
  prop <- (ld - t0[1]) / (t0[Tn] - t0[1])
  rep <- c(rep, "", "## B. 以潛在基底模型為基準的比較", "",
           sprintf("潛在基底的累積變化比例：%s（線性預期 %s）。", paste(sprintf("%.2f", prop), collapse = ", "), paste(sprintf("%.2f", (t0 - t0[1]) / (t0[Tn] - t0[1])), collapse = ", ")), "",
           "| 模型 | df | Δχ² | Δdf | p | ΔAIC | ΔBIC |", "|---|---|---|---|---|---|---|",
           sprintf("| %s | %s | %s | %s | %s | %+.1f | %+.1f |", lab_form(btab$model), format(btab$df),
                   ifelse(is.na(btab$d_chisq), "—", sprintf("%.2f", btab$d_chisq)), ifelse(is.na(btab$d_df), "—", format(btab$d_df)),
                   ifelse(is.na(btab$p), "非巢套", sub("^= ", "", vapply(btab$p, fmt_p, ""))), btab$d_aic, btab$d_bic), "",
           "ΔAIC、ΔBIC = 該模型 − 潛在基底，負值代表比基準好；p 不顯著代表加上該模型的限制後配適沒有明顯變差。")
}

seg_base <- rep
rep <- c(seg_A, seg_base, seg_res, seg_adopt, seg_diag, seg_time)
writeLines(c(paste0("form=", chosen), paste0("residuals=", resid_choice), paste0("times=", paste(times, collapse = ",")), paste0("anchor=", anchor)), file.path(out, "chosen.txt"))
saveRDS(fin, file.path(out, "chosen_fit.rds"))
writeLines(rep, file.path(out, "unconditional_report.md"))
cat(paste(rep, collapse = "\n"), "\n")
