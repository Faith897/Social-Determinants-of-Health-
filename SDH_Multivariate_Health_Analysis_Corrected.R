# Corrected reanalysis of the supplied workbook.
# Binary outcomes are numeric 0=No and 1=Yes. Original supplied scripts are preserved.
# Adjusted associations use mean bias-reduced logistic regression (Firth/Jeffreys).
# LASSO is exploratory/predictive only; nested CV evaluates the complete tuning procedure.

suppressPackageStartupMessages({
  library(readxl)
  library(glmnet)
  library(brglm2)
})

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)][1])
if (is.na(script) || !nzchar(script)) stop("Run this file using Rscript.")
project <- normalizePath(dirname(script), winslash = "/", mustWork = TRUE)
workbook <- file.path(project, "MultivariateData.xlsx")
outdir <- file.path(project, "Reanalysis_Output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(workbook)) stop("Workbook not found: ", workbook)

raw <- as.data.frame(readxl::read_excel(workbook, sheet = "Data"))
codebook <- as.data.frame(readxl::read_excel(workbook, sheet = "CODES"))
original_names <- names(raw)
names(raw) <- make.names(original_names, unique = TRUE)
nm <- setNames(names(raw), original_names)

outcome_labels <- c(
  "Cardiovascular Diseases (CVD)", "Diabetes",
  "Any Chronic Respiratory Diseases", "Obese",
  "Mental Health Disorders", "Chronic Kidney Disease (CKD)"
)
predictor_labels <- c(
  "Age", "Income", "Education Level", "Does regular physical exercise",
  "Smoking", "Informal Housing", "Exposure to air pollution",
  "Any form of food insecurity", "Employment", "Exposure to violent crime"
)
outcomes <- unname(nm[outcome_labels])
idvar <- nm[["ResponseID"]]
inc <- nm[["Income"]]
edu <- nm[["Education Level"]]
bin_labels <- setdiff(predictor_labels, c("Age", "Income", "Education Level"))
binvars <- unname(nm[bin_labels])
if (anyNA(c(outcomes, idvar, inc, edu, binvars))) stop("Required workbook column missing.")

for (lab in outcome_labels) {
  vals <- na.omit(as.character(codebook[[lab]]))
  if (!all(c("1=Yes", "0=No") %in% vals)) stop("Unexpected outcome codebook: ", lab)
}
for (lab in bin_labels) {
  codebook_col <- if (lab == "Exposure to violent crime") {
    "Exposure to violent crime (Crime rating in the neighbourhood)"
  } else lab
  vals <- na.omit(as.character(codebook[[codebook_col]]))
  if (!all(c("1=Yes", "0=No") %in% vals)) {
    stop("Unexpected binary predictor codebook: ", codebook_col)
  }
}
income_codes <- na.omit(as.character(codebook[["Income"]]))
education_codes <- na.omit(as.character(codebook[["Education Level"]]))
if (!all(c("1=Lower Income", "2=Middle Income", "3=High Income") %in% income_codes)) {
  stop("Unexpected income codebook.")
}
if (!all(c("1=Basic education", "2=Secondary education", "3=Tertiary Education") %in% education_codes)) {
  stop("Unexpected education codebook.")
}
dat <- raw
if (anyDuplicated(dat[[idvar]])) stop("ResponseID is not unique.")
dat$Age <- as.numeric(dat[[nm[["Age"]]]])
dat[[inc]] <- factor(as.integer(dat[[inc]]), levels = 1:3,
                     labels = c("Lower", "Middle", "High"))
dat[[edu]] <- factor(as.integer(dat[[edu]]), levels = 1:3,
                     labels = c("Basic", "Secondary", "Tertiary"))
for (v in binvars) {
  z <- as.integer(dat[[v]])
  if (!all(na.omit(z) %in% c(0L, 1L))) stop("Non-binary predictor values in ", v)
  dat[[v]] <- factor(z, levels = 0:1, labels = c("No", "Yes"))
}
for (v in outcomes) {
  z <- as.integer(dat[[v]])
  if (!all(na.omit(z) %in% c(0L, 1L))) stop("Non-binary outcome values in ", v)
  dat[[v]] <- z
}
required <- c(idvar, outcomes, unname(nm[predictor_labels]))
keep <- complete.cases(dat[, required, drop = FALSE])
a <- dat[keep, , drop = FALSE]
n <- nrow(a)
a$Age10 <- (a$Age - mean(a$Age)) / 10
pred <- c("Age10", inc, edu, binvars)
formula_x <- reformulate(pred)
mm <- model.matrix(formula_x, data = a)
x <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
assign <- attr(mm, "assign")[colnames(mm) != "(Intercept)"]
term_names <- attr(terms(formula_x), "term.labels")
col_to_term <- setNames(term_names[assign], colnames(x))

if (!all(vapply(a[outcomes], function(z) all(sort(unique(z)) %in% c(0, 1)), logical(1)))) {
  stop("Outcome recoding check failed.")
}

# Bias-reduced logistic regression with explicit event=1 and stable iteration controls.
fit_br <- function(yname, data = a) {
  f <- reformulate(pred, response = yname)
  glm(f, data = data, family = binomial("logit"),
      method = brglm2::brglmFit, type = "AS_mean",
      control = brglm2::brglmControl(maxit = 1000, epsilon = 1e-8,
                                     slowit = 0.25, response_adjustment = 0.1))
}
models <- setNames(lapply(outcomes, fit_br), outcome_labels)

pretty_term <- function(term) {
  if (term == "Age10") return("Age (per 10 years)")
  specs <- list(
    list(v = inc, label = "Income", ref = c(Middle = "Lower", High = "Lower")),
    list(v = edu, label = "Education", ref = c(Secondary = "Basic", Tertiary = "Basic"))
  )
  for (v in binvars) {
    orig <- names(nm)[match(v, nm)]
    specs[[length(specs) + 1L]] <- list(v = v, label = orig, ref = c(Yes = "No"))
  }
  for (s in specs) {
    if (startsWith(term, s$v)) {
      level <- substring(term, nchar(s$v) + 1L)
      if (level %in% names(s$ref)) return(paste0(s$label, ": ", level, " vs ", s$ref[[level]]))
    }
  }
  term
}
coef_rows <- list()
model_rows <- list()
for (i in seq_along(outcomes)) {
  m <- models[[i]]
  sm <- coef(summary(m))
  sm <- sm[rownames(sm) != "(Intercept)", , drop = FALSE]
  b <- sm[, "Estimate"]
  se <- sm[, "Std. Error"]
  coef_rows[[i]] <- data.frame(
    Outcome = outcome_labels[i], Term = rownames(sm),
    Predictor = vapply(rownames(sm), pretty_term, character(1)),
    Log_OR = b, Adjusted_OR = exp(b),
    CI95_Lower = exp(b - qnorm(.975) * se),
    CI95_Upper = exp(b + qnorm(.975) * se),
    SE = se, P_value = sm[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  model_rows[[i]] <- data.frame(
    Outcome = outcome_labels[i], N = nobs(m),
    Events_Yes = sum(m$y == 1), NonEvents_No = sum(m$y == 0),
    Model_DF = length(coef(m)) - 1L,
    Events_per_Model_DF = sum(m$y == 1) / (length(coef(m)) - 1L),
    Converged = isTRUE(m$converged)
  )
}
or_table <- do.call(rbind, coef_rows)
or_table$BH_q_value <- p.adjust(or_table$P_value, method = "BH")
model_table <- do.call(rbind, model_rows)
write.csv(or_table, file.path(outdir, "adjusted_ORs_Firth.csv"), row.names = FALSE)
write.csv(model_table, file.path(outdir, "model_sample_and_convergence.csv"), row.names = FALSE)

# Generalized VIFs for categorical terms and VIFs for one-df terms.
if (requireNamespace("car", quietly = TRUE)) {
  vifs <- list()
  for (i in seq_along(models)) {
    v <- tryCatch(car::vif(models[[i]]), error = function(e) e)
    if (inherits(v, "error")) {
      vifs[[i]] <- data.frame(Outcome = outcome_labels[i], Term = NA,
                              GVIF = NA, Df = NA, Adjusted_GVIF = NA,
                              Diagnostic = conditionMessage(v))
    } else if (is.matrix(v)) {
      vifs[[i]] <- data.frame(Outcome = outcome_labels[i], Term = rownames(v),
                              GVIF = v[, 1], Df = v[, 2], Adjusted_GVIF = v[, 3],
                              Diagnostic = "GVIF^(1/(2*Df))")
    } else {
      vifs[[i]] <- data.frame(Outcome = outcome_labels[i], Term = names(v),
                              GVIF = as.numeric(v), Df = 1,
                              Adjusted_GVIF = sqrt(as.numeric(v)),
                              Diagnostic = "sqrt(VIF)")
    }
  }
  write.csv(do.call(rbind, vifs), file.path(outdir, "diagnostic_VIF.csv"), row.names = FALSE)
}

# Box-Tidwell age-by-log(age) test, using the Firth fit and BH adjustment.
a$AgeLogBT <- (a$Age * log(a$Age) - mean(a$Age * log(a$Age))) / 100
other_pred <- setdiff(pred, "Age10")
bt_rows <- list()
for (i in seq_along(outcomes)) {
  f <- reformulate(c("Age10", "AgeLogBT", other_pred), response = outcomes[i])
  m <- glm(f, data = a, family = binomial("logit"), method = brglm2::brglmFit,
           type = "AS_mean",
           control = brglm2::brglmControl(maxit = 1000, epsilon = 1e-8,
                                          slowit = 0.25, response_adjustment = 0.1))
  s <- coef(summary(m))
  bt_rows[[i]] <- data.frame(Outcome = outcome_labels[i],
                             Age_nonlinearity_P = s["AgeLogBT", "Pr(>|z|)"])
}
age_bt <- do.call(rbind, bt_rows)
age_bt$BH_q_value <- p.adjust(age_bt$Age_nonlinearity_P, method = "BH")
write.csv(age_bt, file.path(outdir, "diagnostic_age_linearity.csv"), row.names = FALSE)

# Cook's D > 4/n flags for sensitivity-only fits; primary estimates retain all records.
influence <- list()
influence_ids <- list()
sensitivity <- list()
threshold <- 4 / n
for (i in seq_along(models)) {
  m <- models[[i]]
  cd <- cooks.distance(m)
  flagged <- which(is.finite(cd) & cd > threshold)
  influence[[i]] <- data.frame(
    Outcome = outcome_labels[i], N = n, Threshold_4_over_n = threshold,
    Influential_N = length(flagged), Maximum_Cooks_D = max(cd, na.rm = TRUE)
  )
  if (length(flagged)) {
    influence_ids[[i]] <- data.frame(
      Outcome = outcome_labels[i], Row = flagged,
      ResponseID = as.character(a[[idvar]][flagged]), Cooks_D = cd[flagged]
    )
  }
  reduced <- a[-flagged, , drop = FALSE]
  mr <- fit_br(outcomes[i], reduced)
  sr <- coef(summary(mr))
  full <- or_table[or_table$Outcome == outcome_labels[i], ]
  for (term in full$Term) {
    if (term %in% rownames(sr)) {
      b <- sr[term, "Estimate"]
      se <- sr[term, "Std. Error"]
      fr <- full[match(term, full$Term), ]
      sensitivity[[length(sensitivity) + 1L]] <- data.frame(
        Outcome = outcome_labels[i], Predictor = pretty_term(term),
        Full_OR = fr$Adjusted_OR, Full_CI95_Lower = fr$CI95_Lower,
        Full_CI95_Upper = fr$CI95_Upper, Full_P_value = fr$P_value,
        Full_BH_q_value = fr$BH_q_value, OR_without_influential = exp(b),
        CI95_Lower_without_influential = exp(b - qnorm(.975) * se),
        CI95_Upper_without_influential = exp(b + qnorm(.975) * se),
        P_value_without_influential = sr[term, "Pr(>|z|)"],
        OR_ratio_without_vs_full = exp(b) / fr$Adjusted_OR,
        Removed_n = length(flagged)
      )
    }
  }
}
write.csv(do.call(rbind, influence), file.path(outdir, "diagnostic_influence_summary.csv"), row.names = FALSE)
if (length(influence_ids)) {
  write.csv(do.call(rbind, influence_ids), file.path(outdir, "diagnostic_influential_ResponseIDs.csv"), row.names = FALSE)
}
sens <- do.call(rbind, sensitivity)
sens$BH_q_value_without_influential <- p.adjust(sens$P_value_without_influential, method = "BH")
write.csv(sens, file.path(outdir, "diagnostic_influence_sensitivity.csv"), row.names = FALSE)

# Repeated nested CV: repeated stratified outer 5-fold, with 5-fold lambda
# tuning inside each training set. No outer test outcomes select lambda.
folds <- function(y, k) {
  f <- integer(length(y))
  for (cl in sort(unique(y))) {
    ids <- sample(which(y == cl))
    f[ids] <- rep(seq_len(k), length.out = length(ids))
  }
  f
}
auc <- function(y, p) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
calibration <- function(y, p) {
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  z <- tryCatch(suppressWarnings(glm(y ~ lp, family = binomial())), error = function(e) NULL)
  if (is.null(z)) return(c(Intercept = NA_real_, Slope = NA_real_))
  c(Intercept = unname(coef(z)[1]), Slope = unname(coef(z)[2]))
}
set.seed(20261010)
repeats <- 10L
outer_k <- 5L
inner_k <- 5L
repeat_rows <- list()
fold_rows <- list()
selection_rows <- list()
prediction_rows <- list()
for (i in seq_along(outcomes)) {
  y <- a[[outcomes[i]]]
  for (rep in seq_len(repeats)) {
    outer <- folds(y, outer_k)
    p <- rep(NA_real_, n)
    for (k in seq_len(outer_k)) {
      test <- which(outer == k); train <- which(outer != k)
      inner <- folds(y[train], inner_k)
      fit <- cv.glmnet(x[train, , drop = FALSE], y[train],
                       family = "binomial", alpha = 1,
                       nfolds = inner_k, foldid = inner,
                       type.measure = "deviance", standardize = TRUE)
      p[test] <- as.numeric(predict(fit, newx = x[test, , drop = FALSE],
                                    s = "lambda.min", type = "response"))
      cf <- as.matrix(coef(fit, s = "lambda.min"))
      cols <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
      selected <- unique(na.omit(unname(col_to_term[cols])))
      for (term in unique(term_names)) {
        selection_rows[[length(selection_rows) + 1L]] <- data.frame(
          Outcome = outcome_labels[i], Repeat = rep, Fold = k,
          Predictor = term, Selected = term %in% selected
        )
      }
      fold_rows[[length(fold_rows) + 1L]] <- data.frame(
        Outcome = outcome_labels[i], Repeat = rep, Fold = k,
        N_Train = length(train), N_Test = length(test),
        Events_Test = sum(y[test] == 1),
        Lambda_Min = fit$lambda.min, Lambda_1SE = fit$lambda.1se
      )
    }
    cal <- calibration(y, p)
    repeat_rows[[length(repeat_rows) + 1L]] <- data.frame(
      Outcome = outcome_labels[i], Repeat = rep, AUC = auc(y, p),
      Brier = mean((y - p)^2),
      Sensitivity_at_0_5 = mean(p[y == 1] >= .5),
      Specificity_at_0_5 = mean(p[y == 0] < .5),
      Calibration_Intercept = cal[1], Calibration_Slope = cal[2]
    )
    prediction_rows[[length(prediction_rows) + 1L]] <- data.frame(
      Outcome = outcome_labels[i], ID = a[[idvar]], Response = y,
      Repeat = rep, Probability = p
    )
  }
}
cv_rep <- do.call(rbind, repeat_rows)
cv_fold <- do.call(rbind, fold_rows)
sel <- do.call(rbind, selection_rows)
cv_pred <- do.call(rbind, prediction_rows)
metrics <- c("AUC", "Brier", "Sensitivity_at_0_5", "Specificity_at_0_5",
             "Calibration_Intercept", "Calibration_Slope")
cv_summary <- do.call(rbind, lapply(split(cv_rep, cv_rep$Outcome), function(z) {
  out <- data.frame(Outcome = z$Outcome[1], Repeats = nrow(z))
  for (m in metrics) {
    out[[paste0(m, "_Mean")]] <- mean(z[[m]], na.rm = TRUE)
    out[[paste0(m, "_SD")]] <- sd(z[[m]], na.rm = TRUE)
  }
  out
}))
bootstrap_auc_ci <- function(y, predictions, B = 2000L) {
  ids0 <- which(y == 0); ids1 <- which(y == 1)
  if (!length(ids0) || !length(ids1)) return(c(Lower = NA_real_, Upper = NA_real_))
  boot <- replicate(B, {
    sampled <- c(sample(ids0, length(ids0), replace = TRUE),
                 sample(ids1, length(ids1), replace = TRUE))
    mean(vapply(seq_len(ncol(predictions)), function(j) {
      auc(y[sampled], predictions[sampled, j])
    }, numeric(1)))
  })
  setNames(quantile(boot, c(.025, .975), names = FALSE, type = 6),
           c("Lower", "Upper"))
}
set.seed(20261011)
cv_auc_ci <- do.call(rbind, lapply(split(cv_pred, cv_pred$Outcome), function(z) {
  wide <- reshape(z[, c("ID", "Response", "Repeat", "Probability")],
                  idvar = c("ID", "Response"), timevar = "Repeat",
                  direction = "wide")
  pcols <- paste0("Probability.", seq_len(repeats))
  prediction_matrix <- as.matrix(wide[, pcols, drop = FALSE])
  ci <- bootstrap_auc_ci(wide$Response, prediction_matrix)
  data.frame(Outcome = z$Outcome[1], AUC_CI_Lower = ci[1],
             AUC_CI_Upper = ci[2], AUC_CI_Resamples = 2000L)
}))
cv_summary <- merge(cv_summary, cv_auc_ci, by = "Outcome", all.x = TRUE,
                    sort = FALSE)
cv_summary <- cv_summary[match(outcome_labels, cv_summary$Outcome), ]
sel_summary <- aggregate(Selected ~ Outcome + Predictor, sel, mean)
names(sel_summary)[3] <- "Selection_Frequency"
write.csv(cv_summary, file.path(outdir, "nested_CV_performance.csv"), row.names = FALSE)
write.csv(cv_auc_ci, file.path(outdir, "nested_CV_AUC_bootstrap_CI.csv"), row.names = FALSE)
write.csv(cv_rep, file.path(outdir, "nested_CV_repeat_metrics.csv"), row.names = FALSE)
write.csv(cv_fold, file.path(outdir, "nested_CV_fold_tuning.csv"), row.names = FALSE)
write.csv(sel_summary, file.path(outdir, "nested_CV_selection_frequency.csv"), row.names = FALSE)

# Full-sample LASSO penalties and term selection are descriptive only.
lambda_rows <- list()
full_selected <- list()
for (i in seq_along(outcomes)) {
  y <- a[[outcomes[i]]]
  set.seed(4100 + i)
  fid <- folds(y, 10L)
  fit <- cv.glmnet(x, y, family = "binomial", alpha = 1,
                   nfolds = 10, foldid = fid, type.measure = "deviance")
  lambda_rows[[i]] <- data.frame(
    Outcome = outcome_labels[i], Lambda_Min = fit$lambda.min,
    Lambda_1SE = fit$lambda.1se, Deviance_Min = min(fit$cvm, na.rm = TRUE), N = n
  )
  for (choice in c("lambda.min", "lambda.1se")) {
    cf <- as.matrix(coef(fit, s = choice))
    cols <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
    terms <- unique(na.omit(unname(col_to_term[cols])))
    if (length(terms)) full_selected[[length(full_selected) + 1L]] <-
      data.frame(Outcome = outcome_labels[i], Lambda = choice, Predictor = terms)
  }
}
write.csv(do.call(rbind, lambda_rows), file.path(outdir, "full_sample_LASSO_lambdas.csv"), row.names = FALSE)
if (length(full_selected)) {
  write.csv(do.call(rbind, full_selected), file.path(outdir, "full_sample_LASSO_selected_predictors.csv"), row.names = FALSE)
}

# Exploratory CCA: treatment-coded predictors, binary outcomes, standardized inputs.
# In-sample roots/loadings/redundancy are descriptive and not validated inference.
xs <- scale(x)
ys <- scale(as.matrix(a[, outcomes, drop = FALSE]))
storage.mode(ys) <- "numeric"
if (qr(xs)$rank < ncol(xs) || qr(ys)$rank < ncol(ys)) stop("CCA matrix is rank deficient.")
cca <- cancor(xs, ys, xcenter = FALSE, ycenter = FALSE)
nr <- length(cca$cor)
xc <- cca$xcoef[, seq_len(nr), drop = FALSE]
yc <- cca$ycoef[, seq_len(nr), drop = FALSE]
xload <- cor(xs, xs %*% xc)
yload <- cor(ys, ys %*% yc)
xve <- colMeans(xload^2); yve <- colMeans(yload^2)
cca_roots <- data.frame(
  Root = seq_len(nr), Canonical_Correlation = cca$cor,
  X_Variance_Extracted = xve, Y_Variance_Extracted = yve,
  X_Redundancy_Explained_by_Y = xve * cca$cor^2,
  Y_Redundancy_Explained_by_X = yve * cca$cor^2
)
write.csv(cca_roots, file.path(outdir, "CCA_roots_and_redundancy.csv"), row.names = FALSE)
loadings <- data.frame(Variable = c(colnames(xs), outcome_labels),
                       Set = c(rep("Predictor", ncol(xs)), rep("Outcome", ncol(ys))),
                       rbind(xload, yload), check.names = FALSE)
names(loadings)[-(1:2)] <- paste0("Root", seq_len(nr))
write.csv(loadings, file.path(outdir, "CCA_structure_loadings.csv"), row.names = FALSE)
coefs <- data.frame(Variable = c(colnames(xs), outcome_labels),
                    Set = c(rep("Predictor", ncol(xs)), rep("Outcome", ncol(ys))),
                    rbind(xc, yc), check.names = FALSE)
names(coefs)[-(1:2)] <- paste0("Root", seq_len(nr))
write.csv(coefs, file.path(outdir, "CCA_standardized_coefficients.csv"), row.names = FALSE)

# Final machine-readable audit and environment record.
writeLines(c(
  paste("Run date:", format(Sys.Date(), "%Y-%m-%d")),
  paste("R:", R.version.string),
  paste("Workbook MD5:", unname(tools::md5sum(workbook))),
  paste("Complete cases:", n, "of", nrow(dat)),
  paste("Outcome event coding:", "1=Yes; 0=No"),
  paste("References:", "Income Lower; Education Basic; binary predictors No"),
  paste("Codebook mappings verified:", "all binary variables 1=Yes, 0=No"),
  paste("Adjusted models:", "brglm2 mean bias reduction (AS_mean/Firth)"),
  paste("Nested CV:", repeats, "repeats of", outer_k, "-fold outer CV;",
        inner_k, "-fold inner tuning by binomial deviance"),
  paste("AUC intervals:", "95% stratified participant-cluster percentile bootstrap",
        "of mean AUC across repeats; conditional on fitted out-of-fold predictions;",
        2000, "resamples; seed 20261011"),
  paste("Random seed:", 20261010),
  paste("Package versions:",
        paste(c(paste0("readxl ", packageVersion("readxl")),
                paste0("glmnet ", packageVersion("glmnet")),
                paste0("brglm2 ", packageVersion("brglm2")),
                if (requireNamespace("car", quietly = TRUE))
                  paste0("car ", packageVersion("car"))), collapse = "; "))
), file.path(outdir, "run_metadata.txt"))
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
cat("Corrected reanalysis completed: ", outdir, "\n", sep = "")


