library(dplyr)
library(openxlsx)
library(survival)
library(survminer)
library(coxphf)
library(forestplot)

# Load the dataset
data <- read.xlsx("data/Assignment Data Nov 2025.xlsx")

# remove the rows with missing values
df <- na.omit(data)

# draw KM plot to study the difference in PFS between the treatment groups
# corresponding column names: "treatment", "time.pfs", "event.pfs"

# Create a survival object
surv_object <- Surv(time = df$time.pfs, event = df$event.pfs)

# Fit the survival curves
fit <- survfit(surv_object ~ df$treatment)

# Plot the survival curves
p <- ggsurvplot(fit, data = df, pval = TRUE, conf.int = TRUE, risk.table = TRUE)

# Run Cox regression to estimate the hazard ratio (HR) between treatment groups
# Check wether Cox regression analyses have met the assumption of proportional hazards
cox_model <- coxph(surv_object ~ df$treatment)
ph_test <- cox.zph(cox_model)
ph_p_value <- ph_test$table[1, "p"]

# Check the Cox model summary
beta <- coef(cox_model)[1]
se <- sqrt(vcov(cox_model)[1, 1])
cox_treatment_result <- data.frame(
    HR = exp(beta),
    HR_lower = exp(beta - 1.96 * se),
    HR_upper = exp(beta + 1.96 * se),
    p_value = summary(cox_model)$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
)

# List the column names of the 8 gene expressions
cols_to_exclude <- c("treatment", "time.pfs", "event.pfs", "id")
cols_exp <- setdiff(colnames(df), cols_to_exclude)

# Check whether Cox regression analyses have met the assumption of proportional hazards
ph_assumptions <- vector("list", length(cols_exp))
for (i in seq_along(cols_exp)) {
    gene <- cols_exp[i]
    cox_formula <- as.formula(sprintf("Surv(time.pfs, event.pfs) ~ %s + strata(treatment)", gene))
    model <- coxph(cox_formula, data = df)
    ph_test <- cox.zph(model)
    ph_assumptions[[i]] <- data.frame(
        gene = gene,
        p_value = ph_test$table[1, "p"],
        stringsAsFactors = FALSE
    )
}
ph_assumptions <- bind_rows(ph_assumptions)

# Run Cox regression on each gene expression with PFS, stratifying the treatment
cox_gene_results <- vector("list", length(cols_exp))
for (i in seq_along(cols_exp)) {
    gene <- cols_exp[i]
    cox_formula <- as.formula(sprintf("Surv(time.pfs, event.pfs) ~ %s + strata(treatment)", gene))
    model <- coxph(cox_formula, data = df)
    beta <- coef(model)[1]
    se <- sqrt(vcov(model)[1, 1])
    cox_gene_results[[i]] <- data.frame(
        HR = exp(beta),
        HR_lower = exp(beta - 1.96 * se),
        HR_upper = exp(beta + 1.96 * se),
        p_value = summary(model)$coefficients[1, "Pr(>|z|)"],
        stringsAsFactors = FALSE
    )
}
cox_gene_results <- bind_rows(cox_gene_results)

# Still run Cox regression on each gene expression with PFS, stratifying the treatment, but considering the small sample issue
cox_gene_small_sample <- vector("list", length(cols_exp))
for (i in seq_along(cols_exp)) {
    gene <- cols_exp[i]
    cox_formula <- as.formula(sprintf("Surv(time.pfs, event.pfs) ~ %s + strata(treatment)", gene))
    penalized_fit <- coxphf::coxphf(cox_formula, data = df)
    beta_pena <- penalized_fit$coefficients[1]
    se_pena <- sqrt(diag(penalized_fit$var))[1]
    cox_gene_small_sample[[i]] <- data.frame(
        HR = exp(beta_pena),
        HR_lower = exp(beta_pena - 1.96 * se_pena),
        HR_upper = exp(beta_pena + 1.96 * se_pena),
        p_value = penalized_fit$prob[1],
        stringsAsFactors = FALSE
    )
}
cox_gene_small_sample <- bind_rows(cox_gene_small_sample)

# ==============================================================================
# PREDICTIVE BIOMARKER ANALYSIS
# ==============================================================================
# Divide patients into two groups based on median gene expression scores
# and study whether any scores are potential predictive biomarkers

biomarker_results <- vector("list", length(cols_exp))
biomarker_plots <- vector("list", length(cols_exp))
cox_pvalues <- vector("list", length(cols_exp))
for (i in seq_along(cols_exp)) {
    gene <- cols_exp[i]

    # Calculate median for this gene
    gene_median <- median(df[[gene]])

    # Create binary groups: High (above median) and Low (at or below median)
    df[[paste0(gene, "_group")]] <- ifelse(df[[gene]] > gene_median, "High", "Low")

    # Ensure roughly equal group sizes by checking distribution
    group_counts <- table(df[[paste0(gene, "_group")]], useNA = "ifany")

    cat("\n", gene, "median:", round(gene_median, 3), "\n")
    cat("Group sizes:", group_counts, "\n")

    # Test for predictive biomarker effect using treatment interaction
    # A predictive biomarker shows differential treatment effect between high/low groups

    # Create interaction term for Cox regression: gene_group * treatment
    cox_formula_interaction <- as.formula(paste("Surv(time.pfs, event.pfs) ~",
                                              paste0(gene, "_group"), "+ treatment +",
                                              paste0(gene, "_group"), ": treatment"))

    # Fit Cox model with interaction
    cox_interaction <- coxphf::coxphf(cox_formula_interaction, data = df)

    # Extract interaction p-value (key test for predictive biomarker)
    interaction_pvalue <- cox_interaction$prob[grep(":", names(cox_interaction$coefficients))]

    # Calculate treatment effect (HR) in High/Low group
    # High expression group: Treatment B vs A
    df_high_group <- df[df[[paste0(gene, "_group")]] == "High", ]
    if(length(unique(df_high_group$treatment)) == 2) {  # both A and B are there
        cox_high <- coxphf::coxphf(Surv(time.pfs, event.pfs) ~ treatment, data = df_high_group)
        hr_high <- exp(cox_high$coefficients[1])
        hr_high_se <- sqrt(diag(cox_high$var))[1]
        hr_high_ci <- c(exp(cox_high$coefficients[1] - 1.96 * hr_high_se), exp(cox_high$coefficients[1] + 1.96 * hr_high_se))
        high_pvalue <- cox_high$prob[1]
    } else {
        high_pvalue <- NA
    }

    # Low expression group: Treatment B vs A
    df_low_group <- df[df[[paste0(gene, "_group")]] == "Low" & !is.na(df[[paste0(gene, "_group")]]), ]
    if(length(unique(df_low_group$treatment)) == 2) {
        cox_low <- coxphf::coxphf(Surv(time.pfs, event.pfs) ~ treatment, data = df_low_group)
        hr_low <- exp(cox_low$coefficients[1])
        hr_low_se <- sqrt(diag(cox_low$var))[1]
        hr_low_ci <- c(exp(cox_low$coefficients[1] - 1.96 * hr_low_se), exp(cox_low$coefficients[1] + 1.96 * hr_low_se))
        low_pvalue <- cox_low$prob[1]
    } else {
        low_pvalue <- NA
    }

    # Store results
    cox_pvalues[[i]] <- data.frame(
        gene = gene,
        interaction_pvalue = interaction_pvalue,
        high_group_pvalue = high_pvalue,
        low_group_pvalue = low_pvalue,
        stringsAsFactors = FALSE
    )

    biomarker_results[[i]] <- data.frame(
        gene = gene,
        hr_high = hr_high,
        hr_high_lower = hr_high_ci[1],
        hr_high_upper = hr_high_ci[2],
        high_pvalue = high_pvalue,
        hr_low = hr_low,
        hr_low_lower = hr_low_ci[1],
        hr_low_upper = hr_low_ci[2],
        low_pvalue = low_pvalue,
        stringsAsFactors = FALSE
    )

    # Create survival plots for each gene
    df_surv_plot <- df
    df_surv_plot$group_treatment <- paste(df_surv_plot[[paste0(gene, "_group")]],
                                           df_surv_plot$treatment, sep = "_")

    surv_object <- Surv(time = df_surv_plot$time.pfs,
                         event = df_surv_plot$event.pfs)
    fit <- survfit(surv_object ~ df_surv_plot$group_treatment)

    biomarker_plots[[i]] <- ggsurvplot(
        fit,
        data = df_surv_plot,
        pval = FALSE,
        risk.table = TRUE,
        legend.title = paste(gene, "Expression"),
        legend.labs = c("High_A", "High_B", "Low_A", "Low_B"),
        title = paste("Survival by", gene, "Expression and Treatment"),
        subtitle = paste("Interaction p-value:", round(interaction_pvalue, 4)),
        xlab = "Time (PFS)",
        ylab = "Survival Probability",
        palette = c("red", "darkred", "blue", "darkblue")
    )
}
biomarker_results <- bind_rows(biomarker_results)
cox_pvalues <- bind_rows(cox_pvalues)


# Create forest plot showing treatment effects by biomarker groups

# Prepare data for forest plot - showing HR for treatment effect in high vs low groups
df_forest <- data.frame()
for (i in 1:nrow(biomarker_results)) {
    gene <- biomarker_results$gene[i]

    # High group
    df_forest <- rbind(df_forest, data.frame(
        subgroup = paste(gene, "High"),
        hr = biomarker_results$hr_high[i],
        lower = biomarker_results$hr_high_lower[i],
        upper = biomarker_results$hr_high_upper[i],
        pvalue = biomarker_results$high_pvalue[i]
    ))

    # Low group
    df_forest <- rbind(df_forest, data.frame(
        subgroup = paste(gene, "Low"),
        hr = biomarker_results$hr_low[i],
        lower = biomarker_results$hr_low_lower[i],
        upper = biomarker_results$hr_low_upper[i],
        pvalue = biomarker_results$low_pvalue[i]
    ))
}

# Create a forest plot visualization
table_text <- cbind(
    c("Biomarker subgroup", df_forest$subgroup),
    c("HR (95% CI)", sprintf("%.2f (%.2f, %.2f)", df_forest$hr, df_forest$lower, df_forest$upper)),
    c("p-value", signif(df_forest$pvalue, 3))
)

options(repr.plot.width = 10, repr.plot.height = 8, repr.plot.res = 150)

p_forest <- forestplot(
    labeltext = table_text,
    mean = c(NA, df_forest$hr),
    lower = c(NA, df_forest$lower),
    upper = c(NA, df_forest$upper),
    zero = 1,
    xlog = TRUE,
    clip = c(0.05, 20),
    boxsize = 0.2,
    line.margin = 0.2,
    xlab = "Hazard Ratio (log scale)",
    title = "Treatment Effect by Biomarker Subgroups\nTreatment B vs A (HR < 1 favors Treatment B)",
    txt_gp = fpTxtGp(
        label = grid::gpar(fontsize = 12),
        ticks = grid::gpar(fontsize = 18),
        xlab = grid::gpar(fontsize = 18, fontface = "bold"),
        title = grid::gpar(fontsize = 14, fontface = "bold")
    ),
    col = fpColors(box = "royalblue", line = "darkblue", zero = "red")
)


