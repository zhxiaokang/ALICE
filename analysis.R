library(dplyr)
library(openxlsx)
library(survival)
library(survminer)

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

# save the plot and risk table
ggsave("output/survival_plot_treatment.png", plot = print(p), width = 8, height = 6)

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

# Run Cox regression on each gene expression with PFS without stratifying the treatment
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

