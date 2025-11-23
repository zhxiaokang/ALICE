```R
suppressPackageStartupMessages({
    library(dplyr)
    library(openxlsx)
    library(survival)
    library(survminer)
    library(coxphf)
    library(forestplot)
    library(stepp)
})
```


```R
# Load the dataset
data <- read.xlsx("data/Assignment Data Nov 2025.xlsx")
```


```R
# remove the rows with missing values
df <- na.omit(data)
```

## Question 1: Is there a difference in PFS between the treatment groups? Please elaborate.


```R
# draw KM plot to study the difference in PFS between the treatment groups

# Create a survival object
surv_object <- Surv(time = df$time.pfs, event = df$event.pfs)

# Fit the survival curves
fit <- survfit(surv_object ~ df$treatment)

# Plot the survival curves
ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE)
```

    Warning message:
    “[1m[22mUsing `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    [36mℹ[39m Please use `linewidth` instead.
    [36mℹ[39m The deprecated feature was likely used in the [34mggpubr[39m package.
      Please report the issue at [3m[34m<https://github.com/kassambara/ggpubr/issues>[39m[23m.”
    [1m[22mIgnoring unknown labels:
    [36m•[39m [32mcolour[39m : [34m"Strata"[39m



    
![png](analysis_files/analysis_4_1.png)
    



```R
# Run Cox regression to estimate the hazard ratio (HR) between treatment groups
# Check whether Cox regression analyses have met the assumption of proportional hazards
cox_model <- coxph(surv_object ~ df$treatment)
ph_test <- cox.zph(cox_model)
ph_pvalue <- ph_test$table[1, "p"]
if (ph_pvalue < 0.05) {
    cat("Warning: Proportional hazards assumption may be violated for treatment variable\n")
}
ph_pvalue
```


0.530258218382148



```R
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
cox_treatment_result
```


<table class="dataframe">
<caption>A data.frame: 1 × 4</caption>
<thead>
	<tr><th></th><th scope=col>HR</th><th scope=col>HR_lower</th><th scope=col>HR_upper</th><th scope=col>p_value</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>df$treatmentB</th><td>0.6925678</td><td>0.3622747</td><td>1.323996</td><td>0.2665201</td></tr>
</tbody>
</table>




```R
# Considering the sample size is 44, use coxphf to fit Cox Regression with Firth's Penalized Likelihood
penalized_fit <- coxphf::coxphf(Surv(time.pfs, event.pfs) ~ treatment, df)
beta_pena <- penalized_fit$coefficients[1]
se_pena <- sqrt(diag(penalized_fit$var))[1]
cox_treatment_pena <- data.frame(
    HR = exp(beta_pena),
    HR_lower = exp(beta_pena - 1.96 * se_pena),
    HR_upper = exp(beta_pena + 1.96 * se_pena),
    p_value = penalized_fit$prob[1],
    stringsAsFactors = FALSE
)
cox_treatment_pena
```


<table class="dataframe">
<caption>A data.frame: 1 × 4</caption>
<thead>
	<tr><th></th><th scope=col>HR</th><th scope=col>HR_lower</th><th scope=col>HR_upper</th><th scope=col>p_value</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>treatmentB</th><td>0.6900251</td><td>0.3610582</td><td>1.31872</td><td>0.2610211</td></tr>
</tbody>
</table>



## Answer to Question 1 based on the analysis above
In the two-group (A = Chemotherapy + placebo, B = Chemotherapy + atezolizumab (immunotherapy / PD-L1 inhibitor)) comparison of progression-free survival (PFS), the Cox model estimated Hazard Ration (HR) is 0.69 (95% CI 0.36–1.32; P-value = 0.27). The Cox regression analysis has met the assumption of proportional hazards. Considering the sample size (n=44), the small-sample bias was taken care by fitting a Firth-penalized Cox model, and it gave similar results. And the Kaplan-Meier log-rank test was P-value = 0.26.

The consistent results indicate that there is no statistically significant difference in PFS between the two treatment groups. Though the estimated hazard ratio as 0.69 implies a lower risk of progression with atezolizumab, but the wide 95% confidence interval includes 1, ranging from benefit (HR=0.36) to increase in risk (HR=1.32), which indicates uncertainty.

## Question 2: Please assess whether any of the numerical gene expression scores are associated with progression-free survival in this dataset.


```R
# List the column names of the 8 gene expressions
cols_to_exclude <- c("treatment", "time.pfs", "event.pfs", "id")
cols_exp <- setdiff(colnames(df), cols_to_exclude)
cols_exp
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'TIS'</li><li>'APM'</li><li>'Cytotoxic.Cells'</li><li>'Cell.Adhesion'</li><li>'Macrophages'</li><li>'HRD'</li><li>'PD.L2'</li><li>'Claudin.Low'</li></ol>




```R
# Check whether Cox regression analyses have met the assumption of proportional hazards
ph_assumptions <- vector("list", length(cols_exp))
for (i in seq_along(cols_exp)) {
    gene <- cols_exp[i]
    cox_formula <- as.formula(sprintf("Surv(time.pfs, event.pfs) ~ %s + strata(treatment)", gene))
    model <- coxph(cox_formula, data = df)
    ph_test <- cox.zph(model)
    ph_pvalue <- ph_test$table[1, "p"]
    if (ph_pvalue < 0.05) {
        cat("Warning: Proportional hazards assumption may be violated for gene", gene, "\n")
    }
    ph_assumptions[[i]] <- data.frame(
        gene = gene,
        p_value = ph_pvalue,
        stringsAsFactors = FALSE
    )
}
ph_assumptions <- bind_rows(ph_assumptions)
ph_assumptions
```


<table class="dataframe">
<caption>A data.frame: 8 × 2</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>p_value</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>TIS            </td><td>0.9690752</td></tr>
	<tr><td>APM            </td><td>0.6267775</td></tr>
	<tr><td>Cytotoxic.Cells</td><td>0.9436848</td></tr>
	<tr><td>Cell.Adhesion  </td><td>0.3037935</td></tr>
	<tr><td>Macrophages    </td><td>0.9034903</td></tr>
	<tr><td>HRD            </td><td>0.7543653</td></tr>
	<tr><td>PD.L2          </td><td>0.5904097</td></tr>
	<tr><td>Claudin.Low    </td><td>0.3680182</td></tr>
</tbody>
</table>




```R
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
cox_gene_results
```


<table class="dataframe">
<caption>A data.frame: 8 × 4</caption>
<thead>
	<tr><th></th><th scope=col>HR</th><th scope=col>HR_lower</th><th scope=col>HR_upper</th><th scope=col>p_value</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>TIS</th><td>0.7628480</td><td>0.5670273</td><td>1.0262946</td><td>0.07369308</td></tr>
	<tr><th scope=row>APM</th><td>0.7400094</td><td>0.5816535</td><td>0.9414778</td><td>0.01425094</td></tr>
	<tr><th scope=row>Cytotoxic.Cells</th><td>0.6986305</td><td>0.4960194</td><td>0.9840029</td><td>0.04014207</td></tr>
	<tr><th scope=row>Cell.Adhesion</th><td>1.0751782</td><td>0.9310604</td><td>1.2416038</td><td>0.32355132</td></tr>
	<tr><th scope=row>Macrophages</th><td>0.6662358</td><td>0.4045021</td><td>1.0973245</td><td>0.11066955</td></tr>
	<tr><th scope=row>HRD</th><td>1.0050754</td><td>0.6577921</td><td>1.5357078</td><td>0.98132629</td></tr>
	<tr><th scope=row>PD.L2</th><td>0.6466131</td><td>0.4341052</td><td>0.9631500</td><td>0.03197793</td></tr>
	<tr><th scope=row>Claudin.Low</th><td>0.8173325</td><td>0.6223500</td><td>1.0734030</td><td>0.14689231</td></tr>
</tbody>
</table>




```R
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
cox_gene_small_sample
```


<table class="dataframe">
<caption>A data.frame: 8 × 4</caption>
<thead>
	<tr><th></th><th scope=col>HR</th><th scope=col>HR_lower</th><th scope=col>HR_upper</th><th scope=col>p_value</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>TIS</th><td>0.7665290</td><td>0.5701434</td><td>1.0305596</td><td>0.07918711</td></tr>
	<tr><th scope=row>APM</th><td>0.7520218</td><td>0.5956806</td><td>0.9493959</td><td>0.01646902</td></tr>
	<tr><th scope=row>Cytotoxic.Cells</th><td>0.6988466</td><td>0.4976061</td><td>0.9814722</td><td>0.04124580</td></tr>
	<tr><th scope=row>Cell.Adhesion</th><td>1.0697955</td><td>0.9309081</td><td>1.2294042</td><td>0.30978378</td></tr>
	<tr><th scope=row>Macrophages</th><td>0.6551237</td><td>0.4007076</td><td>1.0710730</td><td>0.08984742</td></tr>
	<tr><th scope=row>HRD</th><td>0.9794952</td><td>0.6464268</td><td>1.4841755</td><td>0.92183398</td></tr>
	<tr><th scope=row>PD.L2</th><td>0.6286213</td><td>0.4235690</td><td>0.9329407</td><td>0.02452061</td></tr>
	<tr><th scope=row>Claudin.Low</th><td>0.8174394</td><td>0.6271234</td><td>1.0655116</td><td>0.11523445</td></tr>
</tbody>
</table>




```R
q_BH <- p.adjust(cox_gene_results$p_value, method = "BH")
q_BH
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>0.14738616614732</li><li>0.107045514436072</li><li>0.107045514436072</li><li>0.369772938962337</li><li>0.177071276074732</li><li>0.981326293252811</li><li>0.107045514436072</li><li>0.195856415196937</li></ol>



## Answer to Question 2 based on the analysis above
A Cox model was fit for each gene stratified by treatment (assumption of proportional hazards was met, and small-sample bias made no big difference and therefore ignored).

Three gene expression scores are significantly (P-value < 0.05) associated with PFS:
APM: HR = 0.75 (95% CI 0.58-0.94, P-value = 0.01)
Cytotoxic.Cells: HR = 0.70 (95% CI 0.50-0.98, P-value = 0.04)
PD.L2: HR = 0.65 (95% CI 0.43-0.96, P-value = 0.03).

And their hazard ratios and the 95% confidence intervals (HR,95%CI < 1) indicate that their higher expressions associate with better PFS.

Be noted: since multiple tests were run and the multiple P-values were adjusted by Benjamini-Hochberg, but then none of them passed the significance threshold 0.05.

## Question 3 a): Are any of the scores potential predictive biomarkers for the experimental treatment? 


```R
# PREDICTIVE BIOMARKER ANALYSIS

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
    cox_interaction <- coxph(cox_formula_interaction, data = df)

    # Check whether Cox regression analyses have met the assumption of proportional hazards
    ph_test_interaction <- cox.zph(cox_interaction)
    ph_test_interaction_coef <- ph_test_interaction$table
    ph_test_interaction_pvalue <- ph_test_interaction_coef[grep(":", rownames(ph_test_interaction_coef)), "p"]

    if (ph_test_interaction_pvalue < 0.05) {
        cat("Warning: Proportional hazards assumption may be violated for interaction term in gene", gene, "\n")
    }
    
    # Extract interaction p-value (key test for predictive biomarker)
    interaction_coef <- summary(cox_interaction)$coefficients
    interaction_pvalue <- interaction_coef[grep(":", rownames(interaction_coef)), "Pr(>|z|)"]

    # Calculate treatment effect (HR) in High/Low group
    # High expression group: Treatment B vs A
    df_high_group <- df[df[[paste0(gene, "_group")]] == "High", ]
    if(length(unique(df_high_group$treatment)) == 2) {  # both A and B are there
        cox_high <- coxph(Surv(time.pfs, event.pfs) ~ treatment, data = df_high_group)
        ph_test_high <- cox.zph(cox_high)
        ph_test_high_coef <- ph_test_high$table
        ph_test_high_pvalue <- ph_test_high_coef[1, "p"]
        if (ph_test_high_pvalue < 0.05) {
            cat("Warning: Proportional hazards assumption may be violated for high group in gene", gene, "\n")
        }
        hr_high <- exp(coef(cox_high))
        hr_high_ci <- exp(confint(cox_high))
        high_pvalue <- summary(cox_high)$coefficients[1, "Pr(>|z|)"]
    } else {
        high_pvalue <- NA
    }

    # Low expression group: Treatment B vs A
    df_low_group <- df[df[[paste0(gene, "_group")]] == "Low" & !is.na(df[[paste0(gene, "_group")]]), ]
    if(length(unique(df_low_group$treatment)) == 2) {
        cox_low <- coxph(Surv(time.pfs, event.pfs) ~ treatment, data = df_low_group)
        ph_test_low <- cox.zph(cox_low)
        ph_test_low_coef <- ph_test_low$table
        ph_test_low_pvalue <- ph_test_low_coef[1, "p"]
        if (ph_test_low_pvalue < 0.05) {
            cat("Warning: Proportional hazards assumption may be violated for low group in gene", gene, "\n")
        }
        hr_low <- exp(coef(cox_low))
        hr_low_ci <- exp(confint(cox_low))
        low_pvalue <- summary(cox_low)$coefficients[1, "Pr(>|z|)"]
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
biomarker_results

cox_pvalues <- bind_rows(cox_pvalues)
cox_pvalues

# plot the surv plot for each gene
for (i in seq_along(cols_exp)) {
    print(biomarker_plots[[i]])
}
```

    
     TIS median: 6.698 
    Group sizes: 22 22 
    
     APM median: 12.256 
    Group sizes: 22 22 
    
     Cytotoxic.Cells median: 4.054 
    Group sizes: 22 22 
    
     Cell.Adhesion median: 8.327 
    Group sizes: 22 22 
    
     Macrophages median: 6.599 
    Group sizes: 22 22 
    
     HRD median: 5.33 
    Group sizes: 22 22 
    
     PD.L2 median: 4.47 
    Group sizes: 22 22 
    
     Claudin.Low median: 0.587 
    Group sizes: 22 22 



<table class="dataframe">
<caption>A data.frame: 8 × 9</caption>
<thead>
	<tr><th></th><th scope=col>gene</th><th scope=col>hr_high</th><th scope=col>hr_high_lower</th><th scope=col>hr_high_upper</th><th scope=col>high_pvalue</th><th scope=col>hr_low</th><th scope=col>hr_low_lower</th><th scope=col>hr_low_upper</th><th scope=col>low_pvalue</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>treatmentB...1</th><td>TIS            </td><td>0.6723760</td><td>0.2376642</td><td>1.9022197</td><td>0.45440565</td><td>1.0434967</td><td>0.4223039</td><td>2.578440</td><td>0.92649947</td></tr>
	<tr><th scope=row>treatmentB...2</th><td>APM            </td><td>0.6919174</td><td>0.2474098</td><td>1.9350475</td><td>0.48275114</td><td>0.7849128</td><td>0.3027815</td><td>2.034761</td><td>0.61826654</td></tr>
	<tr><th scope=row>treatmentB...3</th><td>Cytotoxic.Cells</td><td>0.5930615</td><td>0.2216705</td><td>1.5866883</td><td>0.29809130</td><td>1.4419755</td><td>0.5774651</td><td>3.600725</td><td>0.43309189</td></tr>
	<tr><th scope=row>treatmentB...4</th><td>Cell.Adhesion  </td><td>1.5650404</td><td>0.5899939</td><td>4.1514858</td><td>0.36817947</td><td>0.3688455</td><td>0.1319942</td><td>1.030704</td><td>0.05713405</td></tr>
	<tr><th scope=row>treatmentB...5</th><td>Macrophages    </td><td>0.4107398</td><td>0.1558328</td><td>1.0826168</td><td>0.07195040</td><td>1.3493513</td><td>0.5265632</td><td>3.457798</td><td>0.53258338</td></tr>
	<tr><th scope=row>treatmentB...6</th><td>HRD            </td><td>0.5755754</td><td>0.2241522</td><td>1.4779556</td><td>0.25095137</td><td>0.6685407</td><td>0.2518553</td><td>1.774617</td><td>0.41885991</td></tr>
	<tr><th scope=row>treatmentB...7</th><td>PD.L2          </td><td>0.5391296</td><td>0.2104444</td><td>1.3811759</td><td>0.19804280</td><td>1.0683051</td><td>0.4179819</td><td>2.730443</td><td>0.89023744</td></tr>
	<tr><th scope=row>treatmentB...8</th><td>Claudin.Low    </td><td>0.3280494</td><td>0.1183411</td><td>0.9093744</td><td>0.03214688</td><td>1.7697791</td><td>0.7067351</td><td>4.431813</td><td>0.22289849</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 8 × 4</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>interaction_pvalue</th><th scope=col>high_group_pvalue</th><th scope=col>low_group_pvalue</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>TIS            </td><td>0.571131654</td><td>0.45440565</td><td>0.92649947</td></tr>
	<tr><td>APM            </td><td>0.725314287</td><td>0.48275114</td><td>0.61826654</td></tr>
	<tr><td>Cytotoxic.Cells</td><td>0.210064141</td><td>0.29809130</td><td>0.43309189</td></tr>
	<tr><td>Cell.Adhesion  </td><td>0.039144281</td><td>0.36817947</td><td>0.05713405</td></tr>
	<tr><td>Macrophages    </td><td>0.061758082</td><td>0.07195040</td><td>0.53258338</td></tr>
	<tr><td>HRD            </td><td>0.851462920</td><td>0.25095137</td><td>0.41885991</td></tr>
	<tr><td>PD.L2          </td><td>0.363157271</td><td>0.19804280</td><td>0.89023744</td></tr>
	<tr><td>Claudin.Low    </td><td>0.009921877</td><td>0.03214688</td><td>0.22289849</td></tr>
</tbody>
</table>



    [1m[22mIgnoring unknown labels:
    [36m•[39m [32mcolour[39m : [34m"TIS Expression"[39m
    [1m[22mIgnoring unknown labels:
    [36m•[39m [32mcolour[39m : [34m"APM Expression"[39m



    
![png](analysis_files/analysis_17_4.png)
    


    [1m[22mIgnoring unknown labels:
    [36m•[39m [32mcolour[39m : [34m"Cytotoxic.Cells Expression"[39m



    
![png](analysis_files/analysis_17_6.png)
    


    [1m[22mIgnoring unknown labels:
    [36m•[39m [32mcolour[39m : [34m"Cell.Adhesion Expression"[39m



    
![png](analysis_files/analysis_17_8.png)
    


    [1m[22mIgnoring unknown labels:
    [36m•[39m [32mcolour[39m : [34m"Macrophages Expression"[39m



    
![png](analysis_files/analysis_17_10.png)
    


    [1m[22mIgnoring unknown labels:
    [36m•[39m [32mcolour[39m : [34m"HRD Expression"[39m



    
![png](analysis_files/analysis_17_12.png)
    


    [1m[22mIgnoring unknown labels:
    [36m•[39m [32mcolour[39m : [34m"PD.L2 Expression"[39m



    
![png](analysis_files/analysis_17_14.png)
    


    [1m[22mIgnoring unknown labels:
    [36m•[39m [32mcolour[39m : [34m"Claudin.Low Expression"[39m



    
![png](analysis_files/analysis_17_16.png)
    



    
![png](analysis_files/analysis_17_17.png)
    



```R
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

df_forest
```


<table class="dataframe">
<caption>A data.frame: 16 × 5</caption>
<thead>
	<tr><th scope=col>subgroup</th><th scope=col>hr</th><th scope=col>lower</th><th scope=col>upper</th><th scope=col>pvalue</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>TIS High            </td><td>0.6723760</td><td>0.2376642</td><td>1.9022197</td><td>0.45440565</td></tr>
	<tr><td>TIS Low             </td><td>1.0434967</td><td>0.4223039</td><td>2.5784402</td><td>0.92649947</td></tr>
	<tr><td>APM High            </td><td>0.6919174</td><td>0.2474098</td><td>1.9350475</td><td>0.48275114</td></tr>
	<tr><td>APM Low             </td><td>0.7849128</td><td>0.3027815</td><td>2.0347610</td><td>0.61826654</td></tr>
	<tr><td>Cytotoxic.Cells High</td><td>0.5930615</td><td>0.2216705</td><td>1.5866883</td><td>0.29809130</td></tr>
	<tr><td>Cytotoxic.Cells Low </td><td>1.4419755</td><td>0.5774651</td><td>3.6007252</td><td>0.43309189</td></tr>
	<tr><td>Cell.Adhesion High  </td><td>1.5650404</td><td>0.5899939</td><td>4.1514858</td><td>0.36817947</td></tr>
	<tr><td>Cell.Adhesion Low   </td><td>0.3688455</td><td>0.1319942</td><td>1.0307039</td><td>0.05713405</td></tr>
	<tr><td>Macrophages High    </td><td>0.4107398</td><td>0.1558328</td><td>1.0826168</td><td>0.07195040</td></tr>
	<tr><td>Macrophages Low     </td><td>1.3493513</td><td>0.5265632</td><td>3.4577975</td><td>0.53258338</td></tr>
	<tr><td>HRD High            </td><td>0.5755754</td><td>0.2241522</td><td>1.4779556</td><td>0.25095137</td></tr>
	<tr><td>HRD Low             </td><td>0.6685407</td><td>0.2518553</td><td>1.7746167</td><td>0.41885991</td></tr>
	<tr><td>PD.L2 High          </td><td>0.5391296</td><td>0.2104444</td><td>1.3811759</td><td>0.19804280</td></tr>
	<tr><td>PD.L2 Low           </td><td>1.0683051</td><td>0.4179819</td><td>2.7304428</td><td>0.89023744</td></tr>
	<tr><td>Claudin.Low High    </td><td>0.3280494</td><td>0.1183411</td><td>0.9093744</td><td>0.03214688</td></tr>
	<tr><td>Claudin.Low Low     </td><td>1.7697791</td><td>0.7067351</td><td>4.4318135</td><td>0.22289849</td></tr>
</tbody>
</table>




```R
# Create a forest plot visualization
table_text <- cbind(
    c("Biomarker subgroup", df_forest$subgroup),
    c("HR (95% CI)", sprintf("%.2f (%.2f, %.2f)", df_forest$hr, df_forest$lower, df_forest$upper)),
    c("p-value", signif(df_forest$pvalue, 3))
)

options(repr.plot.width = 10, repr.plot.height = 8, repr.plot.res = 150)

forestplot(
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
```


    
![png](analysis_files/analysis_19_0.png)
    


## Answer to Question 3a based on the analysis above
For each gene, the patients were divided into two groups by the median of the gene expression scores. A Cox model with the interaction term of the gene subgroup and treatment was fit (assumption of proportional hazards was met, and small-sample bias made no big difference and therefore ignored) to study whether the gene is a predictive biomarker -- treatment effect differs between High and Low subgroups. A KM plot was also drawn with the interaction P-value indicated. And furthermore, a Cox model was also fit to study how the treatment affects PFS within each subgroup.

In the end, a forest plot collects the key information and visualizes them.

As we go through the analysis for all the genes, we can find that there are two gene expression scores potential predictive biomarkers for the treatment:
Claudin.Low: interaction P-value = 0.01, and in its High subgroup, treatment B (atezolizumab) has a lower risk of progression at HR = 0.33 (95% CI 0.12-0.91, P-value=0.03);
Cell.Adhesion: interaction P-value = 0.04, and in its Low subgroup, treatment B (atezolizumab) has a lower risk of progression at HR = 0.37 (95% CI 0.13-1.03, P-value=0.06).

And one more gene expression score with moderate potential:
Macrophages: interaction P-value = 0.06, and in its High subgroup, treatment B (atezolizumab) has a lower risk of progression at HR = 0.41 (95% CI 0.16-1.08, P-value=0.07).

## Question 3 b): While the median is a commonly used cutoff when dichotomizing a population by a numeric value, it is not always the best cutoff. Suggest strategies for optimizing the cutoff value of a predictive or prognostic biomarker. (A brief discussion of methods is sufficient, no need to perform the calculations).

## Answer to Question 3b
Though dichotomization simplifies the statistical analysis and can lead to a clear interpretation of the results as shown above, but it can also be problematic:
1. Information lost, so the statistical power is reduced, especially in our case, where the patient number is not high.
2. May increase the risk of false positive.
3. Individuals close to the cutoff but on opposite sides are characterized as very different rather than very similar.
4. Hide potential non-linearity in the relation between the variable and treatment.

Some strategies may be considered:
1. Retain the original continuous scores to do the analysis. Then no information is lost, but the fitting requirement may be too high, and a simple standard Cox model may not be enough.
2. Multiple categories based on prior knowledge, e.g. Gleason score in prostate cancer grading system. But a general definition for the gene expression scores' levels may be missing.
3. Sustainable Technology Promotion Platform (STEPP), sliding window to formulate subgroups with overlaps. STEPP requires no pre-defined split, so it can realize the "multiple categories" purpose above.
4. Machine learning or deep learning methods that may be able to discover non-linearity in the relation between gene expression scores and the treatment. They can address the potential issue mentioned in strategy 1.


```R

```
