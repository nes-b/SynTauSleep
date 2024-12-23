# ---
# Title: PSG in aSyn & Tau ####
#   date: "2024-07-25"
# contact: nils.briel@usz.ch
# name: Nils Briel
# ---

# 0 Environment setup ####
# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(readxl)
library(mice)
library(caret)
library(forestmodel)
library(ggpubr)
library(ggstatsplot)
library(gtsummary)
library(GGally)
library(car)
library(multcomp)
library(ggstatsplot) 
library(caret)
library(factoextra)
library(data.table)
library(ggfortify)
library(scales)

# 1 Load dataset ####
# Create sample data (replace with your actual data)
df <- read_excel("../data/psg_raw.xlsx", skip = 1, na = "NaN")
if (!file.exists("../results")) {dir.create('../results')}
setwd("../results")
if (!file.exists("plots")) {dir.create('plots')}
if (!file.exists("figures")) {dir.create('figures')}
  
# select sp and epi params
sel_epi <- c(4,7,8,12,38:40,1,2,52,13:21,23:27,35:37)
df <- df[,sel_epi]
sel_num <- c(2:5,10:24)
sel_fac <- which(!colnames(df) %in% colnames(df)[sel_num])
df[,sel_fac] <- lapply(df[,sel_fac], as.factor) %>% as.data.frame()
df[,sel_num] <- lapply(df[,sel_num], as.numeric) %>% as.data.frame()
x <- df$Dx_label 
df$Dx_label <- ifelse(df$Dx_label == "CBD", "CBS", as.character(df$Dx_label)) %>% as.factor()
df$Dx_label <- factor(df$Dx_label, levels = c("DLB",'PD','MSA','PSP','CBS'))
df$sex <- ifelse(df$sex == "0", "M","F") %>% as.factor()
colnames(df)[9] <- "Proteinopathy"
colnames(df)[25] <- "RWA" 
df$ID <- 1:nrow(df)

# select sleep params
sel_sp <- c(8,9,4,10:ncol(df))
df_sel <- df[,sel_sp]

# define color scales
Dx_cols <- c("darkred","#cc3300","#ff9933","darkblue","blue")
P_cols <- c("#cc3300","blue")

# 2 Descriptive  ####
# Basic stats
tbl_prot <- tbl_summary(
  df,
  include = c(age, disdur, sex, LED, sESS, TST_min, SPT, N2_latency_min, Sl_eff_perc, R_perc_spt, N1_perc_spt, N2_perc_spt, N3_perc_spt, WAKE_perc_spt, RBD, RWA, PLMS_index, Arousal_index, stridor, AHI, ODI, antidepressants, antipsychotics),
  label = list(age ~ "Age (y)", 
               disdur ~ "Disease Duration (m)", 
               sex ~ "Sex (M)", 
               LED ~ "LED (mg/d)",
               sESS ~ "ESS",
               TST_min ~ "TST (min)",
               SPT ~ "SPT (min)",
               N2_latency_min  ~ "Latency to N2 (min)",
               Sl_eff_perc ~ "%Sleep Efficiency",
               R_perc_spt ~ "% REM",
               N1_perc_spt ~ "% N1", 
               N2_perc_spt ~ "% N2", 
               N3_perc_spt ~ "% N3",
               WAKE_perc_spt ~ "% W",
               RBD ~ "RBD", 
               RWA ~ "RWA",
               PLMS_index ~ "PLMS Index", 
               Arousal_index ~"Arousal Index",
               stridor ~ "Stridor",
               AHI ~ "AHI",
               ODI ~ "ODI",
               antidepressants ~"Antidepressants use",
               antipsychotics ~ "Antipsychotics use"),
  by = c(Proteinopathy), 
  missing = "no" 
) %>%
  add_n() %>% 
  add_p() %>% 
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels()

# Across diagnosis labels:
tbl_dg <- tbl_summary(
  df,
  include = c(age, disdur, sex, LED, sESS, TST_min, SPT, N2_latency_min, Sl_eff_perc, R_perc_spt, N1_perc_spt, N2_perc_spt, N3_perc_spt, WAKE_perc_spt, RBD, RWA, PLMS_index, Arousal_index, stridor, AHI, ODI, antidepressants, antipsychotics),
  label = list(age ~ "Age (y)", 
               disdur ~ "Disease Duration (m)", 
               sex ~ "Sex (M)", 
               LED ~ "LED (mg/d)",
               sESS ~ "ESS",
               TST_min ~ "TST (min)",
               SPT ~ "SPT (min)",
               N2_latency_min  ~ "Latency to N2 (min)",
               Sl_eff_perc ~ "%Sleep Efficiency",
               R_perc_spt ~ "% REM",
               N1_perc_spt ~ "% N1", 
               N2_perc_spt ~ "% N2", 
               N3_perc_spt ~ "% N3",
               WAKE_perc_spt ~ "% W",
               RBD ~ "RBD", 
               RWA ~ "RWA",
               PLMS_index ~ "PLMS Index", 
               Arousal_index ~"Arousal Index",
               stridor ~ "Stridor",
               AHI ~ "AHI",
               ODI ~ "ODI",
               antidepressants ~"Antidepressants use",
               antipsychotics ~ "Antipsychotics use"),
  by = c(Dx_label), 
  missing = "no" 
) %>%
  add_n() %>% 
  add_p() %>% 
  modify_header(label = "**Variable**") %>%
  bold_labels()

gtsummary::tbl_merge(tbls = list(tbl_prot, tbl_dg), tab_spanner = c("Proteinopathy", "Clinical Diagnosis")) %>% 
  as_flex_table() %>% flextable::save_as_html(
    path = "tbl_baseline_characteristics.html")


# 3 Testing for Group variance ####
## 3.1 Group variance across proteinopathy ####
# Here, we use R's glm regression function that regresses out the variables "age","sex", "disdur", "antidepressants", "antipsychotics" and leaves the remaining estimate on differences in parameter values between "proteinopathies", and a variable-specific p-value.
# Analysis of covariance (ANCOVA) would be of choice here, because we want to control covariance with potential confounders and if we fulfilled assumptions of linearity, normal distribution and equal variances.
mod_ls <- list()
sum_ls <- list()
p_list <- list()
est_ls <- list()
tbl_ls <- list()
vars <- colnames(df_sel[,3:21])
vars_poisson <- c("stridor","RBD","RWA") # factor vars
vars_gaussian <- vars[!vars %in% vars_poisson] # continuous vars
vars_transf <- c("N1_perc_spt", "N2_latency_min","N3_perc_spt", "WAKE_perc_spt","WASO_min", "Sl_eff_perc","Arousal_index", 
 "PLMS_index","AHI","cAHI","ODI")

for(i in vars){
  subset_data <- df[, c(i, "Proteinopathy", "age","sex", "disdur", "antidepressants", "antipsychotics")] %>% 
    .[complete.cases(.),]
  vec <- as.numeric(as.matrix(subset_data[,c(i)]))
  
  # core code
  formula <- paste0(i," ~ Proteinopathy + sex + age + disdur + antidepressants + antipsychotics")
  if(i %in% vars_gaussian){
    subset_data[,c(1,3,5)] <- sapply(subset_data[,c(1,3,5)], as.numeric) 
    if(i %in% vars_transf){
      subset_data[,c(1)] <- log(subset_data[,1]+1)
    }
    subset_data[,c(1,3,5)] <- sapply(subset_data[,c(1,3,5)], function(x) rescale(x,to=c(-1,1)))
    model <- glm(formula, data = as.data.table(subset_data), family = 'gaussian')
  }else{
    subset_data[,c(1,3,5)] <- sapply(subset_data[,c(1,3,5)], as.numeric) 
    subset_data[,c(1,3,5)] <- sapply(subset_data[,c(1,3,5)], function(x) rescale(x,to=c(-1,1)))
    subset_data[,i] <- as.numeric(as.matrix(subset_data[,c(i)])) %>% as.factor()
    model <- glm(formula, data = as.data.table(subset_data), family = 'binomial')
    print(i)
  }
  mod_ls[[i]] <- model
  sum_ls[[i]] <- summary(model)
  p_list[[i]] <- summary(model)$coefficients[2,4]
  est_ls[[i]] <- summary(model)$coefficients[2,1]
}

# Create a data frame with the variable names and p-values
p_values_df <- data.frame(variable = names(p_list), est = unlist(est_ls), p_value = unlist(p_list))
p_values_df$p_value_adj <- p.adjust(p_values_df$p_value, method = "fdr")
p_values_df$variable <- reorder(p_values_df$variable, -p_values_df$p_value_adj)

# Plot the p-values
l1 <- ggplot(p_values_df, aes(y = variable, fill = est, x = -log(p_value_adj))) +
  geom_point(size = 2, shape=21, color="black") + viridis::scale_fill_viridis() +
  geom_vline(xintercept = -log(1e-1), linetype = "dashed", color = "red") +
  labs(y = "Variable", color = "Estimate / Beta") +
  theme_bw()
p_values_df <- arrange(p_values_df, p_value_adj)

pdf("plots/prot_lollipop.pdf",width = 4, height = 3)
plot(l1)
dev.off()


int <- subset(p_values_df, p_value_adj <0.05)$variable %>% as.character()
for(i in int){
  tbl_ls[[i]] <- tbl_regression(mod_ls[[i]], exponentiate = T,
                                estimate_fun = function(x) style_number(x, digits = 2))
}
m <- tbl_merge(
  tbls = tbl_ls,
  tab_spanner = c(names(tbl_ls))
)
m %>% as_flex_table() %>% flextable::save_as_html(
  path = "tbl_multivar_regression.html")

formula <- paste0(paste(vars, collapse = "+"),"~ Proteinopathy+ sex + age + disdur + antidepressants + antipsychotics")

## Forest plot
library(forestmodel)
pdf("plots/prot_forest_.pdf",width = 10, height = 22)
plot(forest_model(theme = theme_void(),model_list =  mod_ls[int]))
dev.off()
pdf("plots/prot_forest_for_poster.pdf",width = 7, height = 8)
plot(forest_model(theme = theme_void(),model_list =  mod_ls[int]))
dev.off()


## 3.2 Group variance across Diagnosis labels ####
# Here, we apply a similar regression model to the diagnosis labels.
# We correct for multiple testing using FDR method (less conservative).
# Eventually, we plot the most significant features' tables.
mod_ls <- list()
sum_ls <- list()
sum_tbls <- list()
p_list <- list()
est_ls <- list()
tbl_ls <- list()
f_ls <- list()
tuk_ls <- list()

for(i in vars){
  print(i)
  subset_data <- df[, c(i, "Dx_label", "age","sex", "disdur", "antidepressants", "antipsychotics")] %>% 
    .[complete.cases(.),]
  formula <- paste0(i," ~ Dx_label + age + sex + disdur + antidepressants + antipsychotics")
  if(i %in% vars_gaussian){
    subset_data[,c(1,3,5)] <- sapply(subset_data[,c(1,3,5)], as.numeric) 
    if(i %in% vars_transf){
      subset_data[,1] <- log(subset_data[,1]+1)
    }
    subset_data[,c(1,3,5)] <- sapply(subset_data[,c(1,3,5)], function(x) rescale(x,to=c(-1,1)))
    model <- glm(formula, data = subset_data, family = 'gaussian')
  }else{
    subset_data[,c(1,3,5)] <- sapply(subset_data[,c(1,3,5)], as.numeric) 
    subset_data[,c(1,3,5)] <- sapply(subset_data[,c(1,3,5)], function(x) rescale(x,to=c(-1,1)))
    subset_data[,i] <- as.numeric(as.matrix(subset_data[,c(i)])) %>% as.factor()
    model <- glm(formula, data = as.data.table(subset_data), family = 'binomial')
  }  
  mod_ls[[i]] <- model
  sum_ls[[i]] <- Anova(model, type = "III")
  p_list[[i]] <- Anova(model, type = "III")[1,3]
  est_ls[[i]] <- Anova(model, type = "III")[1,1]
}

# Create a data frame with the variable names and p-values
p_values_df <- data.frame(variable = names(p_list), est = unlist(est_ls), p_value = unlist(p_list))
p_values_df$p_value_adj <- p.adjust(p_values_df$p_value, method = "fdr")
p_values_df$variable <- reorder(p_values_df$variable, p_values_df$est)

# lollipop plots with feature effect estimates
l2 <- ggplot(p_values_df, aes(y = variable, fill = est, x = -log(p_value_adj))) +
  geom_point(size = 2, shape=21, color="black") + viridis::scale_fill_viridis() +
  geom_vline(xintercept = -log(1e-1), linetype = "dashed", color = "red") +
  labs(y = "Variable", color = "Estimate / Beta") +
  theme_bw()

pdf("plots/dx_labels_lollipop.pdf",width = 4, height = 3)
plot(l2)
dev.off()

# filter top significant features
p_values_df <- arrange(p_values_df, p_value_adj)
int <- subset(p_values_df, p_value_adj <0.05)$variable %>% as.character() 

write.csv(p_values_df, "../results/anova_p_values_df.csv")

### 3.2.1 Tukey post-hoc test ####
# perform Tukey post-hoc test to find single group differences
for(k in int){
  model <- mod_ls[[k]]
  tuk_ls[[k]] <- tbl_regression(summary(glht(model, linfct = mcp("Dx_label" = "Tukey"))), exponentiate = T)
}

tbl_merge(
  tbls = tuk_ls,
  tab_spanner = names(tuk_ls)
)%>% as_flex_table() %>% flextable::save_as_html(
  path = "tbl_diagnoses.html")

## 3.3 Group-wise comparisons  ####
# ... of 3 significant parameters between diagnosis and proteinopathy groups.
# Y depicts parameters corrected for "sex + disdur + antidepressants + antipsychotics", as included in the ANCOVA model.

### 3.3.1 Box-/Violin Plots ####
int <- c("WAKE_perc_spt","N2_perc_spt", "R_perc_spt")
pl_list_1 <- list()
for(i in 1:3){
  pl_list_1[[i]] <- ggbetweenstats(df, x="Dx_label",
                                   y= !!int[i],
                                   type = "nonparametric",
                                   plot.type = "boxviolin",
                                   centrality.plotting = F,
                                   title = int[i],
                                   pairwise.comparisons = T, 
                                   p.adjust.method = "fdr") + 
    theme_bw()+ scale_color_manual(values = Dx_cols)   
}

pl_list_2 <- list()
for(i in 1:3){
  pl_list_2[[i]] <- ggbetweenstats(df, x="Proteinopathy",
                                   y= !!int[i],
                                   type = "nonparametric",
                                   plot.type = "boxviolin",
                                   centrality.plottin = F,
                                   title = int[i],
                                   p.adjust.method = "fdr") + 
    theme_bw()+ scale_color_manual(values = P_cols)   
}

p1 <- ggarrange(plotlist = pl_list_1, ncol=1, common.legend = TRUE, legend="bottom") 
p2 <- ggarrange(plotlist = pl_list_2, ncol=1, common.legend = TRUE, legend="bottom")

cairo_pdf("plots/violin_combined.pdf",width = 10, height = 12)
plot(ggarrange(plotlist = list(p2,p1), ncol=2, common.legend = TRUE, legend="bottom", widths = c(0.7,1)))
dev.off()

## 3.4 Radar plots ####
# commented out because ggradar only running on my linux.
### 3.4.1 Sleep features across diagnosis and proteinopathy labels ####
# devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
library(ggradar)
df_radar <- df[complete.cases(df), c("Dx_label","Sl_eff_perc", "RBD","R_perc_spt", "N1_perc_spt", "N2_perc_spt", "N3_perc_spt","WAKE_perc_spt", "TST_min")]
df_radar[,2:ncol(df_radar)] <- lapply(df_radar[,2:ncol(df_radar)], as.numeric) %>% as.data.frame()
df_radar$RBD <- rescale(df_radar$RBD, to=c(0, 1))
means <- aggregate(df_radar[,2:ncol(df_radar)],
                     by = list(Dx_label = df_radar$Dx_label),
                     FUN = mean) %>%
  rbind(.,as.numeric(c("NU", rep(0,ncol(df_radar)-1)))) %>%
  mutate_each(funs(rescale(., to = c(0, 1))), -Dx_label)
cairo_pdf("plots/ggradar_dx.pdf",width = 7, height = 5)
plot(ggradar(means,
             values.radar = c("0%","50%","100%"),
             grid.min = 0,
             grid.mid = max(means[,2:ncol(means)])/2,
             grid.max = max(means[,2:ncol(means)]),
             base.size = 10,
             font.radar = "arial",
             axis.labels = c("SE\n(% range)","RBD","REM\n(% spt, range)", "N1\n(% spt, range)", "N2\n(% spt, range)", "N3\n(% spt, range)", "WAKE\n(% spt, range)", "TST\n(min, % range)"),
             group.point.size = 2,
             group.line.width = 1,
             axis.label.size	=4,
             grid.line.width = 0.5,
             gridline.min.colour = "grey",
             gridline.mid.colour = "grey42",
             gridline.max.colour = "grey35",
) + theme_void() +
  scale_color_manual(values=Dx_cols))
dev.off()

### 3.4.2 Sleep features across proteinopathy labels ####
df_radar <- df[complete.cases(df), c("Proteinopathy","Sl_eff_perc", "RBD","R_perc_spt", "N1_perc_spt", "N2_perc_spt", "N3_perc_spt","WAKE_perc_spt", "TST_min")]
df_radar[,2:ncol(df_radar)] <- lapply(df_radar[,2:ncol(df_radar)], as.numeric) %>% as.data.frame()
df_radar$RBD <- rescale(df_radar$RBD, to=c(0, 1))

means <- aggregate(df_radar[,2:ncol(df_radar)],
                   by = list(Proteinopathy = df_radar$Proteinopathy),
                   FUN = mean) %>%
  rbind(.,as.numeric(c("NU", rep(0,ncol(df_radar)-1)))) %>%
  mutate_each(funs(rescale(., to = c(0, 1))), -Proteinopathy)

cairo_pdf("plots/ggradar_prot.pdf",width = 7, height = 5)
plot(ggradar(means[-3,],
             values.radar = c("0%","50%","100%"),
             grid.min = 0,
             grid.mid = max(means[,2:ncol(means)])/2,
             grid.max = max(means[,2:ncol(means)]),
             base.size = 10,
             font.radar = "arial",
             axis.labels = c("SE\n(% range)","RBD","REM\n(% spt, range)", "N1\n(% spt, range)", "N2\n(% spt, range)", "N3\n(% spt, range)", "WAKE\n(% spt, range)", "TST\n(min, % range)"),
             group.point.size = 2,
             group.line.width = 1,
             axis.label.size	=4,
             grid.line.width = 0.5,
             gridline.min.colour = "grey",
             gridline.mid.colour = "grey42",
             gridline.max.colour = "grey35",
) + theme_void() +
  scale_color_manual(values=P_cols))
dev.off()


# 4 LDA ####
df_s <- df_sel %>%.[,3:ncol(.)] %>% na.omit()
lib <- df[,c(1:9,28)]
lib <- lib[which(lib$ID %in% df_s$ID),]
rownames(lib) <- lib$ID
df_s <- tibble::column_to_rownames(df_s, var = 'ID')
vec_s <- data.frame(NP = as.numeric(as.factor(lib$Proteinopathy))) %>% as.matrix()
df_s <- as.data.frame(lapply(df_s[,1:ncol(df_s)], as.numeric))
scaled_data <- scale(df_s)

## 4.1 Sleep only ####
library(caret)
library(MASS)
vec_dx <- data.frame(NP = as.numeric(as.factor(lib$Dx_label))) %>% as.matrix()
theme_set(theme_bw())
set.seed(123)
index <- createDataPartition(vec_dx, p = 0.8, list = FALSE)
lda_df_tr <- scaled_data[index, ]
lda_df_te <- scaled_data[-index, ]

# transform
preproc.param <- lda_df_tr %>% preProcess(method = c("center", "scale"))
train.transformed <- preproc.param %>% predict(lda_df_tr) %>% cbind(NP = vec_dx[index], .) %>% as.data.frame()
test.transformed <- preproc.param %>% predict(lda_df_te) %>% cbind(NP = vec_dx[-index], .) %>% as.data.frame()

# Fit the model
model_s <- lda(NP~., data = train.transformed)
predictions <- predict(model_s, test.transformed)
saveRDS(confusionMatrix(table(list(predicted=predictions$class, observed=test.transformed$NP))), "confmat_lda_sleep.Rds")

# plot
lda.data <- cbind(train.transformed, predict(model_s)$x)
lda.data <- cbind(lda.data, Dx_label=lib$Dx_label[index], Proteinopathy=lib$Proteinopathy[index])

g_style <- theme_test() + 
  theme(legend.position = "bottom")
l1 <- ggplot(lda.data, aes(LD1,LD2,fill = Dx_label), color = "grey42") +
  coord_equal(ratio = 1)+ 
  stat_density2d(geom="density_2d", aes(color = Proteinopathy, alpha=..level..),
                 size=1, contour=T, contour_var = "count") + scale_color_manual(values = P_cols) +
  geom_point(shape = 21, size = 3) + 
  scale_fill_manual(values = Dx_cols) + g_style

l2 <- ggplot(lda.data, aes(LD1,LD2,fill = Proteinopathy), color = "grey42") + 
  scale_fill_manual(values = P_cols) +
  coord_equal(ratio = 1)+ 
  geom_point(shape = 21, size = 3) + 
  scale_fill_manual(values = P_cols) + g_style

## 4.2 Sleep + Demo ####
# without LED
set.seed(123)
df_d <- df[,c(ncol(df),8,9,1:7,10:(ncol(df)-1))] %>% na.omit() %>% dplyr::select(-"LED")
df_d <- tibble::column_to_rownames(df_d, var = 'ID')
vec_d <- data.frame(NP = as.numeric(as.factor(lib$Proteinopathy))) %>% as.matrix()
df_d <- as.data.frame(lapply(df_d[,3:ncol(df_d)], as.numeric))
scaled_data_d <- scale(df_d)

set.seed(123)
lda_df_d_tr <- scaled_data_d[index, ]
lda_df_d_te <- scaled_data_d[-index, ]

#transform
preproc.param_d <- lda_df_d_tr %>% preProcess(method = c("center", "scale"))
train.transformed_d <- preproc.param_d %>% predict(lda_df_d_tr) %>% cbind(NP = vec_dx[index], .) %>% as.data.frame()
test.transformed_d <- preproc.param_d %>% predict(lda_df_d_te) %>% cbind(NP = vec_dx[-index], .) %>% as.data.frame()

# Fit the model
model_d <- lda(NP~., data = train.transformed_d)
predictions_d <- predict(model_d, test.transformed_d)
saveRDS(confusionMatrix(table(list(predicted=predictions_d$class, observed=test.transformed_d$NP))), "confmat_lda_sleep_demo.Rds")

# plot
lda.data_d <- cbind(train.transformed_d, predict(model_d)$x)
lda.data_d <- cbind(lda.data_d, Dx_label=lib$Dx_label[index], Proteinopathy=lib$Proteinopathy[index])
lda.data_d$Dx_label <- factor(lda.data_d$Dx_label, levels = c("DLB",'PD','MSA','PSP','CBS'))

l3 <- ggplot(lda.data_d, aes(LD1,LD2,fill = Dx_label), color = "grey42") +
  coord_equal(ratio = 1)+ 
  stat_density2d(geom="density_2d", aes(color = Proteinopathy, alpha=..level..),
                 size=1, contour=T, contour_var = "count") + scale_color_manual(values = P_cols) +
  geom_point(shape = 21, size = 3) + 
  scale_fill_manual(values = Dx_cols) + g_style

l4 <- ggplot(lda.data_d, aes(LD1,LD2,fill = Proteinopathy), color = "grey42") +
  coord_equal(ratio = 1)+ 
  geom_point(shape = 21, size = 3) + 
  scale_fill_manual(values = P_cols) + g_style


cairo_pdf("plots/lda_1.pdf",width = 7, height = 4)
plot(l1|l2)
dev.off()

cairo_pdf("plots/lda_2.pdf",width = 7, height = 4)
plot(l3|l4)
dev.off()

## 4.3 Discrimination Proteinopathy ####
print("Sleep only")
predictions <- predict(model_s, test.transformed, method = "debiased")
test.transformed$Prot_true <- ifelse(test.transformed$NP >=4, "Tau","Syn") %>% as.factor()
test.transformed$Prot_pred <- ifelse(as.numeric(predictions$class) >=4, "Tau","Syn") %>% as.factor()
conf_matrix <- confusionMatrix(table(list(predicted=test.transformed$Prot_pred, observed=test.transformed$Prot_true)), mode = "sens_spec")
saveRDS(conf_matrix, "../results/confmat_lda_sleep_prot.Rds")

print("Sleep + Demo")
predictions_d <- predict(model_d, test.transformed_d, method = "debiased")
test.transformed_d$Prot_true <- ifelse(test.transformed_d$NP >=4, "Tau","Syn") %>% as.factor()
test.transformed_d$Prot_pred <- ifelse(as.numeric(predictions_d$class) >=4, "Tau","Syn") %>% as.factor()
conf_matrix_d <- confusionMatrix(table(list(predicted=test.transformed_d$Prot_pred, observed=test.transformed_d$Prot_true)), mode = "sens_spec")
saveRDS((conf_matrix_d), "../results/confmat_lda_sleep_demo_prot.Rds")

## LD1 plots
library(ggthemes)
g_style <-  theme_few() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) 

# Sleep model
LD1_s <- data.frame(coef(model_s),Var = rownames(coef(model_s))) %>% 
  arrange(desc(LD1), .by_group = TRUE)
ld1_p1 <- ggplot(LD1_s, aes(x=reorder(Var, LD1), y = 0, fill = LD1, size = abs(LD1))) + 
  geom_point(shape = 22) + coord_equal(ratio = 10) + scale_fill_distiller(palette = "RdBu", limits = c(-3,3)) + scale_x_discrete(guide = guide_axis(angle = 45), position = "top") + g_style
# Demo model
LD1_d <- data.frame(coef(model_d),Var = rownames(coef(model_d))) %>% 
  arrange(desc(LD1), .by_group = TRUE)
ld1_p2 <- ggplot(LD1_d, aes(x=reorder(Var, LD1), y = 0, fill = LD1, size = abs(LD1))) + 
  geom_point(shape = 22) + coord_equal(ratio = 10) + scale_fill_distiller(palette = "RdBu", limits = c(-1.5,1.5)) + scale_x_discrete(guide = guide_axis(angle = 45), position = "top") + g_style

cairo_pdf("plots/ld1_p.pdf",width = 7, height = 7)
plot(ggarrange(plotlist = list(ld1_p1,ld1_p2), ncol=1))
dev.off()

## LD2 plots
g_style <- theme_few() + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank()) 
LD2_s <- data.frame(coef(model_s),Var = rownames(coef(model_s))) %>% 
  arrange(desc(LD1), .by_group = TRUE)
ld2_p1 <- ggplot(LD2_s, aes(y=reorder(Var, LD2), x=0, fill = LD2, size = abs(LD2))) + geom_point(shape = 22) + scale_fill_distiller(palette = "RdBu") + coord_equal(ratio = 0.05)+ scale_y_discrete(guide = guide_axis(position = "right"))  + g_style
LD2_d <- data.frame(coef(model_d),Var = rownames(coef(model_d))) %>% 
  arrange(desc(LD2), .by_group = TRUE)
ld2_p2 <- ggplot(LD2_d, aes(y=reorder(Var, LD2), x=0, fill = LD2, size = abs(LD2))) + geom_point(shape = 22)+ scale_fill_distiller(palette = "RdBu", limits = c(-2,2)) + coord_equal(ratio = 0.05)+ scale_y_discrete(guide = guide_axis(position = "right"))  + g_style

cairo_pdf("plots/ld2_p.pdf",width = 7, height = 5)
plot(ggarrange(plotlist = list(ld2_p1,ld2_p2), ncol=2))
dev.off()

set.seed(123)
index <- createDataPartition(vec_d, p = 0.8, list = FALSE)
lda_df_tr <- scaled_data[index, ]
lda_df_te <- scaled_data[-index, ]

#transform
preproc.param <- lda_df_tr %>% preProcess(method = c("center", "scale"))
train.transformed <- preproc.param %>% predict(lda_df_tr) %>% cbind(NP = vec_d[index], .) %>% as.data.frame()
test.transformed <- preproc.param %>% predict(lda_df_te) %>% cbind(NP = vec_d[-index], .) %>% as.data.frame()

# Fit the model
model_s_p <- lda(NP~., data = train.transformed)
predictions <- model_s_p %>% predict(test.transformed)
confusionMatrix(table(list(predicted=predictions$class, observed=test.transformed$NP)))
plot(model_s_p, col=test.transformed$NP)

# plot
lda.data <- cbind(train.transformed, predict(model_s_p)$x)
lda.data <- cbind(lda.data, Dx_label=lib$Dx_label[index], Proteinopathy=lib$Proteinopathy[index])
l1 <- ggplot(lda.data, aes(LD1,LD2,fill = Dx_label), color = "grey42") +
  geom_point(shape = 21, size = 3) + scale_fill_manual(values = Dx_cols) +
  coord_equal(ratio = 1)+ 
  theme_bw()+ 
  theme(legend.position = "bottom")

ldahist(data = predictions$x[,], g=test.transformed$NP)

# 5 Parsimonous Models ####
## 5.1 Resampling and model training ####
library(pROC)
library(cvms)
df[,sel_num] <- sapply(df[,sel_num], as.numeric) 
df[,sel_num] <- sapply(df[,sel_num], function(x) rescale(x,to=c(0,1)))

#discussion from 20240724: invert comparison or case/control relation
df$Proteinopathy <- factor(df$Proteinopathy, levels = c("Syn","Tau"))

set.seed(123) # Set the random seed for reproducibility
index <- createDataPartition(df$Proteinopathy, p = 0.70, list = FALSE)
resamp <- createResample(df$Proteinopathy,times = 500)

# Create a formula for the logistic regression model
f_rem <- formula("Proteinopathy ~ RBD + RWA")
f_sl <- formula("Proteinopathy ~ RBD + RWA + R_perc_spt + Sl_eff_perc + WAKE_perc_spt + Arousal_index + N2_perc_spt + TST_min")
f_sl_dem <- formula("Proteinopathy ~ RBD + RWA + R_perc_spt + Sl_eff_perc + 
                    WAKE_perc_spt + Arousal_index + N2_perc_spt + TST_min + 
                    age  + sex + disdur")

smpl_ls <- lapply(1:length(resamp), function(i){
  train_data <- df[unique(resamp[[i]]), ]
  test_data <- df[-unique(resamp[[i]]), ]
  
  # Fit the logistic regression models
  m_rem <- glm(f_rem, data = train_data, family = binomial)
  m_sl <- glm(f_sl, data = train_data, family = binomial)
  m_sl_dem <- glm(f_sl_dem, data = train_data, family = binomial)
  
  # Calculate the AUC-ROC for each model
  auc_rem <- roc(test_data$Proteinopathy, predict(m_rem, newdata= test_data, type = "response"))
  auc_sl <- roc(test_data$Proteinopathy, predict(m_sl, newdata= test_data, type = "response"))
  auc_sl_dem <- roc(test_data$Proteinopathy, predict(m_sl_dem, newdata= test_data, type = "response"))
  
  data_rem <- data.frame(
    "target" = as.character(rescale(as.numeric(test_data$Proteinopathy),to=c(0,1))),
    "prediction" = as.character(round(predict(m_rem, newdata= test_data, type = "response"),
                                      digits = 0)), stringsAsFactors = FALSE) %>% drop_na() %>%
    evaluate(data = ., target_col = "target",  prediction_cols = "prediction",  type = 'binomial') %>%
    .[,1:8] %>% t() %>% round(., 3) %>% as.data.frame(.) %>% 
    mutate(mod = "RBD",metric = rownames(.),  sampling_no = i)
  data_sl <- data.frame(
    "target" = as.character(rescale(as.numeric(test_data$Proteinopathy),to=c(0,1))),
    "prediction" = as.character(round(predict(m_sl, newdata= test_data, type = "response"),
                                      digits = 0)), stringsAsFactors = FALSE) %>% drop_na() %>%
    evaluate(data = ., target_col = "target",  prediction_cols = "prediction",  type = 'binomial') %>%
    .[,1:8] %>% t() %>% round(., 3) %>% as.data.frame(.) %>% 
    mutate(mod = "RBD+Sleep",metric = rownames(.),  sampling_no = i)
  data_sl_dem <- data.frame(
    "target" = as.character(rescale(as.numeric(test_data$Proteinopathy),to=c(0,1))),
    "prediction" = as.character(round(predict(m_sl_dem, newdata= test_data, type = "response"),
                                      digits = 0)), stringsAsFactors = FALSE) %>% drop_na() %>%
    evaluate(data = ., target_col = "target",  prediction_cols = "prediction",  type = 'binomial') %>%
    .[,1:8] %>% t() %>% round(., 3) %>% as.data.frame(.)%>%
    mutate(mod = "RBD+Sleep+Demo",metric = rownames(.),  sampling_no = i)
  
  data_rem$V1[8] = auc_rem$auc
  data_sl$V1[8] = auc_sl$auc
  data_sl_dem$V1[8] = auc_sl_dem$auc
    
  df_temp <- rbind(data_rem,data_sl,data_sl_dem)
  colnames(df_temp)[1] <- 'Value'
  return(df_temp)
})

dat <- do.call(rbind,smpl_ls) 
dat$mod <- factor(dat$mod, levels = c("RBD",'RBD+Sleep','RBD+Sleep+Demo'
                                      ))
pdf("plots/auc_violins.pdf", width = 6.5, height = 6) 
ggbetweenstats(subset(dat,metric == "Balanced Accuracy"), x="mod",y= "Value", type = "parametric", 
               plot.type = "boxviolin",centrality.plottin = F, 
               title = "PSG Models: AUC Tau vs. Syn", p.adjust.method = "fdr",results.subtitle = F) + 
  theme_test() + scale_color_manual(values=c("grey","grey30","grey10","grey10"))
dev.off()

## 5.2 Classification metrics ####
set.seed(42)
dat_1 <- aggregate(Value ~ metric, data =  subset(dat, mod == "RBD"), mean, na.rm = TRUE) %>% tibble::column_to_rownames("metric") %>% mutate(Value = round(Value,2))
dat_2 <- aggregate(Value ~ metric, data =  subset(dat, mod == "RBD+Sleep"), mean, na.rm = TRUE)%>% tibble::column_to_rownames("metric") %>% mutate(Value = round(Value,2))
dat_3 <- aggregate(Value ~ metric, data =  subset(dat, mod == "RBD+Sleep+Demo"), mean, na.rm = TRUE)%>% tibble::column_to_rownames("metric") %>% mutate(Value = round(Value,2))

pdf("plots/confmat.pdf") 
gridExtra::grid.table(do.call(cbind, list(dat_1,dat_2,dat_3)) %>% 
                        `colnames<-`(c("RBD", "RBD+Sleep","RBD+Sleep+Demo"))
)
dev.off()

## 5.3 ROC ####
train_data <- df[unique(resamp[[250]]), ]
test_data <- df[-unique(resamp[[250]]), ]

# Fit the logistic regression models
m_rem <- glm(f_rem, data = train_data, family = binomial)
m_sl <- glm(f_sl, data = train_data, family = binomial)
m_dem <- glm(f_dem, data = train_data, family = binomial)
m_sl_dem <- glm(f_dem, data = train_data, family = binomial)
auc_rem <- roc(test_data$Proteinopathy, predict(m_rem, newdata= test_data, type = "response"))
auc_sl <- roc(test_data$Proteinopathy, predict(m_sl, newdata= test_data, type = "response"))
auc_dem <- roc(test_data$Proteinopathy, predict(m_dem, newdata= test_data, type = "response"))
auc_sl_dem <- roc(test_data$Proteinopathy, predict(m_sl_dem, newdata= test_data, type = "response"))

# Plot the ROC curves
cairo_pdf("plots/ROC.pdf",width = 6, height = 6)
plot(auc_rem, col = "grey", main = "ROC Curves", xlab = "False Positive Rate", ylab = "True Positive Rate")
plot(auc_sl, col = "grey20", add = TRUE)
plot(auc_sl_dem, col = "grey10", add = TRUE)
auc_rem$auc
auc_sl$auc
auc_sl_dem$auc
# Add a legend
legend("bottomright", legend = c("RBD", "RBD+Sleep", "RBD+Sleep+Demo"), col = c("grey", "grey20", "grey10"), lty = 1)
dev.off()


# 6. Differentiate MSA from other Syn
df <- subset(df, Proteinopathy == "Syn")
df$MSA_vs_LB <- ifelse(df$Dx_label == "MSA", "MSA", "LBD")
df$MSA_vs_LB <- factor(df$MSA_vs_LB, levels = c("LBD", "MSA"))

f_msa <- formula("MSA_vs_LB ~ RBD + RWA + R_perc_spt + Sl_eff_perc + WAKE_perc_spt + Arousal_index + N2_perc_spt + TST_min")

set.seed(123)
resamp_msa <- createResample(df$MSA_vs_LB, times = 500)

smpl_ls_msa <- lapply(1:length(resamp_msa), function(i) {
  train_data <- df[unique(resamp_msa[[i]]), ]
  test_data <- df[-unique(resamp_msa[[i]]), ]
  
  # Fit the logistic regression model
  m_msa <- glm(f_msa, data = train_data, family = binomial)
  
  # Calculate the AUC-ROC
  auc_msa <- roc(test_data$MSA_vs_LB, predict(m_msa, newdata = test_data, type = "response"))
  
  data_msa <- data.frame(
    "target" = as.character(rescale(as.numeric(test_data$MSA_vs_LB), to = c(0,1))),
    "prediction" = as.character(round(predict(m_msa, newdata = test_data, type = "response"), digits = 0)),
    stringsAsFactors = FALSE
  ) %>%
    drop_na() %>%
    evaluate(data = ., target_col = "target", prediction_cols = "prediction", type = 'binomial') %>%
    .[,1:8] %>%
    t() %>%
    round(., 3) %>%
    as.data.frame(.) %>%
    mutate(mod = "MSA_vs_LB", metric = rownames(.), sampling_no = i)
  
  data_msa$V1[8] = auc_msa$auc
  
  colnames(data_msa)[1] <- 'Value'
  return(data_msa)
})

dat_msa <- do.call(rbind, smpl_ls_msa)

## Plotting

pdf("plots/auc_violins_msa.pdf", width = 6.5, height = 6)
ggbetweenstats(
  subset(dat_msa, metric == "Balanced Accuracy"),
  x = "mod",
  y = "Value",
  type = "parametric",
  plot.type = "boxviolin",
  centrality.plotting = FALSE,
  title = "PSG Model: AUC MSA vs Other Synucleinopathies",
  results.subtitle = FALSE
) +
  theme_test() +
  scale_color_manual(values = c("grey10"))
dev.off()

dat_msa_metrics <- aggregate(Value ~ metric, data = dat_msa, median, na.rm = TRUE) %>%
  tibble::column_to_rownames("metric") %>%
  mutate(Value = round(Value, 2))

pdf("plots/confmat_msa.pdf")
gridExtra::grid.table(dat_msa_metrics %>% `colnames<-`(c("MSA vs LBD")))
dev.off()

# save full model
m_msa <- glm(f_msa, data = df, family = binomial)
saveRDS(m_msa, "../results/logmod_msa.Rds")

sessionInfo()
