library(readxl)
library(dplyr)
library(stringr)
library(gtsummary)
library(gt)
library(survminer)
setwd("/Users/micoria/Documents/work/MedSci/12.17")
Basic_data <- read_excel("/Users/micoria/Documents/work/MedSci/12.17/78例病人完整基础资料(1).xlsx",col_names = FALSE)
colnames(Basic_data) <- Basic_data[2, ]
head(Basic_data)
Basic_data <- Basic_data[-c(1,2), ]
rownames(Basic_data) <- NULL
head(Basic_data)
Basic_data$年龄<-as.numeric(Basic_data$年龄)
summary(Basic_data$年龄)
Basic_data <- Basic_data %>%
  mutate(
    年龄组 = cut(年龄, 
              breaks = c(0, 60, 100),  
              labels = c("<=60",">60"),
              right = TRUE)  
  )
convert_to_halfwidth <- function(x) {
  x <- chartr("（）", "()", x) 
  x <- gsub("　", " ", x)      
  return(x)
}
Basic_data$肿瘤大小 <- convert_to_halfwidth(Basic_data$肿瘤大小)
matches <- grepl("\\)\\s*(\\d+)", Basic_data$肿瘤大小)
print(matches)
Basic_data$size <- str_match(Basic_data$肿瘤大小, "\\)\\s*(\\d+)")[, 2]
print(Basic_data$size)
Basic_data$size[c(2, 4, 9, 13, 16, 20, 24, 25, 26, 30,49,57,61 )] <- c("2", "3", "2", "6", "2", "3","2", "1","1","1","2","2","1" )
Basic_data$size_continuously<-Basic_data$size
Basic_data <- Basic_data[-c(79:84), ]
Basic_data$pN <- substr(Basic_data$淋巴结, 1, 2)
Basic_data$pN[Basic_data$pN == "Nx"] <- "N0"
Basic_data$pN[c(30,36,37)] <- c("N0","N2","N0")
print(Basic_data$pN)
Basic_data$pT <-substr(Basic_data$TNM,1,2)
print(Basic_data$pT)
Basic_data$HER2_status<-Basic_data$HER2阳性与否
Basic_data$HER2_status[Basic_data$HER2_status == "1"]<-"阳性"
Basic_data$HER2_status[Basic_data$HER2_status == "2"]<-"阴性"
Basic_data <- Basic_data %>% rename(age = 年龄组, tumor_size = size, histological_grade = 病理分级,TMN_stage = 分期)
Basic_data$TMN_stage <- gsub("[ABC]", "", Basic_data$TMN_stage)
Basic_data$tumor_size <- as.numeric(Basic_data$tumor_size)
Basic_data <- Basic_data %>%
  mutate(
    tumor_size = cut(tumor_size, 
              breaks = c(0, 3, 10),  
              labels = c("<=3",">3"),
              right = TRUE)  
 )
pacman::p_load(
  rio,          # File import
  here,         # File locator
  skimr,        # get overview of data
  tidyverse,    # data management + ggplot2 graphics 
  gtsummary,    # summary statistics and tests
  rstatix,      # summary statistics and statistical tests
  janitor,      # adding totals and percents to tables
  scales,       # easily convert proportions to percents  
  flextable     # converting tables to pretty images
)
Basic_data$protein<- Basic_data$双变性量化打分
Basic_data$protein<- as.numeric(Basic_data$protein)
Basic_data <- Basic_data %>%
  mutate(
    RNF114 = cut(protein, 
              breaks = c(0, 4, 100),  
              labels = c("低表达(%)","高表达组(%)"),
              right = FALSE)  
  )
#for making a table for descrptive statitics
require(readr)
require(dplyr)
require(gtsummary)
require(gt)
require(webshot2)
require(htmlwidgets)
Basic_data$tumor_size <- as.factor(Basic_data$tumor_size)
Basic_data$pT <- as.factor(Basic_data$pT)
Basic_data$pN <- as.factor(Basic_data$pN)
Basic_data$histological_grade<- as.factor(Basic_data$histological_grade)
Basic_data$TMN_stage<- as.factor(Basic_data$TMN_stage)
Basic_data$HER2_status <- as.factor(Basic_data$HER2_status)
Basic_data$RNF114 <- as.factor(Basic_data$RNF114)
table_summary <- Basic_data %>%
  select(RNF114, age, tumor_size, pT, pN,histological_grade, TMN_stage,HER2_status) %>%
  tbl_summary(by = RNF114) %>%
  add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>%
  add_overall() %>%
  add_n() %>%
  modify_header(label ~ "") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**RNF114**") %>%
  modify_footnote(
    all_stat_cols() ~ "Median (IQR) or Frequency (%)"
  ) %>%
  modify_caption("**表1.各临床病例指标与PIWIL2蛋白表达的关系**") %>%
  as_gt()
gtsave(data = table_summary, 
       filename = "表1.各临床病例指标与RNF114蛋白表达的关系.png")
Basic_data$time <- Basic_data$`生存期（月）`
Basic_data$status <- Basic_data$`末次随访患者生存状况（删失1，死亡2，复发转移3，存活4）`
Basic_data$status[Basic_data$status == "死亡"] <- "1"
Basic_data$status[Basic_data$status == "存活"] <- "0"
print(Basic_data$status)
Basic_data$status[30] <- 0
Basic_data$status <- as.numeric(Basic_data$status)
Basic_data$time <- as.numeric(Basic_data$time)
print(Basic_data$time)
require("survival")

#age
fit_age <- survfit(Surv(time, status) ~ age, data = Basic_data)
print(fit_age)
KM_age <- ggsurvplot(fit_age, data = Basic_data,  pval = TRUE,  xlab = "Time in months")
ggsave("survival_plot_age.png", plot = KM_age$plot, width = 8, height = 6, dpi = 300)
#tumor_size
fit_tumor_size <- survfit(Surv(time, status) ~ tumor_size, data = Basic_data)
KM_tumor_size <- ggsurvplot(fit_tumor_size, data = Basic_data,pval = TRUE,xlab = "Time in months")
ggsave("survival_plot_tumor_size.png", plot = KM_tumor_size$plot, width = 8, height = 6, dpi = 300)
#RNF114
fit_RNF114 <- survfit(Surv(time, status) ~ RNF114, data = Basic_data)
KM_RNF114 <- ggsurvplot(fit_RNF114, data = Basic_data,pval = TRUE,xlab = "Time in months")
ggsave("survival_plot_RNF114.png", plot = KM_RNF114$plot, width = 8, height = 6, dpi = 300)
#pT
fit_pT <- survfit(Surv(time, status) ~ pT, data = Basic_data)
KM_pT <- ggsurvplot(fit_pT, data = Basic_data,pval = TRUE,xlab = "Time in months")
ggsave("survival_plot_pT.png", plot = KM_pT$plot, width = 8, height = 6, dpi = 300)
#pN
fit_pN <- survfit(Surv(time, status) ~ pN, data = Basic_data)
KM_pN <- ggsurvplot(fit_pN, data = Basic_data,pval = TRUE,xlab = "Time in months")
ggsave("survival_plot_pN.png", plot = KM_pN$plot, width = 8, height = 6, dpi = 300)
#histological_grade
fit_histological_grade <- survfit(Surv(time, status) ~ histological_grade, data = Basic_data)
KM_histological_grade <- ggsurvplot(fit_histological_grade, data = Basic_data,pval = TRUE,xlab = "Time in months")
ggsave("survival_plot_histological_grade.png", plot = KM_histological_grade$plot, width = 8, height = 6, dpi = 300)
#TMN_stage
fit_TMN_stage <- survfit(Surv(time, status) ~ TMN_stage, data = Basic_data)
KM_TMN_stage <- ggsurvplot(fit_TMN_stage, data = Basic_data,pval = TRUE,xlab = "Time in months")
ggsave("survival_plot_TMN_stage.png", plot = KM_TMN_stage$plot, width = 8, height = 6, dpi = 300)
#HER2_status
fit_HER2_status <- survfit(Surv(time, status) ~ HER2_status, data = Basic_data)
KM_HER2_status <- ggsurvplot(fit_HER2_status, data = Basic_data,pval = TRUE,xlab = "Time in months")
ggsave("survival_plot_HER2_status.png", plot = KM_HER2_status$plot, width = 8, height = 6, dpi = 300)

#cox model
# 加载必要包
library(survival)
library(dplyr)
install.packages("coxphf")
library(coxphf)

# 数据清理与预处理
Basic_data <- Basic_data %>% 
  mutate(
    pT = factor(pT, levels = c("T1", "T2", "T3", "T4")),
    pN = factor(pN, levels = c("N0", "N1", "N2", "N3")),
    tumor_size = factor(tumor_size, levels = c("<=3", ">3"))
  )

# 创建单因素Cox模型公式
covariates <- c("protein", "age", "tumor_size", "pT", "pN", "histological_grade", "TMN_stage", "HER2_status")
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(time, status) ~', x)))

# 拟合单因素Cox模型
univ_models <- lapply(univ_formulas, function(x) tryCatch(coxph(x, data = Basic_data), error = function(e) NA))

# 提取结果
# Extract results with consistent row format
univ_results <- lapply(univ_models, function(x) {
  if (inherits(x, "coxph")) {
    x <- summary(x)
    
    # Handle multiple rows for categorical variables
    beta <- signif(x$coef[, 1], digits = 2)
    HR <- signif(x$coef[, 2], digits = 2)
    HR.confint.lower <- signif(x$conf.int[, "lower .95"], 2)
    HR.confint.upper <- signif(x$conf.int[, "upper .95"], 2)
    p.value <- signif(x$wald["pvalue"], digits = 2)
    
    # Combine HR and confidence interval
    HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
    
    # Create a data frame for each variable
    res <- data.frame(
      beta = beta,
      HR = HR,
      p.value = p.value,
      row.names = rownames(x$coef)
    )
    return(res)
  } else {
    # Return NA if model fails
    return(data.frame(beta = NA, HR = NA, p.value = NA, row.names = NULL))
  }
})

# Combine all results into a single data frame
# Combine all results into a single data frame
res <- do.call(rbind, univ_results)
res$Variable <- rep(covariates, sapply(univ_results, nrow))
res$Subcategory <- rownames(res)  # 保留子分类名（例如 T1, T2）

# Sort the data by Variable and Subcategory
res <- res %>%
  arrange(
    factor(Variable, levels = c("pT", "pN", "tumor_size", "protein", "age", "histological_grade", "TMN_stage", "HER2_status")), # 控制主分类顺序
    Subcategory  # 控制子分类顺序
  )
res$Subcategory <- rownames(res)  # Ensure `Subcategory` is created from row names

# Rearrange columns
library(gt)

# 添加分组列（pT, pN, 等）
res <- res %>%
  arrange(
    factor(Variable, levels = c("pT", "pN", "tumor_size", "protein", "age", "histological_grade", "TMN_stage", "HER2_status"))
  )

# 格式化表格
res_table <- gt(res) %>%
  tab_header(
    title = "Univariate Cox Regression Results"  # 表格标题
  ) %>%
  tab_row_group(
    group = "pT", rows = res$Variable == "pT"
  ) %>%
  tab_row_group(
    group = "pN", rows = res$Variable == "pN"
  ) %>%
  fmt_number(
    columns = c("beta", "p.value"),
    decimals = 2  # 控制小数点位数
  ) %>%
  fmt_number(
    columns = "HR",
    pattern = "{x}"
  ) %>%
  cols_label(
    Variable = "Variable",
    Subcategory = "Subcategory",
    beta = "β",
    HR = "HR (95% CI)",
    p.value = "p-value"
  ) %>%
  opt_row_striping() %>%  # 添加条纹行
  tab_style(
    style = cell_borders(
      sides = "all",
      color = "gray",
      weight = px(1)
    ),
    locations = cells_body()
  )

# 保存表格为图片
gtsave(res_table, "Optimized_Univariate_Cox_Analysis_Results.png")
library(writexl)
write_xlsx(res, "Univariate_Cox_Analysis_Results.xlsx")
#spearman
Basic_data$年龄 <- as.numeric(Basic_data$年龄)
Basic_data$双变性量化打分<-as.numeric(Basic_data$双变性量化打分)
result = cor.test(Basic_data$年龄, Basic_data$双变性量化打分, method = "spearman")
print(result)
Basic_data$size_continuously <- as.numeric(Basic_data$size_continuously)
result = cor.test(Basic_data$size_continuously, Basic_data$双变性量化打分, method = "spearman")
print(result)
Basic_data$HER2阳性与否<- as.numeric(Basic_data$HER2阳性与否)
result = cor.test(Basic_data$HER2阳性与否, Basic_data$双变性量化打分, method = "spearman")
print(result)
Basic_data$histological_grade<- as.numeric(Basic_data$histological_grade)
result = cor.test(Basic_data$histological_grade, Basic_data$双变性量化打分, method = "spearman")
print(result)
summary(Basic_data$pN)
Basic_data$pN_n <- as.numeric(factor(Basic_data$pN, 
                                         levels = c("N0", "N1", "N2", "N3")))
result = cor.test(Basic_data$pN_n, Basic_data$双变性量化打分, method = "spearman")
print(result)
summary(Basic_data$pT)
Basic_data$pT_n <- as.numeric(factor(Basic_data$pT, 
                                     levels = c("T1", "T2", "T3", "T4")))
result = cor.test(Basic_data$pT_n, Basic_data$双变性量化打分, method = "spearman")
print(result)
summary(Basic_data$TMN_stage)
Basic_data$TMN_stage_n <- as.numeric(factor(Basic_data$TMN_stage, 
                                     levels = c("I", "II", "III", "IV")))
result = cor.test(Basic_data$TMN_stage_n, Basic_data$双变性量化打分, method = "spearman")
print(result)


table(Basic_data$pT, Basic_data$status)
table(Basic_data$pN, Basic_data$status)
summary(Basic_data$tumor_size)



