# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Read the hap.py summary data
happy_data <- read.csv("/home/zhonggr/projects/250224_DFSP_WES/data/benchmark/NA12878/hpy/NA12878_happy.summary.csv")

# View the structure
str(happy_data)
head(happy_data)

# 1. Performance Metrics Bar Plot
performance_metrics <- happy_data |>
    filter(Filter == "ALL") |>
    select(Type, METRIC.Recall, METRIC.Precision, METRIC.F1_Score) |>
    pivot_longer(
        cols = starts_with("METRIC"),
        names_to = "Metric",
        values_to = "Value"
    ) |>
    mutate(Metric = gsub("METRIC\\.", "", Metric))

p1 <- ggplot(performance_metrics, aes(x = Type, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = paste0(round(Value * 100, 1), "%")), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  scale_y_continuous(labels = scales::percent_format(), 
                     limits = c(0, max(performance_metrics$Value) * 1.1)) +
  labs(title = "Performance Metrics by Variant Type",
       x = "Variant Type", y = "Performance (%)",
       fill = "Metric") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)

# 2. Truth vs Query Comparison
confusion_data <- happy_data |>
  filter(Filter == "ALL") |>
  select(Type, TRUTH.TOTAL, TRUTH.TP, TRUTH.FN, QUERY.TOTAL, QUERY.FP) |>
  mutate(
    True_Positives = TRUTH.TP,
    False_Negatives = TRUTH.FN,
    False_Positives = QUERY.FP,
    True_Negatives = NA  # Not directly available from hap.py
  ) |>
  select(Type, True_Positives, False_Negatives, False_Positives) |>
  pivot_longer(cols = -Type, names_to = "Category", values_to = "Count")

p2 <- ggplot(confusion_data, aes(x = Type, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = scales::comma(Count)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  scale_y_log10(labels = scales::comma_format(),
                limits = c(1, max(confusion_data$Count, na.rm = TRUE) * 2)) +
  labs(title = "Variant Classification Results",
       x = "Variant Type", y = "Count (log scale)",
       fill = "Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)

# 3. Detailed Performance Analysis
performance_detailed <- happy_data |>
  filter(Filter == "ALL") |>
  mutate(
    Sensitivity = METRIC.Recall,
    Specificity = 1 - (QUERY.FP / (QUERY.FP + TRUTH.TP)), # Approximate
    Missing_Rate = 1 - METRIC.Recall,
    False_Discovery_Rate = 1 - METRIC.Precision
  ) |>
  select(Type, Sensitivity, Missing_Rate, METRIC.Precision, False_Discovery_Rate) |>
  pivot_longer(cols = -Type, names_to = "Metric", values_to = "Value")

p3 <- ggplot(performance_detailed, aes(x = Metric, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = paste0(round(Value * 100, 1), "%")), 
            position = position_dodge(width = 0.9), 
            hjust = -0.1, size = 3) +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0, max(performance_detailed$Value, na.rm = TRUE) * 1.1)) +
  labs(title = "Detailed Performance Analysis",
       x = "Performance Metric", y = "Rate (%)",
       fill = "Variant Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

print(p3)

# 4. Create a summary table
summary_table <- happy_data |>
  filter(Filter == "ALL") |>
  mutate(
    `Total Truth` = scales::comma(TRUTH.TOTAL),
    `Found (TP)` = scales::comma(TRUTH.TP),
    `Missed (FN)` = scales::comma(TRUTH.FN),
    `False Pos` = scales::comma(QUERY.FP),
    `Sensitivity %` = round(METRIC.Recall * 100, 3),
    `Precision %` = round(METRIC.Precision * 100, 1),
    `F1-Score` = round(METRIC.F1_Score, 4)
  ) |>
  select(Type, `Total Truth`, `Found (TP)`, `Missed (FN)`, 
         `False Pos`, `Sensitivity %`, `Precision %`, `F1-Score`)

print(summary_table)

# Save as a formatted table
library(knitr)
library(kableExtra)

kable(summary_table, format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))

  # 5. Create a combined dashboard
library(patchwork)

# Combine plots
combined_plot <- (p1 | p2) / p3

combined_plot <- combined_plot + 
  plot_annotation(
    title = "NA12878 Variant Calling Performance Analysis",
    subtitle = "Mutect2 Pipeline vs GIAB Truth Set",
    caption = "Data: hap.py benchmarking results"
  )

print(combined_plot)

# Save the plot
ggsave("NA12878_performance_analysis.png", 
       plot = combined_plot, 
       width = 12, height = 8, 
       dpi = 300)