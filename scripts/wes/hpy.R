# Load required libraries

source(here::here("bin/R/lib/study_lib.R"))

# Read the hap.py summary data
# hap_dir <- "/home/zhonggr/projects/250224_DFSP_WES/data/benchmark/NA12878/hpy"
hap_dir <- "/home/zhonggr/projects/250224_DFSP_WES/data/benchmark/NA12878/hpy_default"
summary_file <- dir_ls(here(hap_dir), glob = "*.summary.csv")
happy_data <- read_csv(summary_file)

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
        vjust = -0.5, size = 3
    ) +
    scale_y_continuous(
        labels = scales::percent_format(),
        limits = c(0, max(performance_metrics$Value) * 1.1)
    ) +
    labs(
        x = "Variant Type", y = "Performance (%)",
        fill = "Metric"
    ) +
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
        True_Negatives = NA # Not directly available from hap.py
    ) |>
    select(Type, True_Positives, False_Negatives, False_Positives) |>
    pivot_longer(cols = -Type, names_to = "Category", values_to = "Count")

p2 <- ggplot(confusion_data, aes(x = Type, y = Count, fill = Category)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = scales::comma(Count)),
        position = position_dodge(width = 0.9),
        vjust = -0.5, size = 3
    ) +
    scale_y_log10(
        labels = scales::comma_format(),
        limits = c(1, max(confusion_data$Count, na.rm = TRUE) * 2)
    ) +
    labs(
        x = "Variant Type", y = "Count (log scale)",
        fill = "Classification"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)

# 3. Create a summary table (removed redundant p3 plot)
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
    select(
        Type, `Total Truth`, `Found (TP)`, `Missed (FN)`,
        `False Pos`, `Sensitivity %`, `Precision %`, `F1-Score`
    )

print(summary_table)

# Save as a formatted table
library(knitr)
library(kableExtra)

kable(summary_table, format = "html") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"))

# 5. Create a combined dashboard
library(patchwork)

# Combine plots - simplified to avoid redundancy
combined_plot <- p1 / p2

combined_plot <- combined_plot +
    plot_annotation(
        title = "Mutect2 Pipeline vs GIAB Truth Set"
    )


print(combined_plot)

# Save the plot
SavePlot(
    plot= combined_plot & plot_theme(),
    width = 5, 
    height = 6,
    dir = hap_dir, 
    filename = "NA12878_performance_analysis.png"
)
