# Wczytanie danych i statystyki opisowe

library(moments)
library(knitr)

data_dir <- "data"
report_dir <- "report"
arff_file <- file.path(data_dir, "Dry_Bean_Dataset.arff")

if (!dir.exists(data_dir)) {
  stop("Nie znaleziono folderu data.")
}

if (!file.exists(arff_file)) {
  stop("Nie znaleziono pliku: ", arff_file)
}

if (!dir.exists(report_dir)) {
  dir.create(report_dir, recursive = TRUE)
}

read_arff_simple <- function(path) {
  lines <- readLines(path, warn = FALSE)

  attribute_lines <- grep("^@ATTRIBUTE", lines, ignore.case = TRUE, value = TRUE)
  column_names <- sub("^@ATTRIBUTE[[:space:]]+([^[:space:]]+).*", "\\1", attribute_lines, ignore.case = TRUE)

  data_start <- grep("^@DATA", lines, ignore.case = TRUE)
  if (length(data_start) == 0) {
    stop("Plik ARFF nie zawiera sekcji @DATA.")
  }

  data_lines <- lines[(data_start[1] + 1):length(lines)]
  data_lines <- data_lines[nzchar(trimws(data_lines))]

  dataset <- read.csv(
    text = paste(data_lines, collapse = "\n"),
    header = FALSE,
    stringsAsFactors = FALSE
  )

  names(dataset) <- column_names
  dataset
}

bean_data <- read_arff_simple(arff_file)

numeric_columns <- names(bean_data)[sapply(bean_data, is.numeric)]
categorical_columns <- names(bean_data)[!names(bean_data) %in% numeric_columns]

numeric_summary <- data.frame(
  Zmienna = numeric_columns,
  Srednia = sapply(bean_data[numeric_columns], mean, na.rm = TRUE),
  Mediana = sapply(bean_data[numeric_columns], median, na.rm = TRUE),
  Minimum = sapply(bean_data[numeric_columns], min, na.rm = TRUE),
  Maksimum = sapply(bean_data[numeric_columns], max, na.rm = TRUE),
  OdchylenieStandardowe = sapply(bean_data[numeric_columns], sd, na.rm = TRUE),
  Skosnosc = sapply(bean_data[numeric_columns], moments::skewness, na.rm = TRUE),
  row.names = NULL
)

numeric_summary[-1] <- lapply(numeric_summary[-1], function(x) round(x, 4))

latex_tables <- list(
  knitr::kable(
    numeric_summary,
    format = "latex",
    booktabs = TRUE,
    position = "H",
    caption = "Statystyki opisowe zmiennych numerycznych",
    label = "statystyki-opisowe-numeryczne"
  )
)

if (length(categorical_columns) > 0) {
  for (column_name in categorical_columns) {
    frequency_table <- as.data.frame(table(bean_data[[column_name]], useNA = "ifany"))
    names(frequency_table) <- c(column_name, "liczebnosc")

    latex_tables[[length(latex_tables) + 1]] <- knitr::kable(
      frequency_table,
      format = "latex",
      booktabs = TRUE,
      position = "H",
      caption = paste("Liczebnosci dla zmiennej", column_name),
      label = paste0("liczebnosci-", column_name)
    )
  }
}

latex_output <- unlist(lapply(latex_tables, function(table) c(table, "")))
latex_output_file <- file.path(report_dir, "tabele_statystyki_opisowe.tex")
writeLines(latex_output, latex_output_file)

cat("Wczytano obserwacji:", nrow(bean_data), "\n")
cat("Wczytano zmiennych:", ncol(bean_data), "\n\n")
cat("Zapisano tabele LaTeX w pliku:", latex_output_file, "\n")



# 2.2.3 Transformacje danych

library(dplyr)
library(ggplot2)

selected_variables <- c(
  "Area",
  "Perimeter",
  "MajorAxisLength",
  "MinorAxisLength",
  "ConvexArea",
  "EquivDiameter",
  "Eccentricity",
  "roundness",
  "Compactness",
  "Solidity",
  "AspectRation",
  "Extent"
)

bean_selected <- bean_data[, c(selected_variables, "Class")]

# Analiza skośności
skewness_values <- sapply(
  bean_selected[selected_variables],
  moments::skewness,
  na.rm = TRUE
)

skewness_table <- data.frame(
  Zmienna = names(skewness_values),
  Skosnosc = round(skewness_values, 4)
)

print(skewness_table)

# Transformacja logarytmiczna
high_skew_vars <- names(skewness_values[abs(skewness_values) > 1])

bean_transformed <- bean_selected

for (var in high_skew_vars) {
  bean_transformed[[var]] <- log(bean_transformed[[var]])
}

cat("\nZmienne po log-transformacji:\n")
print(high_skew_vars)

bean_processed <- bean_transformed


# 2.2.5 Obserwacje odstające

outlier_indices_iqr <- c()

outlier_indices_iqr <- c()

for (class_name in unique(bean_processed$Class)) {
  
  class_subset <- bean_processed[
    bean_processed$Class == class_name,
  ]
  
  class_indices <- which(
    bean_processed$Class == class_name
  )
  
  for (var in selected_variables) {
    
    q1 <- quantile(
      class_subset[[var]],
      0.25,
      na.rm = TRUE
    )
    
    q3 <- quantile(
      class_subset[[var]],
      0.75,
      na.rm = TRUE
    )
    
    iqr_value <- q3 - q1
    
    lower_bound <- q1 - 1.5 * iqr_value
    upper_bound <- q3 + 1.5 * iqr_value
    
    local_outliers <- which(
      class_subset[[var]] < lower_bound |
        class_subset[[var]] > upper_bound
    )
    
    global_outliers <- class_indices[local_outliers]
    
    outlier_indices_iqr <- union(
      outlier_indices_iqr,
      global_outliers
    )
  }
}

cat("\nLiczba outlierów (IQR):", length(outlier_indices_iqr), "\n")


# Mahalanobis

numeric_data <- bean_processed[, selected_variables]

cov_matrix <- cov(numeric_data) + diag(1e-6, ncol(numeric_data))
mean_vector <- colMeans(numeric_data)

mahalanobis_distance <- mahalanobis(
  numeric_data,
  center = mean_vector,
  cov = cov_matrix
)

threshold <- qchisq(0.999, df = length(selected_variables))

outlier_indices_mahal <- which(mahalanobis_distance > threshold)

cat("\nLiczba outlierów (Mahalanobis):", length(outlier_indices_mahal), "\n")


# Łączenie outlierów
all_outliers <- union(outlier_indices_iqr, outlier_indices_mahal)
#all_outliers <- outlier_indices_iqr

cat("\nŁączna liczba outlierów:", length(all_outliers), "\n")


# Usuwanie outlierów
bean_clean <- bean_processed[-all_outliers, ]

cat("\nLiczba obserwacji po czyszczeniu:", nrow(bean_clean), "\n")


# FINALNE SKALOWANIE (PO CZYSZCZENIU)

bean_scaled <- data.frame(
  scale(bean_clean[, selected_variables])
)

bean_scaled$Class <- bean_clean$Class

cat("\nDane zostały przeskalowane po usunięciu obserwacji odstających.\n")


# STATYSTYKI DLA OCZYSZCZONYCH DANYCH

numeric_columns_clean <- names(bean_scaled)[
  sapply(bean_scaled, is.numeric)
]

categorical_columns_clean <- names(bean_scaled)[
  !names(bean_scaled) %in% numeric_columns_clean
]

numeric_summary_clean <- data.frame(
  Zmienna = numeric_columns_clean,
  Srednia = sapply(bean_scaled[numeric_columns_clean], mean, na.rm = TRUE),
  Mediana = sapply(bean_scaled[numeric_columns_clean], median, na.rm = TRUE),
  Minimum = sapply(bean_scaled[numeric_columns_clean], min, na.rm = TRUE),
  Maksimum = sapply(bean_scaled[numeric_columns_clean], max, na.rm = TRUE),
  OdchylenieStandardowe = sapply(bean_scaled[numeric_columns_clean], sd, na.rm = TRUE),
  Skosnosc = sapply(bean_scaled[numeric_columns_clean], moments::skewness, na.rm = TRUE),
  row.names = NULL
)

numeric_summary_clean[-1] <- lapply(numeric_summary_clean[-1], round, 4)

latex_tables_clean <- list(
  knitr::kable(
    numeric_summary_clean,
    format = "latex",
    booktabs = TRUE,
    caption = "Statystyki opisowe po preprocessingu (po czyszczeniu i standaryzacji)"
  )
)

for (column_name in categorical_columns_clean) {
  
  new_frequency_table <- as.data.frame(
    table(bean_scaled[[column_name]], useNA = "ifany")
  )
  
  names(new_frequency_table) <- c(column_name, "liczebnosc")
  
  latex_tables_clean[[length(latex_tables_clean) + 1]] <-
    knitr::kable(
      new_frequency_table,
      format = "latex",
      booktabs = TRUE,
      caption = paste("Liczebnosci dla", column_name)
    )
}

latex_output_clean <- unlist(lapply(latex_tables_clean, function(x) c(x, "")))

writeLines(
  latex_output_clean,
  file.path(report_dir, "tabele_statystyki_opisowe_clean.tex")
)