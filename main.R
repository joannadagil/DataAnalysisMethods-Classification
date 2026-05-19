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



