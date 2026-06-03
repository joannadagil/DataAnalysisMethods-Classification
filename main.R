# ----------------------------------------------------------
## Wczytanie danych i statystyki opisowe
# ----------------------------------------------------------

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


# Usuwanie outlierów (tu kończy się stary kod przed zmianami)
bean_clean <- bean_processed[-all_outliers, ]
cat("\nLiczba obserwacji po czyszczeniu:", nrow(bean_clean), "\n")

# --------------------------------------------------------------
# PODZIAŁ I SKALOWANIE CECH
# --------------------------------------------------------------

library(nnet)

class_column <- "Class"
# ZMIANA 1: Pracujemy na WYCZYSZCZONYCH danych (bean_clean)
bean_clean[[class_column]] <- as.factor(bean_clean[[class_column]])

# Sprawdzenie rozkładu klas po czyszczeniu
cat("\nRozkład klas po czyszczeniu:\n")
print(table(bean_clean[[class_column]]))

# ----------------------------------------------------------
## PODZIAŁ NA ZBIÓR UCZĄCY I TESTOWY (ze stratąfikacją)
# ----------------------------------------------------------
set.seed(42)

train_indices <- unlist(
  tapply(
    seq_len(nrow(bean_clean)),
    bean_clean[[class_column]],
    function(indices) {
      sample(indices, size = floor(0.7 * length(indices)))
    }
  )
)

train_data <- bean_clean[train_indices, ]
test_data  <- bean_clean[-train_indices, ]

cat("\nLiczba obserwacji w zbiorze uczącym:", nrow(train_data), "\n")
cat("Liczba obserwacji w zbiorze testowym:", nrow(test_data), "\n\n")

# --- WYKRESY PROPORCJI KLAS ---
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Brakuje pakietu ggplot2.")
}
library(ggplot2)

class_levels <- levels(bean_clean[[class_column]])

get_class_proportions <- function(data, set_name) {
  counts <- table(factor(data[[class_column]], levels = class_levels))
  data.frame(
    Zbior = set_name,
    Class = factor(names(counts), levels = class_levels),
    Liczebnosc = as.integer(counts),
    Proportion = as.numeric(counts) / sum(counts)
  )
}

class_proportions <- rbind(
  get_class_proportions(bean_clean, "Cały zbiór (oczyszczony)"),
  get_class_proportions(train_data, "Zbiór uczący"),
  get_class_proportions(test_data, "Zbiór testowy")
)

class_proportions$Procent <- paste0(round(100 * class_proportions$Proportion, 1), "%")

plot_class_proportions <- ggplot(class_proportions, aes(x = Class, y = Proportion)) +
  geom_col() +
  geom_text(aes(label = Procent), vjust = -0.3, size = 3) +
  facet_wrap(~ Zbior, ncol = 1) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     limits = c(0, max(class_proportions$Proportion) * 1.15)) +
  labs(title = "", x = "Klasa fasoli", y = "Udział obserwacji") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

class_proportions_pdf <- file.path(report_dir, "proporcje_klas_original_train_test.pdf")
ggsave(filename = class_proportions_pdf, plot = plot_class_proportions, width = 10, height = 9)
cat("Zapisano wykres proporcji klas w pliku:", class_proportions_pdf, "\n")

# ----------------------------------------------------------
## SKALOWANIE ZMIENNYCH NUMERYCZNYCH 
## ZMIANA 2: Parametry liczone TYLKO ze zbioru uczącego!
# ----------------------------------------------------------

predictor_columns <- setdiff(names(bean_clean), class_column)
numeric_predictors <- predictor_columns[sapply(bean_clean[predictor_columns], is.numeric)]

train_scaled <- train_data
test_scaled <- test_data

# Złota zasada: wyciągamy średnią i odchylenie tylko ze zbioru treningowego
train_means <- sapply(train_data[numeric_predictors], mean, na.rm = TRUE)
train_sds <- sapply(train_data[numeric_predictors], sd, na.rm = TRUE)

# Skalujemy oba zbiory używając parametrów treningowych
train_scaled[numeric_predictors] <- scale(
  train_data[numeric_predictors],
  center = train_means,
  scale = train_sds
)

test_scaled[numeric_predictors] <- scale(
  test_data[numeric_predictors],
  center = train_means,
  scale = train_sds
)

cat("\nDane zostały poprawnie podzielone i przeskalowane (bez wycieku danych).\n")

# --- STATYSTYKI OPISOWE DLA ZBIORU UCZĄCEGO (po skalowaniu) ---
numeric_summary <- data.frame(
  Zmienna = numeric_predictors,
  Srednia = sapply(train_scaled[numeric_predictors], mean, na.rm = TRUE),
  Mediana = sapply(train_scaled[numeric_predictors], median, na.rm = TRUE),
  Minimum = sapply(train_scaled[numeric_predictors], min, na.rm = TRUE),
  Maksimum = sapply(train_scaled[numeric_predictors], max, na.rm = TRUE),
  OdchylenieStandardowe = sapply(train_scaled[numeric_predictors], sd, na.rm = TRUE),
  Skosnosc = sapply(train_scaled[numeric_predictors], moments::skewness, na.rm = TRUE),
  row.names = NULL
)

numeric_summary[-1] <- lapply(numeric_summary[-1], function(x) round(x, 4))

latex_tables <- list(
  knitr::kable(
    numeric_summary,
    format = "latex",
    booktabs = TRUE,
    position = "H",
    caption = "Statystyki opisowe zmiennych numerycznych dla zbioru uczącego (po skalowaniu)",
    label = "statystyki-opisowe-numeryczne-scaled"
  )
)

latex_output <- unlist(lapply(latex_tables, function(table) c(table, "")))
latex_output_file <- file.path(report_dir, "tabele_statystyki_opisowe_przeskalowane.tex")
writeLines(latex_output, latex_output_file)

cat("Zapisano tabele LaTeX w pliku:", latex_output_file, "\n")


# ----------------------------------------------------------
## 3.1 WIELOMIANOWA REGRESJA LOGISTYCZNA
# ----------------------------------------------------------

# dopasowywanie modelu
model_multinom <- multinom(
  Class ~ .,
  data = train_scaled,
  trace = FALSE,
  maxit = 1000,
  MaxNWts = 10000
)

summary(model_multinom)

# predykcja
pred_multinom_class <- predict(
  model_multinom,
  newdata = test_scaled,
  type = "class"
)

pred_multinom_prob <- predict(
  model_multinom,
  newdata = test_scaled,
  type = "probs"
)

# ocena
conf_matrix_multinom <- table(
  Rzeczywista = test_scaled$Class,
  Przewidziana = pred_multinom_class
)

print(conf_matrix_multinom)

accuracy_multinom <- mean(pred_multinom_class == test_scaled$Class)

cat("Accuracy dla wielomianowej regresji logistycznej:", accuracy_multinom, "\n")

calculate_class_metrics <- function(conf_matrix) {
  classes <- rownames(conf_matrix)

  metrics <- data.frame(
    Klasa = classes,
    Precision = NA,
    Recall = NA,
    F1 = NA
  )

  for (class_name in classes) {
    TP <- conf_matrix[class_name, class_name]
    FP <- sum(conf_matrix[, class_name]) - TP
    FN <- sum(conf_matrix[class_name, ]) - TP

    precision <- ifelse(TP + FP == 0, NA, TP / (TP + FP))
    recall <- ifelse(TP + FN == 0, NA, TP / (TP + FN))
    f1 <- ifelse(
      is.na(precision) | is.na(recall) | precision + recall == 0,
      NA,
      2 * precision * recall / (precision + recall)
    )

    metrics[metrics$Klasa == class_name, "Precision"] <- precision
    metrics[metrics$Klasa == class_name, "Recall"] <- recall
    metrics[metrics$Klasa == class_name, "F1"] <- f1
  }

  metrics
}

metrics_multinom <- calculate_class_metrics(conf_matrix_multinom)

metrics_multinom[, c("Precision", "Recall", "F1")] <- round(
  metrics_multinom[, c("Precision", "Recall", "F1")],
  4
)

print(metrics_multinom)


# latex
predictions_multinom <- data.frame(
  Rzeczywista = test_scaled$Class,
  Przewidziana = pred_multinom_class,
  pred_multinom_prob
)


results_multinom <- data.frame(
  Metoda = "Wielomianowa regresja logistyczna",
  Accuracy = round(accuracy_multinom, 4)
)

latex_multinom_tables <- list(
  knitr::kable(
    as.data.frame.matrix(conf_matrix_multinom),
    format = "latex",
    booktabs = TRUE,
    position = "H",
    caption = "Macierz pomyłek dla wielomianowej regresji logistycznej",
    label = "conf-matrix-multinom"
  ),
  knitr::kable(
    results_multinom,
    format = "latex",
    booktabs = TRUE,
    position = "H",
    caption = "Dokładność klasyfikacji dla wielomianowej regresji logistycznej",
    label = "accuracy-multinom"
  ),
  knitr::kable(
    metrics_multinom,
    format = "latex",
    booktabs = TRUE,
    position = "H",
    caption = "Miary jakości klasyfikacji dla wielomianowej regresji logistycznej",
    label = "metrics-multinom"
  )
)

latex_multinom_output <- unlist(lapply(latex_multinom_tables, function(table) c(table, "")))

latex_multinom_output_file <- file.path(report_dir, "wyniki_regresja_logistyczna.tex")

writeLines(latex_multinom_output, latex_multinom_output_file)

cat("Zapisano wyniki regresji logistycznej w pliku:", latex_multinom_output_file, "\n")


# ----------------------------------------------------------
## 3.2. RANDOM FOREST
# ----------------------------------------------------------
library(randomForest)

# Ustawienie ziarna dla powtarzalności wyników
set.seed(42)

# Budowa modelu lasu losowego (500 drzew, obliczanie istotności cech)
model_rf <- randomForest(
  Class ~ .,
  data = train_scaled,
  ntree = 500,
  importance = TRUE
)

# Wyświetlenie podsumowania modelu (w tym błędu Out-of-Bag)
print(model_rf)

# Predykcja klas na poprawnie wydzielonym i przeskalowanym zbiorze testowym
pred_rf_class <- predict(
  model_rf,
  newdata = test_scaled,
  type = "response"
)

# Predykcja prawdopodobieństw przynależności do klas (opcjonalnie)
pred_rf_prob <- predict(
  model_rf,
  newdata = test_scaled,
  type = "prob"
)

# Generowanie macierzy pomyłek
conf_matrix_rf <- table(
  Rzeczywista = test_scaled$Class,
  Przewidziana = pred_rf_class
)
print("Macierz pomyłek dla Lasu Losowego:")
print(conf_matrix_rf)

# Obliczenie ogólnej dokładności (Accuracy)
accuracy_rf <- mean(pred_rf_class == test_scaled$Class)
cat("Accuracy dla Random Forest:", round(accuracy_rf, 4), "\n")

# Wyznaczenie miar jakości dla poszczególnych klas (Precision, Recall, F1)
metrics_rf <- calculate_class_metrics(conf_matrix_rf)
metrics_rf[, c("Precision", "Recall", "F1")] <- round(
  metrics_rf[, c("Precision", "Recall", "F1")],
  4
)
print("Metryki dla poszczególnych klas:")
print(metrics_rf)

# Wyciągnięcie i zaokrąglenie ważności cech (Feature Importance)
importance_rf <- as.data.frame(importance(model_rf))
importance_rf$Zmienna <- rownames(importance_rf)
importance_rf <- importance_rf[, c("Zmienna", "MeanDecreaseAccuracy", "MeanDecreaseGini")]
importance_rf[, 2:3] <- round(importance_rf[, 2:3], 4)

# Przygotowanie ramki danych z ogólnym wynikiem dokładności
results_rf <- data.frame(
  Metoda = "Random Forest (Las Losowy)",
  Accuracy = round(accuracy_rf, 4)
)

# Generowanie kodu tabel LaTeX przy użyciu pakietu knitr
latex_rf_tables <- list(
  knitr::kable(
    as.data.frame.matrix(conf_matrix_rf),
    format = "latex",
    booktabs = TRUE,
    position = "H",
    caption = "Macierz pomyłek dla lasu losowego (Random Forest)",
    label = "conf-matrix-rf"
  ),
  knitr::kable(
    results_rf,
    format = "latex",
    booktabs = TRUE,
    position = "H",
    caption = "Dokładność klasyfikacji dla lasu losowego",
    label = "accuracy-rf"
  ),
  knitr::kable(
    metrics_rf,
    format = "latex",
    booktabs = TRUE,
    position = "H",
    caption = "Miary jakości klasyfikacji dla lasu losowego dla poszczególnych odmian fasoli",
    label = "metrics-rf"
  ),
  knitr::kable(
    importance_rf,
    format = "latex",
    booktabs = TRUE,
    position = "H",
    caption = "Istotność zmiennych w modelu Random Forest",
    label = "importance-rf"
  )
)

# Połączenie tabel i zapis do pliku wynikowego .tex dla raportu
latex_rf_output <- unlist(lapply(latex_rf_tables, function(table) c(table, "")))
latex_rf_output_file <- file.path(report_dir, "wyniki_random_forest.tex")
writeLines(latex_rf_output, latex_rf_output_file)

cat("Wyniki metody Random Forest zostały pomyślnie zapisane w:", latex_rf_output_file, "\n")


# ----------------------------------------------------------
## 3.3. SIEĆ NEURONOWA
# ----------------------------------------------------------
# klasyczny klasyfikator wieloklasowy
# warstwa wejściowa z tyloma wejściami ile jest cech
# jedna warstwa ukryta z 8 neuronami
# warstwa wyjściowa ma tyle neuronów ile klas - gatunków fasoli
# dla klasyfikacji wieloklasowej zastosowano funkcję softmax.

set.seed(42)

model_neural_net <- nnet(
  Class ~ .,         # wykorzystujemy wszystkie cechy
  data = train_scaled,
  size = 8,          # liczba neuronów w warstwie ukrytej
  decay = 0.001,     # regularyzacja wag
  #softmax = TRUE,    # interpretacja wyniku jako prawdopodobieństwo
  maxit = 1000,      # maksymalna liczba iteracji uczenia
  MaxNWts = 10000,   # zwiększenie limitu wag w sieci
  trace = TRUE
)

summary(model_neural_net)


# predykcja klas i prawdopodobieństw

pred_nn_class <- predict(
  model_neural_net,
  newdata = test_scaled,
  type = "class"
)

pred_nn_class <- factor(pred_nn_class, levels = levels(test_scaled$Class))

pred_nn_prob <- predict(
  model_neural_net,
  newdata = test_scaled,
  type = "raw"
)


# ocena jakości klasyfikacji

conf_matrix_nn <- table(
  Rzeczywista = test_scaled$Class,
  Przewidziana = pred_nn_class
)

print(conf_matrix_nn)

accuracy_nn <- mean(pred_nn_class == test_scaled$Class)

cat("Accuracy dla sieci neuronowej:", accuracy_nn, "\n")

metrics_nn <- calculate_class_metrics(conf_matrix_nn)

metrics_nn[, c("Precision", "Recall", "F1")] <- round(
  metrics_nn[, c("Precision", "Recall", "F1")],
  4
)

print(metrics_nn)



# Zapis predykcji i wyników

predictions_nn <- data.frame(
  Rzeczywista = test_scaled$Class,
  Przewidziana = pred_nn_class,
  pred_nn_prob
)

# Opcjonalnie: czytelniejsze nazwy kolumn z prawdopodobieństwami
names(predictions_nn)[3:ncol(predictions_nn)] <- paste0(
  "P_",
  make.names(colnames(pred_nn_prob))
)

write.csv(
  predictions_nn,
  file = file.path(report_dir, "predykcje_siec_neuronowa.csv"),
  row.names = FALSE
)

results_nn <- data.frame(
  Metoda = "Sieć neuronowa",
  Accuracy = round(accuracy_nn, 4)
)

latex_nn_tables <- list(
  knitr::kable(
    as.data.frame.matrix(conf_matrix_nn),
    format = "latex",
    booktabs = TRUE,
    position = "H",
    caption = "Macierz pomyłek dla sieci neuronowej",
    label = "conf-matrix-nn"
  ),
  knitr::kable(
    results_nn,
    format = "latex",
    booktabs = TRUE,
    position = "H",
    caption = "Dokładność klasyfikacji dla sieci neuronowej",
    label = "accuracy-nn"
  ),
  knitr::kable(
    metrics_nn,
    format = "latex",
    booktabs = TRUE,
    position = "H",
    caption = "Miary jakości klasyfikacji dla sieci neuronowej",
    label = "metrics-nn"
  )
)

latex_nn_output <- unlist(lapply(latex_nn_tables, function(table) c(table, "")))

latex_nn_output_file <- file.path(report_dir, "wyniki_siec_neuronowa.tex")

writeLines(latex_nn_output, latex_nn_output_file)

cat("Zapisano wyniki sieci neuronowej w pliku:", latex_nn_output_file, "\n")
cat("Zapisano predykcje sieci neuronowej w pliku:", file.path(report_dir, "predykcje_siec_neuronowa.csv"), "\n")




# ----------------------------------------------------------
## 3.4. METODA HYBRYDOWA
# ----------------------------------------------------------




# ----------------------------------------------------------
## 5. SZTUCZNE DANE
# ----------------------------------------------------------

