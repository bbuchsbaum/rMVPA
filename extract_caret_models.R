library(caret)

# Extract model info for random forest and sparse PLS
rf_info <- getModelInfo("rf")$rf
spls_info <- getModelInfo("spls")$spls

# Print structure for random forest
cat("=== Random Forest (rf) Model Info ===\n")
cat("\nModel Type:\n")
print(rf_info$type)

cat("\nParameters:\n")
print(rf_info$parameters)

cat("\nGrid Function:\n")
print(rf_info$grid)

cat("\nFit Function:\n")
print(rf_info$fit)

cat("\nPredict Function:\n")
print(rf_info$predict)

cat("\nProb Function (if exists):\n")
if (!is.null(rf_info$prob)) {
  print(rf_info$prob)
}

cat("\n\n=== Sparse PLS (spls) Model Info ===\n")
cat("\nModel Type:\n")
print(spls_info$type)

cat("\nParameters:\n")
print(spls_info$parameters)

cat("\nGrid Function:\n")
print(spls_info$grid)

cat("\nFit Function:\n")
print(spls_info$fit)

cat("\nPredict Function:\n")
print(spls_info$predict)

cat("\nProb Function (if exists):\n")
if (!is.null(spls_info$prob)) {
  print(spls_info$prob)
}

# Also check what libraries are required
cat("\n\nRequired Libraries:\n")
cat("RF:", rf_info$library, "\n")
cat("SPLS:", spls_info$library, "\n")