# ============================================================
# CHEM 124 - Lab 3: PCA of coffee NMR data
# Each ROW = one coffee sample
# Each COLUMN = one metabolite peak area (divided by trigonelline)
# ============================================================

# ---- STEP 1: Load the data ----
data <- read.csv(file.choose())

print(data)

# ---- STEP 2: Keep only the metabolite numbers ----
numbers <- data[, 4:ncol(data)]

# ---- STEP 3: Run the PCA ----
pca <- prcomp(numbers, scale. = TRUE)

percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
print(percent)

# ---- STEP 4: Scores plot (one dot per sample) ----
dot_colors <- c("darkgreen", "red", "blue")   # Robusta, Arabica, Unknown

plot(pca$x[, 1], pca$x[, 2],
     pch = 19, cex = 2, col = dot_colors,
     xlim = c(-3, 3), ylim = c(-3, 3),
     xlab = paste0("PC1 (", percent[1], "%)"),
     ylab = paste0("PC2 (", percent[2], "%)"),
     main = "PCA scores of coffee extracts")
abline(h = 0, v = 0, col = "grey")
text(pca$x[, 1], pca$x[, 2], labels = data$group, pos = 3)

# ---- STEP 5: Loadings plot (which metabolites drive each PC) ----
plot(0, 0, type = "n",
     xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = paste0("PC1 loading (", percent[1], "%)"),
     ylab = paste0("PC2 loading (", percent[2], "%)"),
     main = "PCA loadings by metabolite")
abline(h = 0, v = 0, col = "grey")
arrows(0, 0, pca$rotation[, 1], pca$rotation[, 2], length = 0.1, col = "navy")
text(pca$rotation[, 1], pca$rotation[, 2], labels = colnames(numbers), pos = 4, cex = 0.8)


# ---- STEP 6: Print the numbers for your report ----
round(pca$x[, 1:2], 2)         # scores (where each sample sits)
round(pca$rotation[, 1:2], 2)  # loadings (how much each metabolite matters)