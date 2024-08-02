# Load necessary libraries
library(tidyverse)
library(ggvis)
library(gmodels)
library(GGally)
library(gridExtra)
library(cluster)
library(car)
library(factoextra)

liver <- read.csv("C:/Users/Raeds/Documents/Code/r-scripts/DSFSFinalPaper/raw-cirrhosis-data.csv", header = TRUE)

liver$ID <- NULL

liver <- liver %>%
  mutate(Status = case_when(
    Status == "C" ~ "Censored",
    Status == "CL" ~ "Censored (transplant)",
    Status == "D" ~ "Death",
    TRUE ~ Status
  ))

# Delete rows 313 to 418, missing crucial data
liver <- liver[-c(313:418), ]

# Convert Age in days to Age in years and round to the nearest decimal
liver <- liver %>%
  mutate(Age = (round(Age / 365.25, 1)))

# Identify categorical and continuous variables
categorical_vars <- c("Status", "Drug", "Sex", "Ascites", "Hepatomegaly", "Spiders", "Edema", "Stage")
continuous_vars <- c("N_Days", "Age", "Bilirubin", "Cholesterol", "Albumin", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin")

# Handling missing values for categorical variables by deleting the row
liver <- liver %>%
  filter(if_any(all_of(categorical_vars), ~ !is.na(.)))

# Handling missing values for continuous variables by replacing with the mean
liver <- liver %>%
  mutate(across(all_of(continuous_vars), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))

liver <- liver %>%
  mutate(
    Ascites = ifelse(Ascites == "Y", 1, 0),
    Hepatomegaly = ifelse(Hepatomegaly == "Y", 1, 0),
    Spiders = ifelse(Spiders == "Y", 1, 0),
    Edema = ifelse(Edema == "Y", 1, 0),
    Sex = ifelse(Sex == "M", 1, 0)
  )


summary(liver)


# EDA: Histograms
par(mfrow=c(2,3))
hist(liver$Age, main="Age Distribution", xlab="Age (Years)", col="skyblue")
hist(liver$Bilirubin, main="Bilirubin Distribution", xlab="Bilirubin (mg/dl)", col="salmon")
hist(liver$Cholesterol, main="Cholesterol Distribution", xlab="Cholesterol (mg/dl)", col="lightgreen")
hist(liver$Albumin, main="Albumin Distribution", xlab="Albumin (gm/dl)", col="coral")
hist(liver$SGOT, main="SGOT Distribution", xlab="SGOT (U/ml)", col="lightblue")
hist(liver$Prothrombin, main="Prothrombin Distribution", xlab="Prothrombin Time (s)", col="orchid")

# EDA: Bar Graphs
par(mfrow = c(3, 3))

# Plot 1
barplot(table(liver$Status), main="Distribution of Status", col="lightgreen", ylab="Count")

# Plot 2
barplot(table(liver$Drug), main="Distribution of Drug", col="lightblue", ylab="Count")

# Plot 3
barplot(table(liver$Sex), main="Distribution of Sex", col="lightcoral", ylab="Count")

# Plot 4
barplot(table(liver$Ascites), main="Distribution of Ascites", col="lightsalmon", ylab="Count")

# Plot 5
barplot(table(liver$Hepatomegaly), main="Distribution of Hepatomegaly", col="lightsteelblue", ylab="Count")

# Plot 6
barplot(table(liver$Spiders), main="Distribution of Spiders", col="lightpink", ylab="Count")

# Plot 7
barplot(table(liver$Edema), main="Distribution of Edema", col="lightgoldenrod", ylab="Count")

# Plot 8
barplot(table(liver$Stage), main="Distribution of Stage", col="lightseagreen", ylab="Count")

# Reset the layout to default
par(mfrow = c(1, 1))

# EDA: Correlation Heatmap
cor_matrix <- cor(liver[, continuous_vars])
heatmap(cor_matrix, col = colorRampPalette(c("blue", "white", "red"))(30))

results <- prcomp(liver[, continuous_vars], scale = TRUE)

pca_df <- as.data.frame(results$x)
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = liver$Unnamed.0), alpha = 0.6) +
  labs(title = "PCA: First Two Principal Components") +
  theme_minimal()

loadings <- results$rotation

# Screeplot for PCA
screeplot(results, type = "lines", main = "Scree Plot for PCA")

# Rank variables based on their influence in PC1 and PC2
ranked_vars_PC1 <- order(abs(loadings[,"PC1"]), decreasing = TRUE)
ranked_vars_PC2 <- order(abs(loadings[, "PC2"]), decreasing = TRUE)
# Get the names of the ranked variables for PC1 and PC2
PC1_ranked_vars <- rownames(loadings)[ranked_vars_PC1]
PC2_ranked_vars <- rownames(loadings)[ranked_vars_PC2]
# Print out the loadings from the PCA Results
print("Variables ranked by influence in PC1:")
print(PC1_ranked_vars)
print("Variables ranked by influence in PC2:")
print(PC2_ranked_vars)

#scatterplot based on pca
ggplot(liver, aes(x = Bilirubin, y = Copper, color = Status)) +
  geom_point(aes(size = Albumin), alpha = 0.7) +
  labs(title = "Scatter Plot of Bilirubin, Copper, Albumin, and Status",
       x = "Bilirubin",
       y = "Copper",
       color = "Status",
       size = "Albumin") +
  theme_minimal()

liver.cont = liver[continuous_vars]
# Scale continuous data
liver.scaled = scale(liver.cont)

head(liver.scaled)

elbow <- fviz_nbclust(liver.scaled, kmeans, method = "wss") +
   labs(subtitle = "Elbow Method")

sil <- fviz_nbclust(liver.scaled, kmeans, method = "silhouette") +
   labs(subtitle = "Silhouette Method")

gap <- fviz_nbclust(liver.scaled, kmeans, nstart = 25,
                  nboot= 50, method = "gap_stat") +
   labs(subtitle = "Gap Statistic Method")


grid.arrange(elbow, sil, gap, ncol = 2)

#grouping
set.seed(123456789)
group <- kmeans(liver.scaled[,c("Bilirubin", "Copper", "Albumin", "N_Days")],
                centers = 4, nstart = 25)

ggplot(liver, aes(x = Bilirubin, y = Copper, color = factor(group$cluster))) +
  geom_point(aes(size = Albumin), alpha = 0.8) +
  scale_size_continuous(range = c(2, 12)) +
  labs(
    title = "Bilirubin vs Copper",
    x = "Bilirubin",
    y = "Copper",
    size = "Albumin"
  )

# Create stacked bar plots for gender against death
plot1 <- ggplot(liver, aes(x = Ascites, fill = Status)) +
  geom_bar(position = "stack") +
  labs(title = "Stacked Bar Plot of Ascites and Status",
       x = "Ascites",
       fill = "Status") +
  theme_minimal()

plot2 <- ggplot(liver, aes(x = Hepatomegaly, fill = Status)) +
  geom_bar(position = "stack") +
  labs(title = "Stacked Bar Plot of Hepatomegaly and Status",
       x = "Hepatomegaly",
       fill = "Status") +
  theme_minimal()

plot3 <- ggplot(liver, aes(x = Spiders, fill = Status)) +
  geom_bar(position = "stack") +
  labs(title = "Stacked Bar Plot of Spiders and Status",
       x = "Spiders",
       fill = "Status") +
  theme_minimal()

plot4 <- ggplot(liver, aes(x = Edema, fill = Status)) +
  geom_bar(position = "stack") +
  labs(title = "Stacked Bar Plot of Edema and Status",
       x = "Edema",
       fill = "Status") +
  theme_minimal()

plot5 <- ggplot(liver, aes(x = Sex, fill = Status)) +
  geom_bar(position = "stack") +
  labs(title = "Stacked Bar Plot of Gender and Status",
       x = "Gender",
       fill = "Status") +
  theme_minimal()

# Modify Status to binary (1 for Death, 0 for others)
liver$Status <- ifelse(liver$Status == "Death", 1, 0)

# Fit the logistic regression model
model <- glm(Status ~ Bilirubin + Copper + Bilirubin + as.factor(Sex) + as.factor(Ascites) + as.factor(Spiders) + as.factor(Hepatomegaly),
             family = binomial, data = liver)

# Summary of the model to view coefficients and statistics
summary(model)



