# Read in the data from the txt file
data <- read.table("babies.txt", header = TRUE)

# View the first 5 rows
head(data, 5)

# 1. Data Processing and Summaries
# Print a description table to summarize each variable
bwt_description <- summary(data$bwt)
bwt_description

gestation_description <- summary(data$gestation)
gestation_description

parity_description <- summary(data$parity)
parity_description

age_description <- summary(data$age)
age_description

height_description <- summary(data$height)
height_description

weight_description <- summary(data$weight)
weight_description

smoke_description <- summary(data$smoke)
smoke_description

# View data types of variables
str(data)

# View histograms of the variables
hist(data$bwt, main = "Birthweight", xlab = "Birthweight (oz)", col = "lightblue")
hist(data$gestation, main = "Gestation", xlab = "Gestation (days)", breaks=50, col = "red")
hist(data$parity, main = "Parity", xlab = "Parity (0 = First Pregnancy)", col = "lightgreen")
hist(data$age, main = "Mother's Age", xlab = "Age (years)", col = "lightgrey")
hist(data$height, main = "Mother's Height", xlab = "Height (inches)", breaks=50, col = "lightpink")
hist(data$weight, main = "Mother's Weight", xlab = "Weight (lbs)", breaks=50, col = "lightgreen")
hist(data$smoke, main = "Smoking Status", xlab = "Smoker (0 = No)", col = "lightblue")

# Run a hypothesis test to compare the observed skew in gestation to simulated skews
skew = function(X) {
  Xbar = mean(X); S = sd(X)
  Z = (X - Xbar)/S
  mean(Z^3)
}

mu = 5; sigma = 2
n = nrow(data)
sample_skews = replicate(1000, {
  X = rnorm(n, mean = mu, sd = sigma)
  skew(X)
})

gestation_obs_skew = skew(data$gestation)
gestation_p_value = mean(abs(sample_skews) >= abs(gestation_obs_skew))
cat("Observed skew for gestation: ", gestation_obs_skew, "\n",
    "The p-value for the gestation hypothesis test is: ", gestation_p_value, "\n", sep="")

# Repeat hypothesis test for age
age_obs_skew = skew(data$age)
age_p_value = mean(abs(sample_skews) >= abs(age_obs_skew))
cat("Observed skew for age: ", age_obs_skew, "\n", 
    "The p-value for the age hypothesis test is: ", age_p_value, "\n", sep = "")

# Repeat hypothesis test for height
height_obs_skew = skew(data$height)
height_p_value = mean(abs(sample_skews) >= abs(height_obs_skew))
cat("Observed skew for height: ", height_obs_skew, "\n",
    "The p-value for the height hypothesis test is: ", height_p_value, "\n", sep = "")

# Repeat hypothesis test for weight
weight_obs_skew = skew(data$weight)
weight_p_value = mean(abs(sample_skews) >= abs(weight_obs_skew))
cat("Observed skew for weight: ", weight_obs_skew, "\n",
    "The p-value for the weight hypothesis test is: ", weight_p_value, "\n", sep = "")

# clean the dataframe after determining outliers
cleaned_df <- data[data$gestation < 500 & data$age < 50 & data$height < 80 & data$weight < 500 & data$smoke <= 1, ]
head(cleaned_df, 5)

# 2. Numeric Summary of Birthweights Between Smokers and Non-Smokers
# Get smoker and non-smoker df's
smokers <- cleaned_df[cleaned_df$smoke == 1, ]
non_smokers <- cleaned_df[cleaned_df$smoke == 0, ]

# Print the min, max, mean, median, and std. dev. for both groups
stats <- data.frame(
  Statistic = c("Minimum", "Maximum", "Mean", "Median", "Standard Deviation"),
  Smokers = round(c(min(smokers$bwt), max(smokers$bwt), mean(smokers$bwt), median(smokers$bwt), sd(smokers$bwt)), 2),
  NonSmokers = round(c(min(non_smokers$bwt), max(non_smokers$bwt), mean(non_smokers$bwt), median(non_smokers$bwt), sd(non_smokers$bwt)), 2)
)

print(stats, row.names = FALSE)

# Print the quartiles for both groups
smoker_quartiles <- quantile(smokers$bwt)
nonsmoker_quartiles <- quantile(non_smokers$bwt)

smoker_q1_bwt <- smoker_quartiles[2]
smoker_q2_bwt <- smoker_quartiles[3]
smoker_q3_bwt <- smoker_quartiles[4]
print("Smoker Mothers:")
cat("Q1:", smoker_quartiles[2], "Q2:", smoker_quartiles[3], "Q3:", smoker_quartiles[4], "\n\n")

nonsmoker_q1_bwt <- nonsmoker_quartiles[2]
nonsmoker_q2_bwt <- nonsmoker_quartiles[3]
nonsmoker_q3_bwt <- nonsmoker_quartiles[4]

print("Non-Smoker Mothers:")
cat("Q1:", nonsmoker_quartiles[2], "Q2:", nonsmoker_quartiles[3], "Q3:", nonsmoker_quartiles[4], "\n")

# 3. Graphically Summarize Birthweights Between Smokers and Non-Smokers
# Produce boxplot for both distributions
boxplot(bwt ~ smoke, 
        data = cleaned_df, 
        main = "Birthweight by Smoking Status",
        xlab = "Smoking Status (0 = Non-Smoker, 1 = Smoker)",
        ylab = "Birthweight (oz)",
        col = c("lightblue", "lightcoral"),
        names = c("Non-Smoker", "Smoker"))

# Produce density plot for both distributions
library(lattice)
densityplot(~ bwt, group = smoke, data = cleaned_df, plot.points = FALSE,
            auto.key = list(
              columns = 2, 
              title = "Smokes?", 
              space = "right",
              text = c("No","Yes")),
            main = "Birthweight by Smoking Status",
            xlab = "Birthweight (oz)",
            col = c("blue", "red"))

# 4. Incidence Comparison of Low-Birth-Weight Between Smokers and Non-Smokers
# Get percentage of birthweights under 100 oz for both groups
low_bwt_smoker <- nrow(smokers[smokers$bwt < 100, ] ) / nrow(smokers) * 100
cat('low_bwt_smoker:', low_bwt_smoker, "%")
low_bwt_nonsmoker <- nrow(non_smokers[non_smokers$bwt < 100, ] ) / nrow(non_smokers) * 100
cat('low_bwt_nonsmoker:', low_bwt_nonsmoker, "%")

# Generate histograms for both groups to inspect how proportions change with threshold
hist(non_smokers$bwt, main = "Non-Smoker Birthweights", xlab = "Birthweight (oz)", breaks=20)
hist(smokers$bwt, main = "Smoker Birthweights", xlab = "Birthweight (oz)", breaks=20)

# Test different thresholds to see how percentages change
cat("Babies under 120 ounces:\n",
    "Smokers: ", formatC(nrow(smokers[smokers$bwt < 120, ]) / nrow(smokers), format="f", digits=4), "\n",
    "Non-smokers: ", formatC(nrow(non_smokers[non_smokers$bwt < 120, ]) / nrow(non_smokers), format="f", digits=4), "\n\n",
    "Babies under 80 ounces:\n",
    "Smokers: ", formatC(nrow(smokers[smokers$bwt < 80, ]) / nrow(smokers), format="f", digits=4), "\n",
    "Non-smokers: ", formatC(nrow(non_smokers[non_smokers$bwt < 80, ]) / nrow(non_smokers), format="f", digits=4), "\n\n",
    "Babies under 60 ounces:\n",
    "Smokers: ", formatC(nrow(smokers[smokers$bwt < 60, ]) / nrow(smokers), format="f", digits=4), "\n",
    "Non-smokers: ", formatC(nrow(non_smokers[non_smokers$bwt < 60, ]) / nrow(non_smokers), format="f", digits=4), "\n\n",
    "Babies under 58 ounces:\n",
    "Smokers: ", formatC(nrow(smokers[smokers$bwt < 58, ]) / nrow(smokers), format="f", digits=4), "\n",
    "Non-smokers: ", formatC(nrow(non_smokers[non_smokers$bwt < 58, ]) / nrow(non_smokers), format="f", digits=4),
    sep="")

# Advanced Analysis: Correlation between birthweight and gestation period
# Plot scatterplot and regression line to show relationship
# between birthweight and gestation period
plot(cleaned_df$gestation, cleaned_df$bwt, 
     col = ifelse(cleaned_df$smoke == 1, "red", "blue"),
     pch = 19, 
     xlab = "Gestation (days)", 
     ylab = "Birth Weight (oz)",
     main = "Birth Weight vs. Gestation by Smoking Status",
     cex = 0.5)

smoke_1 <- subset(cleaned_df, smoke == 1)
smoke_0 <- subset(cleaned_df, smoke == 0)

abline(lm(bwt ~ gestation, data = smoke_1), col = "red", lwd = 2)
abline(lm(bwt ~ gestation, data = smoke_0), col = "blue", lwd = 2)


legend("topleft", legend = c("Smokers", "Non-smokers"), col = c("red", "blue"), pch = 19)
