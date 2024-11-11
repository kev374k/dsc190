# Read in the data
data <- read.table("hcmv.txt", header=TRUE)

# Question 1
# Use a random sample from a number range representing the DNA sequence
# to generate simulated palindromes
set.seed(236)

generate_palindrome_sites <- function() {
  sort(sample(1:229354, 296, replace=FALSE))
}

simulated_palindrome_data <- replicate(1000, generate_palindrome_sites())

calc_summary_stats <- function(locations) {
  list(
    mean = mean(locations),
    median = median(locations),
    variance = var(locations),
    std_dev = sd(locations)
  )
}

simulated_stats <- apply(simulated_palindrome_data, 2, calc_summary_stats)
real_stats <- calc_summary_stats(data$location)

mean_simulated_means <- mean(sapply(simulated_stats, function(x) x$mean))
mean_simulated_medians <- mean(sapply(simulated_stats, function(x) x$median))
mean_simulated_variances <- mean(sapply(simulated_stats, function(x) x$variance))
mean_simulated_std_devs <- mean(sapply(simulated_stats, function(x) x$std_dev))

comparison_table <- data.frame(
  Statistic = c("Mean", "Median", "Variance", "Standard Deviation"),
  Real_Locations = c(real_stats$mean, real_stats$median, real_stats$variance, real_stats$std_dev),
  Simulated_Locations = c(mean_simulated_means, mean_simulated_medians, mean_simulated_variances, mean_simulated_std_devs)
)

print(comparison_table)


# Plot the distribution of simulated means and std devs and compare to the real one
par(mfrow = c(1, 2))

simulated_means <- sapply(simulated_stats, function(x) x$mean)
hist(simulated_means, breaks = 30, col = "#CEE6F2", main = "Simulated Mean Locations", xlab = "Mean of Palindrome Locations", cex.main = 0.8, cex.lab = 0.8)
abline(v = real_stats$mean, col = "red", lwd = 2)
legend("topright", legend = "Real Data Mean", col = "red", lwd = 2, cex=0.5)

simulated_std_devs <- sapply(simulated_stats, function(x) x$std_dev)
hist(simulated_std_devs, breaks = 30, col = "#A1BE95", main = "Simulated Standard Deviations", xlab = "Standard Deviation of Palindrome Locations", cex.main = 0.8, cex.lab = 0.8)
abline(v = real_stats$std_dev, col = "red", lwd = 2)
legend("topright", legend = "Real Data Std Dev", col = "red", lwd = 2, cex=0.5)

par(mfrow = c(1, 1))

# Qustion 2
# Plot spacing distribution for both simulated and real data
options(scipen = 900)
par(mfrow = c(1, 2))

differences <- diff(data$location)
spacing_data <- apply(replicate(100, generate_palindrome_sites())
                      , 2, function(x) diff(x))

xlim <- range(c(differences, spacing_data), na.rm = TRUE)
ymax <- max(c(
  hist(differences, breaks = 30, plot = FALSE)$density,
  hist(spacing_data, breaks = 30, plot = FALSE)$density
), na.rm = TRUE)

hist(differences, breaks = seq(xlim[1], xlim[2], length.out = 30),
     main = "Spacing Between Palindromes in CMV DNA",
     xlab = "Spacing (in base pairs)",
     ylab = "Density",
     col = "#89ABE3",
     freq = FALSE,
     xlim = xlim,
     ylim = c(0, ymax),
     cex.main = 0.7,
     cex.axis = 0.6)

hist(spacing_data, breaks = seq(xlim[1], xlim[2], length.out = 30),
     main = "Simulated Spacing Between Palindromes", 
     xlab = "Spacing (in base pairs)", 
     ylab = "Density", 
     col = "#FFA07A", 
     freq = FALSE, 
     xlim = xlim, 
     ylim = c(0, ymax), 
     cex.main = 0.7, 
     cex.axis = 0.6)

par(mfrow = c(1, 1)) 

# Plot distribution of pair sums for both real and simulate data
par(mfrow = c(1, 2))  

pair_sums <- (data[-length(data$location)] + data$location[-1])$location

simulated_pair_sums <- apply(replicate(100, generate_palindrome_sites()), 2, function(x) x[-length(x)] + x[-1])  # Adjust based on the simulation structure

xlim <- range(c(pair_sums, simulated_pair_sums), na.rm = TRUE)
ymax <- max(
  hist(pair_sums, breaks = 30, plot = FALSE)$density
)

hist(pair_sums, breaks = seq(xlim[1], xlim[2], length.out = 30), 
     main = "Pair Sums of Palindromes in CMV DNA", 
     xlab = "Pair Sum Total", 
     ylab = "Density", 
     col = "#B85042", 
     freq = FALSE, 
     xlim = xlim, 
     ylim = c(0, ymax), 
     cex.main = 0.8, 
     cex.axis = 0.6)

hist(simulated_pair_sums, breaks = seq(xlim[1], xlim[2], length.out = 30), 
     main = "Simulated Pair Sums of Palindromes", 
     xlab = "Pair Sum Total", 
     ylab = "Density", 
     col = "#A7BEAE", 
     freq = FALSE, 
     xlim = xlim, 
     ylim = c(0, ymax), 
     cex.main = 0.8, 
     cex.axis = 0.6)

par(mfrow = c(1, 1))

# Plot triplet sum distributions for both real and simulated data

options(scipen = 900) 
par(mfrow = c(1, 2))  

# Calculate the triplet sums
palindrome_length <- length(data$location)
triplet_sums <- data$location[1:(palindrome_length - 2)] +
  data$location[2:(palindrome_length - 1)] +
  data$location[3:palindrome_length]

# Simulated triplet sums
simulated_triplet_sums <- apply(replicate(100, generate_palindrome_sites()), 2, 
                                function(x) {
                                  return(x[1:(length(x) - 2)] + 
                                           x[2:(length(x) - 1)] + 
                                           x[3:length(x)])
                                })

# Define limits for x-axis and y-axis
xlim <- range(c(triplet_sums, na.omit(simulated_triplet_sums)), na.rm = TRUE)
ymax <- max(hist(triplet_sums, breaks = 30, plot = FALSE)$density)


# Plot the histogram for triplet sums
hist(triplet_sums, breaks = seq(xlim[1], xlim[2], length.out = 30), 
     main = "Triplet Sums of Palindromes in DNA", 
     xlab = "Triplet Sum Total", 
     ylab = "Density", 
     freq = FALSE, 
     col = "#735DA5", 
     xlim = xlim, 
     ylim = c(0, ymax),
     cex.main = 0.8, 
     cex.axis = 0.5)

# Plot the histogram for simulated triplet sums
hist(simulated_triplet_sums, breaks = seq(xlim[1], xlim[2], length.out = 30), 
     main = "Simulated Triplet Sums of Palindromes", 
     xlab = "Triplet Sum Total", 
     ylab = "Density", 
     freq = FALSE, 
     col = "#D3C5E5", 
     xlim = xlim, 
     ylim = c(0, ymax), 
     cex.main = 0.8, 
     cex.axis = 0.5)

# Reset layout
par(mfrow = c(1, 1)) 


# Question 3
# Split the data into 57 equally sized intervals and count the number of intervals
# that exhibit each palindrome count
# Additionally, calculate the expected counts from the Poisson distribution

interval_size <- floor(229354 / 57)
interval_counts <- hist(data$location, breaks = seq(0, 229354, by = interval_size), plot = FALSE)$counts

interval_counts_df <- data.frame(
  Interval = 1:57,
  Count = interval_counts
)

lambda_rate <- 296 / 57

count_frequency_df <- as.data.frame(table(interval_counts_df$Count))
colnames(count_frequency_df) <- c("Palindrome_Count", "Intervals_Observed")
count_frequency_df$Palindrome_Count <- as.integer(as.character(count_frequency_df$Palindrome_Count))
count_frequency_df$Intervals_Expected <- dpois(count_frequency_df$Palindrome_Count, lambda_rate) * 57
count_frequency_df$Intervals_Expected <- count_frequency_df$Intervals_Expected * (57 / sum(count_frequency_df$Intervals_Expected))

count_below_10 <- count_frequency_df[count_frequency_df$Palindrome_Count < 10, ]
count_10_or_more <- count_frequency_df[count_frequency_df$Palindrome_Count >= 10, ]

count_10_or_more_agg <- data.frame(
  Palindrome_Count = "10+",
  Intervals_Observed = sum(count_10_or_more$Intervals_Observed),
  Intervals_Expected = sum(count_10_or_more$Intervals_Expected)
)

count_frequency_df <- rbind(count_below_10, count_10_or_more_agg)
count_frequency_df$Palindrome_Count <- as.character(count_frequency_df$Palindrome_Count)

# Plot the distribution of palindrome counts in intervals for both real and expected
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
barplot(
  count_frequency_df$Intervals_Observed,
  names.arg = count_frequency_df$Palindrome_Count,
  col = "blue",
  xlab = "Palindrome Count",
  ylab = "Number of Intervals",
  main = "Observed Frequency of Palindrome Counts",
  cex.main = 0.8
)

barplot(
  count_frequency_df$Intervals_Expected,
  names.arg = count_frequency_df$Palindrome_Count,
  col = "red",
  xlab = "Palindrome Count",
  ylab = "Number of Intervals",
  main = "Expected Frequency of Palindrome Counts (Poisson)",
  cex.main = 0.7
)

mtext("Interval Size = 4,023", outer = TRUE, cex = 1.5)

par(mfrow = c(1, 1))

# Conduct a chi-squared goodness-of-fit test to assess whether the expected and real
# data come from different distributions

chi_square_test <- chisq.test(
  x = count_frequency_df$Intervals_Observed,
  p = count_frequency_df$Intervals_Expected / sum(count_frequency_df$Intervals_Expected),
  rescale.p = TRUE
)

chi_square_test

# Visualize palindrome count distributions when using much larger intervals
larger_interval_size <- floor(229354 / 10)
larger_interval_counts <- hist(data$location, breaks = seq(0, 229354, by = larger_interval_size), plot = FALSE)$counts

larger_interval_counts_df <- data.frame(
  Interval = 1:10,
  Count = larger_interval_counts
)

lambda_rate <- 296 / 10

larger_count_frequency_df <- as.data.frame(table(larger_interval_counts_df$Count))
colnames(larger_count_frequency_df) <- c("Palindrome_Count", "Intervals_Observed")
larger_count_frequency_df$Palindrome_Count <- as.integer(as.character(larger_count_frequency_df$Palindrome_Count))
larger_count_frequency_df$Intervals_Expected <- dpois(larger_count_frequency_df$Palindrome_Count, lambda_rate) * 10
larger_count_frequency_df$Intervals_Expected <- larger_count_frequency_df$Intervals_Expected * (10 / sum(larger_count_frequency_df$Intervals_Expected))


par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
barplot(
  larger_count_frequency_df$Intervals_Observed,
  names.arg = larger_count_frequency_df$Palindrome_Count,
  col = "blue",
  xlab = "Palindrome Count",
  ylab = "Number of Intervals",
  main = "Observed Frequency of Palindrome Counts",
  cex.main = 0.8
)

barplot(
  larger_count_frequency_df$Intervals_Expected,
  names.arg = larger_count_frequency_df$Palindrome_Count,
  col = "red",
  xlab = "Palindrome Count",
  ylab = "Number of Intervals",
  main = "Expected Frequency of Palindrome Counts (Poisson)",
  cex.main = 0.7
)

mtext("Interval Size = 22,935", outer = TRUE, cex = 1.5)

par(mfrow = c(1, 1))

# Visualize distributions of palindrome counts with much smaller intervals
smaller_interval_size <- floor(229354 / 573)
smaller_interval_counts <- hist(data$location, breaks = seq(0, 229354, by = smaller_interval_size), plot = FALSE)$counts

smaller_interval_counts_df <- data.frame(
  Interval = 1:573,
  Count = smaller_interval_counts
)

lambda_rate <- 296 / 573

smaller_count_frequency_df <- as.data.frame(table(smaller_interval_counts_df$Count))
colnames(smaller_count_frequency_df) <- c("Palindrome_Count", "Intervals_Observed")
smaller_count_frequency_df$Palindrome_Count <- as.integer(as.character(smaller_count_frequency_df$Palindrome_Count))
smaller_count_frequency_df$Intervals_Expected <- dpois(smaller_count_frequency_df$Palindrome_Count, lambda_rate) * 573
smaller_count_frequency_df$Intervals_Expected <- smaller_count_frequency_df$Intervals_Expected * (573 / sum(smaller_count_frequency_df$Intervals_Expected))


par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
barplot(
  smaller_count_frequency_df$Intervals_Observed,
  names.arg = smaller_count_frequency_df$Palindrome_Count,
  col = "blue",
  xlab = "Palindrome Count",
  ylab = "Number of Intervals",
  main = "Observed Frequency of Palindrome Counts",
  cex.main = 0.8
)

barplot(
  smaller_count_frequency_df$Intervals_Expected,
  names.arg = smaller_count_frequency_df$Palindrome_Count,
  col = "red",
  xlab = "Palindrome Count",
  ylab = "Number of Intervals",
  main = "Expected Frequency of Palindrome Counts (Poisson)",
  cex.main = 0.7
)

mtext("Interval Size = 400", outer = TRUE, cex = 1.5)

par(mfrow = c(1, 1))

# Conduct chi-squared test again, but with smaller intervals to see if results are similar and valid
smaller_chi_square_test <- chisq.test(
  x = smaller_count_frequency_df$Intervals_Observed,
  p = smaller_count_frequency_df$Intervals_Expected / sum(smaller_count_frequency_df$Intervals_Expected),
  rescale.p = TRUE
)

smaller_chi_square_test


# Question 4
# Determine ideal interval size and plot each interval's palindrome counts to
# find the largest palindrome clusters
library(ggplot2)

# Set parameters
total_length <- 229354
interval_size <- 2500  # Adjust interval size as needed

# Create intervals
intervals <- seq(1, total_length, by = interval_size)
interval_counts <- table(cut(data$location, breaks = intervals))

# Convert to a data frame for plotting
interval_data <- data.frame(
  interval = intervals[-length(intervals)],
  counts = as.numeric(interval_counts)
)

# Create a ggplot with updated aesthetics
ggplot(interval_data, aes(x = interval, y = counts)) +
  geom_bar(stat = "identity", fill = "#408EC6", color = "white", linewidth = 0.2) +  # Use linewidth instead of size
  labs(title = "Distribution of Palindromic Sequences in CMV DNA",
       x = "Interval (Base Pairs)",
       y = "Number of Palindromic Sequences") +
  theme_minimal(base_size = 15) +  # Use a minimal theme with larger base font size
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Center title and make it bold
    axis.title.x = element_text(size = 14, face = "bold"),  # Bold x-axis title
    axis.title.y = element_text(size = 14, face = "bold"),  # Bold y-axis title
    axis.text = element_text(size = 12)  # Increase axis text size
  )

# Conduct Poisson tests to see if identified clusters have significantly more palindromes
total_intervals <- length(intervals)  
total_palindromes <- sum(interval_data$counts)  

lambda <- total_palindromes / total_intervals

observed_count_92500_95000 <- 13  

p_value_92500_95000 <- ppois(observed_count_92500_95000 - 1, lambda, lower.tail = FALSE)

observed_count_195000_197500 <- interval_data$counts[interval_data$interval == 195001]  

p_value_195000_197500 <- ppois(observed_count_195000_197500 - 1, lambda, lower.tail = FALSE)

cat("P-value for interval 92500–95000:", p_value_92500_95000, "\n")
cat("P-value for interval 195000–197500:", p_value_195000_197500, "\n")


# Advanced Analysis
# Conduct Regression with a quadratic function predicting palindrome spacing from location
spacing <- diff(data$location)

location_midpoints <- (data$location[-1] + data$location[-length(data$location)]) / 2

spacing_df <- data.frame(Location = location_midpoints, Spacing = spacing)
spacing_df$Location_Squared <- spacing_df$Location^2

quadratic_model <- lm(Spacing ~ Location + Location_Squared, data = spacing_df)
summary(quadratic_model)

plot(spacing_df$Location, spacing_df$Spacing, type = "p", col = "blue",
     xlab = "Location on DNA Sequence", ylab = "Spacing Between Consecutive Palindromes",
     main = "Spacing Between Consecutive Palindromes Across DNA Sequence")

location_seq <- seq(min(spacing_df$Location), max(spacing_df$Location), length.out = 100)
predicted_spacing <- predict(quadratic_model, newdata = data.frame(Location = location_seq, Location_Squared = location_seq^2))
lines(location_seq, predicted_spacing, col = "red", lwd = 2)