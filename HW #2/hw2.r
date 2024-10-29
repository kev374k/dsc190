# Read in the data tables
survey_p1 <- read.table("videodata.txt", header = TRUE)
survey_p2 <- read.table("videoMultiple.txt", header = TRUE)

# Remove extra columns and check if numeric
survey_p2 <- survey_p2[, !(names(survey_p2) %in% c("other", "other2"))]
all(sapply(survey_p2, is.numeric))

# Combine the two data tables and remove na rows
combined_table <- cbind(survey_p1, survey_p2)
cleaned_df <- na.omit(combined_table)
dim(cleaned_df)

# Question 1
# Calculate point estimate of how many people played games in the week
n_played <- nrow(cleaned_df[cleaned_df$time > 0 | cleaned_df$time == 99, ])
n_total <- nrow(cleaned_df)

point_estimate_fraction <- n_played / n_total
print(paste("Point Estimate: ", round(point_estimate_fraction, 4)))

# Calculate the 95% confidence interval for the number of people playing games
confidence_level = 0.95
z <- qnorm((1 + confidence_level) / 2)

se <- sqrt(point_estimate_fraction * (1 - point_estimate_fraction) / n_total)

lower_interval_estimate_fraction <- point_estimate_fraction - z * se
upper_interval_estimate_fraction <- point_estimate_fraction + z * se
print(paste("95% Confidence Interval: [", round(lower_interval_estimate_fraction, 4), ", ", round(upper_interval_estimate_fraction, 4), "]", sep = ""))

# Question 2
# Create summary of the time column in the df
actual_values = cleaned_df[cleaned_df$time != 99 & cleaned_df$freq != 99, ]
avg_time <- tapply(actual_values$time, actual_values$freq, mean, na.rm = TRUE)
median_time <- tapply(actual_values$time, actual_values$freq, median, na.rm = TRUE)
count <- table(actual_values$freq)

summary <- data.frame(
  frequency = names(avg_time),
  avg_time = as.vector(avg_time),
  median_time = as.vector(median_time),
  count = as.vector(count)
)
summary

# Filter out students who didn't play during the exam week and recreate summary
exam_adjusted_df <- actual_values[actual_values$time != 0, ]
avg_time_adj <- tapply(exam_adjusted_df$time, exam_adjusted_df$freq, mean, na.rm = TRUE)
median_time_adj <- tapply(exam_adjusted_df$time, exam_adjusted_df$freq, median, na.rm = TRUE)
count_adj <- table(exam_adjusted_df$freq)

exam_adj_summary <- data.frame(
  frequency = names(avg_time_adj),
  avg_time = as.vector(avg_time_adj),
  median_time = as.vector(median_time_adj),
  count = as.vector(count_adj)
)

exam_adj_summary

# Create box plots of hours played during exam week by how often students normally play
par(mfrow = c(1, 2), # 1 row, 2 columns
    mar = c(4, 4, 3, 1)) # Adjust margins (bottom, left, top, right)

boxplot(time ~ freq,
        data = actual_values,
        xlab = "",
        ylab = "Hours Played",
        main = "Hours During Exam Week",
        col = "#7A2048",
        cex.main = 0.9,    # Reduce title size
        cex.lab = 0.8,     # Reduce label size
        cex.axis = 0.8)    # Reduce axis text size

boxplot(time ~ freq, 
        data = exam_adjusted_df,
        xlab = "",
        ylab = "Hours Played",
        main = "Overall Hours Played",
        col = "#408EC6",
        cex.main = 0.9,
        cex.lab = 0.8,
        cex.axis = 0.8)

mtext("Frequency: 1=Daily, 2=Weekly, 3=Monthly, 4=Semesterly", 
      side = 1, 
      line = -1, 
      outer = TRUE,
      cex = 0.7)

par(mfrow = c(1, 1))

# Question 3
# Get the point estimate by finding average time spent playing and calculate SE
point_estimate_average <- mean(cleaned_df$time, na.rm=TRUE)
N <- 314
n <- length(cleaned_df$time)

sample_sd <- sd(cleaned_df$time, na.rm = TRUE)
pop_correction_factor = sqrt((N - n) / (N - 1))
standard_error <- (sample_sd / sqrt(n)) * pop_correction_factor

print(paste("point_estimate_average:", point_estimate_average))
print(paste("Standard Error:", standard_error))

# Construct 95% Confidence Interval using point estimate and SE
alpha <- 0.05
z_value <- qnorm(1 - alpha / 2)

lower_interval_estimate_average <- point_estimate_average - z_value * standard_error
upper_interval_estimate_average <- point_estimate_average + z_value * standard_error

print(paste("lower_interval_estimate_average:", lower_interval_estimate_average, "hours"))
print(paste("upper_interval_estimate_average:", upper_interval_estimate_average, "hours"))

# Perform 1000 boostrapped samples to estimate sampling distribution
library(boot)

# Calculate the sample mean for each bootstrapped sample
bootstrap_mean <- function(data, indices) {
  return(mean(data[indices], na.rm = TRUE))
}

# Perform bootstrapping with 1000 iterations
bootstrap_results <- boot(cleaned_df$time, bootstrap_mean, R = 1000)
bootstrap_means <- bootstrap_results$t
bootstrap_mean <- mean(bootstrap_means)
bootstrap_sd <- sd(bootstrap_means)

print(paste("Mean of bootstrap distribution:", bootstrap_mean, "hours"))
print(paste("Standard Deviation of bootstrap distribution:", bootstrap_sd, "hours"))

# Plot histogram of bootstrap sampling distribution and plot CIs
bootstrap_ci <- boot.ci(bootstrap_results, type = "perc")

# Extract lower and upper bounds of the 95% Bootstrap CI
lower_bound_bootstrap <- bootstrap_ci$percent[4]
upper_bound_bootstrap <- bootstrap_ci$percent[5]

hist(bootstrap_means, main = "Bootstrap Sample Means with 95% CIs", 
     xlab = "Sample Mean Time Spent", 
     col = "lightblue", 
     border = "black", 
     probability = TRUE)

# Add vertical lines for the confidence intervals
abline(v = lower_bound_bootstrap, col = "red", lwd = 2, lty = 2)  # Lower CI bound
abline(v = upper_bound_bootstrap, col = "red", lwd = 2, lty = 2)  # Upper CI bound
abline(v = lower_interval_estimate_average, col = "blue", lwd = 2, lty = 1)
abline(v = upper_interval_estimate_average, col = "blue", lwd = 2, lty = 1)


legend("topright", legend = c("Lower Bootstrap CI", "Upper Bootstrap CI", "Lower 95% CI", "Upper 95% CI"),
       col = c("red", "red", "blue", "blue"), lty = c(2, 2, 1, 1), lwd = 2)

#Question 4
#Create bar plot of like variable
like_counts <- table(subset(cleaned_df, like != 99)$like)

# Plot the bar plot
like_counts_plot <- barplot(like_counts, 
                            main = "Counts of Responses to 'Like to Play' Video Games", 
                            xlab = "Response Category", 
                            ylab = "Count of Students", 
                            col = "lightblue", 
                            border = "black", 
                            names.arg = c("Never Played", "Very Much", "Somewhat", "Not Really", "Not At All"))

text(x = like_counts_plot, y = like_counts, label = like_counts, pos = 1, cex = 0.8, col = "red")

# Plot frequency of reasons for liking/disliking video games
reason_frequencies <- colSums(cleaned_df[, !names(cleaned_df) %in% c("action","adv","sim","sport","strategy","time","like","where","freq","busy","educ","sex","math","age","home","work","own","cdrom","email","grade", "time.1")])

bar_colors <- c("relax" = "lightblue",
                "coord" = "lightblue",
                "challenge" = "lightblue",
                "master" = "lightblue",
                "bored" = "lightblue",
                "graphic" = "lightblue",
                "time" = "red",
                "frust" = "red",
                "lonely" = "red",
                "rules" = "red",
                "cost" = "red",
                "boring" = "red",
                "friends" = "red",
                "point" = "red")

reason_plot <- barplot(reason_frequencies, 
                       main = "Frequencies of Reasons for Liking/Disliking Video Games",
                       ylab = "Frequency", 
                       col = bar_colors[names(reason_frequencies)], 
                       border = "black", 
                       las = 2)
text(x = reason_plot, y = reason_frequencies, label = reason_frequencies, pos = 1, cex = 0.8, col = "black")

# Question 5
# Plot the most relevant grouped bar plots showing differences in preferences for video games
clean_p2 <- cleaned_df

clean_p2$Entertainment <- clean_p2$graphic | clean_p2$relax | clean_p2$bored
clean_p2$Mental_Engagement <- clean_p2$coord | clean_p2$challenge | clean_p2$master
clean_p2$Difficulty <- clean_p2$frust | clean_p2$rules
clean_p2$Low_Priority <- clean_p2$time | clean_p2$cost | clean_p2$boring | clean_p2$point
clean_p2$Isolation <- clean_p2$lonely | clean_p2$friends

clean_p2 <- clean_p2[, c("sex", "work", "own", "Entertainment", "Mental_Engagement", "Difficulty", "Low_Priority", "Isolation")]
clean_p2[clean_p2 == 99] <- NA
clean_p2$work <- clean_p2$work > 0
clean_p2 <- na.omit(clean_p2)


# Create cross-tabs
sex_cross_tab <- table(clean_p2$Entertainment, clean_p2$sex)
dimnames(sex_cross_tab) <- list(
  "Entertainment" = c("No", "Yes"),
  "Sex" = c("F", "M")
)

work_cross_tab <- table(clean_p2$Mental_Engagement, clean_p2$work)
dimnames(work_cross_tab) <- list(
  "Mental Engagement" = c("No", "Yes"),
  "Works" = c("Doesn't Work", "Works")
)

own_cross_tab <- table(clean_p2$Entertainment, clean_p2$own)
dimnames(own_cross_tab) <- list(
  "Entertainment" = c("No", "Yes"),
  "Owns PC" = c("Lacks", "Owns")
)

# Set up plotting layout for side-by-side graphs
par(mfrow = c(1, 3), mar = c(4, 4, 4, 2) + 0.1)

# Color palette
cols <- c("#FF6B6B", "#4ECDC4")

# Plot 1: Sex vs Entertainment
bp1 <- barplot(sex_cross_tab, beside = TRUE,
               col = cols,
               main = "Gender vs Entertainment",
               xlab = "Gender",
               ylab = "Count",
               ylim = c(0, max(sex_cross_tab) * 1.2))
legend("topright", 
       legend = rownames(sex_cross_tab),
       fill = cols,
       title = "Entertainment",
       cex = 0.8)

# Plot 2: Work vs Mental Engagement
bp2 <- barplot(work_cross_tab, beside = TRUE,
               col = cols,
               main = "Employment vs\nMental Engagement",
               xlab = "Employment Status",
               ylab = "Count",
               ylim = c(0, max(work_cross_tab) * 1.2))
legend("topright",
       legend = rownames(work_cross_tab),
       fill = cols,
       title = "Mental\nEngagement",
       cex = 0.8)

# Plot 3: PC Ownership vs Entertainment
bp3 <- barplot(own_cross_tab, beside = TRUE,
               col = cols,
               main = "PC Ownership vs\nEntertainment",
               xlab = "PC Ownership",
               ylab = "Count",
               ylim = c(0, max(own_cross_tab) * 1.2))
legend("topright",
       legend = rownames(own_cross_tab),
       fill = cols,
       title = "Entertainment",
       cex = 0.8)

# Reset layout
par(mfrow = c(1, 1))

# Question 6
# Perform chi-squared test to compare student expected grade distribution with target distribution
grade_factors <- factor(cleaned_df$grade, levels=0:4)
observed_counts <- table(grade_factors)
names(observed_counts) <- c("F", "D", "C", "B", "A")
observed_counts_combined <- c(
  D_or_lower = observed_counts["F"] + observed_counts["D"],
  C = observed_counts["C"],
  B = observed_counts["B"],
  A = observed_counts["A"]
)
print(observed_counts_combined)

expected_proportions <- c(0.1, 0.4, 0.3, 0.2)

chisq_test <- chisq.test(x = observed_counts_combined, p = expected_proportions)

chisq_test

# Perform chi-squared test again after imputing failing grades for nonrespondents
combined_table$grade <- ifelse(is.na(combined_table$relax), 0, combined_table$grade)

grade_factors <- factor(combined_table$grade, levels = 0:4)
observed_counts <- table(grade_factors)
names(observed_counts) <- c("F", "D", "C", "B", "A")

observed_counts_combined <- c(
  D_or_lower = observed_counts["F"] + observed_counts["D"],
  C = observed_counts["C"],
  B = observed_counts["B"],
  A = observed_counts["A"]
)

print(observed_counts_combined)

expected_proportions <- c(0.1, 0.4, 0.3, 0.2)

chisq_test <- chisq.test(x = observed_counts_combined, p = expected_proportions)
chisq_test

# Advanced Analysis
# Filter df based on type of video game preferred
video_games_data <- cleaned_df[!(cleaned_df$action == 0 & cleaned_df$adv == 0 & cleaned_df$sim == 0 & cleaned_df$sport == 0 & cleaned_df$strategy == 0), ]
dim(video_games_data)

freq_table <- colSums(video_games_data[, c("action", "adv", "sim", "sport", "strategy")] == 1)
action_time <- sum(cleaned_df$time[cleaned_df$action == 1], na.rm = TRUE)
adv_time <- sum(cleaned_df$time[cleaned_df$adv == 1], na.rm = TRUE)
sim_time <- sum(cleaned_df$time[cleaned_df$sim == 1], na.rm = TRUE)
sport_time <- sum(cleaned_df$time[cleaned_df$sport == 1], na.rm = TRUE)
strategy_time <- sum(cleaned_df$time[cleaned_df$strategy == 1], na.rm = TRUE)

summed_times <- c(
  action_time,
  adv_time,
  sim_time,
  sport_time,
  strategy_time
)
names(summed_times) <- c("action", "adv", "sim", "sport", "strategy")

# Set up the plotting layout for side by side graphs
par(mfrow = c(1, 2), 
    mar = c(6, 4, 3, 1))  # Increased bottom margin for rotated labels

# First barplot
barplot(freq_table, 
        main = "Genre Popularity",    # Shortened title
        xlab = "Genre",                    # Removed xlabel since labels are rotated
        ylab = "Frequency", 
        col = "#B85042",
        las = 2,                      # Rotate x-axis labels
        cex.main = 0.9,              # Reduce title size
        cex.axis = 0.8,              # Reduce axis text size
        cex.names = 0.8)             # Reduce bar labels size

# Second barplot
barplot(summed_times,
        main = "Time by Category",    # Shortened title
        xlab = "Genre",                    # Removed xlabel since labels are rotated
        ylab = "Total Time", 
        col = "#A7BEAE",
        las = 2,
        cex.main = 0.9,
        cex.axis = 0.8,
        cex.names = 0.8)

# Reset par to default after plotting
par(mfrow = c(1, 1))