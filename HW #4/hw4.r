# Read in the data
data <- read.table("gauge.txt", header=TRUE)

# Question 1: Fit the Linear Model from Data

par(mfrow = c(1, 2))

model <- lm(data$gain ~ data$density)
plot(data$density, data$gain, 
     main = "Snow Gain vs Density", 
     xlab = "Snow Density (g/cm³)", 
     ylab = "Gamma Ray Intensity", 
     pch = 19, 
     col = "#89ABE3")

abline(model, col = "#EA738D", lwd = 2)
predicted_gain <- predict(model)
segments(data$density, data$gain, data$density, predicted_gain, col = "purple", lwd = 1, lty = 2)

residuals <- resid(model)
fitted_values <- fitted(model)

plot(fitted_values, residuals, 
     main = "Residual Plot", 
     xlab = "Fitted Values", 
     ylab = "Residuals", 
     pch = 19, 
     col = "#89ABE3")

# Add a horizontal line at zero to indicate no residuals
abline(h = 0, col = "#EA738D", lwd = 2, lty = 2)

r2 <- round(summary(model)$r.squared, 3)

# Present it in text
sprintf("The R² value is %0.3f", r2)

# Question 2: Transform Data and Fit
# Use the log method to apply natural logarithm to gain column before fitting model
par(mfrow = c(1, 2))

data$log_gain <- log(data$gain)

log_model <- lm(log_gain ~ density, data=data)

plot(data$density, data$log_gain,
     main = "ln(Gain) vs Density", 
     xlab = "Snow Density (g/cm³)", 
     ylab = "Natural Log of Gain", 
     pch = 19, 
     col = "#89ABE3")

abline(log_model, col = "#EA738D", lwd = 2)
predicted_log_gain <- predict(log_model)
segments(data$density, data$log_gain, data$density, predicted_log_gain, col = "purple", lwd = 1, lty = 2)

residuals <- resid(log_model)
fitted_values <- fitted(log_model)

plot(fitted_values, residuals, 
     main = "Residual Plot", 
     xlab = "Fitted Values", 
     ylab = "Residuals", 
     pch = 19, 
     col = "#89ABE3")

# Add a horizontal line at zero to indicate no residuals
abline(h = 0, col = "#EA738D", lwd = 2, lty = 2)

# Calculate R^2 of log regression model
cat("Logarithm Regression Results")
r2 <- round(summary(log_model)$r.squared, 3)

# Present it in text
sprintf("The R² value for the log_model is %0.3f", r2)

# Compare to previous regression R^2
cat("Regression Results Before Transformation")
r2 <- round(summary(model)$r.squared, 3)

sprintf("The R² value for our regreesion model is %0.3f", r2)

# Question 3: Effect of Random Reporting Error
# Use a normal distribution and sample errors from it to alter the density in order to
# Simulate reporting errors
# Use the expected/predicted gain of densities before they are altered

set.seed(49)

n_sims <- 500
n <- nrow(data)

observed_mse <- mean(resid(log_model)^2)
small_error_sim_mse <- numeric(n_sims)

for (i in 1:n_sims) {
  sim_density <- data$density + rnorm(n, mean=0, sd=0.008)
  sim_data <- data.frame(log_gain=predicted_log_gain, density=sim_density)
  sim_model <- lm(log_gain ~ density, data=sim_data)
  small_error_sim_mse[i] = mean(resid(sim_model)^2)
}

hist(small_error_sim_mse, breaks = 30, col = "skyblue", border = "white",
     main = "Histogram of Simulated Regression MSEs\nAverage Reporting Error = 0.008",
     xlab = "MSE", xlim = range(c(small_error_sim_mse, observed_mse)))
abline(v = observed_mse, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Actual Regression's MSE"), col = c("red"), lwd = 2, lty = 2)

sim_mse <- numeric(n_sims)

for (i in 1:n_sims) {
  sim_density <- data$density + rnorm(n, mean=0, sd=0.015)
  sim_data <- data.frame(log_gain=predicted_log_gain, density=sim_density)
  sim_model <- lm(log_gain ~ density, data=sim_data)
  sim_mse[i] = mean(resid(sim_model)^2)
}

# Produce histograms of resulting regression MSEs after simulating error

hist(sim_mse, breaks = 30, col = "skyblue", border = "white",
     main = "Histogram of Simulated Regression MSEs\nAverage Reporting Error = 0.015",
     xlab = "MSE", xlim = range(c(sim_mse, observed_mse)))
abline(v = observed_mse, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Actual Regression's MSE"), col = c("red"), lwd = 2, lty = 2)

large_error_sim_mse <- numeric(n_sims)

for (i in 1:n_sims) {
  sim_density <- data$density + rnorm(n, mean=0, sd=0.023)
  sim_data <- data.frame(log_gain=predicted_log_gain, density=sim_density)
  sim_model <- lm(log_gain ~ density, data=sim_data)
  large_error_sim_mse[i] = mean(resid(sim_model)^2)
}

hist(large_error_sim_mse, breaks = 30, col = "skyblue", border = "white",
     main = "Histogram of Simulated Regression MSEs\nAverage Reporting Error = 0.023",
     xlab = "MSE", xlim = range(c(large_error_sim_mse, observed_mse)))
abline(v = observed_mse, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Actual Regression's MSE"), col = c("red"), lwd = 2, lty = 2)


# Question 4: Use the log model and un-transform the predicted log(gain) into normal,
# unscaled gain. Then produce interval estimates and plot these estimates against
# the data.
library(ggplot2)
library(knitr)

density_grid <- seq(min(data$density), max(data$density), length.out = 100)

pred <- predict(log_model, 
                newdata = data.frame(density = density_grid), 
                se.fit = TRUE)

# Calculate prediction interval on log scale using 95% CI
sigma <- summary(log_model)$sigma
t_value <- qt(0.975, df = length(data$gain) - 2)
pred_se <- sqrt(pred$se.fit^2 + sigma^2)

# Create prediction intervals on log scale
log_lower <- pred$fit - t_value * pred_se
log_upper <- pred$fit + t_value * pred_se

# Back-transform to original scale (use e^x since ln was used earlier)
pred_df <- data.frame(
  density = density_grid,
  fit = exp(pred$fit),
  lower = exp(log_lower),
  upper = exp(log_upper)
)

# creates the 95% CI for the uncertainty bands for each point
get_uncertainty_band <- function(density_value) {
  pred_specific <- predict(log_model, 
                           newdata = data.frame(density = density_value), 
                           se.fit = TRUE)
  pred_se_specific <- sqrt(pred_specific$se.fit^2 + sigma^2)
  log_interval <- c(
    pred_specific$fit - t_value * pred_se_specific,
    pred_specific$fit,
    pred_specific$fit + t_value * pred_se_specific
  )
  interval <- exp(log_interval)
  return(interval)
}

# Create visualization
ggplot() +
  geom_point(data = data, aes(x = density, y = gain), alpha = 0.5) +
  geom_line(data = pred_df, aes(x = density, y = fit), color = "blue") +
  geom_ribbon(data = pred_df, 
              aes(x = density, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "blue") +
  labs(x = "Density (g/cm³)", 
       y = "Gain (Original Scale)",
       title = "Gain vs Density with 95% Uncertainty Bands") +
  theme_minimal()

table <- data.frame(
  Density = c(0.001, 0.08, 0.148, 0.223, 0.318, 0.412, 0.508, 0.686),
  Predicted_Value = c(
    round(get_uncertainty_band(0.001)[2], 2),
    round(get_uncertainty_band(0.08)[2], 2),
    round(get_uncertainty_band(0.148)[2], 2),
    round(get_uncertainty_band(0.223)[2], 2),
    round(get_uncertainty_band(0.318)[2], 2),
    round(get_uncertainty_band(0.412)[2], 2),
    round(get_uncertainty_band(0.508)[2], 2),
    round(get_uncertainty_band(0.686)[2], 2)
  ),
  `95%_Prediction_Interval` = c(
    paste(round(get_uncertainty_band(0.001)[1], 2), "to", round(get_uncertainty_band(0.001)[3], 2)),
    paste(round(get_uncertainty_band(0.08)[1], 2), "to", round(get_uncertainty_band(0.08)[3], 2)),
    paste(round(get_uncertainty_band(0.148)[1], 2), "to", round(get_uncertainty_band(0.148)[3], 2)),
    paste(round(get_uncertainty_band(0.223)[1], 2), "to", round(get_uncertainty_band(0.223)[3], 2)),
    paste(round(get_uncertainty_band(0.318)[1], 2), "to", round(get_uncertainty_band(0.318)[3], 2)),
    paste(round(get_uncertainty_band(0.412)[1], 2), "to", round(get_uncertainty_band(0.412)[3], 2)),
    paste(round(get_uncertainty_band(0.508)[1], 2), "to", round(get_uncertainty_band(0.508)[3], 2)),
    paste(round(get_uncertainty_band(0.686)[1], 2), "to", round(get_uncertainty_band(0.686)[3], 2))
  ),
  Actual_Range = sapply(
    c(0.001, 0.08, 0.148, 0.223, 0.318, 0.412, 0.508, 0.686),
    function(d) paste(range(data$gain[data$density == d])[1], "to", range(data$gain[data$density == d])[2])
  )
)

kable(table, col.names = c("Density", "Predicted Value", "95% Uncertainty Bands", "Actual Range"),
      caption = "Prediction Intervals (Original Scale)")

# Question 5: Reverse Prediction
# Reverse the prediction function to predict density from gain, and similar to
# last question, produce interval estimates for density and plot

beta_0 <- coef(log_model)[1]
beta_1 <- coef(log_model)[2]
n <- nrow(data)

get_se <- function(gain_val) {
  return(sigma * sqrt(1 + 1/n + (log(gain_val) - mean(data$log_gain))^2 / sum((data$log_gain - mean(data$log_gain))^2)))
}

# Reverse prediction: Invert the model to find density
reverse_predict_density <- function(gain_value, log_model) {
  beta_0 <- coef(log_model)[1]
  beta_1 <- coef(log_model)[2]
  
  # Calculate point estimate for density
  ln_gain <- log(gain_value)
  density_estimate <- (ln_gain - beta_0) / beta_1
  
  # Get standard error for the log prediction
  se <- get_se(gain_value)
  
  # Calculate t-value for desired confidence level
  t_value <- qt(1 - 0.05/2, df = length(data$log_gain) - 2)
  
  # Calculate confidence intervals on log scale
  margin <- t_value * se
  
  # Transform confidence intervals back to density scale
  lower_density <- (ln_gain - (beta_0 + margin)) / beta_1
  upper_density <- (ln_gain - (beta_0 - margin)) / beta_1
  
  return(data.frame(
    gain = gain_value,
    predicted_density = density_estimate,
    lower_ci = lower_density,
    upper_ci = upper_density
  ))
}

# Reverse predict densities for all gain values in your dataset
gain_values <- data$gain
results <- do.call(rbind, lapply(gain_values, function(g) {
  reverse_predict_density(g, log_model)
}))

# get specifics for table
breaks = c(38.6, 426.7)
specific_results <- do.call(rbind, lapply(breaks, function(g) {
  # Get predicted values
  pred_result <- reverse_predict_density(g, log_model)
  
  actual_density <- data$density[which.min(abs(data$gain - g))]
  
  # Add actual density column to the results
  pred_result$actual_density <- actual_density
  
  return(pred_result)
}))

ggplot() +
  geom_ribbon(data = results, 
              aes(x = gain, ymin = lower_ci, ymax = upper_ci),
              fill = "#FB6542", alpha = 0.3) +
  geom_point(data = results,
             aes(x = gain, y = predicted_density),
             # Add prediction intervals
             color = "#375E97", size = 3) +
  geom_point(data = data.frame(gain = c(38.6, 426.7),
                               density = c(0.508, 0.001)),
             aes(x = gain, y = density),
             color = "#6AB187", size = 3) +
  scale_x_continuous(breaks = c(38.6, 426.7)) +
  labs(x = "Gain", y = "Density",
       title = "Reverse Predictions with Confidence Intervals",
       subtitle = "Blue: Predicted Values, Green: True Values") +
  theme_minimal()

specific_results$difference <- abs(specific_results$predicted_density - specific_results$actual_density)
rownames(specific_results) <- NULL

kable(specific_results, digits = 4,
      col.names = c("Gain", "Predicted Density", "Lower CI", "Upper CI", 
                    "Actual Density", "Prediction Error"),
      caption = "Predictions vs Actual Densities")

# Question 6: Cross-Validation
# Remove values with density 0.508 from the training set and validate the model
# on the unseen data to see how well it generalizes
cross_val_data <- data[data$density != 0.508, ]

cv_model <- lm(log_gain ~ density, data = cross_val_data)
cv_beta_0 <- coef(cv_model)[1]
cv_beta_1 <- coef(cv_model)[2]
cv_predict_density <- function(gain_value) {
  ln_gain <- log(gain_value)
  density_value <- (ln_gain - cv_beta_0) / cv_beta_1
  return(density_value)
}

pred_densities <- cv_predict_density(cross_val_data$gain)
n <- nrow(cross_val_data)
cv_residuals <- cross_val_data$density - pred_densities
sigma <- sqrt(sum(cv_residuals^2)/(n-2))
cv_t_crit <- qt(1 - 0.05/2, n-2)


plot(data$density ~ data$gain, ylab="Density", xlab="Gain", main="Reverse Prediction of Density from Gain with Cross-Validated Model\nMissing Density From Training: 0.508")

gain_seq <- seq(0, 500, length.out=500)
density_seq <- cv_predict_density(gain_seq)
se_fit_seq <- sigma * sqrt(1 + 1/n + (log(gain_seq) - mean(cross_val_data$log_gain))^2 / sum((cross_val_data$log_gain - mean(cross_val_data$log_gain))^2))
lwr_ci_seq <- cv_predict_density(gain_seq) - cv_t_crit * se_fit_seq
upr_ci_seq <- cv_predict_density(gain_seq) + cv_t_crit * se_fit_seq

lines(gain_seq, density_seq, type="l", col = "blue", lwd=2)
lines(gain_seq, lwr_ci_seq, type="l", col = "green", lwd=2)
lines(gain_seq, upr_ci_seq, type="l", col = "green", lwd=2)
legend("topright", 
       legend = c("Regression Line", "Prediction Interval Bounds"), 
       col = c("blue", "green"), 
       lty = c(1, 2), 
       lwd = c(2, 2))

get_cv_se <- function(gain_val) {
  return(sigma * sqrt(1 + 1/n + (log(gain_val) - mean(cross_val_data$log_gain))^2 / sum((cross_val_data$log_gain - mean(cross_val_data$log_gain))^2)))
}

cv_test_estimate <- cv_predict_density(38.6)
cv_test_lower <- cv_test_estimate - cv_t_crit * get_cv_se(38.6)
cv_test_upper <- cv_test_estimate + cv_t_crit * get_cv_se(38.6)

# Do the same process of removing density values of 0.001 from the training data
# and use the reverse function to predict density from a gain value of 426.7 to see
# how well model generalizes

cv_2_data <- data[data$density != 0.001, ]

cv_2_model <- lm(log_gain ~ density, data = cv_2_data)
cv_2_beta_0 <- coef(cv_model)[1]
cv_2_beta_1 <- coef(cv_model)[2]
cv_2_predict_density <- function(gain_value) {
  ln_gain <- log(gain_value)
  density_value <- (ln_gain - cv_2_beta_0) / cv_2_beta_1
  return(density_value)
}

cv_2_pred_densities <- cv_2_predict_density(cv_2_data$gain)
n_2 <- nrow(cv_2_data)
cv_2_residuals <- cv_2_data$density - cv_2_pred_densities
sigma_2 <- sqrt(sum(cv_residuals^2)/(n-2))

get_cv_2_se <- function(gain_val) {
  return(sigma_2 * sqrt(1 + 1/n_2 + (log(gain_val) - mean(cv_2_data$log_gain))^2 / sum((cv_2_data$log_gain - mean(cv_2_data$log_gain))^2)))
}

average_gain_reading_2 <- mean(data[data$density == 0.001, ]$gain)

cv_2_test_estimate <- cv_2_predict_density(average_gain_reading_2)
cv_2_test_lower <- cv_2_test_estimate - cv_t_crit * get_cv_2_se(average_gain_reading_2)
cv_2_test_upper <- cv_2_test_estimate + cv_t_crit * get_cv_2_se(average_gain_reading_2)

results <- data.frame(
  Average_Gain = c(38.6, 426.7),
  Point_Estimate = c(0.5092, -0.0128),
  Lower_Bound = c(0.4780, -0.0446),
  Upper_Bound = c(0.5404, 0.0190)
)

# Print the data frame as a table
print(results)


# Advanced Analysis: Exploring Transformation Choices

# Train a polynomial model with squared and linear features and plot regression results
poly_model <- lm(gain ~ poly(density, 2), data=data)

density_seq <- seq(min(data$density), max(data$density), length.out = 500)
predicted_gain <- predict(poly_model, newdata = data.frame(density = density_seq))

plot(data$density, data$gain, 
     pch = 19, col = "red", 
     xlab = "Density", ylab = "Gain", 
     main = "Polynomial Regression: Gain vs Density")

lines(density_seq, predicted_gain, col = "blue", lwd = 2)

# Train a model using log-10 transformation on gain
data$log_10_gain <- log10(data$gain)

log_10_model <- lm(log_10_gain ~ density, data=data)

plot(data$density, data$log_10_gain,
     main = "log_10(Gain) vs Density", 
     xlab = "Snow Density (g/cm³)", 
     ylab = "Base-10 Log of Gain", 
     pch = 19, 
     col = "#89ABE3")

abline(log_10_model, col = "#EA738D", lwd = 2)
predicted_log_10_gain <- predict(log_10_model)
segments(data$density, data$log_10_gain, data$density, predicted_log_10_gain, col = "purple", lwd = 1, lty = 2)

# Train a model using the square-root transformation on gain and plot regression
data$sqrt_gain <- sqrt(data$gain)
sqrt_model <- lm(sqrt_gain ~ density, data=data)
plot(data$density, data$sqrt_gain,
     main = "Square-root Gain vs Density", 
     xlab = "Snow Density (g/cm³)", 
     ylab = "sqrt(Gain)", 
     pch = 19, 
     col = "#89ABE3")

abline(sqrt_model, col = "#EA738D", lwd = 2)
predicted_sqrt_gain <- predict(sqrt_model)
segments(data$density, data$sqrt_gain, data$density, predicted_sqrt_gain, col = "purple", lwd = 1, lty = 2)
