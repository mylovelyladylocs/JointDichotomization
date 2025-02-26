
#' @title Compute Optimal Joint Cutpoints
#'
#' @description Identifies optimal cutpoints for joint dichotomization of two continuous predictors.
#' to classify a binary outcome using different classification metrics.
#'
#' @name Joint
#'
#' @param mydata A data frame containing variables X1, X2, and Y.
#' @param py Prevalence of Y in the population (default = 0.5).
#'
#' @return A data frame with the optimal cutpoints for each classification metric.


#' @examples
#' data <- data.frame(X1 = rnorm(500), X2 = rnorm(500), Y = sample(0:1, 500, replace = TRUE))
#' Joint(data, py = 0.5)


# Load required packages
library(mvtnorm)  # for rmvnorm and pmvnorm
library(mgcv)     # for gam smoothing

# Read in your data
my_data <- data.frame(
  X1 = rnorm(500),
  X2 = rnorm(500),
  Y = sample(0:1, 500, replace = TRUE)
)
mydata<-my_data
# Inspect the data to verify it loaded correctly
print(head(my_data))
print(summary(my_data))

# (Optional) Rename columns if your CSV's columns differ.
# For example, if your file has columns "Outcome", "Predictor1", "Predictor2":
# names(my_data) <- c("Y", "X1", "X2")

# Function to compute optimal cutpoints using joint dichotomization
#' @export
Joint <- function(mydata, py) {
  library(mgcv)
  # Order and trim the most extreme 5% from the predictors
  sorted_X1 <- sort(mydata$X1)
  sorted_X2 <- sort(mydata$X2)
  trim_lower <- floor(0.05 * length(sorted_X1)) + 1
  trim_upper <- ceiling(0.95 * length(sorted_X1))
  oX1 <- sorted_X1[trim_lower:trim_upper]
  oX2 <- sorted_X2[trim_lower:trim_upper]
  n <- length(oX1)

  # Initialize contingency tables
  a <- matrix(0, nrow = n, ncol = n)
  b <- matrix(0, nrow = n, ncol = n)
  c <- matrix(0, nrow = n, ncol = n)
  d <- matrix(0, nrow = n, ncol = n)

  # Small error term to prevent division by zero
  epsilon <- 1e-6

  # Populate contingency tables with added error term
  for (i in 1:n) {
    for (j in 1:n) {
      a[i, j] <- sum(mydata$X1 >= oX1[i] & mydata$X2 >= oX2[j] & mydata$Y == 1) + epsilon
      b[i, j] <- sum(mydata$X1 >= oX1[i] & mydata$X2 >= oX2[j] & mydata$Y == 0) + epsilon
      c[i, j] <- sum((mydata$X1 < oX1[i] | mydata$X2 < oX2[j]) & mydata$Y == 1) + epsilon
      d[i, j] <- sum((mydata$X1 < oX1[i] | mydata$X2 < oX2[j]) & mydata$Y == 0) + epsilon
    }
  }


  # Compute performance metrics
  ORmat    <- (a * d) / (b * c)
  youden   <- (a / (a + c)) + (d / (b + d)) - 1
  gini     <- 2 * py * (1 - py) - (2 * (a * b / (a + b) + c * d / (c + d)))
  kappa    <- ((a + d) / (a + b + c + d)) - (((a + b) * (a + c) + (c + d) * (b + d)) / ((a + b + c + d)^2))
  chisq    <- ((a + b) * (a + c) * (b + d) * (c + d)) / (a + b + c + d)

  # Create a grid for smoothing
  grid_df <- expand.grid(x = 1:n, y = 1:n)

  mod.OR    <- gam(as.vector(ORmat) ~ te(x, y), data = data.frame(grid_df, z = as.vector(ORmat)))
  OR_smoothed    <- matrix(fitted(mod.OR), ncol = n)

  mod.youden   <- gam(as.vector(youden) ~ te(x, y), data = data.frame(grid_df, z = as.vector(youden)))
  youden_smoothed   <- matrix(fitted(mod.youden), ncol = n)

  mod.gini   <- gam(as.vector(gini) ~ te(x, y), data = data.frame(grid_df, z = as.vector(gini)))
  gini_smoothed   <- matrix(fitted(mod.gini), ncol = n)

  mod.kappa   <- gam(as.vector(kappa) ~ te(x, y), data = data.frame(grid_df, z = as.vector(kappa)))
  kappa_smoothed   <- matrix(fitted(mod.kappa), ncol = n)

  mod.chisq   <- gam(as.vector(chisq) ~ te(x, y), data = data.frame(grid_df, z = as.vector(chisq)))
  chisq_smoothed   <- matrix(fitted(mod.chisq), ncol = n)

  # Find indices of maximum values for each metric
  or.ans     <- which(OR_smoothed == max(OR_smoothed), arr.ind = TRUE)
  youden.ans <- which(youden_smoothed == max(youden_smoothed), arr.ind = TRUE)
  gini.ans   <- which(gini_smoothed == max(gini_smoothed), arr.ind = TRUE)
  kappa.ans  <- which(kappa_smoothed == max(kappa_smoothed), arr.ind = TRUE)
  chisq.ans  <- which(chisq_smoothed == max(chisq_smoothed), arr.ind = TRUE)

  # Construct a table with the optimal cutpoints
  cutpoints_table <- data.frame(
    Method = c("Odds Ratio", "Youden's Index", "Gini", "Kappa", "Chi-Square"),
    X1 = c(oX1[or.ans[1]], oX1[youden.ans[1]], oX1[gini.ans[1]], oX1[kappa.ans[1]], oX1[chisq.ans[1]]),
    X2 = c(oX2[or.ans[2]], oX2[youden.ans[2]], oX2[gini.ans[2]], oX2[kappa.ans[2]], oX2[chisq.ans[2]])
  )

  return(cutpoints_table)
}

# Set the overall prevalence (py) for outcome Y (adjust as necessary)
py <- 0.5

# Compute the joint cutpoints table using your CSV data
joint_results <- Joint(my_data, py)
print(joint_results)
