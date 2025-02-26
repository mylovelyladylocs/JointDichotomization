# JointDichotomization  

### 📊 Identify Optimal Joint Dichotomization Cutpoints  

**JointDichotomization** is an R package for identifying the **optimal cutpoints** for dichotomizing two continuous predictors (`X1` and `X2`) when classifying a binary outcome (`Y`).  
It implements five classification metrics:  
- ✅ **Odds Ratio**  
- ✅ **Youden's Index**  
- ✅ **Gini**  
- ✅ **Kappa**  
- ✅ **Chi-Square**  

---

## **📥 Installation**  

You can install the latest version directly from GitHub:  
```r
# Install from GitHub
devtools::install_github("mylovelyladylocs/JointDichotomization")

# Load the package
library(JointDichotomization)

# Load the package
library(JointDichotomization)

# Simulate some sample data
set.seed(123)
data <- data.frame(
  X1 = rnorm(500),
  X2 = rnorm(500),
  Y = sample(0:1, 500, replace = TRUE)
)

# Compute optimal cutpoints
results <- Joint(data, py = 0.5)

# Print results
print(results)
📖 How It Works
1️⃣ The function trims the extreme 5% of predictor values to reduce noise.
2️⃣ It then constructs contingency tables (a, b, c, d) for all possible cutpoints.
3️⃣ It computes the five classification metrics for all cutpoints.
4️⃣ It applies Generalized Additive Models (GAMs) to smooth the results.
5️⃣ It identifies the cutpoints that maximize each metric.

vignette("Introduction", package = "JointDichotomization")

