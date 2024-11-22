# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 03:42:27 2024

@author: Junhyung Park
"""

# Import necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns

# Load the merged data from Homework #1
merged_data_path = 'output_combined_vwt_data.xlsx'
merged_data = pd.read_excel(merged_data_path)

# --- Task 1: Scatterplot to Examine Correlation ---
# Variables selection for correlation
# x: Trade volume (Value column), y: Virtual Water Trade (VWT column)
x = merged_data['Value']
y = merged_data['VWT']

# Scatter plot of Trade Volume vs Virtual Water Trade
plt.figure(figsize=(12, 8))  # Slightly larger figure size
sns.scatterplot(x=x, y=y, color='blue')
plt.title('Scatterplot of Trade Volume vs Virtual Water Trade (VWT)', fontsize=20, fontweight='bold')
plt.xlabel('Trade Volume', fontsize=18, fontweight='bold')
plt.ylabel('Virtual Water Trade (VWT)', fontsize=18, fontweight='bold')
plt.xticks(fontsize=14, fontweight='bold')  # X-axis tick labels
plt.yticks(fontsize=14, fontweight='bold')  # Y-axis tick labels
plt.grid(True)
plt.show()
#%%
# --- Task 2: Ordinary Least Squares (OLS) Regression ---
# Add a constant to the independent variable for OLS regression
x_with_const = sm.add_constant(x)
model = sm.OLS(y, x_with_const).fit()

# Print the summary of the OLS regression
print(model.summary())

# Convert the regression results to a DataFrame
summary_df = pd.DataFrame({
    "Coefficients": model.params,
    "Standard Errors": model.bse,
    "t-values": model.tvalues,
    "P-values": model.pvalues
})

# Add key metrics like R-squared and F-statistic
additional_metrics = pd.DataFrame({
    "Metrics": ["R-squared", "Adj. R-squared", "F-statistic", "Prob (F-statistic)"],
    "Values": [model.rsquared, model.rsquared_adj, model.fvalue, model.f_pvalue]
})

# Save to Excel using ExcelWriter
with pd.ExcelWriter("ols_regression_results.xlsx") as writer:
    summary_df.to_excel(writer, sheet_name="Regression Coefficients", index=True)
    additional_metrics.to_excel(writer, sheet_name="Metrics", index=False)

print("OLS regression results have been saved to 'ols_regression_results.xlsx'.")