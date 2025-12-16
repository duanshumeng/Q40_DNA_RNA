import pandas as pd
import numpy as np

# Generate sequential index column
column1 = list(range(1, 11))

# Generate random float numbers between 0 and 1, 4 columns and 10 rows
data = np.random.rand(10, 4)

# Create dataframe
df = pd.DataFrame(data, columns=['Column2', 'Column3', 'Column4', 'Column5'])
df.insert(0, 'Column1', column1) # Insert the sequential index column as the first column

# Calculate Pearson correlation coefficient matrix
correlation_matrix = df.iloc[:, 1:].corr()

# Extract pairwise Pearson correlation coefficients between columns 2 to 5
correlation_pairs = correlation_matrix.unstack().sort_values(ascending=False)
correlation_pairs = correlation_pairs[correlation_pairs != 1] # Remove self-correlations

# Print dataframe
print(df)

# Print pairwise Pearson correlation coefficients
print("Pearson Correlation Coefficients:")
print(correlation_pairs)