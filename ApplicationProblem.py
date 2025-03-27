import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

# Example set
# (Replace this with actual data from the machine)
success_failures = np.array([0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1])

# Convert to number of attempts before success
num_attempts = np.diff(np.where(np.concatenate(([1], success_failures == 1, [1])))[0])

# Estimate probability of success (p)
p = 1 / np.mean(num_attempts)  # p = 1/E[X]

# Print estimated probability and expected number ofimport scipy.stats as stats tries
print(f"Sample Estimated Probability (p): {p:.4f}")
print(f"Sample Expected Number of Tries (1/p): {1/p:.2f}")

# Plot the histogram of observed data
sns.histplot(num_attempts, bins=range(1, max(num_attempts) + 2), kde=False, stat="probability", label="Observed Data", color='blue')
plt.xlabel("Number of Attempts Before Success")
plt.ylabel("Probability")
plt.title("Observed vs. Geometric Distribution")

# Overlay the theoretical geometric PMF
x = np.arange(1, max(num_attempts) + 1)
y = stats.geom.pmf(x, p)
plt.bar(x, y, alpha=0.6, color='red', width=0.4, label="Geometric PMF")

# Add legend and display the plot
plt.legend()
plt.show()
