import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

# Simulated dataset: Success (1) or Failure (0) recorded by the machine
# (Replace this with actual data from the machine)
success_failures = np.array([0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1])

# Convert to number of attempts before success
num_attempts = np.diff(np.where(np.concatenate(([1], success_failures == 1, [1])))[0])

# Estimate probability of success (p)
p_hat = 1 / np.mean(num_attempts)  # p = 1/E[X]

# Print estimated probability and expected number ofimport scipy.stats as stats tries
print(f"Estimated Probability of Success (p): {p_hat:.4f}")
print(f"Expected Number of Tries (1/p): {1/p_hat:.2f}")

# Plot the histogram of observed data
sns.histplot(num_attempts, bins=range(1, max(num_attempts) + 2), kde=False, stat="probability")
plt.xlabel("Number of Attempts Before Success")
plt.ylabel("Probability")
plt.title("Observed Distribution of Attempts Before Winning")
plt.xticks(range(1, max(num_attempts) + 1))
plt.show()

# Overlay the theoretical geometric PMF
x = np.arange(1, max(num_attempts) + 1)
y = stats.geom.pmf(x, p_hat)
plt.bar(x, y, alpha=0.6, color='red', label="Geometric PMF")
sns.histplot(num_attempts, bins=range(1, max(num_attempts) + 2), kde=False, stat="probability", label="Observed Data")
plt.xlabel("Number of Attempts Before Success")
plt.ylabel("Probability")
plt.title("Comparison: Observed vs. Theoretical Geometric Distribution")
plt.legend()
plt.show()# Plot the histogram of observed data (no overlay)
sns.histplot(num_attempts, bins=range(1, max(num_attempts) + 2), kde=False, stat="probability")
plt.xlabel("Number of Attempts Before Success")
plt.ylabel("Probability")
plt.title("Observed Distribution of Attempts Before Winning")
plt.xticks(range(1, max(num_attempts) + 1))
plt.show()

# Plot the histogram of observed data (no overlay)
sns.histplot(num_attempts, bins=range(1, max(num_attempts) + 2), kde=False, stat="probability")
plt.xlabel("Number of Attempts Before Success")
plt.ylabel("Probability")
plt.title("Observed Distribution of Attempts Before Winning")
plt.xticks(range(1, max(num_attempts) + 1))
plt.show()# Plot the histogram of observed data (no overlay)
sns.histplot(num_attempts, bins=range(1, max(num_attempts) + 2), kde=False, stat="probability")
plt.xlabel("Number of Attempts Before Success")
plt.ylabel("Probability")
plt.title("Observed Distribution of Attempts Before Winning")
plt.xticks(range(1, max(num_attempts) + 1))
plt.show()
