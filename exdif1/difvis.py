import numpy as np
import matplotlib.pyplot as plt

# Load data from the output file
data = np.loadtxt("radiation_output.txt")

# Extract position (x) and radiation intensity
x = data[:, 0]
radiation = data[:, 1]

# Plot the results
plt.figure(figsize=(8, 5))
plt.plot(x, radiation, label="Radiation Intensity", color="blue", linestyle="-", marker="o")

# Labels and title
plt.xlabel("Position (x)")
plt.ylabel("Radiation Intensity")
plt.title("Radiation Intensity Across the 1D Domain")
plt.legend()
plt.grid(True)

# Show plot
plt.show()