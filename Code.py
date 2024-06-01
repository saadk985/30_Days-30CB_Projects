#Importing relevent Libraries
import matplotlib.pyplot as plt
import numpy as np

#Defining Function
def hill_equation(n, kd, L):
  """
  This function calculates the fractional occupancy (Y) for a given oxygen concentration ([L])
  based on the Hill equation.

  Args:
      n (float): Hill coefficient
      kd (float): Dissociation constant
      L (float): Oxygen concentration

  Returns:
      float: Fractional occupancy (Y)
  """
  return L**n / (kd + L**n)

#Defining values
n = 1.5  # Hill coefficient (typical value for hemoglobin-oxygen binding)
kd = 20  # Dissociation constant (micromolar)
L = np.linspace(0.5e-6, 50e-6, 1000)  # Oxygen concentration range (micromolar)

Y = hill_equation(n, kd, L) * 1000  # Convert to percentage saturation

# Create the plot figure
plt.figure(figsize=(8, 6))  # Set figure size

# Plot the fractional occupancy vs oxygen concentration
plt.plot(L, Y, label="Hemoglobin Saturation (%)")

# Add labels and title
plt.xlabel("Oxygen Concentration (uM)")
plt.ylabel("Hemoglobin Saturation (%)")
plt.title("Hill Plot for Hemoglobin-Oxygen Binding (n={}, Kd={} uM)".format(n, kd))

# Add legend
plt.legend()

# Show the plot
plt.grid(True)  # Add grid lines for better readability
plt.show()  # Display the plot