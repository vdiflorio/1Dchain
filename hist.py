import numpy as np
import matplotlib.pyplot as plt



# Parameters
a = 1.0  # Coefficient for x^3 term
b = 1.0  # Coefficient for x^4 term
beta = 1.0  # Inverse temperature

# Step 1: Load the data and process it as before
pressione_vect = []
differenze_vect =[]
f_l_vect =[]
f_r_vect =[]


with open('condizioni_iniziali.dat', 'r') as f:  # Replace 'your_file.txt' with your actual filename
    for line in f:
        line = line.strip()
        if line and '--' not in line:  # Ignore separator lines
            position = float(line.split()[0])  # Extract the position (first column)
            pressione_tmp = float(line.split()[1])  # Extract the position (first column)
            f_l_tmp = float(line.split()[2])
            f_r_tmp = float(line.split()[3])
            differenze_vect.append(position)
            pressione_vect.append(pressione_tmp)
            f_l_vect.append(f_l_tmp)
            f_r_vect.append(f_r_tmp)

# Step 2: Calculate differences within each block
# position_diffs = []
# for block in all_positions:
#     if len(block) > 1:
#         diffs = np.diff(block)
#         position_diffs.extend(diffs)

pressione = np.mean(f_r_vect)

print(np.mean(differenze_vect))
print(np.mean(f_l_vect))
print(np.mean(pressione))
print(np.mean(pressione_vect))


# Step 3: Plot histogram of position differences
# plt.hist(position_diffs, bins=20, density=True, alpha=0.6, edgecolor='black', label='Position Differences')
plt.hist(differenze_vect, bins=50, density=True, alpha=0.6, edgecolor='black', label='Position Differences')

# Step 4: Define the energy function and density based on Boltzmann distribution
def energy(x, a, b, press):
    return 0.5*(x-1)**2 + a * (x-1)**3*1/3 + b * 0.25*(x-1)**4 -press*(x)

def density(x, a, b, beta, press):
    return np.exp(-beta * energy(x, a, b, press))

# Generate x values for plotting the density
x_vals = np.linspace(min(differenze_vect), max(differenze_vect), 500)
# x_vals = np.linspace(min(position_diffs), max(position_diffs), 500)
density_vals = density(x_vals, a, b, beta, pressione)
density_vals /= np.trapz(density_vals, x_vals)  # Normalize the density

# Step 5: Plot the density function
plt.plot(x_vals, density_vals, color='red', linewidth=2, label='Energy-based Density')

# Add labels and legend
plt.xlabel('Difference between consecutive positions')
plt.ylabel('Probability Density')
plt.title('Histogram of Position Differences with Energy-based Density')
plt.legend()
plt.show()
