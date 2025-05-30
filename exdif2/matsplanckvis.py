import numpy as np
import matplotlib.pyplot as plt

# Function to load data from 3D text file
def load_temperature_data(filename, N):
    data = np.loadtxt(filename)
    temperature = np.zeros((N, N, N))

    for line in data:
        i, j, k, value = int(line[0]), int(line[1]), int(line[2]), line[3]
        temperature[i, j, k] = value
    
    return temperature

# Function to plot selected slices
def plot_slice(temp_data, slice_index, axis='xy', colormap='hot'):
    plt.figure(figsize=(6,6))
    
    if axis == 'xy':
        slice_data = temp_data[slice_index, :, :]
        plt.title(f'XY Plane at Z={slice_index}')
    elif axis == 'xz':
        slice_data = temp_data[:, slice_index, :]
        plt.title(f'XZ Plane at Y={slice_index}')
    elif axis == 'yz':
        slice_data = temp_data[:, :, slice_index]
        plt.title(f'YZ Plane at X={slice_index}')
    
    plt.imshow(slice_data, cmap=colormap, origin='lower')
    plt.colorbar(label='Temperature')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

# Load data from text file
N = 100  # Grid resolution
filename = "temperature_solution.txt"
temperature_data = load_temperature_data(filename, N)

# Plot example slices
plot_slice(temperature_data, slice_index=50, axis='xy')
plot_slice(temperature_data, slice_index=75, axis='xz')
plot_slice(temperature_data, slice_index=25, axis='yz')