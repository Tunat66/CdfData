import matplotlib.pyplot as plt

# Read data from the file
def read_data(file_path):
    with open(file_path, 'r') as file:
        data = [float(line.strip()) for line in file if line.strip()]
    return data

# Plot histogram
def plot_histogram(data, bins, title, xlabel, ylabel, output_file):
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=bins, color='blue', edgecolor='black', alpha=0.7, range=(0.4, 0.6))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(0.4, 0.6)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.savefig(output_file)  # Save the histogram as an image file
    plt.show()

# Main function
if __name__ == "__main__":
    # File path to the data file
    file_path = "massArray.dat"
    # Read the data
    data = read_data(file_path)
    # Plot the histogram
    plot_histogram(
        data=data,
        bins=50,  # Number of bins
        title="Mass Array Histogram",
        xlabel="Mass (GeV/c^2)",
        ylabel="Frequency",
        output_file="mass_histogram.png"  # Output image file
    )