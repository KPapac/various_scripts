#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import numpy as np


def filter_outliers(data):
    """
    Filter outliers from data based on the interquartile range (IQR).
    """
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    return [x for x in data if lower_bound <= x <= upper_bound]


def plot_histogram(data):
    """
    Plot histogram for the filtered data.
    """
    plt.hist(data, bins=30, color='skyblue', edgecolor='black')
    plt.title('Histogram of Filtered Data')
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.savefig(f'{sys.argv[1]}_hist.pdf')
    plt.close()


def main():
    # Read list of numbers from standard input
    numbers = []
    for line in sys.stdin:
        numbers.extend(map(float, line.strip().split()))
    # Filter outliers
    filtered_numbers = filter_outliers(numbers)
    # Plot histogram
    plot_histogram(filtered_numbers)


if __name__ == "__main__":
    main()
