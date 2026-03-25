import os
import sys
import time
import numpy as np
from PIL import Image
from multiprocessing import Pool

def _interpolate_row(args):
    row, new_width = args
    return cubic_spline_1D(row, new_width)

def _interpolate_col(args):
    col, new_height = args
    return cubic_spline_1D(col, new_height)

"""
Use of Thomas algorithm to solve the tridiagonal
system
"""
def Thomas(n, a, b, c, d):
    a, b, c, d = a.copy(), b.copy(), c.copy(), d.copy()
    M = np.zeros(n + 2)

    for i in range(1, n):
        m = a[i] / b[i - 1]
        b[i] -= m * c[i - 1]
        d[i] -= m * d[i - 1]

    M[n] = d[n - 1] / b[n - 1]
    for i in range(n - 2, -1, -1):
        M[i + 1] = (d[i] - c[i] * M[i + 2]) / b[i]

    return M

"""
Resampling is given by the equation:
Δx = width/new_width
x[i] = x[i-1] + Δχ

But then we don't have a good distribution on the polynomial

So softwares use:
x[i] = (n + 0.5) * Δχ - 0.5, where n is the index.

In the new distribution some values appear that are smaller than
the smallest value in x and some larger than the largest value in x.

So the question arises, which polynomial are we using to interpolate them?
"""
def ResampleNoArtifacts(new_size, size):
    return np.linspace(0, size - 1, new_size)

def ReconstructY(x_resampled, x, a, b, c, d):
    size = len(x)
    size_of_resampled = len(x_resampled)
    y_reconstructed = np.zeros(size_of_resampled)

    for i in range(1, size):
        mask = (x_resampled >= x[i - 1]) & (x_resampled <= x[i])
        dx = x_resampled[mask] - x[i - 1]
        y_reconstructed[mask] = (
            a[i - 1] * dx**3
            + b[i - 1] * dx**2
            + c[i - 1] * dx
            + d[i - 1]
        )

    return y_reconstructed

"""
Spline coefficients calculated based on the second derivative M,
which was calculated using Thomas algorithm for tridiagonal systems
"""
def compute_spline_coefficients(x, y, M):
    n = len(x)
    h = np.diff(x).astype(float)
    a = (M[1:n] - M[:n-1]) / (6 * h)
    b = M[:n-1] / 2.0
    c = (y[1:] - y[:-1]) / h - h * (M[1:n] + 2 * M[:n-1]) / 6.0
    d = y[:-1].astype(float)
    return a, b, c, d

"""
cubic_spline_1D(data, new_size):
   - Performs cubic spline interpolation on 1D data (e.g., a single row or column of an image).
   - It calculates second derivatives (`M`) using the Thomas algorithm for solving a tridiagonal system.
   - Computes the spline coefficients (`a`, `b`, `c`, `d`) for cubic interpolation.
   - Resamples the input data to the desired size using cubic interpolation.
"""
def cubic_spline_1D(data, new_size):
    data = np.asarray(data, dtype=float)
    size = len(data)
    h = 1.0

    a_tri = np.full(size - 2, h / 6.0)
    b_tri = np.full(size - 2, 2 * h / 3.0)
    c_tri = np.full(size - 2, h / 6.0)
    rhs   = (data[2:] - 2 * data[1:-1] + data[:-2]) / h

    M = Thomas(size - 2, a_tri, b_tri, c_tri, rhs)

    x = np.arange(size, dtype=float)
    a, b, c, d = compute_spline_coefficients(x, data, M)

    x_resampled = ResampleNoArtifacts(new_size, size)
    return ReconstructY(x_resampled, x, a, b, c, d)

"""
cubic_spline_2D(image_data, new_width, new_height):
   - Handles cubic spline interpolation for 2D images.
   - **Step 1: Interpolating Rows:**
     - Each row in the image data is resampled to match the new width using `cubic_spline_1D`.
   - **Step 2: Interpolating Columns:**
     - Columns of the row-resampled image are extracted and resampled to match the new height using `cubic_spline_1D`.
   - Combines the resampled columns into the final 2D interpolated image by transposing the column data back into rows.
"""
def cubic_spline_2D(image_data, new_width, new_height):
    image_data = np.asarray(image_data, dtype=float)
    height, width = image_data.shape
    workers = os.cpu_count()

    row_args = [(image_data[r], new_width) for r in range(height)]
    with Pool(processes=workers) as pool:
        row_interpolated = np.array(pool.map(_interpolate_row, row_args))
    # row_interpolated shape: (height, new_width)

    col_args = [(row_interpolated[:, c], new_height) for c in range(new_width)]
    with Pool(processes=workers) as pool:
        col_interpolated = np.array(pool.map(_interpolate_col, col_args))
    # col_interpolated shape: (new_width, new_height) → transpose
    return col_interpolated.T

"""
Method to open the image. Returns the data as a numpy array (height, width, 3).
"""
def ImageOpen(path):
    image = Image.open(path).convert("RGB")
    width, height = image.size
    image_data = np.array(image)  # shape: (height, width, 3)
    return image_data, width, height


def createNewImage(data_R, data_G, data_B, width, height, output_path):
    rgb = np.stack([
        np.clip(data_R, 0, 255).astype(np.uint8),
        np.clip(data_G, 0, 255).astype(np.uint8),
        np.clip(data_B, 0, 255).astype(np.uint8),
    ], axis=-1)  # shape: (height, width, 3)

    image = Image.fromarray(rgb, mode="RGB")
    image.save(output_path)


def main():
    if len(sys.argv) != 4:
        print("Usage: python3 bsisequential.py <original_file> <output_file> <multiplicative_factor>")
        sys.exit(1)
    try:
        maximizing_factor = int(sys.argv[3])
    except ValueError:
        print("Error: <multiplicative_factor> must be an integer.")
        sys.exit(1)

    path = "Images/" + sys.argv[1]
    output_path = "Images/" + sys.argv[2]

    image_data, width, height = ImageOpen(path)
    new_width  = width  * maximizing_factor
    new_height = height * maximizing_factor

    start = time.time()

    # Process each color channel separately
    red_channel   = image_data[:, :, 0]
    green_channel = image_data[:, :, 1]
    blue_channel  = image_data[:, :, 2]

    resampled_red   = cubic_spline_2D(red_channel,   new_width, new_height)
    resampled_green = cubic_spline_2D(green_channel, new_width, new_height)
    resampled_blue  = cubic_spline_2D(blue_channel,  new_width, new_height)

    createNewImage(resampled_red, resampled_green, resampled_blue, new_width, new_height, output_path)
    end = time.time()

    print(f"Resampled image saved to {output_path}\nTime taken: {end - start:.2f}s")


if __name__ == "__main__":
    main()