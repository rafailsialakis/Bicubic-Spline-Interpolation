
import sys
import time

from PIL import Image
"""
Use of Thomas algorithm to solve the tridiagonal
system
"""
def Thomas(n, a, b, c, d, M):
    for i in range(1, n):
        m = a[i] / b[i - 1]
        b[i] = b[i] - m * c[i - 1]
        d[i] = d[i] - m * d[i - 1]
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

def Resample(new_size, size, x):
    x_resampled = [0] * new_size
    Dx = size / new_size
    for i in range(new_size):
        x_resampled[i] = (i + 0.5) * Dx - 0.5
        x_resampled[i] = max(0, min(size - 1, x_resampled[i]))
    return x_resampled

def ResampleNoArtifacts(new_size, size, x):
    x_resampled = [0] * new_size
    for i in range(new_size):
        x_resampled[i] = i * ((size-1)/(new_size-1))
    return x_resampled

def ReconstructY(x_resampled, x, a, b, c, d):
    size = len(x)
    size_of_resampled = len(x_resampled)
    y_reconstructed = [0] * size_of_resampled
    for i in range(1, size):
        for j in range(size_of_resampled):
            if x[i - 1] <= x_resampled[j] <= x[i]:
                y_reconstructed[j] = (
                    a[i - 1] * (x_resampled[j] - x[i - 1]) ** 3
                    + b[i - 1] * (x_resampled[j] - x[i - 1]) ** 2
                    + c[i - 1] * (x_resampled[j] - x[i - 1])
                    + d[i - 1]
                )
    return y_reconstructed

"""
Spline coefficients calculated based on the second derivative M,
which was calculated using Thomas algorithm for tridiagonal systems
"""
def compute_spline_coefficients(n, x, y, M, a, b, c, d):
    for i in range(n - 1):
        h = x[i + 1] - x[i]
        a[i] = (M[i + 1] - M[i]) / (6 * h)
        b[i] = M[i] / 2.0
        c[i] = (y[i + 1] - y[i]) / h - h * (M[i + 1] + 2 * M[i]) / 6.0
        d[i] = y[i]
    return a, b, c, d

"""
cubic_spline_1D(data, new_size):
   - Performs cubic spline interpolation on 1D data (e.g., a single row or column of an image).
   - It calculates second derivatives (`M`) using the Thomas algorithm for solving a tridiagonal system.
   - Computes the spline coefficients (`a`, `b`, `c`, `d`) for cubic interpolation.
   - Resamples the input data to the desired size using cubic interpolation.
"""
def cubic_spline_1D(data, new_size):
    size = len(data)
    h = 1
    rhs = [0] * (size - 2)
    M = [0] * size
    a = [0] * (size - 1)
    b = [0] * (size - 1)
    c = [0] * (size - 1)
    d = [0] * (size - 1)
    a_tridiagonal = [h / 6.0] * (size - 2)
    b_tridiagonal = [2 * h / 3.0] * (size - 2)
    c_tridiagonal = [h / 6.0] * (size - 2)

    for i in range(1, size - 1):
        rhs[i - 1] = (data[i + 1] - 2 * data[i] + data[i - 1]) / h

    M = Thomas(size - 2, a_tridiagonal, b_tridiagonal, c_tridiagonal, rhs, M)
    a, b, c, d = compute_spline_coefficients(size, list(range(size)), data, M, a, b, c, d)
    x_resampled = ResampleNoArtifacts(new_size, size, list(range(size)))
    return ReconstructY(x_resampled, list(range(size)), a, b, c, d)

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
    #Interpolation of rows
    row_interpolated = []
    for row in image_data:
        row_interpolated.append(cubic_spline_1D(row, new_width))

    #Interpolation of columns
    column_interpolated = []
    for col_idx in range(new_width):
        column = [row[col_idx] for row in row_interpolated]
        column_resampled = cubic_spline_1D(column, new_height)
        column_interpolated.append(column_resampled)

    #Transpose back to rows for final 2D output.
    #Unpack first the column interpolated matrix and use
    #zip to groups corresponding elements from each input iterable into tuples.

    final_image = list(zip(*column_interpolated))
    return final_image

"""
Method to open the image. Returns the data in a matrix.
Currently is converted to "L"
TODO Convert to RGB and handle all the colours
"""

# Other functions (Thomas, Resample, ReconstructY, compute_spline_coefficients, cubic_spline_1D, cubic_spline_2D) remain unchanged.

def ImageOpen(path):
    image = Image.open(path).convert("RGB")  # Convert to RGB
    width, height = image.size
    image_data = list(image.getdata())
    return [image_data[i * width: (i + 1) * width] for i in range(height)], width, height


def createNewImage(data_R, data_G, data_B, width, height, output_path):
    # Combine R, G, B data into a single list of tuples
    flat_data = [
        (
            int(max(0, min(255, data_R[row][col]))),
            int(max(0, min(255, data_G[row][col]))),
            int(max(0, min(255, data_B[row][col])))
        )
        for row in range(height)
        for col in range(width)
    ]
    image = Image.new("RGB", (width, height))
    image.putdata(flat_data)
    image.save(output_path)



def extract_color_channel(image_data, color_index):
    # Extracts a specific color channel (R=0, G=1, B=2) from the image data
    return [[pixel[color_index] for pixel in row] for row in image_data]


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
    new_width = width * maximizing_factor
    new_height = height * maximizing_factor

    start = time.time()

    # Process each color channel separately
    red_channel = extract_color_channel(image_data, 0)
    green_channel = extract_color_channel(image_data, 1)
    blue_channel = extract_color_channel(image_data, 2)

    resampled_red = cubic_spline_2D(red_channel, new_width, new_height)
    resampled_green = cubic_spline_2D(green_channel, new_width, new_height)
    resampled_blue = cubic_spline_2D(blue_channel, new_width, new_height)

    # Create the final resampled RGB image
    createNewImage(resampled_red, resampled_green, resampled_blue, new_width, new_height, output_path)
    end = time.time()

    print(f"Resampled image saved to {output_path}\nTime taken: {end-start:.2f}s")



if __name__ == "__main__":
    main()
