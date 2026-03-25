/*
 * bsi_parallel.cpp — Parallel Bicubic Spline Interpolation in C++
 *
 * Build:
 *   g++ -O2 -fopenmp -o bsi_parallel bsi_parallel.cpp -lstdc++ -lm
 *
 * Requires stb_image / stb_image_write (header-only, drop alongside this file):
 *   https://github.com/nothings/stb
 *
 * Usage:
 *   ./bsi_parallel <input_file> <output_file> <scale_factor>
 *   e.g.: ./bsi_parallel photo.png out.png 2
 *
 * Images are expected / written relative to an "Images/" subdirectory,
 * matching the Python original's convention.
 */

#define STB_IMAGE_IMPLEMENTATION
#include "../includes/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../includes/stb_image_write.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <omp.h>

// ---------------------------------------------------------------------------
// Thomas algorithm — solves the tridiagonal system in-place.
// Returns M of length n+2 (M[0] = M[n+1] = 0, natural spline BCs).
// ---------------------------------------------------------------------------
static std::vector<double> Thomas(int n,
                                   std::vector<double> a,
                                   std::vector<double> b,
                                   std::vector<double> c,
                                   std::vector<double> d)
{
    // Forward sweep
    for (int i = 1; i < n; ++i) {
        double m = a[i] / b[i - 1];
        b[i] -= m * c[i - 1];
        d[i] -= m * d[i - 1];
    }

    std::vector<double> M(n + 2, 0.0); // M[0] = M[n+1] = 0 (natural BC)

    M[n] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; --i)
        M[i + 1] = (d[i] - c[i] * M[i + 2]) / b[i];

    return M;
}

// ---------------------------------------------------------------------------
// Evenly-spaced resampling grid: linspace(0, size-1, new_size)
// ---------------------------------------------------------------------------
static std::vector<double> ResampleNoArtifacts(int new_size, int size)
{
    std::vector<double> xs(new_size);
    if (new_size == 1) { xs[0] = 0.0; return xs; }
    double step = (double)(size - 1) / (double)(new_size - 1);
    for (int i = 0; i < new_size; ++i)
        xs[i] = i * step;
    return xs;
}

// ---------------------------------------------------------------------------
// Compute cubic spline coefficients a,b,c,d from x, y, and second-derivatives M.
// ---------------------------------------------------------------------------
struct SplineCoeffs {
    std::vector<double> a, b, c, d;
};

static SplineCoeffs compute_spline_coefficients(const std::vector<double>& x,
                                                 const std::vector<double>& y,
                                                 const std::vector<double>& M)
{
    int n = (int)x.size();
    SplineCoeffs sc;
    sc.a.resize(n - 1); sc.b.resize(n - 1);
    sc.c.resize(n - 1); sc.d.resize(n - 1);

    for (int i = 0; i < n - 1; ++i) {
        double h = x[i + 1] - x[i];
        sc.a[i] = (M[i + 1] - M[i]) / (6.0 * h);
        sc.b[i] =  M[i] / 2.0;
        sc.c[i] = (y[i + 1] - y[i]) / h - h * (M[i + 1] + 2.0 * M[i]) / 6.0;
        sc.d[i] =  y[i];
    }
    return sc;
}

// ---------------------------------------------------------------------------
// Reconstruct y values at x_resampled using piecewise cubic polynomials.
// ---------------------------------------------------------------------------
static std::vector<double> ReconstructY(const std::vector<double>& x_resampled,
                                         const std::vector<double>& x,
                                         const SplineCoeffs& sc)
{
    int size        = (int)x.size();
    int new_size    = (int)x_resampled.size();
    std::vector<double> y_out(new_size, 0.0);

    // For each interval [x[i-1], x[i]], fill matching resampled points.
    // Because x_resampled is sorted we could binary-search, but a simple
    // linear scan over intervals is fine and cache-friendly for 1-D data.
    for (int i = 1; i < size; ++i) {
        double x0 = x[i - 1], x1 = x[i];
        for (int j = 0; j < new_size; ++j) {
            double xr = x_resampled[j];
            if (xr >= x0 && xr <= x1) {
                double dx = xr - x0;
                y_out[j] = ((sc.a[i-1]*dx + sc.b[i-1])*dx + sc.c[i-1])*dx + sc.d[i-1];
            }
        }
    }
    return y_out;
}

// ---------------------------------------------------------------------------
// 1-D cubic spline interpolation: resample `data` to `new_size` points.
// ---------------------------------------------------------------------------
static std::vector<double> cubic_spline_1D(const std::vector<double>& data, int new_size)
{
    int size = (int)data.size();
    // const h = 1 (uniform grid)
    int n = size - 2;                              // interior knots
    std::vector<double> a_tri(n, 1.0 / 6.0);
    std::vector<double> b_tri(n, 2.0 / 3.0);
    std::vector<double> c_tri(n, 1.0 / 6.0);
    std::vector<double> rhs(n);
    for (int i = 0; i < n; ++i)
        rhs[i] = data[i + 2] - 2.0 * data[i + 1] + data[i]; // h=1 → /h

    std::vector<double> M = Thomas(n, a_tri, b_tri, c_tri, rhs);

    std::vector<double> x(size);
    for (int i = 0; i < size; ++i) x[i] = (double)i;

    SplineCoeffs sc = compute_spline_coefficients(x, data, M);

    std::vector<double> x_resampled = ResampleNoArtifacts(new_size, size);
    return ReconstructY(x_resampled, x, sc);
}

// ---------------------------------------------------------------------------
// 2-D bicubic spline: row-pass then column-pass, both parallelised with OpenMP.
// image_in  : row-major, height × width
// image_out : row-major, new_height × new_width   (caller allocates)
// ---------------------------------------------------------------------------
static void cubic_spline_2D(const std::vector<double>& image_in,
                              int width, int height,
                              std::vector<double>& image_out,
                              int new_width, int new_height)
{
    // ---- Step 1: interpolate each row from width → new_width ----
    // row_interp : height × new_width
    std::vector<double> row_interp((size_t)height * new_width);

    #pragma omp parallel for schedule(dynamic)
    for (int r = 0; r < height; ++r) {
        std::vector<double> row(width);
        for (int c = 0; c < width; ++c)
            row[c] = image_in[(size_t)r * width + c];

        std::vector<double> out_row = cubic_spline_1D(row, new_width);
        for (int c = 0; c < new_width; ++c)
            row_interp[(size_t)r * new_width + c] = out_row[c];
    }

    // ---- Step 2: interpolate each column from height → new_height ----
    // image_out : new_height × new_width
    image_out.assign((size_t)new_height * new_width, 0.0);

    #pragma omp parallel for schedule(dynamic)
    for (int c = 0; c < new_width; ++c) {
        std::vector<double> col(height);
        for (int r = 0; r < height; ++r)
            col[r] = row_interp[(size_t)r * new_width + c];

        std::vector<double> out_col = cubic_spline_1D(col, new_height);
        for (int r = 0; r < new_height; ++r)
            image_out[(size_t)r * new_width + c] = out_col[r];
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <original_file> <output_file> <scale_factor>\n";
        return 1;
    }

    int scale = std::atoi(argv[3]);
    if (scale <= 0) {
        std::cerr << "Error: scale_factor must be a positive integer.\n";
        return 1;
    }

    std::string input_path  = std::string("Images/") + argv[1];
    std::string output_path = std::string("Images/") + argv[2];

    // ---- Load image ----
    int width, height, channels;
    unsigned char* raw = stbi_load(input_path.c_str(), &width, &height, &channels, 3);
    if (!raw) {
        std::cerr << "Error: could not open " << input_path << "\n";
        return 1;
    }
    int new_width  = width  * scale;
    int new_height = height * scale;

    std::cout << "Input : " << width << "x" << height
              << "  →  Output: " << new_width << "x" << new_height
              << "  (scale=" << scale << ")\n";
    std::cout << "Threads available: " << omp_get_max_threads() << "\n";

    // ---- Split into channels (double arrays) ----
    size_t npix = (size_t)width * height;
    std::vector<double> ch_R(npix), ch_G(npix), ch_B(npix);
    for (size_t i = 0; i < npix; ++i) {
        ch_R[i] = raw[3 * i + 0];
        ch_G[i] = raw[3 * i + 1];
        ch_B[i] = raw[3 * i + 2];
    }
    stbi_image_free(raw);

    // ---- Interpolate each channel ----
    auto t0 = std::chrono::high_resolution_clock::now();

    size_t new_npix = (size_t)new_width * new_height;
    std::vector<double> out_R(new_npix), out_G(new_npix), out_B(new_npix);

    // The three channels are independent — run them concurrently via sections.
    #pragma omp parallel sections
    {
        #pragma omp section
        cubic_spline_2D(ch_R, width, height, out_R, new_width, new_height);

        #pragma omp section
        cubic_spline_2D(ch_G, width, height, out_G, new_width, new_height);

        #pragma omp section
        cubic_spline_2D(ch_B, width, height, out_B, new_width, new_height);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();

    // ---- Pack output image ----
    std::vector<unsigned char> out_raw(new_npix * 3);
    for (size_t i = 0; i < new_npix; ++i) {
        auto clamp = [](double v) -> unsigned char {
            return (unsigned char)std::max(0.0, std::min(255.0, v));
        };
        out_raw[3 * i + 0] = clamp(out_R[i]);
        out_raw[3 * i + 1] = clamp(out_G[i]);
        out_raw[3 * i + 2] = clamp(out_B[i]);
    }

    // ---- Save ----
    // Choose writer based on extension
    std::string ext = output_path.substr(output_path.find_last_of('.') + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    int ok = 0;
    if (ext == "png")
        ok = stbi_write_png(output_path.c_str(), new_width, new_height, 3,
                            out_raw.data(), new_width * 3);
    else if (ext == "jpg" || ext == "jpeg")
        ok = stbi_write_jpg(output_path.c_str(), new_width, new_height, 3,
                            out_raw.data(), 95);
    else if (ext == "bmp")
        ok = stbi_write_bmp(output_path.c_str(), new_width, new_height, 3,
                            out_raw.data());
    else {
        // Default to PNG and warn
        std::cerr << "Warning: unknown extension '" << ext
                  << "', writing as PNG.\n";
        ok = stbi_write_png(output_path.c_str(), new_width, new_height, 3,
                            out_raw.data(), new_width * 3);
    }

    if (!ok) {
        std::cerr << "Error: could not write " << output_path << "\n";
        return 1;
    }

    std::cout << "Resampled image saved to " << output_path << "\n";
    std::cout << "Time taken: " << elapsed << "s\n";
    return 0;
}