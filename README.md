# Bicubic Spline Interpolation

**Bicubic spline interpolation** is a mathematical method commonly used for resizing images and performing smooth interpolation of pixel values. This project provides both a Python and a high-performance parallel C++ implementation, benchmarked side by side.

---

## Table of Contents

1. [Introduction](#introduction)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Performance Benchmarks](#performance-benchmarks)

---

## Introduction

Bicubic interpolation considers the values of nearby pixels to compute new pixel values during image scaling. Compared to simpler methods such as bilinear interpolation, bicubic produces noticeably smoother results, making it the standard choice for high-quality image resizing.

This project includes:

- A **sequential Python** implementation (`bsisequential.py`)
- A **parallel Python** implementation (`bsiparallel.py`) — processes each RGB channel concurrently
- A **parallel C++** implementation (`Cpp/bsi_parallel.cpp`) — uses OpenMP to parallelise both the row-pass and column-pass of the 2D spline, as well as the three colour channels

---

## Prerequisites

**Python:**
- Python 3.6 or newer
- `Pillow` and `numpy`

**C++:**
- GCC with OpenMP support (`g++ -fopenmp`)
- `stb_image` / `stb_image_write` headers (see `make stb`)

---

## Installation

**Python dependencies:**

```bash
pip install pillow numpy
```

**C++ build:**

```bash
make stb     # download stb headers into includes/
make         # compile with -O2 -fopenmp
```

---

## Usage

**Python — sequential:**

```bash
python3 bsisequential.py <input_image> <output_image> <scale_factor>
```

**Python — parallel:**

```bash
python3 bsiparallel.py <input_image> <output_image> <scale_factor>
```

**C++:**

```bash
./bsi_parallel <input_image> <output_image> <scale_factor>
```

All implementations read from and write to the `Images/` subdirectory.

**Benchmark against each other:**

```bash
chmod +x benchmark.sh
./benchmark.sh PIXIL.png
```

---

## Performance Benchmarks

Benchmarks were run on a `8×8` pixel input image (`PIXIL.png`) scaled to various sizes. Each configuration was repeated **3 times**; the table shows the **best** (lowest) elapsed time.

| Scale | Output size | Python (s) | C++ (s)   | Speedup   |
|------:|------------:|-----------:|----------:|----------:|
|     2 |     16×16   |      0.58  | 0.000506  | **1147×** |
|     5 |     40×40   |      0.60  | 0.000664  |  **904×** |
|    10 |     80×80   |      0.63  | 0.004778  |  **132×** |
|    50 |   400×400   |      0.93  | 0.010418  |   **89×** |
|   100 |   800×800   |      1.57  | 0.046915  |   **33×** |
|   200 | 1600×1600   |      3.34  | 0.125016  |   **27×** |

**Overall (all scales combined):**

| | Python | C++ |
|---|---|---|
| Total time | 7.65 s | 0.188 s |
| **Average speedup** | — | **40.6×** |

> The C++ implementation is consistently **25–1100× faster** than the Python equivalent across all tested scales. The extreme speedup at small scales is dominated by Python's interpreter and import overhead. At larger scales — where computation dominates — the C++ advantage stabilises around **27–33×**, reflecting the efficiency of OpenMP parallelism and compiled native code.
