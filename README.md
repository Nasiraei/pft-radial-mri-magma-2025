# Polar Fourier Transform (PFT) for 2D Radial MRI Reconstruction

This repository contains the reference code and documentation for:

> **"Polar Fourier Transform in Practice: Its Efficiency and Characteristics in Reconstructing Radially Acquired MRI Images"**  
> Accepted in *MAGMA* (Magnetic Resonance Materials in Physics, Biology and Medicine), 2025.  
> Authors: Fatemeh Rastegar Jooybari, Ali Aghaeifar, Elham Mohammadi, Klaus Scheffler, Abbas Nasiraei Moghaddam  
> DOI: [10.1007/s10334-025-01284-w](https://doi.org/10.1007/s10334-025-01284-w)

---

## 🔍 Overview

This code implements the **recursive numerical Hankel transform**, which forms the central mathematical component of our **Polar Fourier Transform (PFT)**-based reconstruction algorithm for radially acquired 2D MRI data.

The PFT approach avoids interpolation in k-space by working in native polar coordinates. This module is designed for integration into a larger image reconstruction pipeline (e.g., Siemens ICE), but it can also be adapted for other platforms.

---

## 📁 Repository Contents

- `pft_core.cpp` — Core C++ implementation of the recursive Hankel transform used for radial k-space to image-space mapping
- `README.md` — This file
- `LICENSE` — MIT License for open use

> ⚠️ **Note**: This repository does *not* include Siemens ICE wrapper or integration code. To use this algorithm in a clinical scanner environment, access to Siemens ICE and a C2P agreement is required. The code here is platform-agnostic and intended for academic and research use.

---

## 🛠 Requirements

- C++17 or later
- Standard libraries only (no external dependencies)
- A C++ build environment (e.g., `g++`, `clang`, or MSVC)

### 🔧 Example compilation:

```bash
g++ -std=c++17 -O3 pft_core.cpp -o pft_core
