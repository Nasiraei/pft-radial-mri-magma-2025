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


## 🚀 Usage
This standalone module implements the Hankel transform component of PFT reconstruction. It does not perform full image reconstruction on its own. Instead, it can be:

Reviewed to understand the Hankel transform algorithm

Integrated into your own reconstruction pipeline

Adapted for use in other non-Cartesian MRI applications

✍️ Note: The current version is lightly documented. Additional comments and examples will be added in future updates.

## 📜 Citation
If you use this code in your research or build upon it, please cite the following:

Fatemeh Rastegar Jooybari, Ali Aghaeifar, Elham Mohammadi, Klaus Scheffler, Abbas Nasiraei Moghaddam
Polar Fourier Transform in Practice: Its Efficiency and Characteristics in Reconstructing Radially Acquired MRI Images.
MAGMA – Magnetic Resonance Materials in Physics, Biology and Medicine, 2025.
DOI: 10.1007/s10334-025-01284-w

## 🛡 License

This project is licensed under the MIT License  **with the restriction** that you are free to use, modify, and distribute this code for **academic and non-commercial research purposes**, with proper attribution.

> **Note**: This license does **not** grant any rights under any patent.  
> The underlying method is protected by **US Patent No. US9835708B2 **.  
> For patent licensing (especially for commercial use), please contact the corresponding author.

## 🔏 Patent Notice

The reconstruction algorithm implemented in this repository is based on subject matter covered by **US Patent No. US9835708B2**, owned by Abbas Nasiraei Moghaddam.  
The code provided here is intended for demonstration and academic use only. **No rights under the patent are granted through this repository.**  
For **commercial use** or the **distribution of derivative implementations**, a separate license may be required.  
For licensing inquiries, please contact the corresponding author.

## 📬 Contact
For questions or collaboration inquiries, please contact:
nasiraei@aut.ac.ir or abbas.nasiraei@gmail.com

