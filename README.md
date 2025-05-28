# Orbit Recovery for Spherical Functions

This repository contains MATLAB code for reconstructing 3D volumes from spherical harmonic expansions using frequency marching, with applications to biological structures such as TRPV1 and the ribosome.

The repository is used to reproduce the numerical experiments in the paper "Orbit recovery for spherical functions".

## 🧪 Overview

The pipeline includes:
- Volume loading and normalization
- Projection onto spherical harmonics + Fourier-Bessel basis
- Frequency marching reconstruction
- Quantitative error evaluation and visualization

## 📁 Structure
Main scripts for paper's code reproducibility.


## 📦 Dependencies

This codebase requires the following toolboxes and libraries:

### MATLAB Toolboxes
- **Optimization Toolbox**

### External Libraries
1. **[ASPIRE](https://github.com/ComputationalCryoEM/ASPIRE-Python) - MATLAB version**
   - Used for `cryo_downsample` and volume preprocessing
   - Add to path with:  
     ```matlab
     addpath(genpath('/path/to/ASPIRE'))
     ```

2. **[EasySpin](https://easyspin.org/)**
   - Required for advanced matrix operations
   - Add to path with:  
     ```matlab
     addpath('/path/to/easyspin')
     ```

3. **[Spherical Harmonics Toolbox](https://github.com/polarch/Spherical-Harmonic-Transform/)**
   - Custom or third-party toolbox for spherical harmonic evaluations
   - Add to path with:  
     ```matlab
     addpath('/path/to/SphericalHarmonics')
     ```

## ▶️ Running the Code

Run any of the main scripts directly in MATLAB:

```matlab
main_TRPV1_shells

