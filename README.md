# Time Series Covariance Matrix Analysis

## Project Overview

This project focuses on the analysis of time series data, specifically looking at the covariance matrices of NOAA temperature data. The goal is to study the spectral properties of these matrices, compare them with random matrix theory (Marchenko-Pastur distribution), and investigate the exponents of the return distributions for different regions.

## Physical Systems

The project also extends its analysis to the covariance matrices derived from time series of various physical systems:

- **Spin Systems**: Analysis of **Ising** and **Blume-Capel** models.
- **Contact Process**: Investigation of the **Contact Process** (specifically 1D with diffusion).

## Features

- **NOAA Data Processing**: Tools to handle and process NOAA temperature data (`DataIO.jl`).
- **Time Series Analysis**: Functions for normalizing time series, calculating covariance matrices, and computing eigenvalues (`TimeSeries.jl`).
- **Statistical Analysis**: Comparison with Marchenko-Pastur distribution and bootstrapping methods.
- **Visualization**: Scripts to generate plots (e.g., eigenvalue spectra) are located in the `plots` directory (implied by file list).

## Installation

This project is written in Julia. To set up the environment:

1. Clone the repository.
2. Open a Julia REPL in the project directory.
3. Enter Pkg mode by pressing `]`.
4. Activate the project environment:

    ```julia
    activate .
    ```

5. Instantiate the dependencies:

    ```julia
    instantiate
    ```

## Usage

The project contains several scripts and notebooks for analysis.

### Core Modules

- **`DataIO.jl`**: Handles file I/O, filename parsing, and data loading (support for `.cbor.gz`, `.jld2`, etc.).
- **`TimeSeries.jl`**: Contains core functions for:
  - `normalize_ts`: Normalize time series.
  - `covariance_matrix`: Compute covariance matrices.
  - `covariance_matrix_eigvals`: Compute eigenvalues of covariance matrices.
  - `marchenko_pastur`: Generate Marchenko-Pastur distribution.
