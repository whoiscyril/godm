# Global Optimisation of Defects with Mott-Littleton

## Description

A package for predicting the structures of point defects in crystals using global optimisation algorithms and Mott-Littleton techniques.

Highlights:

1. Prediction of tentative global and local minima for point defects at true infinite dilution.
2. An end-to-end pipeline for structure generation, energy evaluation, and results parsing.
3. Support for random quench and basin-hopping algorithms.
4. Integration with the General Lattice Utility Program (GULP).

This package is a refactored version of MLGOptimiser, intended to make the code more usable for a wider audience. It is designed to have a single entry point through the main function. Use of other methods is at the user's own risk.

## Getting Started

### Dependencies

- General Lattice Utility Program (GULP).
- Optional: [KLMC3](https://github.com/sweetmixture/KLMC3) and [libgulp](https://github.com/sweetmixture/KLMC3-libgulp-6.1.2).
- For Python dependencies, see [pyproject.toml](pyproject.toml).

### Installing

- Recommended: `uv add godm`
- With pip: `pip install godm`
- For an editable installation: `uv pip install -e .`

### Running the program

To run the program, you need a directory named `inputs` containing `CONTROL` and `CONFIG` files.

## Help

For help, run `godm -h`.

## Authors

Cyril Xu — University College London — uccazxu@ucl.ac.uk

## Version History

- 0.1: Initial release.

## License

This project is licensed under the MIT License. Dependencies are subject to their own license requirements.

