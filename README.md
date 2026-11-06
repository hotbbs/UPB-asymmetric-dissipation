# UPB with Asymmetric Dissipation

This repository contains the data and the numerical simulation codes associated with our manuscript:

**"Unconventional Photon Blockade in Second-Order Nonlinearity Coupled Cavities with Asymmetric Dissipation"**

## 📂 Repository Contents

*   **`/figures`**: Contains the final high-resolution image files for all four figures presented in the manuscript (Fig. 1 to Fig. 4).
*   
*   **`/codes`**: Contains the MATLAB scripts used for numerical simulations.
    *   `Fig2a_FB_vs_g.m`: Calculates the logarithmic second-order correlation function $\log_{10}[g^{(2)}(0)]$ in the $F_B-g$ parameter space.
    *   `Fig2b_kA_vs_kB.m`: Calculates the $\log_{10}[g^{(2)}(0)]$ map in the $\kappa_A-\kappa_B$ parameter space.
    *   `Fig3_detuning_dependence.m`: Calculates and plots the dependence of $\log_{10}[g^{(2)}(0)]$ and the mean photon number $\log_{10}\langle N_a \rangle$ on the normalized detuning $\Delta/\kappa_A$ for both Family 1 and Family 2.
    *   *(Note: Fig. 1 is a schematic diagram and Fig. 4 is an energy-level diagram, both generated using illustration software, and thus do not require numerical computation codes.)*

## 💻 Computational Environment

The numerical simulations were carried out using MATLAB. The master equation for the open quantum system is solved using the open-source **Quantum Optics Toolbox for MATLAB** (S. M. Tan, *J. Opt. B*, **1**, 424 (1999)).

*   **Required Software**: MATLAB
*   **Required Toolbox**: Quantum Optics Toolbox (by Sze M. Tan)

## 🚀 Usage

1. Ensure that the Quantum Optics Toolbox is installed and added to your MATLAB path.
2. Run the scripts in the `/codes` directory. The physical parameters (e.g., $F_A$, $\Delta$, $\kappa_A$, $\kappa_B$) are pre-defined at the beginning of each script and can be easily modified to reproduce the results.
3. The scripts will output the numerical data and generate the corresponding 2D density plots or 1D line graphs.

## 📝 Citation

If you use these codes or data in your work, please cite our paper:

> *[在此处填入你论文的完整引用格式，例如：]*
> Z. Q. Zhang and Y. H. Zhou, "Unconventional Photon Blockade in Second-Order Nonlinearity Coupled Cavities with Asymmetric Dissipation", Phys. Rev. A (2026).
