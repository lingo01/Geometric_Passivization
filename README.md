
# Geometric_Passivization

Data and code of the article "Positive Damping Region: A Graphic Tool for Passivization Analysis with Passivity Index"

This source code repository accompanies the following paper:

> X. Peng, X. Ru, Z. Li, J. Zhang, X. Chen, and F. Liu, "Positive Damping Region: A Graphic Tool for Passivization Analysis with Passivity Index," 2025.

The full paper and the source code can be found at: https://github.com/lingo01/Geometric_Passivization

For any questions or uses of the source codes, please feel free to contact the first author, Xiaoyu Peng (pengxy19@tsinghua.org.cn), and the corresponding author, Feng Liu (lfeng@mail.tsinghua.edu.cn).

**CITATION**: If you use this code in your work, whether directly or indirectly, please cite the above paper.

**LICENSE**: This work is licensed under the MIT License. See the [LICENSE](LICENSE) file in the repository for details.


# Introduction to the source code

All code is written in `MATLAB` and can be run directly in `MATLAB 2018a` or later versions.

This repository provides the implementation for generating all main figures in the paper, enabling full reproduction of the results. Each script corresponds to a specific figure as detailed below.


## Code-Figure Correspondence

| Figure in Paper | MATLAB Script           | Description                                      |
|-----------------|------------------------|--------------------------------------------------|
| Fig.3         | `main_Illu_2_Ctrl.m`   | Nichols chart with positive damping region overlay |
| Fig.4        | `main_Simu_1_SISO.m`   | SISO system Nyquist plots with positive damping region overlay |
| Fig.5        | `main_Simu_1_MIMO.m`   | MIMO system Rayleigh quotient region  with positive damping region overlay |
| Fig.6        | `main_Simu_2.m`        | SISO system Nyquist plots with positive damping region (defined by differential passivity, or equivalently, negative imaginariness) overlay |
| Fig.7        | `main_Simu_3_EX2.m`    | 3D positive damping region of the definition Example 2 and Nyquist overlay |
| Fig.8        | `main_Simu_3_EX3.m`    | 3D positive damping regions of the definitions Example 2 and Example 3 and Nyquist overlay |


## How to Generate Each Figure

To reproduce a figure from the paper, simply run the corresponding `.m` script in MATLAB. For example, to generate Figure 3, run:

```matlab
main_Illu_2_Ctrl
```

Each script will generate the figure as shown in the paper. Some scripts will also save the figure as a PDF or PNG file in the current directory (see script comments for details).

### Supporting Functions

The following supporting function files are required and should be present in the same directory:

- `func_PRPregion_EX2.m`
- `func_PRPregion_EX3.m`
- `func_Nyquist.m`
- `func_checkNyquist.m`

These functions are called by the main scripts to generate regions, perform Nyquist calculations, and check passivity.


## Requirements

-  `MATLAB 2018a`  or newer version
- Control System Toolbox  of MATLAB software for calling `tf`, `bode`, `nyquist`


## Notes

- The code is organized for clarity and reproducibility of the paper's results.
- All figures are generated using default parameters as in the paper. You may adjust parameters in the scripts for further exploration.


## Reproducing the Article

To reproduce all main results and figures:

1. Ensure all `.m` files are in the same directory.
2. Open MATLAB and set the working directory to this folder.
3. Run each script as listed in the table above to generate the corresponding figure.

