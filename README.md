# Ice-GasIce-Crack-ElastoPlastic-Model
Code for elasto-plastic initiation model at unsaturated inclined open crack tip under ice-gas-ice pressure
This project provides a complete set of numerical simulation tools, including COMSOL multiphysics simulation models and MATLAB script files, specifically designed to analyze the elasto-plastic initiation criteria, stress field distribution, and crack evolution law at the tip of unsaturated inclined open cracks under the coupling effect of ice-gas-ice three-phase pressure. All models have been validated and can directly reproduce relevant research results.
File Structure:
COMSOL Model Files	2	Used to build the unsaturated crack geometric model, define ice-gas-ice pressure boundary conditions, and solve the elasto-plastic stress field and crack initiation criteria
MATLAB Script Files	10	Used for post-processing COMSOL calculation results, extracting mechanical parameters at crack tips, plotting evolution curves, and quantitatively analyzing initiation thresholds and sensitivity analysis
Software Requirements:
COMSOL Multiphysics: 6.3 64-bit (Core calculation software, requires installation of the Structural Mechanics Module)
MATLAB: R2020b or later (Compatible with COMSOL data export format, it is recommended to install the Optimization Toolbox)

This project code is open source under the MIT License, see the LICENSE file in the root directory of the repository for details.
