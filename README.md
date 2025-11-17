# Time-Division Multiplexing for Parallel Transmission at Ultra-High Field with Limited RF Channels
**Felix Glang<sup>1,2*</sup>, Georgiy A. Solomakha<sup>1</sup>, Dario Bosch<sup>1,3,4</sup>, Klaus Scheffler<sup>1,3</sup>, Nikolai I. Avdievich<sup>1</sup>**\
<sup>1</sup>Magnetic Resonance Center, Max Planck Institute for Biological Cybernetics, Tübingen, Germany \
<sup>2</sup>Institute of Biomedical Imaging, Graz University of Technology, Graz, Austria \
<sup>3</sup>Department of Biomedical Magnetic Resonance, Eberhard Karls University Tübingen, Tübingen, Germany \
<sup>4</sup>MRI Core Facility of the Medical Faculty, University of Tübingen, Otfried-Müller-Straße 51, Tübingen, 72076, Germany 

![Schematic of the proposed multiplexing approach, which includes 8 high-power single pole double throw (SPDT) RF switches that are controlled by the trigger output of the scanner to route 8 transmit (Tx) channels alternately to each row of a double-row 16ch Tx array during the sequence.](tdm_ptx.png)

This repository contains MATLAB code for optimizing **time-division multiplexing parallel transmit (pTx)** pulses.
Time-division multiplexing allows driving a larger number of transmit elements (e.g. 16) with a smaller number of RF channels (e.g. 8), resulting in improved pTx performance.
This is enabled on the hardware level by high-power absorptive single pole double throw (SPDT) RF switches based on a set of lumped-element λ/4 transformers and PIN-diodes. These can alternately route 8 RFPAs to each row of a 16-element double-row Tx coil array, controlled by the optical trigger output of the scanner during the sequence.
The demo code optimizes and simulates both conventional and multiplexed kT points pTx pulses.

## Repository Structure

```
├── tdm_ptx_demo.m         % Demo script for basic use
├── data.mat               % Example dataset (sensitivity maps, masks, etc.)
├── utils/                 % Helper functions
│   └── [...]
├── loss_functions/        % Loss functions used in optimization
    └── [...]
```
The demo dataset is simulated for a double-row 16-element folded-end dipole transceiver array developed for human brain imaging at 9.4T in the Duke voxel model using CST Studio Suite 2021 (Dassault Systèmes, Vélizy-Villacoublay, France).

## Optimization Approach
Pulse design is formulated as a magnitude least-squares optimization based on the spatial domain method with VOP-based SAR constraints and optional joint optimization of k-space locations.
Solutions are obtained using the interior-point method implemented in Matlab's fmincon with user-supplied analytical Jacobians for the cost and constraint functions.

## Getting Started

### Prerequisites
* MATLAB (R2022b or newer recommended)
* Optimization Toolbox 

### Running the Demo
To run the demo, simply run ```tdm_ptx_demo.m```

This will:
* Load the sample dataset from `data.mat`
* Compute a standard CP mode excitation as a baseline
* Optimize a "full16" kT points pTx pulse, i.e., using all 16 coil elements simultaneously
* Re-arrange the obtained pulse into multiplexing shape without re-optimization ("multi8_conv")
* Optimize a multiplexed kT points pulse
* Compare flip angle maps, NRMSE and SAR values for the obtained pulses


## License

This project is licensed under the MIT License — see the LICENSE file for details.

## Citation

If you use this code in your research, please cite:

```
@article{Glang2025,
  author    = {Felix Glang and Georgiy A. Solomakha and Dario Bosch and Klaus Scheffler and Nikolai I. Avdievich},
  title     = {Time-Division Multiplexing for Parallel Transmission at Ultra-High Field with Limited RF Channels},
  year      = {2025},
  journal   = {Magnetic Resonance in Medicine},
  note      = {Under review},
}
```

## Contact
Feel free to open an issue or reach out for questions or suggestions!

