### Date created
2010-05-27

### Changelog

#### 1.0.0 - 2020-05-27

Added original source code

## Extended incompressible Large-Eddy Simulation (LES) library for OpenFOAM

## Description

The extended library is based on my research work on LES in 2005-2010. I have modified some existing models (ending with "2" or "3") so that the standard Smagorinsky's coefficient (Cs) is calculated instead of using ck and ce model coefficients.

### Included Models

- *Smagorinsky3* - Smagorinsky-Germano model (constant Cs)
- *dynSmagorinsky2* - Germano-Lilly dynamic model (globally averaged Cs)
- *locDynSmagorinsky* - Germano-Lilly localized dynamic model (local Cs, with or without backscatter)
- *oneEqEddy2* - ksgs equation model
- *dynOneEqEddy2* - Dynamic ksgs equation model
- *locDynOneEqEddy2* - Localized dynamic ksgs equation model
- *PiomelliLiu* - Piomelli & Liu (1995) localized dynamic model (with or without backscatter)
- *dynYouMoin* - You& Moin (2006) model (global average)
- *lagrangianScaleInvariant* - Lagrangian scale-invariant model, Meneveau et al. (1996)
- *lagrangianScaleInvariantAlgebraic* - Algegraic Lagrangian scale-invariant model, Meneveau et al. (1996)

### Requirements

The library was developed for [OpenFOAM 1.7.x](https://github.com/OpenCFD/OpenFOAM-1.7.x). Small modifications would be for OpenFOAM 2.x and newer versions.

### Additional Sources

1. The OpenFOAM Foundation https://openfoam.org
