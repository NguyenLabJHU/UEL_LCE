# LCE User Element (UEL) for Abaqus

A user element (UEL) implementation of a liquid crystal elastomer (LCE) model for Abaqus.  
This code provides a finite element framework to simulate the time-dependent viscoelastic response of LCEs.

---

## Overview
This repository contains the UEL subroutine and supporting files for simulating LCEs within Abaqus.

The implementation is based on the nonlinear viscoelasticity theory developed in:

> Wang, Zheliang, et al. *"A nonlinear viscoelasticity theory for nematic liquid crystal elastomers."*  
> Journal of the Mechanics and Physics of Solids 163 (2022): 104829.  
> https://doi.org/10.1016/j.jmps.2022.104829

---

## Features
- Finite strain viscoelastic formulation  
- Director–network coupled constitutive model  
- Ready-to-use UEL code compatible with Abaqus  

---

## Usage
1. Compile the UEL subroutine with Abaqus Fortran compiler (Intel/ifort recommended).  
2. Include the compiled UEL when submitting Abaqus jobs, e.g.:

   ```bash
   abaqus job=example inp=example.inp user=UEL_LCE.for
