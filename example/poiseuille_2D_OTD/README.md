# Linear simulation of plane Poiseuille flow to compute the optimally time
dependent (OTD) modes

This example presents the general concept of the OTD method to find an
orthonormal basis aligned with the most unstable direction in the flow.
The example also includes the possibility to compute reduced-order FTLEs.

Plane poiseuille flow at Re = 5000 (subcritical)
The initial conditions for the perturbations is chosen to be the leading 
optimal initial conditions for maximum transient growth at t=25.

# Tool list:
- otd.f		: Main routines to compute the OTD modes
- otd_tools.f	: Auxiliary routines for OTD computation
- linalg.f	: Linear algebra routines using LAPACK
- inc_src/OTD	: unique include file for OTD module

# Requirements
This module is developped for Nek5000 v17

# Dependencies
This module requires the KTH framework to run
