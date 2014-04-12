Free Energy calculation for binding free energy
========================================

This code automates the free energy calculations for binding between
two organic molecules at different concentrations, components and
state points.


Programs
----------------------------------------

*packmol
*gromacs-4.6.1 
*plumed-1.3
*python


Workflow
----------------------------------------

1. Combine structures with packmol 
2. Solvate with gen_box and balance charge (3 nm)
3. NPT simulation to achieve desired pressure, 10ns
4. Well-tempered metadynamics in NVT with Bussi-Donadio-Parrinello thermostat


More details
----------------------------------------

TIP4P2005-EW, OPLS-AA, 100mM salt, 200ps Equil, 300K, tau=0.5ps, 8A
cutoff, LINCS (bonds+angles), 2fs timestep, 2D WMetaD C-C-(center)
orientation, 1D min dist for multiple ions. Try 40 ns to begin with.

