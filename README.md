# watersheds
Codes from "A universal approach for drainage basins".

# 1. Fixed boundary conditions (FBC)

## 1.1 Fractional Brownian Motion (fBm) code

fbm2d.sh - Compile and run the FORTRAN code. If necessary, edit the input parameters.  
fbm2d.f90 - FORTRAN code for fBm landscape.  

How should you run the code?  

./fbm2d.sh

## 1.2 Input file layout

lx ly  
i j h_ij  
.  
.  
.  
  
(See fbm2d_l0032_h0.7_s1231.dat)  

## 1.3 Invasion Percolation Based Algorithm (IPBA)

watershed_FBC.sh - Compile and run the codes. If necessary, edit the input parameters.  
watershed_FBC.c - C code for IPBA with FBC.  
plot.py - Python code for ploting.  

How should you run the code?  

./watershed_FBC.sh

# 2. Periodic boundary conditions (PBC)

## 2.1 Input file layout

lx ly lon lat delta R  
i j h_ij  
.  
.  
.  
  
(See GEBCO_2014_2D_90N_90S_f8.dat)

## 2.2 Invasion Percolation Based Algorithm (IPBA)

watershed_PBC.sh - Compile and run the codes. If necessary, edit the input parameters.  
watershed_PBC.c - C code for IPBA with PBC.  
plot.py - Python code for ploting.  
  
How should you run the code?  

./watershed_PBC.sh

# 3. Watershed lines

The watershed lines are easily obtained from the "sigma_*.dat" files.
