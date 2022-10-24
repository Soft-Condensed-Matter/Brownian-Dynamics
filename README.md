
<h1 align="center">
  <br>
  <a href="https://github.com/Soft-Condensed-Matter/Monte-Carlo/blob/master/MC.png"><img src="https://github.com/Soft-Condensed-Matter/Monte-Carlo/blob/master/MC.png" alt="Markdownify" width="200"></a>
  <br>
  Brownian Dynamics
  <br>
</h1>

<h3 align="center">Monte Carlo homogeneous fluid simulation  </h3>
<h4 align="center">FORTRAN code to simulate two-dimensional Lennard-Jones colloids   </h4>

<p align="center">
  <a href="#code">Code</a> •
  <a href="#compile">Compile</a> •
  <a href="#execute">Execute</a> •
    <a href="#visualization">Visualization</a> •
  <a href="#observables">Observables</a>
</p>

## Code
* BD_2D.f90
  - Disk colloidal particles interacting with the Lennard-Jones potential using the Ermak and MaCammon alghoritm. The code is able to compute structural and dynamical properties such as the radial distribution function and the mean square displacement, respectively.The code also makes a file with frames that are used to build a the simulation movie with vmd.

## Compile
```bash
Serial:
   # GNU compiler
   $ gfortran -O3 BD_2D.f90 -o bdlj
   
   # Intel oneAPI
   $ ifort -O3 BD_2D.f90 -o bdlj
   
   # Nvidia HPC SDK
   $ nvfortran -O3 BD_2D.f90 -o bdlj

Parallel (share memory):
   # GNU compiler

```   

## Execute
* Simulation parameters
Simulation conditions as the number of particles, the volume fraction, temperature and number of simulation steps are specified in the BD.inp file. 
<i>Each time the code is executed files with results are rewritten</i>

* Run simulation
```bash
# Run the code
$./bdlj
```

## Visualization
* Simulation movie file
By default the simulation movie file creation is not activated, to active it uncomment line 401 in file BD_2D.f90. Larger simulations and/or great number of particles will create larger files (>MB) that will be hard to manipulate in personal computers

* Simulation movie
The file is created in <i> *.xyz</i> format that could be visualizated in xmakemol or vmd
```bash
# Load movie file
$ vmd -f MCMovie.xyz
```

* Initial configuration
The code also creates a snapshot of the initial and final configuration configuration, to see them run
```bash
# Load initial config snapshot file
$ vmd -f BDSnap.xyz

# Load final config snapshot file
$ vmd -f BDPic.xyz
```

## Observables
* Dynamic properties
The mean square displacement is saved in file MSD.dat. This file can be plotted with any software for graphics, like gnuplot. The file saves the cartersian contributions of W(t) and the total contribution. To see the total MSD using gnuplot run
```bash
gnuplot> set logscale
gnuplot> p "MSD.dat" u 1:4  w l
```
	
* Microscopic properties
The radial distribution function is saved in file BDGr.dat and can be plotted with any graphical software like gnuplot
```bash
gnuplot> p "BDGr.dat" title '{/Symbol f}=0.4' with line line style 1 line width 2
```


## License

GNU LGPL-v3


## Cite this code
[![DOI](https://zenodo.org/badge/542237113.svg)](https://zenodo.org/badge/latestdoi/542237113)

* BibTex
```
@software{MC_HS2022,
  author = {Alexis Torres-Carbajal},
  doi = {10.5281/zenodo.7117628,
  month = {10},
  title = {{Hard-Sphere Monte Carlo}},
  url = {https://github.com/Soft-Condensed-Matter/Monte-Carlo},
  version = {1.0},
  year = {2022}
}
```

---

> Twitter [@Alpixels](https://twitter.com/Alpixels)
