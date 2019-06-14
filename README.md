HERCULES: A *H*igh-order Finit*E*-difference Solve*R* for In*C*ompressible Bo*U*ndary *L*ay*E*r Flow*S*
-------------------------------------------------------------------------------------------------------

![](doc/source/images/image.png)

HERCULES is an open-source computational fluid dynamics (CFD) code for simulating incompressible boundary layer flows. HERCULES is developed for high-performance turbulence simulations, and it can be used to conduct direct numerical simulation (DNS) of neutrally and stably stratified turbulent open-/closed-channel flows, as well as Ekman layer flows. HERCULES is written in Fortran 90. It has been tested on a number of HPC systems, e.g., ARL HPC Excalibur, AFRL HPC Lightening, and TACC Stampede, and is shown to have excellent parallel efficiency with up to 10,000 CPU cores. 

HERCULES is configured for turbulent channel flow simulations in a rectangular wall-bounded domain with periodic boundaries in the horizontal directions. It solves the Navier-Stokes equations and the temperature equation using a high-order finite-difference approach. Spectral discretization can also be used for horizontal derivatives. 

Documentation
-------------

Refer to https://herculescode.rtfd.io for HERCULES installation and tutorials.

To build the documentation locally, go to the **doc** folder and run:

`./Allwmake`

The built documentation is located at **doc/HERCULESDoc.html**

Citation
--------

Ping He (2016), A high order finite difference solver for massively parallel simulations of stably stratified turbulent channel flows, *Computers and Fluids*, 127, 161-173. https://doi.org/10.1016/j.compfluid.2015.12.012

```
@article{HERCULESPaper,
	Author = {Ping He},
	Doi = {10.1016/j.compfluid.2015.12.012},
	Journal = {Computers \& Fluids},
	Pages = {161--173},
	Title = {A high order finite difference solver for massively parallel simulations of stably stratified turbulent channel flows},
	Volume = {127},
	Year = {2015}}
```

License
-------

GNU General Public License (GPL), version 3; see the files in **license** for details.

Contact and Feedback
--------------------
If you have questions, or want to contribute to the code, please contact: Ping He (friedenhe@gmail.com)

