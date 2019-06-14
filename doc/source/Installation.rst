.. _Installation:

Installation 
------------

To install HERCULES, you need a Linux environment and these software packages: 

1. A GNU or Intel C compiler, e.g., gcc or icc. 
2. A GNU or Intel Fortran compiler, e.g., gfortran or ifort 
3. A MPI software, e.g., openmpi, mvapich2, or mpich 

**NOTE**: if you use the GNU C and Fortran compilers, your MPI software should be also compiled by GNU.

The installation of HERCULES is straightforward. If you use **GNU compilers**, run::

    sh install_GNU.sh

The installation will be automatically done. 

**NOTE**: HERCULES depends on two external libs: FFTW and 2DECOMP_FFT. So when you run install_GNU.sh, these two libs will be (automatically) compiled first. 

If you use **Intel compilers**, run this instead::

    sh install_Intel.sh

Similarly, you can install HERCULES on **Cray** by running::

    sh install_Cray_Intel.sh

After the installation is done, an executive named **hercules.exe** should be generated. This is the main program you will use for DNS simulations. In addition, you should see a file named **parameters.input** which defines the input parameters for the simulations. HERCULES comes with a default **parameters.input** file for DNS of plane closed-channel flows at Re_tau=180.