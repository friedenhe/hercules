!*************************************************************************
! This file is part of HERCULES.
!
! HERCULES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! HERCULES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2015, Ping He, MEAS, NCSU
!*************************************************************************



module module_variables

    use module_parameters

    implicit none

    ! By default the variables are on z pencil

    ! u, v, w: velocity in the x, y, and z directions, p: pressure.
    ! phi: pseudo pressure, t: temperature
    real(mytype),dimension(:,:,:),allocatable,save :: u, v, w, p, phi, t

    ! interpolated variables using either 2nd or 4th order interpolation scheme
    ! e.g., uhx is the half-way interpolated value of u along the x direction,
    ! i.e., uhx(i,j,k) is located at the mid point between u(i,j,k) and u(i-1,j,k)
    real(mytype),dimension(:,:,:),allocatable,save :: uhx,uhy,uhz,vhx,vhy,vhz
    real(mytype),dimension(:,:,:),allocatable,save :: whx,why,thx,thy,thz

    ! mean variables and turbulence statistics.
    real(mytype),dimension(:,:,:),allocatable,save :: mean_u,mean_v,mean_w,mean_t
    real(mytype),dimension(:,:,:),allocatable,save :: mean_uu,mean_vv,mean_ww
    real(mytype),dimension(:,:,:),allocatable,save :: mean_uw,mean_vw,mean_tw,mean_tt

    ! ff, gg, hh, ss: differential convective terms used in
    ! momentum and temperature equations for RK3 and AB2 schemes
    ! 0 means value at the previous time step (n-1)
    real(mytype),dimension(:,:,:),allocatable,save :: ff,ff0,gg,gg0,hh,hh0,ss,ss0

    ! right hand side terms for momentum equations, 
    ! for pressure and temperature equations,
    ! we reuse these variables to save memory
    real(mytype),dimension(:,:,:),allocatable,save :: qu,qv,qw

    ! coefficients for Thomas algorithm, the 1st and 2nd dimensions is for l and m
    ! note: a,b,c have ghost cells in k. aa,bb,cc are temporary arrays
    ! the a,b,c values for each equations are stored in aa_u,bb_u,cc_u, etc.
    real(mytype),   dimension(:,:,:),allocatable,save :: aa,cc
    complex(mytype),dimension(:,:,:),allocatable,save :: bb
    complex(mytype),dimension(:,:,:),allocatable,save :: bb_u,bb_w,bb_p,bb_t
    real(mytype),   dimension(:),    allocatable,save :: aa_u,cc_u,aa_w,cc_w, &
                                                         aa_p,cc_p,aa_t,cc_t

    ! temporary variables (z pencil), note that they do not have any ghost cells 
    complex(mytype),dimension(:,:,:),allocatable,save :: c_var,c_src
    ! temporary array for mpi_2dfft, they are on x and y pencil
    complex(mytype),dimension(:,:,:),allocatable,save :: c_tmp1,c_tmp2

    !!!! variables only used for spectral method !!!
    ! Fourier coefficients of the fundamental variables (z pencil) 
    complex(mytype),dimension(:,:,:),allocatable,save :: c_u,c_v,c_w,c_p,c_phi,c_t
    ! Fourier coefficients of the right hand side terms and ff, gg, hh,ss
    complex(mytype),dimension(:,:,:),allocatable,save :: c_qu,c_qv,c_qw
    complex(mytype),dimension(:,:,:),allocatable,save :: c_ff,c_gg,c_hh,c_ss

end module
