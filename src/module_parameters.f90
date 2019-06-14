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



module module_parameters

    ! enable C standard for fftw
    use, intrinsic :: iso_c_binding
    ! 2d decomposition module
    use decomp_2d

    implicit none

    ! define pi
    real(mytype),parameter                 :: pi=3.1415926535897932_mytype

    !**************************************************************************
    !*********************    User Defined Parameters     *********************
    !**************************************************************************

    ! p_row/p_col: domain decomposition along the x/y directions
    ! make sure to use p_row*p_col processors when running the solver
    integer,protected,save                 :: p_row
    integer,protected,save                 :: p_col

    ! Lx, Ly, Lz: domain sizes in x, y, and z directions.
    ! x: streamwise, y: spanwise, z: vertical
    real(mytype),protected,save            :: lx
    real(mytype),protected,save            :: ly
    real(mytype),protected,save            :: lz
    ! grids in x,y,z direction
    integer,protected,save                 :: nx
    integer,protected,save                 :: ny
    integer,protected,save                 :: nz

    ! zstretch: this save control the slope in the boundary layer
    ! the higher this value the finer grids near the walls
    real(mytype),protected,save            :: zstretch

    ! channel: 0: closed channel, 1: open channel
    integer,protected,save                 :: ochannel

    ! dp_opt: options to drive the flow. dp_opt=3 is for Ekman layer
    ! dp_opt=1: constant friction. 
    ! dp_opt=2: constant bulk velocity. 
    ! dp_opt=3: constant geostrophic wind.
    integer,protected,save                 :: dp_opt

    ! issk=1: nonlinear term is casted in skew-symmetric form (only for cds=3)
    integer,protected,save                 :: issk

    ! cds: differential scheme for spatial derivatives 
    ! 1: 2nd order CDS
    ! 2: 4th order CDS
    ! 3: Spectral method
    integer,protected,save                 :: cds

    ! dts: time advancement scheme for the convective terms. 1: AB2; 2: RK3
    ! for the viscous terms, Crank-Nicolson is used
    integer,protected,save                 :: dts

    ! dt: time step; imax: how many steps to run
    real(mytype),protected,save            :: dt
    integer,protected,save                 :: imax

    ! nu: kinematic viscosity; alpha: thermal diffusivity; got: g/T
    ! fc: Coriolis parameter; ug, vg: geostrophic wind speeds
    real(mytype),protected,save            :: nu
    real(mytype),protected,save            :: alpha
    real(mytype),protected,save            :: got
    real(mytype),protected,save            :: fc
    real(mytype),protected,save            :: ug
    real(mytype),protected,save            :: vg

    ! isscalar: include temperature field? t_ref: a reference temperature
    integer,protected,save                 :: isscalar
    real(mytype),protected,save            :: t_ref

    ! tbot, ttop: temperature boundary conditions at the bottom and top walls
    real(mytype),protected,save            :: tbot
    real(mytype),protected,save            :: ttop

    ! is_ri_var: varying Richardson number? ri_str, ri_end: ri at the start
    ! and end of the simulation (only for is_ri_var=1)
    integer,protected,save                 :: is_ri_var
    real(mytype),protected,save            :: ri_str
    real(mytype),protected,save            :: ri_end

    ! isnoise: add noise to the initial fields? if yes, what is the mag? noise_mag
    integer,protected,save                 :: isnoise
    real(mytype),protected,save            :: noise_mag

    ! isdamp: Rayleigh damping? nzdamp: how many points? 
    ! cdamp: damping coefficient
    integer,protected,save                 :: isdamp
    integer,protected,save                 :: nzdamp
    real(mytype),protected,save            :: cdamp

    ! u_mrf: moving speed of the reference frame, set it to 0 for stationary frame
    real(mytype),protected,save            :: u_mrf

    ! ibackup/iinstfl/imeanfl:
    ! how many steps to save backup/instantaneous/mean fields
    integer,protected,save                 :: ibackup
    integer,protected,save                 :: iinstfl
    integer,protected,save                 :: imeanfl

    ! restart: read initial fields from files? istmsr: output time-series?
    integer,protected,save                 :: restart
    integer,protected,save                 :: istmsr

    ! isxy2d: whether to output 2d x-y slice data. 
    ! xy2d_id: where (k indices) to output
    integer,protected,save                 :: isxy2d
    integer,dimension(50),protected,save   :: xy2d_id
    ! isxz2d: whether to output 2d x-z slice data. 
    ! xz2d_id: where (j indices) to output
    integer,protected,save                 :: isxz2d
    integer,dimension(50),protected,save   :: xz2d_id
    ! isyz2d: whether to output 2d y-z slice data. 
    ! yz2d_id: where (i indices) to output
    integer,protected,save                 :: isyz2d
    integer,dimension(50),protected,save   :: yz2d_id
    ! intv_2d: output frequency (time step)
    integer,protected,save                 :: intv_2d

    ! these are namelist parameters read from "parameters.input"
    namelist /domain/ &
    p_row, p_col, lx, ly, lz, nx, ny, nz, zstretch,  ochannel

    namelist /modeling/ &
    dp_opt, cds, dts,isnoise, noise_mag, isdamp, nzdamp, cdamp,issk

    namelist /constants/ &
    nu, alpha, got, fc, ug, vg, dt, imax, u_mrf, isscalar, t_ref, &
           tbot, ttop, is_ri_var, ri_str, ri_end

    namelist /io/ &
    restart, istmsr, ibackup, iinstfl, imeanfl, &
    isxy2d, xy2d_id, isxz2d, xz2d_id, isyz2d, yz2d_id, intv_2d



    !**************************************************************************
    !************************   Derived Parameters   **************************
    !**************************************************************************

    ! number of processors, nprc=p_row*p_col
    integer,protected,save                 :: nprc

    ! iii: sqrt(-1)
    complex(mytype),parameter              :: iii=cmplx(0.d0,1.d0,kind=mytype)

    ! dx, dy: grid spacing in the x and y directions, dx2=dx*dx
    real(mytype),protected,save            :: dx, dy
    real(mytype),protected,save            :: dx2,dy2

    ! ghst: ghost cell number
    integer,parameter                      :: ghst=2

    ! inner domain start and end indices on z pencil
    integer,protected,save                 :: istr3,iend3,jstr3,jend3,kstr3,kend3

    ! nia, nja, nka: the dimensions of 3D array in i, j, and k indices
    ! 1, 2, and 3 means x, y, and z pencil
    integer,protected,save                 :: nia1, nja1, nka1
    integer,protected,save                 :: nia2, nja2, nka2
    integer,protected,save                 :: nia3, nja3, nka3

    ! i and j offset indices, used in mpi
    integer,protected,save                 :: i_offset,j_offset

    ! div_hist_max: historically max divergence errors
    ! time_sim: simulation time
    ! cfl_max: maximal CFL number in the domain
    ! these parameters are non-protected
    real(mytype),save                      :: div_hist_max=0.d0
    real(mytype),save                      :: time_sim=0.d0
    real(mytype),save                      :: cfl_max=0.d0

    ! dpx_drive: mean pressure difference to drive the flow
    ! if dp_opt=1, dpx_drive=1; if dp_opt=2 or 3, dpx_drive varies to 
    ! ensure constant flow rate or geostrophic wind. 
    ! ubar: mean velocity in the channel
    ! dpx_coef: coefficient used for constant mass flow rate case
    real(mytype),protected,save            :: dpx_drive, dpx_coef, ubar

    ! wx1=e^(i*2pi/nx), wy1=e^(i*2pi/ny), used in fft Poisson solver
    complex(mytype),protected,save         :: wx1
    complex(mytype),protected,save         :: wy1

    ! coefficients for RK3 scheme, see Orlandi 2000
    real(mytype),dimension(3),parameter    ::rkc1=(/8.d0/15.d0,5.d0/12.d0,3.d0/4.d0/)
    real(mytype),dimension(3),parameter    ::rkc2=(/0.d0,-17.d0/60.d0,-5.d0/12.d0/)
    real(mytype),dimension(3),parameter    ::rkc3=(/8.d0/15.d0,2.d0/15.d0,1.d0/3.d0/)
    real(mytype),protected,save            ::rk_a,rk_b,rk_c

    ! dz: cell spacing in z direction. these are derived parameters
    ! dz_t, dz_b: top and bottom spacing between cell vertex
    ! l_t, l_b : linear interpolation coefficient for top and bottom
    real(mytype),dimension(:),allocatable,protected,save :: &
                                           dz,l_t,l_b,dz_t,dz_b

    ! myid: the id of the current processor, ierr: error flag
    ! myid_row/colindx: i/j index of pij for the current processor
    integer,protected,save                 :: myid,myid_rowindx,myid_colindx
    integer,save                           :: ierr

    ! pij: processor id for given row and col indices
    integer,dimension(:,:),allocatable,protected,save    :: pij

    ! fft and ifft plans for fftw. 1 and 2 means x and y pencil
    type(C_PTR),protected,save             :: fft_plan1,ifft_plan1
    type(C_PTR),protected,save             :: fft_plan2,ifft_plan2




contains



    ! this function reads input parameters for the simulation,
    ! allocates and initiates parameters, calculates derived parameters, etc
    subroutine initiate_parameters

        implicit none

        ! first read parameters from the file "parameters.input"
        open(10,file='parameters.input')
        read(10,nml=domain)
        read(10,nml=modeling)
        read(10,nml=constants)
        read(10,nml=io)
        close(10)

        ! test if nx and ny are divisible by p_row, and if ny and nz are
        ! divisible by p_col
        if (mod(nx,p_row) .ne. 0 .or. mod(ny,p_row) .ne. 0) then
            print *,'nx and ny are not divisible by p_row! Quit!'
            stop
        endif

        if (mod(ny,p_col) .ne. 0 .or. mod(nz,p_col) .ne. 0) then
            print *,'ny and nz are not divisible by p_col! Quit!'
            stop
        endif

        ! nx and ny must be even numbers
        if (mod(nx,2) .ne. 0 .or. mod(ny,2) .ne. 0) then
            print *,'Please use even numbers for nx and ny! Quit!'
            stop
        endif

        ! calculate some derived parameters
        ! dx, dy: grid spacing in the x and y directions, dx2=dx*dx
        dx=lx/nx
        dy=ly/ny
        dx2=dx*dx
        dy2=dy*dy

        ! dpx_coef: coefficient used for constant mass flow rate case
        dpx_coef=nu*1.d0

        ! wx1=e^(i*2pi/nx), wy1=e^(i*2pi/ny), used in fft Poisson solver
        wx1=exp(iii*2.d0*pi/nx)
        wy1=exp(iii*2.d0*pi/ny)

        ! allocate and initiate some arrays
        ! allocate pij
        allocate( pij(p_row,p_col) )
        ! allocate grid parameters
        allocate( dz(nz+2*ghst) )
        allocate( l_t(nz+2*ghst) )
        allocate( l_b(nz+2*ghst) )
        allocate( dz_t(nz+2*ghst) )
        allocate( dz_b(nz+2*ghst) )

        ! initiate them
        ! grid parameters
        dz(:)=0.d0
        l_t(:)=0.d0
        l_b(:)=0.d0
        dz_t(:)=0.d0
        dz_b(:)=0.d0
        ! pij
        pij(:,:)=0

        ! initiate the mean pressure difference for driving the flow,
        ! if constant friction is used, dpx_drive=1
        ! if constant mass flow rate is used, dpx_drive varies during the simulation
        if (dp_opt .eq. 1) then
            dpx_drive = 1.d0
        else
            dpx_drive = 0.d0 ! give inital values for dpx_drive and ubar
            ubar      = 0.d0
        endif

        return

    end subroutine



    ! this function initiates mpi and tests if nprc=p_row*p_col
    subroutine initiate_mpi

        use mpi

        implicit none

        ! initialize mpi
        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, nprc,ierr )

        ! before we go, do some tests
        ! test if nprc=p_row*p_col
        if (nprc .ne. p_row*p_col) then
            if (myid .eq. 0) then
                print *,'nprc error'
                print *,'please set nprc=p_row*p_col!'
                stop
            endif
        endif

        return

    end subroutine



    ! this function finalizes mpi
    subroutine finalize_mpi

        use mpi

        implicit none

        call MPI_FINALIZE( ierr )

        return

    end subroutine



    ! this function initiates i, j, k, indices, e.g., istr3, iend3, i/j_offset, etc.
    subroutine initiate_indices

        implicit none

        integer                           :: i,j,indx

        ! initiate indices for z pencil
        istr3=1+ghst
        iend3=zsize(1)+ghst
        jstr3=1+ghst
        jend3=zsize(2)+ghst
        kstr3=1+ghst
        kend3=zsize(3)+ghst

        ! initiate array size for x pencil
        nia1=xsize(1)+2*ghst
        nja1=xsize(2)+2*ghst
        nka1=xsize(3)+2*ghst

        ! initiate array size for y pencil
        nia2=ysize(1)+2*ghst
        nja2=ysize(2)+2*ghst
        nka2=ysize(3)+2*ghst

        ! initiate array size for z pencil
        nia3=zsize(1)+2*ghst
        nja3=zsize(2)+2*ghst
        nka3=zsize(3)+2*ghst

        !calculate i and j indices offset based on myid
        i_offset=int(myid/p_col)*nx/p_row
        j_offset=mod(myid,p_col)*ny/p_col

        ! calculate the correspondence of processor id and row-column indices
        ! For example, if p_row=4, p_col=3, and myid=5, then
        ! myid_rowindx=2,myid_colindx=3, that being said, this processor
        ! is assigned to the decomposed domain at 2nd row and 3rd column.
        indx=0
        do i=1,p_row
        do j=1,p_col

            pij(i,j)=indx

            if ( myid .eq. pij(i,j) ) then

                myid_rowindx=i
                myid_colindx=j

            endif

            indx=indx+1

        enddo
        enddo

        return
    end subroutine



    ! this function calculates mesh parameters, e.g., dz,l_t,dz_t,etc.
    subroutine get_mesh_param

        implicit none

        integer                      :: k
        real(mytype)                 :: zmin,zmax,zmean,stretch_sysm,coeff_a
        real(mytype),dimension(nka3) :: zz

        ! calculate z
        zmin=0.d0
        zmax=lz

        ! is it an open-channel or closed-channel?
        if (ochannel .eq. 0) then

            ! closed channel
            zmean=0.5d0*(zmin+zmax)
            stretch_sysm=0.5d0
            coeff_a=(2.d0+nz)*stretch_sysm ! this is a parameter

        elseif (ochannel .eq. 1) then

            ! open channel
            zmean=(zmin+zmax)
            stretch_sysm=1.d0
            coeff_a=(1.d0+nz)*stretch_sysm ! this is a parameter

        endif

        ! calculate z
        zz(:)=0.d0
        do k=1,nz+1,1
            zz(k+ghst)=zmean+(zmin-zmean) / tanh(2.d0*zstretch) * &
                          tanh(2.d0*zstretch*(k-coeff_a)/(1.d0-coeff_a))
        enddo

        ! calculate dz
        dz(:)=0.d0
        do k=kstr3,kend3,1
            dz(k)=zz(k+1)-zz(k)
        enddo
        dz(kstr3-1)=dz(kstr3)
        dz(kend3+1)=dz(kend3)
        dz(kstr3-ghst)=dz(kstr3)
        dz(kend3+ghst)=dz(kend3)

        ! calculate interpolation coefficients l_t,l_b
        l_t(:)=0.d0
        l_b(:)=0.d0
        do k=kstr3-1,kend3+1,1    ! we can extend lt and lb
            l_t(k)=dz(k) / (dz(k)+dz(k+1))
            l_b(k)=dz(k) / (dz(k)+dz(k-1))
        enddo

        ! calculate dz_t and dz_b
        dz_t(:)=0.d0
        dz_b(:)=0.d0
        do k=kstr3-1,kend3+1,1 ! we can extend dz_t and dz_b
            dz_t(k)=0.5d0 * ( dz(k)+dz(k+1) )
            dz_b(k)=0.5d0 * ( dz(k)+dz(k-1) )
        enddo

        ! output z to files
        if (myid .eq. 0) then

            write(*,*) 'Writing Geo Info to Files.'
            open(11,file='results/geo.info')
            do k=kstr3,kend3,1
                write(11,'(4e14.6)') (zz(k+1)+zz(k))/2.d0,dz(k),l_t(k),l_b(k)
            enddo
            close(11)

        endif

        return

    end subroutine



    ! calculate bulk velocity
    subroutine calculate_ubar(u)

        use mpi

        implicit none

        real(mytype),dimension(:,:,:),intent(in) :: u
        real(mytype)                             :: ubar_local
        integer                                  :: i,j,k

        ! here we calculate ubar
        ubar_local=0.d0
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            ubar_local=ubar_local+u(i,j,k)*dz(k)

        enddo
        enddo
        enddo

        ubar_local=ubar_local/lz/zsize(1)/zsize(2)

        ! exchange ubar_local between processors
        ! sum up all the ubar_local between all processors
        ! and get ubar, this will be used in constant mass flow rate control
        call MPI_ALLREDUCE(ubar_local,ubar,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
        ubar=ubar/nprc ! we need to do averaging


        return

    end subroutine



    ! calculate fft plan for fftw
    subroutine get_fft_plan

        ! for using fftw
        use, intrinsic :: iso_c_binding

        implicit none

        ! for using fftw
        include 'include/fftw3.f03'

        complex(mytype),dimension(nx) :: c_tmp1
        complex(mytype),dimension(ny) :: c_tmp2

        ! create plan (fft) for x direction
        fft_plan1 =fftw_plan_dft_1d(nx,c_tmp1,c_tmp1,FFTW_FORWARD,FFTW_MEASURE)
        ! create plan (ifft) for x direction
        ifft_plan1=fftw_plan_dft_1d(nx,c_tmp1,c_tmp1,FFTW_BACKWARD,FFTW_MEASURE)

        ! create plan (fft) for y direction
        fft_plan2 =fftw_plan_dft_1d(ny,c_tmp2,c_tmp2,FFTW_FORWARD,FFTW_MEASURE)
        ! create plan (ifft) for y direction
        ifft_plan2=fftw_plan_dft_1d(ny,c_tmp2,c_tmp2,FFTW_BACKWARD,FFTW_MEASURE)

        return

    end subroutine


    ! this function updates ri if is_ri_var=1
    subroutine update_ri(stps)

        implicit none

        integer,intent(in)                 :: stps

        got=ri_str+(ri_end-ri_str)*1.d0/imax*stps

        return

    end subroutine


    ! this function updates dpx_drive
    subroutine update_dpx_drive

        implicit none

        ! check which option to drive the flow
        if     (dp_opt .eq. 1) then

            ! constant friction
            dpx_drive=1.d0

        elseif (dp_opt .eq. 2) then

            ! constant mass flow rate
            dpx_drive=dpx_drive - (ubar+u_mrf-1.d0)*dpx_coef

        endif

        return

    end subroutine



    ! assign coefficient for the RK3 or AB2 scheme
    subroutine assign_rk_coeff(a_in,b_in,c_in)

        implicit none

        real(mytype), intent(in)  :: a_in,b_in,c_in

        if (dts .eq. 2) then

            rk_a=a_in
            rk_b=b_in
            rk_c=c_in

        else

            ! AB2 scheme
            rk_a= 1.5d0
            rk_b=-0.5d0
            rk_c= 1.d0

        endif

        return

    end subroutine



end module



