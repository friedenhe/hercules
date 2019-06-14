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



module module_navier_stokes

    use module_parameters
    implicit none

contains


    ! the main function for HERCULES, it solves the Navier-Stokes for incompressible flows
    subroutine drive

        use module_variables
        use module_io
        use module_tools
        use module_boundary
        use module_poisson_solver
        use mpi
        use module_spectral

        implicit none

        ! for using fftw
        include 'include/fftw3.f03'

        ! indx_bk/m/ins are indices for output 
        ! backup/mean/instantaneous variables
        ! indx_slc: indices for output slice data
        integer   :: indx_bk,indx_m,indx_ins,indx_slc
        integer   :: stps

        ! calculate mesh parameters, e.g., dz
        call get_mesh_param

        ! allocate arrays for all variables
        ! initiate them using default values
        ! or read them from files (e.g., init_u.dat)
        call initiate_fields

        ! assign initial ri if the varying ri is used, i.e., is_ri_var=1
        if (is_ri_var .eq. 1) then
            call update_ri(0)
        endif

        ! calculate bulk velocity if dp_opt=2
        if (dp_opt .eq. 2) then
            call calculate_ubar(u)
        endif

        ! create fftw plans, this needs to be done only once
        call get_fft_plan

        ! get the coefficient a,b,c for the u, v, w, t, and p equations
        ! we store them in aa_u,bb_u, etc., 
        ! and then assign them to aa,bb,cc using subroutine assign_abc
        call get_u_eqn_coeff(aa_u,bb_u,cc_u)
        call get_w_eqn_coeff(aa_w,bb_w,cc_w)
        call get_p_eqn_coeff(aa_p,bb_p,cc_p)
        call get_t_eqn_coeff(aa_t,bb_t,cc_t)

        ! open new files for outputing time-series
        if (myid .eq. 0 .and. istmsr .eq. 1) then
            open(99,file='results/tmsr_u.dat',status='replace')
            close(99)
            open(98,file='results/tmsr_v.dat',status='replace')
            close(98)
            open(97,file='results/tmsr_w.dat',status='replace')
            close(97)
            open(96,file='results/tmsr_t.dat',status='replace')
            close(96)
            open(95,file='results/tmsr_p.dat',status='replace')
            close(95)
        endif
        
        ! for spectral method, we need to calculate the Fourier coefficients
        ! and do truncation, and update fields in physical space 
        ! note that for spectral method, only the 2/3 truncation dealiasing
        ! is avaiable right now        
        if (cds .eq. 3) then
            
            ! convert variables from physical space to phase space
            call real_2_cmplx(u,c_u)
            call real_2_cmplx(v,c_v)
            call real_2_cmplx(w,c_w)
            call real_2_cmplx(p,c_p)
            call real_2_cmplx(t,c_t)
            call real_2_cmplx(ff,c_ff)
            call real_2_cmplx(gg,c_gg)
            call real_2_cmplx(hh,c_hh)
            call real_2_cmplx(ss,c_ss)
            !!!!!!!!!!! do spectral truncation !!!!!!!!!!!!
            call spectral_truncation(c_u)
            call spectral_truncation(c_v)
            call spectral_truncation(c_w)
            call spectral_truncation(c_p)
            call spectral_truncation(c_t)
            !!!!!!!!!! get updated variables after truncation !!!!!!!!
            call cmplx_2_real(c_u,u)
            call cmplx_2_real(c_v,v)
            call cmplx_2_real(c_w,w)
            call cmplx_2_real(c_p,p)
            call cmplx_2_real(c_t,t)

        endif

        ! update boundary values for u,v,w, t, and p
        call update_boundary_uvw(u,v,w)
        call update_boundary_p(p)
        call update_boundary_t(t)
        

        ! indices to save backup/mean/inst fields, or slice data
        indx_bk=1 
        indx_m=1
        indx_ins=1
        indx_slc=1

        ! main loop, we use fractional time step approach for time advancement
        do stps=1,imax,1

            ! print cpu time and CFL number
            call screen_cpu_time(u,v,w)

            ! calculate new ri if is_ri_var=1
            if (is_ri_var .eq. 1) then
                call update_ri(stps)
            endif

            ! do the time advancement
            if (cds .eq. 3) then
                
                ! subroutine for spectral method
                call time_advancement_sp
                
            else
                
                ! subroutine for finite-difference method
                call time_advancement_fd
                
            endif

            ! update bulk velocity if dp_opt=2
            if (dp_opt .eq. 2) then
                call calculate_ubar(u)
            endif

            !!!!!!!!!!!!!!!!! output data !!!!!!!!!!!!!!!!
            ! output time-series
            if (istmsr .eq. 1) then
                ! output time-series at the channel center
                if (myid_rowindx .eq. p_row/2+1 .and. &
                    myid_colindx .eq. p_col/2+1 ) then
                    call output_time_series(u,v,w,p,t)
                endif
            endif

            ! write backup files
            if (indx_bk/ibackup .eq. 1 ) then
                call output_backup(u,v,w,p,t,ff,gg,hh,ss)
                indx_bk=0
            endif

            ! write mean fields
            if (indx_m/imeanfl .eq. 1) then
                call output_mean_fields(mean_u,mean_v,mean_w,mean_t,mean_uu, &
                                        mean_vv,mean_ww,mean_uw,mean_vw,mean_tt, &
                                        mean_tw,stps)
                indx_m=0
            endif

            ! write inst fields
            if (indx_ins/iinstfl .eq. 1) then
                call output_inst_fields(u,v,w,p,t,stps)
                indx_ins=0
            endif

            ! output 2d slices
            if (isxy2d .eq. 1 .or. isxz2d .eq. 1 .or. isyz2d .eq. 1) then
                if (indx_slc/intv_2d .eq. 1) then
                    call output_2d_slices(u,v,w,p,t,stps)
                    indx_slc=0
                endif
            endif

            indx_bk=indx_bk+1
            indx_m=indx_m+1
            indx_ins=indx_ins+1
            indx_slc=indx_slc+1

        enddo

        ! destroy fftw plans
        call fftw_destroy_plan(fft_plan1)
        call fftw_destroy_plan(ifft_plan1)
        call fftw_destroy_plan(fft_plan2)
        call fftw_destroy_plan(ifft_plan2)

    end subroutine



    ! this function does the time-advancement using finite-difference method
    subroutine time_advancement_fd

        use module_variables
        use module_tools
        use module_boundary
        use module_cnvdiff
        use module_poisson_solver
        use module_io
        use mpi

        implicit none

        integer   :: t_its,t_its_max

        ! if it is RK3 then we need to do 3 intermediate steps
        ! for fractional time step. While for AB2, 1 step is enough
        if (dts .eq. 2) then
            t_its_max=3
        else
            t_its_max=1
        endif

        ! for RK3 we need to do 3 intermediate sub-steps for each time step
        ! for AB2, t_its_max=1
        do t_its=1,t_its_max,1
        
            ! if it is RK3 scheme, write some info to screen
            if (dts .eq. 2 .and. myid .eq. 0) then

                write (*,218) t_its
                218 format(' RK3, Step: ',i2)

            endif
            
            ! check if we need to do rayleigh damping at the top
            if (isdamp .eq. 1) then

                call rayleigh_damping(u)
                call rayleigh_damping(v)
                call rayleigh_damping(w)
                call rayleigh_damping(p)
                call rayleigh_damping(t)

            endif

            ! before calculating ff, gg, hh, assign them to the previous time step
            ! ff,gg,hh are used in get_momentum_cnvdiff
            call assign_fgh_old(ff,gg,hh,ff0,gg0,hh0)

            ! calculate interpolated fields,
            ! which will be used in cnvdiff calculation
            call get_interp_fields_uvwhxyz(u,v,w,uhx,uhy,uhz,vhx,vhy,vhz,whx,why)
            call update_boundary_uvwhxyz(u,v,uhx,uhy,uhz,vhx,vhy,vhz,whx,why)

            ! assign coefficients for RK3 and AB2 scheme
            ! these coefficients, i.e., rk_a,rk_b,rk_c, will be used
            ! in get_momentum_cnvdiff and assign_abc
            call assign_rk_coeff(rkc1(t_its),rkc2(t_its),rkc3(t_its))

            ! calculate the spatial derivatives 
            ! and get the right hand side for u, v, and w eqns
            call get_momentum_cnvdiff &
           (u,v,w,p,t,uhx,uhy,uhz,vhx,vhy,vhz,whx,why,ff,ff0,gg,gg0,hh,hh0,qu,qv,qw)

            ! fractional time step
            ! first step, calculating u* using p^n

            ! assign the previously calculated aa_u,bb_u,cc_u, to aa,bb,cc
            call assign_abc(aa_u,bb_u,cc_u,aa,bb,cc)
            ! solve the u equation using fft solver
            call poisson_solver_fft(aa,bb,cc,qu,u,0,0,0,1)

            ! assign the previously calculated aa_v,bb_v,cc_v, to aa,bb,cc
            ! coeffs a,b,c for v eqn is same with u eqn,
            call assign_abc(aa_u,bb_u,cc_u,aa,bb,cc)
            ! solve the v equation using fft solver
            call poisson_solver_fft(aa,bb,cc,qv,v,0,0,0,0)

            ! assign the previously calculated aa_w,bb_w,cc_w, to aa,bb,cc
            call assign_abc(aa_w,bb_w,cc_w,aa,bb,cc)
            ! solve the w equation using fft solver
            call poisson_solver_fft(aa,bb,cc,qw,w,1,0,0,0)
            ! update the boundary values
            call update_boundary_uvw(u,v,w)

            ! based on u*, calculate the pseudo pressure phi
            ! here we reuse qu for the source term of p eqn to save memory usage
            call get_p_eqn_src(u,v,w,qu)
            ! assign the previously calculated aa_p,bb_p,cc_p, to aa,bb,cc
            call assign_abc_p(aa_p,bb_p,cc_p,aa,bb,cc)
            ! solve the p equation using fft solver
            call poisson_solver_fft(aa,bb,cc,qu,phi,0,1,0,0)
            ! update boundary values
            call update_boundary_p(phi)

            ! calculate u^(n+1) using phi
            call calculate_uvw(u,v,w,phi)
            call update_boundary_uvw(u,v,w)

            ! calculate temperature
            if (isscalar .eq. 1) then

                ! assign ss to previous steps
                call assign_s_old(ss,ss0)
                
                ! calculate interpolated fields
                call get_interp_fields_thxyz(t,thx,thy,thz)
                call update_boundary_thxyz(thx,thy,thz)

                ! here we reuse qu for the source term of t eqn to save memory
                call get_temperature_cnvdiff(u,v,w,t,thx,thy,thz,ss,ss0,qu)
                ! assign the previously calculated aa_t,bb_t,cc_t, to aa,bb,cc
                call assign_abc(aa_t,bb_t,cc_t,aa,bb,cc)
                ! solve the t equation using fft solver
                call poisson_solver_fft(aa,bb,cc,qu,t,0,0,1,0)
                ! update boundary values
                call update_boundary_t(t)

            endif

            ! calculate the final DIV error
            if (dts .ne. 2) then
                call screen_div_error(u,v,w,1,1)
            elseif (t_its .ne. 3) then
                call screen_div_error(u,v,w,1,0)
            else
                call screen_div_error(u,v,w,1,1)
            endif

            ! calculate p^{n+1} using phi
            call calculate_new_p(phi,p)
            call update_boundary_p(p)

        enddo

        ! calculate mean fields
        call calculate_mean_fields &
            (u,v,w,t,uhx,vhy,mean_u,mean_v,mean_w,mean_t, mean_uu,mean_vv, &
             mean_ww,mean_uw,mean_vw,mean_tt,mean_tw)

        return

    end subroutine



    ! this function does the time-advancement 
    ! using spectral method in horizontal directions
    subroutine time_advancement_sp

        use module_spectral
        use module_variables
        use mpi
        use module_tools
        use module_boundary
        use module_io
        use module_poisson_solver

        implicit none

        integer   :: t_its,t_its_max

        ! if it is RK3 then we need to do 3 intermediate steps
        ! for fractional time step. While for AB2, 1 step is enough
        if (dts .eq. 2) then
            t_its_max=3
        else
            t_its_max=1
        endif

        ! for RK3 we need to do 3 intermediate sub-steps for each time step
        ! for AB2, t_its_max=1
        do t_its=1,t_its_max,1
        
            ! if it is RK3 scheme, write some info to screen
            if (dts .eq. 2 .and. myid .eq. 0) then

                write (*,218) t_its
                218 format(' RK3, Step: ',i2)

            endif
            
            ! check if we need to do rayleigh damping
            if (isdamp .eq. 1) then

                call rayleigh_damping_spec(c_u)
                call rayleigh_damping_spec(c_v)
                call rayleigh_damping_spec(c_w)
                call rayleigh_damping_spec(c_p)
                call rayleigh_damping_spec(c_t)

            endif

            ! assign the coefficient for RK3 or AB2
            call assign_rk_coeff(rkc1(t_its),rkc2(t_its),rkc3(t_its))
            
            ! calculate the right-hand side terms
            call get_momentum_cnvdiff_spec &
                (u,v,w,p,t,c_u,c_v,c_w,c_p,c_ff,c_gg,c_hh,c_qu,c_qv,c_qw)
                                         
            ! fractional time step
            ! first step, calculating u* using p^n

            ! assign the previously calculated aa_u,bb_u,cc_u, to aa,bb,cc
            call assign_abc(aa_u,bb_u,cc_u,aa,bb,cc)
            ! solve the u equation using fft solver
            call poisson_solver_fft_spec(aa,bb,cc,c_qu,c_u,0,0,0,1)

            ! assign the previously calculated aa_v,bb_v,cc_v, to aa,bb,cc
            ! coeffs a,b,c for v eqn is same with u eqn,
            call assign_abc(aa_u,bb_u,cc_u,aa,bb,cc)
            ! solve the v equation using fft solver
            call poisson_solver_fft_spec(aa,bb,cc,c_qv,c_v,0,0,0,0)

            ! assign the previously calculated aa_w,bb_w,cc_w, to aa,bb,cc
            call assign_abc(aa_w,bb_w,cc_w,aa,bb,cc)
            ! solve the w equation using fft solver
            call poisson_solver_fft_spec(aa,bb,cc,c_qw,c_w,1,0,0,0)
            
            ! we only need to update w in physical space, since we need to
            ! calculate dwdz for the source term of p equation
            call cmplx_2_real(c_w,w)
            call update_boundary_uvw(u,v,w)
            
            ! based on u*, calculate the pseudo pressure phi
            ! here we reused qu for the source term of p eqn to save memory usage
            call get_div_spec(c_u,c_v,w,c_qu)
            c_qu(:,:,:)=c_qu(:,:,:)/dt
            ! assign the previously calculated aa_p,bb_p,cc_p, to aa,bb,cc
            call assign_abc_p(aa_p,bb_p,cc_p,aa,bb,cc)
            ! solve the p equation using fft solver
            call poisson_solver_fft_spec(aa,bb,cc,c_qu,c_phi,0,1,0,0)

            ! we need to update phi in physical space for calculating new p
            call cmplx_2_real(c_phi,phi)
            call update_boundary_p(phi)
            
            ! calculate p^{n+1} using phi
            call calculate_new_p_spec(phi,c_phi,c_p)
            
            ! calculate u^(n+1) using phi
            call calculate_uvw_spec(phi,c_phi,c_u,c_v,c_w)

            !!!!!!!!!!! do spectral truncation !!!!!!!!!!!!
            call spectral_truncation(c_u)
            call spectral_truncation(c_v)
            call spectral_truncation(c_w)
            call spectral_truncation(c_p)
            
            ! update u,v,w,p in physical space
            call cmplx_2_real(c_u,u)
            call cmplx_2_real(c_v,v)
            call cmplx_2_real(c_w,w)
            call cmplx_2_real(c_p,p)
            call update_boundary_p(p)
            call update_boundary_uvw(u,v,w)

            ! calculate temperature
            if (isscalar .eq. 1) then

                ! here we reuse qu for the source term of t eqn to save memory
                call get_temperature_cnvdiff_spec &
                     (u,v,w,t,c_t,c_ss,c_qu)
                ! assign the previously calculated aa_t,bb_t,cc_t, to aa,bb,cc
                call assign_abc(aa_t,bb_t,cc_t,aa,bb,cc)
                ! solve the t equation using fft solver
                call poisson_solver_fft_spec(aa,bb,cc,c_qu,c_t,0,0,1,0)
                
                ! update t in physical space
                call spectral_truncation(c_t)
                call cmplx_2_real(c_t,t)
                call update_boundary_t(t)

            endif
            
            ! calculate the final DIV error
            if (dts .ne. 2) then
                call screen_div_error_spec(u,v,w,c_u,c_v,1,1)
            elseif (t_its .ne. 3) then
                call screen_div_error_spec(u,v,w,c_u,c_v,1,0)
            else
                call screen_div_error_spec(u,v,w,c_u,c_v,1,1)
            endif

            
        enddo

        ! calculate mean fields
        call calculate_mean_fields_spec &
                (u,v,w,t,mean_u,mean_v,mean_w,mean_t,mean_uu,mean_vv,mean_ww, &
                 mean_uw,mean_vw,mean_tt,mean_tw)

        return

    end subroutine



    ! this function assigns the previously
    ! calculated aa_u, bb_u, cc_u, etc. to aa,bb,cc
    subroutine assign_abc(aa_in,bb_in,cc_in,aa,bb,cc)

        use decomp_2d

        implicit none

        real(mytype),   dimension(:),    intent(in)  :: aa_in,cc_in
        complex(mytype),dimension(:,:,:),intent(in)  :: bb_in
        real(mytype),   dimension(:,:,:),intent(out) :: aa,cc
        complex(mytype),dimension(:,:,:),intent(out) :: bb

        integer              :: i,j,k

        ! note that in get_u_eqn_coeff, we did not add 1.0
        ! to bb, here we add it. Also, if we use RK3, we need to
        ! multiply rk_c to a,b,c, the rk_c for the current RK sub-step
        ! has been assign previously by calling assign_rk_coeff
        if (dts .eq. 2) then

            do k=1,zsize(3)+2*ghst,1
            do j=1,zsize(2),1
            do i=1,zsize(1),1
                
                aa(i,j,k)=rk_c*aa_in(k)
                cc(i,j,k)=rk_c*cc_in(k) 
                bb(i,j,k)=1.d0+rk_c*bb_in(i,j,k)

            enddo
            enddo
            enddo

        else

            do k=1,zsize(3)+2*ghst,1
            do j=1,zsize(2),1
            do i=1,zsize(1),1

                aa(i,j,k)=aa_in(k)
                cc(i,j,k)=cc_in(k)
                bb(i,j,k)=1.d0+bb_in(i,j,k)

            enddo
            enddo
            enddo

        endif

        return

    end subroutine



    ! this function assigns the previously calculated
    ! aa_p, bb_p, cc_p, etc. to aa,bb,cc for p equation
    subroutine assign_abc_p(aa_in,bb_in,cc_in,aa,bb,cc)

        use decomp_2d

        implicit none

        real(mytype),   dimension(:),    intent(in)  :: aa_in,cc_in
        complex(mytype),dimension(:,:,:),intent(in)  :: bb_in
        real(mytype),   dimension(:,:,:),intent(out) :: aa,cc
        complex(mytype),dimension(:,:,:),intent(out) :: bb

        integer              :: i,j,k

        if (dts .eq. 2) then

            do k=1,zsize(3)+2*ghst,1
            do j=1,zsize(2),1
            do i=1,zsize(1),1
                aa(i,j,k)=rk_c*aa_in(k)
                cc(i,j,k)=rk_c*cc_in(k)
                bb(i,j,k)=rk_c*bb_in(i,j,k)
            enddo
            enddo
            enddo

        else

            do k=1,zsize(3)+2*ghst,1
            do j=1,zsize(2),1
            do i=1,zsize(1),1
                aa(i,j,k)=aa_in(k)
                cc(i,j,k)=cc_in(k)
                bb(i,j,k)=bb_in(i,j,k)
            enddo
            enddo
            enddo

        endif

        return

    end subroutine



    ! this function assigns ff to the previous time steps, e.g., ff0
    subroutine assign_fgh_old(ff,gg,hh,ff0,gg0,hh0)

        implicit none

        real(mytype),dimension(:,:,:),intent(in)    :: ff,gg,hh
        real(mytype),dimension(:,:,:),intent(inout) :: ff0,gg0,hh0

        ! assign ff etc to ff0 etc.
        ff0(:,:,:)=ff(:,:,:)
        gg0(:,:,:)=gg(:,:,:)
        hh0(:,:,:)=hh(:,:,:)

        return

    end subroutine



    ! this function assigns ss to the previous time steps, e.g., ss0
    subroutine assign_s_old(ss,ss0)

        implicit none

        real(mytype),dimension(:,:,:),intent(in)  :: ss
        real(mytype),dimension(:,:,:),intent(out) :: ss0

        ss0(:,:,:)=ss(:,:,:)

        return

    end subroutine



    ! calculate u^(n+1) using phi
    subroutine calculate_uvw(u,v,w,phi)
        
        implicit none

        real(mytype),dimension(:,:,:),intent(in)     :: phi
        real(mytype),dimension(:,:,:),intent(inout)  :: u,v,w

        integer                               :: i,j,k
        
        ! check to use 2nd or 4th schemes
        if ( cds .eq. 1) then

            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1

                u(i,j,k)= u(i,j,k) - rk_c*dt* &
                           ( phi(i+1,j,k) - phi(i,j,k) ) / dx

            enddo
            enddo
            enddo

            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1

                v(i,j,k)= v(i,j,k) - rk_c*dt* &
                           ( phi(i,j+1,k) - phi(i,j,k) ) / dy

            enddo
            enddo
            enddo

        elseif ( cds .eq. 2) then

            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1

                u(i,j,k)= u(i,j,k) - rk_c*dt*( -phi(i+2,j,k) + &
                                          27.d0*phi(i+1,j,k) - &
                                          27.d0*phi(i,j,k)   + &
                                                phi(i-1,j,k) ) / 24.d0 / dx
            enddo
            enddo
            enddo

            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1

                v(i,j,k)= v(i,j,k) - rk_c*dt*( -phi(i,j+2,k) + &
                                          27.d0*phi(i,j+1,k) - &
                                          27.d0*phi(i,j,k)   + &
                                                phi(i,j-1,k) ) / 24.d0 / dy
            enddo
            enddo
            enddo

        endif

        do k=kstr3,kend3-1,1 ! note this range is different
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            w(i,j,k)= w(i,j,k) -  rk_c*dt* &
                         ( phi(i,j,k+1) - phi(i,j,k) ) / dz_t(k)

        enddo
        enddo
        enddo

        return

    end subroutine



    ! calculate p^(n+1) using p^n and the pseudo pressure phi
    subroutine calculate_new_p(phi,p)

        implicit none

        real(mytype),dimension(:,:,:),intent(in)  :: phi
        real(mytype),dimension(:,:,:),intent(out)  :: p

        integer              :: i,j,k
        real(mytype)         :: phixx,phiyy,phizz

        ! check to use 2nd or 4th schemes
        if ( cds .eq. 1) then

            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1

                phixx = ( phi(i+1,j,k) - &
                     2.d0*phi(i,j,k)   + &
                          phi(i-1,j,k) ) /dx2

                phiyy = ( phi(i,j+1,k) - &
                     2.d0*phi(i,j,k)   + &
                          phi(i,j-1,k) ) /dy2

                phizz= ( (phi(i,j,k+1)-phi(i,j,k))/dz_t(k) -  &
                         (phi(i,j,k)-phi(i,j,k-1))/dz_b(k) ) / dz(k)

                p(i,j,k)=p(i,j,k)+phi(i,j,k) - &
                     rk_c*0.5d0*dt*nu*(phixx+phiyy+phizz)

            enddo
            enddo
            enddo

        elseif ( cds .eq. 2) then

            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1

                phixx= ( -phi(i+2,j,k) + &
                    16.d0*phi(i+1,j,k) - &
                    30.d0*phi(i,j,k)   + &
                    16.d0*phi(i-1,j,k) - &
                          phi(i-2,j,k) ) / 12.d0 / dx2

                phiyy= ( -phi(i,j+2,k) + &
                    16.d0*phi(i,j+1,k) - &
                    30.d0*phi(i,j,k)   + &
                    16.d0*phi(i,j-1,k) - &
                          phi(i,j-2,k) ) / 12.d0 / dy2

                phizz= ( (phi(i,j,k+1)-phi(i,j,k))/dz_t(k) -  &
                         (phi(i,j,k)-phi(i,j,k-1))/dz_b(k) ) / dz(k)

                p(i,j,k)=p(i,j,k)+phi(i,j,k) - &
                     rk_c*0.5d0*dt*nu*(phixx+phiyy+phizz)

            enddo
            enddo
            enddo

        endif



        return

    end subroutine



end module
