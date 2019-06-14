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



module module_poisson_solver

    use module_parameters
    implicit none

contains



    ! this function calculates a,b,c for p equation
    subroutine get_p_eqn_coeff(aa,bb,cc)

        use mpi

        implicit none

        real(mytype),   dimension(:),    intent(out)  :: aa,cc
        complex(mytype),dimension(:,:,:),intent(out)  :: bb

        integer              :: i,j,k,l,m
        complex(mytype)      :: kx,ky

        aa(:)=0.d0
        bb(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        cc(:)=0.d0

        ! calculate coefficients a,b,c for the pressure equation
        do k=kstr3,kend3,1

            aa(k)=1.d0/dz(k)/dz_b(k)
            cc(k)=1.d0/dz(k)/dz_t(k)

            do j=1,zsize(2),1
            do i=1,zsize(1),1

                ! deal with index
                ! convert i to l
                if ( (i+i_offset) .le. nx/2) then
                    l=i+i_offset-1
                else
                    l=i+i_offset-nx-1
                endif
                ! convert j to m
                if ( (j+j_offset) .le. ny/2) then
                    m=j+j_offset-1
                else
                    m=j+j_offset-ny-1
                endif

                ! horizontal modified wavenumber,
                ! check to use 2nd or 4th schemes
                if ( cds .eq. 1) then

                    kx = ( wx1**l + wx1**(-l) -2.d0 ) /dx2

                    ky = ( wy1**m + wy1**(-m) -2.d0 ) /dy2

                elseif ( cds .eq. 2) then

                    kx =  (  wx1**(3.d0*l)  - &
                       54.d0*wx1**(2.d0*l)  + &
                      783.d0*wx1**l         - &
                     1460.d0                + &
                      783.d0*wx1**(-l)      - &
                       54.d0*wx1**(-2.d0*l) + &
                             wx1**(-3.d0*l) ) /576.d0/dx2

                    ky =  (  wy1**(3.d0*m)  - &
                       54.d0*wy1**(2.d0*m)  + &
                      783.d0*wy1**m         - &
                     1460.d0                + &
                      783.d0*wy1**(-m)      - &
                       54.d0*wy1**(-2.d0*m) + &
                             wy1**(-3.d0*m) ) /576.d0/dy2

                elseif ( cds .eq. 3) then

                    kx =  -(2.d0*pi/lx*l)**2.d0
                    ky =  -(2.d0*pi/ly*m)**2.d0

                endif

                bb(i,j,k)= kx + ky - aa(k) - cc(k)

                ! Neumann boundary condition
                if (k .eq. kstr3) then
                    bb(i,j,k)=bb(i,j,k)+aa(k)
                    aa(k)=0.d0
                endif
                if (k .eq. kend3) then
                    bb(i,j,k)=bb(i,j,k)+cc(k)
                    cc(k)=0.d0
                endif

            enddo
            enddo
        enddo

        return

    end subroutine



    ! this function calculates the right hand side for p equation
    subroutine get_p_eqn_src(u,v,w,qp)

        use mpi
        use module_tools

        implicit none

        real(mytype),dimension(:,:,:),intent(in)    :: u,v,w
        real(mytype),dimension(:,:,:),intent(out)   :: qp

        integer              :: i,j,k
        real(mytype)         :: dudx,dvdy,dwdz

        ! calculate qp
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            ! check to use 2nd or 4th schemes
            if ( cds .eq. 1) then

                dudx = ( u(i,j,k)-u(i-1,j,k) )/dx
                dvdy = ( v(i,j,k)-v(i,j-1,k) )/dy

            elseif ( cds .eq. 2) then

                dudx = (  -u(i+1,j,k) + &
                     27.d0*u(i,j,k)   - &
                     27.d0*u(i-1,j,k) + &
                           u(i-2,j,k) ) / 24.d0 / dx

                dvdy = (  -v(i,j+1,k) + &
                     27.d0*v(i,j,k)   - &
                     27.d0*v(i,j-1,k) + &
                           v(i,j-2,k) ) / 24.d0 / dy

            endif

            ! z direction is 2nd scheme
            dwdz = ( w(i,j,k)-w(i,j,k-1) )/dz(k)

            qp(i,j,k)=(dudx+dvdy+dwdz)/dt

        enddo
        enddo
        enddo

        return

    end subroutine



    ! this function calculates a,b,c for u equation
    subroutine get_u_eqn_coeff(aa,bb,cc)

        use mpi

        implicit none

        real(mytype),   dimension(:),    intent(out) :: aa,cc
        complex(mytype),dimension(:,:,:),intent(out) :: bb

        integer              :: i,j,k,l,m
        complex(mytype)      :: kx,ky

        aa(:)=0.d0
        bb(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        cc(:)=0.d0

        ! calculate coefficient for the u equation
        do k=kstr3,kend3,1

            aa(k)=-0.5d0*dt*nu/dz(k)/dz_b(k)
            cc(k)=-0.5d0*dt*nu/dz(k)/dz_t(k)

            do j=1,zsize(2),1
            do i=1,zsize(1),1

                ! deal with index
                ! convert i to l
                if ( (i+i_offset) .le. nx/2) then
                    l=i+i_offset-1
                else
                    l=i+i_offset-nx-1
                endif
                ! convert j to m
                if ( (j+j_offset) .le. ny/2) then
                    m=j+j_offset-1
                else
                    m=j+j_offset-ny-1
                endif

                ! horizontal modified wavenumber,
                ! check to use 2nd or 4th schemes
                if ( cds .eq. 1) then

                    kx = 0.5d0 * dt * ( wx1**l + wx1**(-l) - 2.d0 ) *nu/dx2

                    ky = 0.5d0 * dt * ( wy1**m + wy1**(-m) - 2.d0 ) *nu/dy2

                elseif ( cds .eq. 2) then

                    kx = dt * ( -wx1**(2.d0*l) + &
                           16.d0*wx1**l        - &
                           30.d0               + &
                           16.d0*wx1**(-l)     - &
                                 wx1**(-2.d0*l)  ) /24.d0*nu/dx2

                    ky = dt * ( -wy1**(2.d0*m) + &
                           16.d0*wy1**m        - &
                           30.d0               + &
                           16.d0*wy1**(-m)     - &
                                 wy1**(-2.d0*m)  ) /24.d0*nu/dy2

                elseif ( cds .eq. 3) then

                    kx = 0.5d0*dt*nu*( -(2.d0*pi/lx*l)**2.d0 ) 
                    ky = 0.5d0*dt*nu*( -(2.d0*pi/ly*m)**2.d0 ) 

                endif

                ! note that here we should do 1-kx-ky-aa-cc for bb
                ! however, to be convenient to implement RK3, we
                ! add the missing 1 in the subroutine "assign_abc"
                bb(i,j,k)= - kx - ky - aa(k) - cc(k)

                ! boundary condition
                ! we need to add the missing terms to qu 
                ! in Poisson equation solver
                if (k .eq. kstr3) then
                    bb(i,j,k)=bb(i,j,k)-aa(k) ! non-slip
                    aa(k)=0.d0
                endif

                ! we need to add the missing terms to qu 
                ! in Poisson equation solver
                if (k .eq. kend3) then

                    if (ochannel .eq. 0) then
                        bb(i,j,k)=bb(i,j,k)-cc(k)  ! non-slip
                        cc(k)=0.d0
                    elseif (ochannel .eq. 1) then
                        bb(i,j,k)=bb(i,j,k)+cc(k)  ! slip
                        cc(k)=0.d0
                    endif

                endif

            enddo
            enddo
        enddo

        return

    end subroutine



    ! this function calculate a,b,c for w equation
    subroutine get_w_eqn_coeff(aa,bb,cc)

        use mpi

        implicit none

        real(mytype),   dimension(:),    intent(out) :: aa,cc
        complex(mytype),dimension(:,:,:),intent(out) :: bb

        integer              :: i,j,k,l,m
        complex(mytype)      :: kx,ky

        aa(:)=0.d0
        bb(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        cc(:)=0.d0

       ! calculate coefficient for the u equation
        do k=kstr3,kend3-1,1  ! note this index is different

            aa(k)=-0.5d0*dt*nu/dz(k)/dz_t(k)
            cc(k)=-0.5d0*dt*nu/dz(k+1)/dz_t(k)

            do j=1,zsize(2),1
            do i=1,zsize(1),1

                ! deal with index
                ! convert i to l
                if ( (i+i_offset) .le. nx/2) then
                    l=i+i_offset-1
                else
                    l=i+i_offset-nx-1
                endif
                ! convert j to m
                if ( (j+j_offset) .le. ny/2) then
                    m=j+j_offset-1
                else
                    m=j+j_offset-ny-1
                endif

                ! horizontal modified wavenumber,
                ! check to use 2nd or 4th schemes
                if ( cds .eq. 1) then

                    kx = 0.5d0 * dt * ( wx1**l + wx1**(-l) - 2.d0 ) *nu/dx2

                    ky = 0.5d0 * dt * ( wy1**m + wy1**(-m) - 2.d0 ) *nu/dy2

                elseif ( cds .eq. 2) then

                    kx = dt * ( -wx1**(2.d0*l) + &
                           16.d0*wx1**l        - &
                           30.d0               + &
                           16.d0*wx1**(-l)     - &
                                 wx1**(-2.d0*l)  ) /24.d0*nu/dx2

                    ky = dt * ( -wy1**(2.d0*m) + &
                           16.d0*wy1**m        - &
                           30.d0               + &
                           16.d0*wy1**(-m)     - &
                                 wy1**(-2.d0*m)  ) /24.d0*nu/dy2

                elseif ( cds .eq. 3) then

                    kx = 0.5d0*dt*nu*( -(2.d0*pi/lx*l)**2.d0 ) 
                    ky = 0.5d0*dt*nu*( -(2.d0*pi/ly*m)**2.d0 ) 

                endif

                ! note that here we should do 1-kx-ky-aa-cc for bb
                ! however, to be convenient to implement RK3, we
                ! add the missing 1 in the subroutine "assign_abc"
                bb(i,j,k)= - kx - ky - aa(k) - cc(k)

                ! Non-slip condition
                if (k .eq. kstr3) then
                    aa(k)=0.d0
                endif
                if (k .eq. kend3-1) then
                    cc(k)=0.d0
                endif

            enddo
            enddo
        enddo

        return

    end subroutine



    ! this function calculates a,b,c for t equation
    subroutine get_t_eqn_coeff(aa,bb,cc)

        use mpi

        implicit none

        real(mytype),   dimension(:),     intent(out) :: aa,cc
        complex(mytype),dimension(:,:,:) ,intent(out) :: bb

        integer              :: i,j,k,l,m
        complex(mytype)      :: kx,ky

        aa(:)=0.d0
        bb(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        cc(:)=0.d0

        ! calculate coefficients a,b,c for the pressure equation
        do k=kstr3,kend3,1

            aa(k)=-0.5d0*dt*alpha/dz(k)/dz_b(k)
            cc(k)=-0.5d0*dt*alpha/dz(k)/dz_t(k)

            do j=1,zsize(2),1
            do i=1,zsize(1),1

                ! deal with index
                ! convert i to l
                if ( (i+i_offset) .le. nx/2) then
                    l=i+i_offset-1
                else
                    l=i+i_offset-nx-1
                endif
                ! convert j to m
                if ( (j+j_offset) .le. ny/2) then
                    m=j+j_offset-1
                else
                    m=j+j_offset-ny-1
                endif

                ! horizontal modified wavenumber,
                ! check to use 2nd or 4th schemes
                if ( cds .eq. 1) then

                    kx = 0.5d0*dt*( wx1**l + wx1**(-l) - 2.d0 )*alpha/dx2

                    ky = 0.5d0*dt*( wy1**m + wy1**(-m) - 2.d0 )*alpha/dy2

                elseif ( cds .eq. 2) then

                    kx = dt * ( -wx1**(2.d0*l) + &
                           16.d0*wx1**l        - &
                           30.d0               + &
                           16.d0*wx1**(-l)     - &
                                 wx1**(-2.d0*l)  ) /24.d0*alpha/dx2

                    ky = dt * ( -wy1**(2.d0*m) + &
                           16.d0*wy1**m        - &
                           30.d0               + &
                           16.d0*wy1**(-m)     - &
                                 wy1**(-2.d0*m)  ) /24.d0*alpha/dy2

                elseif ( cds .eq. 3) then

                    kx = 0.5d0*dt*alpha*( -(2.d0*pi/lx*l)**2.d0 ) 
                    ky = 0.5d0*dt*alpha*( -(2.d0*pi/ly*m)**2.d0 ) 

                endif

                ! note that here we should do 1-kx-ky-aa-cc for bb
                ! however, to be convenient to implement RK3, we
                ! add the missing 1 in the subroutine "assign_abc"
                bb(i,j,k)= - kx - ky - aa(k) - cc(k)

                ! Dirichlet condition
                ! we need to add the missing terms to qt 
                ! in Poisson equation solver
                if (k .eq. kstr3) then

                    bb(i,j,k)=bb(i,j,k)-aa(k)
                    aa(k)=0.d0

                endif

                ! we need to add the missing terms to qt 
                ! in Poisson equation solver
                if (k .eq. kend3) then

                    bb(i,j,k)=bb(i,j,k)-cc(k)
                    cc(k)=0.d0

                endif

            enddo
            enddo
        enddo

        return

    end subroutine


    ! Poisson equation fft solver,
    ! we call it Poisson equation solver,
    ! although we also use it to solve Helmholtz equation for U,V,W,T
    subroutine poisson_solver_fft(aa,bb,cc,r_src,r_var,zstagger,delnull,tbnd,ubnd)

        ! C standard, for using fftw
        use, intrinsic :: iso_c_binding
        use module_tools
        ! these are temporary arrays, their dimensions are: nx/p_row,ny/p_col,nz
        ! they are on z pencil
        use module_variables, only : c_src,c_var

        implicit none

        ! for using fftw
        include 'include/fftw3.f03'

        ! these are Thomas algorithm coefficients
        real(mytype),   dimension(:,:,:),intent(inout) :: aa,cc
        complex(mytype),dimension(:,:,:),intent(inout) :: bb

        ! 3d real array, input from source term
        real(mytype),   dimension(:,:,:),   intent(in) :: r_src

        ! these are flags
        ! zstagger: if is staggered grid in z direction.
        ! set it to 1 for w equation only
        ! delnull: whether to delete null space, set it to 1 for p equation only
        ! tbnd: whether to treat the Dirichlet boundary condition,
        ! set it to 1 for temperature equation only
        ! ubnd: special treatment for u boundary condition (for u_mrf),
        ! this is because in a moving coordinate, the bottom wall u is not zero!!
        ! set it to 1 for u equation only
        integer,                            intent(in) :: zstagger,delnull,tbnd,ubnd

        ! 3d real array, output for the Poisson fft solution
        real(mytype),dimension(:,:,:),     intent(out) :: r_var

        ! temporary variable used in Thomas algorithm
        complex(mytype)                                :: tmp
        ! temporary variable used in T and u eqn boundary condition
        real(mytype)                                   :: aa_t_zstr,cc_t_zend
        real(mytype)                                   :: aa_u_zstr,cc_u_zend
        integer   :: i,j,k,kmax3


        ! initiate variables
        r_var(:,:,:)=0.d0
        c_var(:,:,:)=cmplx(0.d0,0.d0 , kind=mytype )
        c_src(:,:,:)=cmplx(0.d0,0.d0 , kind=mytype )

        ! assign the real source term r_src to complex source term c_src
        ! note that by default, we are on z pencil
        do k=1,zsize(3),1
        do j=1,zsize(2),1
        do i=1,zsize(1),1
            c_src(i,j,k)=cmplx(r_src(i+ghst,j+ghst,k+ghst),0.d0, &
                                                kind=mytype)
        enddo
        enddo
        enddo

        ! mpi 2d fft for c_src
        call mpi_2d_fft(c_src,-1)

        ! check some flags
        ! check if it is w equation, the z range will be different
        if (zstagger .eq. 1) then
            kmax3=kend3-1
        elseif (zstagger .eq. 0) then
            kmax3=kend3
        else
            print *, 'Invalid zstagger flag!'
            stop
        endif
        ! check if we need to remove the null space for pressure Poisson equation
        ! we set the mean component of aa,cc,c_src to 0,
        ! and therefore the mean pressure at the 1st grid level is 0
        if (delnull .eq. 1) then
            ! remove null space
            if (myid .eq. 0) then
                aa(1,1,kstr3)=0.d0
                cc(1,1,kstr3)=0.d0
                c_src(1,1,1)=cmplx(0.d0,0.d0 , kind=mytype )
            endif
        elseif (delnull .eq. 0) then

        else
            print *, 'delnull flag error!'
            stop
        endif
        ! check if it is t equation.
        if (tbnd .eq. 1) then
            ! add the mising terms to qt for the Dirichlet condition
            if (myid .eq. 0) then
                ! applied boundary condition for temperature
                aa_t_zstr=-0.5d0*rk_c*dt*alpha/dz(kstr3)/dz_b(kstr3)
                c_src(1,1,1)=c_src(1,1,1)-2.d0*aa_t_zstr*tbot
                cc_t_zend=-0.5d0*rk_c*dt*alpha/dz(kend3)/dz_t(kend3)
                c_src(1,1,nz)=c_src(1,1,nz)-2.d0*cc_t_zend*ttop
            endif

        elseif (tbnd .eq. 0) then

        else
            print *, 'tbnd flag error!'
            stop
        endif
        ! check if it is u equation. We need some treatment for 
        ! the moving reference frame
        if (ubnd .eq. 1) then
            ! add the missing terms to qu for the Dirichlet condition
            if (myid .eq. 0) then
                ! applied boundary condition for u_mrf
                aa_u_zstr=-0.5d0*rk_c*dt*nu/dz(kstr3)/dz_b(kstr3)
                c_src(1,1,1)=c_src(1,1,1)+2.d0*aa_u_zstr*u_mrf
                if (ochannel .eq. 0) then
                    cc_u_zend=-0.5d0*rk_c*dt*nu/dz(kend3)/dz_t(kend3)
                    c_src(1,1,nz)=c_src(1,1,nz)+2.d0*cc_u_zend*u_mrf
                endif
            endif

        elseif (ubnd .eq. 0) then

        else
            print *, 'ubnd flag error!'
            stop
        endif

        ! Thomas algorithm to calculate c_var
        ! forward loop
        do k=kstr3+1,kmax3,1
        do j=1,zsize(2),1
        do i=1,zsize(1),1

            tmp=aa(i,j,k)/bb(i,j,k-1)
            bb(i,j,k)=bb(i,j,k)-cc(i,j,k-1)*tmp
            c_src(i,j,k-ghst)=c_src(i,j,k-ghst)-c_src(i,j,k-ghst-1)*tmp

        enddo
        enddo
        enddo
        ! backward loop
        do j=1,zsize(2),1
        do i=1,zsize(1),1

            c_var(i,j,kmax3-ghst)=c_src(i,j,kmax3-ghst)/bb(i,j,kmax3)

        enddo
        enddo
        
        do k=kmax3-1,kstr3,-1
        do j=1,zsize(2),1
        do i=1,zsize(1),1

            c_var(i,j,k-ghst)=( c_src(i,j,k-ghst)- &
                    cc(i,j,k)*c_var(i,j,k-ghst+1) ) / bb(i,j,k)

        enddo
        enddo
        enddo

        ! mpi 2d ifft, transfer c_var to r_var
        call mpi_2d_fft(c_var,1)

        ! assign c_var to r_var
        do k=kstr3,kmax3,1
        do j=1,zsize(2),1
        do i=1,zsize(1),1

            r_var(i+ghst,j+ghst,k)=real(c_var(i,j,k-ghst), kind=mytype )

        enddo
        enddo
        enddo

        return

    end subroutine



end module
