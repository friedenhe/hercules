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



module module_tools

    use module_parameters
    implicit none


contains



    ! this function allocates arrays, initiates fields etc
    subroutine initiate_fields

        use module_variables

        implicit none

        ! allocate arrays, by default, all the variables are on z pencil
        ! field variables
        allocate(   u(nia3,nja3,nka3) )
        allocate(   v(nia3,nja3,nka3) )
        allocate(   w(nia3,nja3,nka3) )
        allocate(   p(nia3,nja3,nka3) )
        allocate( phi(nia3,nja3,nka3) )
        allocate(   t(nia3,nja3,nka3) )

        ! convective terms
        allocate(  ff(nia3,nja3,nka3) )
        allocate(  gg(nia3,nja3,nka3) )
        allocate(  hh(nia3,nja3,nka3) )
        allocate(  ss(nia3,nja3,nka3) )

        if (cds .ne. 3) then ! for spectral method, we don't need these vars
            ! source terms
            allocate(  qu(nia3,nja3,nka3) )
            allocate(  qv(nia3,nja3,nka3) )
            allocate(  qw(nia3,nja3,nka3) )
            ! convection terms at t=n-1
            allocate( ff0(nia3,nja3,nka3) )
            allocate( gg0(nia3,nja3,nka3) )
            allocate( hh0(nia3,nja3,nka3) )
            allocate( ss0(nia3,nja3,nka3) )
            ! interpolated variables
            allocate( uhx(nia3,nja3,nka3) )
            allocate( uhy(nia3,nja3,nka3) )
            allocate( uhz(nia3,nja3,nka3) )
            allocate( vhx(nia3,nja3,nka3) )
            allocate( vhy(nia3,nja3,nka3) )
            allocate( vhz(nia3,nja3,nka3) )
            allocate( whx(nia3,nja3,nka3) )
            allocate( why(nia3,nja3,nka3) )
            allocate( thx(nia3,nja3,nka3) )
            allocate( thy(nia3,nja3,nka3) )
            allocate( thz(nia3,nja3,nka3) )
        endif

        ! mean and turbulence statistics
        allocate(  mean_u(nia3,nja3,nka3) )
        allocate(  mean_v(nia3,nja3,nka3) )
        allocate(  mean_w(nia3,nja3,nka3) )
        allocate(  mean_t(nia3,nja3,nka3) )
        allocate( mean_uu(nia3,nja3,nka3) )
        allocate( mean_vv(nia3,nja3,nka3) )
        allocate( mean_ww(nia3,nja3,nka3) )
        allocate( mean_uw(nia3,nja3,nka3) )
        allocate( mean_vw(nia3,nja3,nka3) )
        allocate( mean_tt(nia3,nja3,nka3) )
        allocate( mean_tw(nia3,nja3,nka3) )

        ! coefficients for Thomas algorithm,
        ! aa, bb, cc are also on z pencil
        allocate(   aa(zsize(1),zsize(2),zsize(3)+2*ghst) )
        allocate(   bb(zsize(1),zsize(2),zsize(3)+2*ghst) )
        allocate(   cc(zsize(1),zsize(2),zsize(3)+2*ghst) )
        allocate( aa_u(                  zsize(3)+2*ghst) )
        allocate( bb_u(zsize(1),zsize(2),zsize(3)+2*ghst) )
        allocate( cc_u(                  zsize(3)+2*ghst) )
        allocate( aa_w(                  zsize(3)+2*ghst) )
        allocate( bb_w(zsize(1),zsize(2),zsize(3)+2*ghst) )
        allocate( cc_w(                  zsize(3)+2*ghst) )
        allocate( aa_p(                  zsize(3)+2*ghst) )
        allocate( bb_p(zsize(1),zsize(2),zsize(3)+2*ghst) )
        allocate( cc_p(                  zsize(3)+2*ghst) )
        allocate( aa_t(                  zsize(3)+2*ghst) )
        allocate( bb_t(zsize(1),zsize(2),zsize(3)+2*ghst) )
        allocate( cc_t(                  zsize(3)+2*ghst) )

        ! temporary arrays
        ! c_src, and c_var are on z pencil
        allocate( c_src(zsize(1),zsize(2),zsize(3)) )
        allocate( c_var(zsize(1),zsize(2),zsize(3)) )
        ! c_tmp1, and c_tmp2 are on x and y pencil
        allocate( c_tmp1(xsize(1),xsize(2),xsize(3)) )
        allocate( c_tmp2(ysize(1),ysize(2),ysize(3)) )
        
        if (cds .eq. 3) then
            !!! variables for spectral method
            ! Fourier coefficients of the fundamental variables
            allocate(  c_u(zsize(1),zsize(2),zsize(3)) )
            allocate(  c_v(zsize(1),zsize(2),zsize(3)) )
            allocate(  c_w(zsize(1),zsize(2),zsize(3)) )
            allocate(  c_p(zsize(1),zsize(2),zsize(3)) )
            allocate(c_phi(zsize(1),zsize(2),zsize(3)) )
            allocate(  c_t(zsize(1),zsize(2),zsize(3)) )
            
            allocate(  c_qu(zsize(1),zsize(2),zsize(3)) )
            allocate(  c_qv(zsize(1),zsize(2),zsize(3)) )
            allocate(  c_qw(zsize(1),zsize(2),zsize(3)) )
            allocate(  c_ff(zsize(1),zsize(2),zsize(3)) )
            allocate(  c_gg(zsize(1),zsize(2),zsize(3)) )
            allocate(  c_hh(zsize(1),zsize(2),zsize(3)) )
            allocate(  c_ss(zsize(1),zsize(2),zsize(3)) )
        endif

        ! initiate the fields with zeros
        ! field variables
        u(:,:,:)=0.d0
        v(:,:,:)=0.d0
        w(:,:,:)=0.d0
        t(:,:,:)=0.d0
        p(:,:,:)=0.d0
        phi(:,:,:)=0.d0
        
        ! convective terms
        ff(:,:,:)=0.d0
        gg(:,:,:)=0.d0
        hh(:,:,:)=0.d0
        ss(:,:,:)=0.d0
        
        if (cds .ne. 3) then
            ! source terms
            qu(:,:,:)=0.d0
            qv(:,:,:)=0.d0
            qw(:,:,:)=0.d0
            ! convection terms at t=n-1
            ff0(:,:,:)=0.d0
            gg0(:,:,:)=0.d0
            hh0(:,:,:)=0.d0
            ss0(:,:,:)=0.d0
            ! these are interpolated variables for u,v,w,t
            uhx(:,:,:)=0.d0
            uhy(:,:,:)=0.d0
            uhz(:,:,:)=0.d0
            vhx(:,:,:)=0.d0
            vhy(:,:,:)=0.d0
            vhz(:,:,:)=0.d0
            whx(:,:,:)=0.d0
            why(:,:,:)=0.d0
            thx(:,:,:)=0.d0
            thy(:,:,:)=0.d0
            thz(:,:,:)=0.d0
        endif
        
        ! Thomas algorithm coefficients
        aa(:,:,:)=0.d0
        bb(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        cc(:,:,:)=0.d0
        aa_u(:)=0.d0
        bb_u(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        cc_u(:)=0.d0
        aa_w(:)=0.d0
        bb_w(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        cc_w(:)=0.d0
        aa_p(:)=0.d0
        bb_p(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        cc_p(:)=0.d0
        aa_t(:)=0.d0
        bb_t(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        cc_t(:)=0.d0
        
        ! mean and turbulence statistics
        mean_u(:,:,:)=0.d0
        mean_v(:,:,:)=0.d0
        mean_w(:,:,:)=0.d0
        mean_t(:,:,:)=0.d0
        mean_uu(:,:,:)=0.d0
        mean_vv(:,:,:)=0.d0
        mean_ww(:,:,:)=0.d0
        mean_uw(:,:,:)=0.d0
        mean_vw(:,:,:)=0.d0
        mean_tt(:,:,:)=0.d0
        mean_tw(:,:,:)=0.d0

        ! these are temporary arrays
        c_src(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        c_var(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        c_tmp1(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        c_tmp2(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        
        if (cds .eq. 3) then
            ! Fourier coefficients
            c_u(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
            c_v(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
            c_w(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
            c_t(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
            c_p(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
            c_phi(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
            
            c_qu(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
            c_qv(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
            c_qw(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
            c_ff(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
            c_gg(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
            c_hh(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
            c_ss(:,:,:)=cmplx(0.d0,0.d0,kind=mytype)
        endif

        if (restart .eq. 1) then

            ! read initial fields from files
            call read_init_fields(u,v,w,p,t,ff,gg,hh,ss)

            ! add some noise to the initial fields?
            if (isnoise .eq. 1) then
                call add_noise2fields(u,v,w,p,t)
            endif

        elseif (restart .eq. 0) then

            ! use default initial values (all zeros)

            ! add some noise to the initial fields?
            if (isnoise .eq. 1) then
                call add_noise2fields(u,v,w,p,t)
            endif

        else

            print *, 'restart error!'

        endif

        return
    end subroutine



    ! this function reads initial fields from files
    subroutine read_init_fields(u,v,w,p,t,ff,gg,hh,ss)

        use decomp_2d_io

        implicit none

        real(mytype),dimension(:,:,:),intent(out) :: u, v, w, t, p
        real(mytype),dimension(:,:,:),intent(out) :: ff,gg,hh,ss

        if (myid .eq. 0) then
            print *, 'Reading Initial Fields.'
        endif

        call decomp_2d_read_one(3,u(istr3:iend3,jstr3:jend3,kstr3:kend3), &
            'results/init_u.dat')

        call decomp_2d_read_one(3,v(istr3:iend3,jstr3:jend3,kstr3:kend3), &
            'results/init_v.dat')

        call decomp_2d_read_one(3,w(istr3:iend3,jstr3:jend3,kstr3:kend3), &
            'results/init_w.dat')

        call decomp_2d_read_one(3,p(istr3:iend3,jstr3:jend3,kstr3:kend3), &
            'results/init_p.dat')

        call decomp_2d_read_one(3,t(istr3:iend3,jstr3:jend3,kstr3:kend3), &
            'results/init_t.dat')

        ! subtract the moving frame velocity u_mrf from u
        u(istr3:iend3,jstr3:jend3,kstr3:kend3)= &
                              u(istr3:iend3,jstr3:jend3,kstr3:kend3)-u_mrf

        ! for RK3, we don't need to read f,g,h,s from the previous steps
        if (dts .ne. 2) then

            call decomp_2d_read_one(3,ff(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                'results/init_ff.dat')

            call decomp_2d_read_one(3,gg(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                'results/init_gg.dat')

            call decomp_2d_read_one(3,hh(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                'results/init_hh.dat')

            call decomp_2d_read_one(3,ss(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                'results/init_ss.dat')

        endif


        if (myid .eq. 0) then
            print *, 'Reading Initial Fields Finished!'
        endif

        return

    end subroutine



    ! this function adds some noise to the fields
    subroutine add_noise2fields(u,v,w,p,t)

        implicit none

        real(mytype),dimension(:,:,:),intent(inout) :: u, v, w, t, p

        integer                      :: i,j,k
        real(mytype)                 :: tmp

        if (myid .eq. 0) then
            print *,'Adding some noise to the initial fields'
        endif

        call random_seed

        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            call random_number(tmp)
            u(i,j,k)=u(i,j,k)+u(i,j,k)*tmp*noise_mag
            call random_number(tmp)
            v(i,j,k)=v(i,j,k)+v(i,j,k)*tmp*noise_mag
            call random_number(tmp)
            w(i,j,k)=w(i,j,k)+w(i,j,k)*tmp*noise_mag
            call random_number(tmp)
            p(i,j,k)=p(i,j,k)+p(i,j,k)*tmp*noise_mag
            call random_number(tmp)
            t(i,j,k)=t(i,j,k)+t(i,j,k)*tmp*noise_mag

        enddo
        enddo
        enddo

        return

    end subroutine



    ! this function calculates mean fields and turbulence statistics
    subroutine calculate_mean_fields &
    (u,v,w,t,uhx,vhy,mean_u,mean_v,mean_w,mean_t, mean_uu,mean_vv, &
     mean_ww,mean_uw,mean_vw,mean_tt,mean_tw)

        implicit none

        real(mytype),dimension(:,:,:),intent(in) :: u,v,w,t,uhx, vhy
        real(mytype),dimension(:,:,:),intent(out):: &
         mean_u,mean_v,mean_w,mean_t, mean_uu,mean_vv,mean_ww,mean_uw, mean_vw, &
         mean_tt,mean_tw

        integer                           :: i,j,k

        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            ! note that here we need to add the speed of moving reference frame
            mean_uu(i,j,k)=mean_uu(i,j,k)+(u(i,j,k)+u_mrf)**2.d0
            mean_vv(i,j,k)=mean_vv(i,j,k)+v(i,j,k)**2.d0
            mean_ww(i,j,k)=mean_ww(i,j,k)+w(i,j,k)**2.d0
            mean_uw(i,j,k)=mean_uw(i,j,k)+0.5d0*( w(i,j,k)+w(i,j,k-1) ) * &
                                                ( uhx(i,j,k)+u_mrf )
            mean_vw(i,j,k)=mean_vw(i,j,k)+0.5d0*( w(i,j,k)+w(i,j,k-1) ) * vhy(i,j,k) 
                                                 
            mean_tt(i,j,k)=mean_tt(i,j,k)+t(i,j,k)**2.d0
            mean_tw(i,j,k)=mean_tw(i,j,k)+t(i,j,k) * 0.5d0 * &
                                           ( w(i,j,k)+w(i,j,k-1) )

            ! note that here we need to add the speed of moving reference frame
            mean_u(i,j,k)=mean_u(i,j,k) + u(i,j,k)+u_mrf
            mean_v(i,j,k)=mean_v(i,j,k) + v(i,j,k)
            mean_w(i,j,k)=mean_w(i,j,k) + w(i,j,k)
            mean_t(i,j,k)=mean_t(i,j,k) + t(i,j,k)

        enddo
        enddo
        enddo


        return

    end subroutine



    ! this function calculates the interpolated fields for u,v,w
    ! we don't need to call this function for spectral method
    subroutine get_interp_fields_uvwhxyz(u,v,w,uhx,uhy,uhz,vhx,vhy,vhz,whx,why)

        implicit none

        real(mytype),dimension(:,:,:),intent(in) :: u,v,w
        real(mytype),dimension(:,:,:),intent(out):: uhx,uhy,uhz,vhx,vhy,vhz,&
                                                          whx,why

        integer                           :: i,j,k

        ! calculate interpolation variables
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1
        
            if (cds .eq. 1) then ! linear interpolation
            
                uhx(i,j,k)=0.5d0*( u(i,j,k)+u(i-1,j,k) )
                uhy(i,j,k)=0.5d0*( u(i,j+1,k)+u(i,j,k) )
                vhx(i,j,k)=0.5d0*( v(i+1,j,k)+v(i,j,k) )
                vhy(i,j,k)=0.5d0*( v(i,j,k)+v(i,j-1,k) )
                whx(i,j,k)=0.5d0*( w(i+1,j,k)+w(i,j,k) )
                why(i,j,k)=0.5d0*( w(i,j+1,k)+w(i,j,k) )
            
            elseif (cds .eq. 2) then ! 4th order interpolation
            
                uhx(i,j,k)=( -u(i+1,j,k)+9.d0*u(i,j,k)+ &
                    9.d0*u(i-1,j,k)-u(i-2,j,k) ) / 16.d0
            
                uhy(i,j,k)=( -u(i,j+2,k)+9.d0*u(i,j+1,k)+ &
                    9.d0*u(i,j,k)-u(i,j-1,k) ) / 16.d0
            
                vhx(i,j,k)=( -v(i+2,j,k)+9.d0*v(i+1,j,k)+ &
                    9.d0*v(i,j,k)-v(i-1,j,k) ) / 16.d0
            
                vhy(i,j,k)=( -v(i,j+1,k)+9.d0*v(i,j,k)+ &
                    9.d0*v(i,j-1,k)-v(i,j-2,k) ) / 16.d0
            
                whx(i,j,k)=( -w(i+2,j,k)+9.d0*w(i+1,j,k)+ &
                    9.d0*w(i,j,k)-w(i-1,j,k) ) / 16.d0
            
                why(i,j,k)=( -w(i,j+2,k)+9.d0*w(i,j+1,k)+ &
                    9.d0*w(i,j,k)-w(i,j-1,k) ) / 16.d0
            
            endif

            ! this is always 2nd order interpolation
            uhz(i,j,k)=u(i,j,k+1)*l_t(k)+u(i,j,k)*(1.d0-l_t(k))
            ! this is always 2nd order interpolation
            vhz(i,j,k)=v(i,j,k+1)*l_t(k)+v(i,j,k)*(1.d0-l_t(k))

        enddo
        enddo
        enddo


        return

    end subroutine


    ! this function calculates the interpolated fields for T
    ! we don't need to call this function for spectral method
    subroutine get_interp_fields_thxyz(t,thx,thy,thz)

        implicit none

        real(mytype),dimension(:,:,:),intent(in) :: t
        real(mytype),dimension(:,:,:),intent(out):: thx,thy,thz

        integer                           :: i,j,k

        ! calculate interpolation variables
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1
        
            if (cds .eq. 1) then
            
                thx(i,j,k)=0.5d0*( t(i+1,j,k)+t(i,j,k) )
                thy(i,j,k)=0.5d0*( t(i,j+1,k)+t(i,j,k) )
            
            elseif (cds .eq. 2) then
            
                thx(i,j,k)=( -t(i+2,j,k)+9.d0*t(i+1,j,k)+ &
                         9.d0*t(i,j,k)-t(i-1,j,k) ) / 16.d0
            
                thy(i,j,k)=( -t(i,j+2,k)+9.d0*t(i,j+1,k)+ &
                         9.d0*t(i,j,k)-t(i,j-1,k) ) / 16.d0
            
            endif

            ! this is always 2nd order interpolation
            thz(i,j,k)=t(i,j,k+1)*l_t(k)+t(i,j,k)*(1.d0-l_t(k))

        enddo
        enddo
        enddo


        return

    end subroutine



    ! this function conducts 2d fft or ifft
    ! we use 2d decomposition lib for array transposition between processors
    subroutine mpi_2d_fft(c_inout3,direction)

        ! C standard for using fftw
        use, intrinsic :: iso_c_binding
        ! this is a temporary array, it is on x/y pencil, it is used for 2d fft
        use module_variables, only : c_tmp1,c_tmp2

        implicit none

        ! for using fftw
        include 'include/fftw3.f03'

        ! flag: -1=fft; 1=ifft
        integer,                         intent(in)    :: direction
        ! note that c_inout3 is on z pencil
        complex(mytype),dimension(:,:,:),intent(inout) :: c_inout3

        integer   :: i,j,k,l,m


        ! we always remove the Nyquist wavenumber
        if (direction .eq. 1) then
        
            do k=1,zsize(3),1
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
            
                ! remove Nyquist wavenumber
                if (l .eq. -nx/2) then
                    c_inout3(i,j,k)=cmplx(0.d0,0.d0,kind=mytype)
                endif
            
                if (m .eq. -ny/2) then
                    c_inout3(i,j,k)=cmplx(0.d0,0.d0,kind=mytype)
                endif
                
            enddo
            enddo
            enddo
            
        endif


        ! mpi 2d fft
        ! first transpose from z to y
        call transpose_z_to_y(c_inout3,c_tmp2)

        ! check it is a forward (direction=-1) or a backward (direction=-1) fft
        ! note that we are on y pencil now!!
        if (direction .eq. -1) then

            ! do 1d fft for c_tmp2, we are on y pencil
            do k=1,ysize(3),1
            do i=1,ysize(1),1 ! y pencil index

                call fftw_execute_dft(fft_plan2,c_tmp2(i,:,k), &
                                                c_tmp2(i,:,k))

            enddo
            enddo
            ! ***NOTE*** since fftw is not normalized,
            ! we need to normalize it after the forward fft
            c_tmp2(:,:,:)=c_tmp2(:,:,:)/ny

        elseif (direction .eq. 1) then

            ! do 1d ifft for c_tmp2, we are on y pencil
            do k=1,ysize(3),1
            do i=1,ysize(1),1 ! y pencil index

                call fftw_execute_dft(ifft_plan2,c_tmp2(i,:,k), &
                                                 c_tmp2(i,:,k))

            enddo
            enddo

        else

            print *, "direction is either -1 (forward) or 1 (backward)"
            stop

        endif

        ! transpose from y to x
        call transpose_y_to_x(c_tmp2,c_tmp1)

        ! check if it is a forward (direction=-1) or a backward (direction=-1) fft
        ! NOTE we are on x pencil now!!!
        if (direction .eq. -1) then

            ! do 1d fft for c_tmp1, we are on x pencil
            do k=1,xsize(3),1
            do j=1,xsize(2),1 ! x pencil index

                call fftw_execute_dft(fft_plan1,c_tmp1(:,j,k), &
                                                c_tmp1(:,j,k))

            enddo
            enddo
            ! ***NOTE*** since fftw is not normalized,
            ! we need to normalize it after the forward fft
            c_tmp1(:,:,:)=c_tmp1(:,:,:)/nx

        elseif (direction .eq. 1) then

            ! do 1d ifft for c_tmp1, we are on x pencil
            do k=1,xsize(3),1
            do j=1,xsize(2),1 ! x pencil index

                call fftw_execute_dft(ifft_plan1,c_tmp1(:,j,k), &
                                                 c_tmp1(:,j,k))

            enddo
            enddo

        else

            print *, "direction is either -1 (forward) or 1 (backward)"
            stop

        endif

        ! transpose back from x pencil to z pencil
        ! the output is c_inout3
        call transpose_x_to_y(c_tmp1,c_tmp2)
        call transpose_y_to_z(c_tmp2,c_inout3)

        ! we always remove the Nyquist wave number
        if (direction .eq. -1) then
        
            do k=1,zsize(3),1
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
            
                ! remove Nyquist wavenumber
                if (l .eq. -nx/2) then
                    c_inout3(i,j,k)=cmplx(0.d0,0.d0,kind=mytype)
                endif
            
                if (m .eq. -ny/2) then
                    c_inout3(i,j,k)=cmplx(0.d0,0.d0,kind=mytype)
                endif
                
            enddo
            enddo
            enddo
            
        endif


        ! mpi 2d fft finished

        return

    end subroutine



    ! Rayleigh damping layer at the top of the boundary layer
    subroutine rayleigh_damping(varin)

        use mpi

        implicit none

        real(mytype),dimension(:,:,:),intent(inout) :: varin

        real(mytype)                      :: zd,sigz,mean_local,mean_all
        integer                           :: i,j,k


        ! do the damping for the topmost nzdamp points
        do k=kend3-nzdamp,kend3,1
            
            ! first calculate the mean of the variable on the local processor
            mean_local=0.d0
            do j=jstr3,jend3,1
            do i=istr3,iend3,1
                mean_local=mean_local+varin(i,j,k)
            enddo
            enddo
            mean_local=mean_local/zsize(1)/zsize(2)

            ! exchange mean_local between processors
            ! sum up all the mean_local between all processors
            call MPI_ALLREDUCE(mean_local,mean_all,1,real_type,MPI_SUM, &
                                                   MPI_COMM_WORLD,ierr)
            mean_all=mean_all/nprc ! calculate mean all
            
            ! do the damping
            zd=(k-kend3+nzdamp)*1.d0/nzdamp
            sigz=cdamp*0.5d0*( 1.d0-cos(pi*zd) )
            do j=jstr3,jend3,1
            do i=istr3,iend3,1
                varin(i,j,k)=varin(i,j,k)-sigz*( varin(i,j,k)-mean_all )
            enddo
            enddo

        enddo


        return

    end subroutine
    


    


end module



