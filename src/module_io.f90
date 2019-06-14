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



module module_io

    use module_parameters
    implicit none

contains



    ! this function prints cpu time, CFL number etc to screen
    subroutine screen_cpu_time(u,v,w)

        use mpi

        implicit none

        real(mytype),dimension(:,:,:),intent(in) :: u, v, w

        real(mytype)                      :: cpu_t,cfl
        integer                           :: i,j,k

        ! calculate cfl number
        cfl=0.d0
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            if ( ( dt*abs(u(i,j,k))/dx + &
                   dt*abs(v(i,j,k))/dy + &
                   dt*abs(w(i,j,k))/dz(k) ) .gt. cfl) then

                cfl=dt*abs(u(i,j,k))/dx + &
                    dt*abs(v(i,j,k))/dy + &
                    dt*abs(w(i,j,k))/dz(k)

             endif

        enddo
        enddo
        enddo

        ! exchange CFL between processors
        ! and find the max CFL number in the domain
        call MPI_REDUCE(cfl,cfl_max,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierr)

        ! calculate the simulation time
        time_sim = time_sim+dt

        ! print the CFL max to screen
        if (myid .eq. 0) then
            call cpu_time(cpu_t)
            print *, " "
            write (*,184) time_sim,cpu_t,cfl_max
            184 format(' Time = ',es12.5,'; CPU Time = ',es11.4,'; CFL max = ',f5.2)
        endif

        return

    end subroutine



    ! this function calculates the divergence of the flow fields and prints to screen
    subroutine screen_div_error(u,v,w,isfinal,ishistout)

        use mpi
        use module_tools

        implicit none

        real(mytype),dimension(:,:,:),intent(in)  :: u,v,w
        integer,intent(in)                        :: isfinal,ishistout

        integer                           :: i,j,k
        real(mytype)                      :: dudx,dvdy,dwdz
        real(mytype)                      :: div_avg, div_max
        real(mytype)                      :: div_avg_all,div_max_all
        real(mytype)                      :: div_hist_max_all

        ! calculate the divergence
        div_avg=0.d0
        div_max=0.d0
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            ! check to use 2nd or 4th schemes
            if ( cds .eq. 1) then

                dudx = ( u(i,j,k)-u(i-1,j,k) ) / dx
                dvdy = ( v(i,j,k)-v(i,j-1,k) ) / dy

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

            dwdz = ( w(i,j,k)-w(i,j,k-1) ) / dz(k)
            div_avg=div_avg+abs(dudx+dvdy+dwdz)

            if ( abs(dudx+dvdy+dwdz) .gt. div_max ) then
                div_max=(dudx+dvdy+dwdz)
            endif

            if (    u(i,j,k) .ne. u(i,j,k) .or. &
                    v(i,j,k) .ne. v(i,j,k) .or. &
                    w(i,j,k) .ne. w(i,j,k) ) then
                print *, 'NaN error! ',i,' ',j,' ',k
                stop
            endif

            if (isfinal .eq. 1) then
                if ( abs(dudx+dvdy+dwdz) .gt. div_hist_max ) then
                    div_hist_max=(dudx+dvdy+dwdz)
                endif
            endif

        enddo
        enddo
        enddo
        div_avg=div_avg/nx/ny/nz

        ! exchange div_avg div_max between processors
        ! sum up all the div_avg between all processors
        call MPI_REDUCE(div_avg,div_avg_all,1,real_type,MPI_SUM,0, &
                                                   MPI_COMM_WORLD,ierr)
        ! find the max div_max between all processors
        call MPI_REDUCE(div_max,div_max_all,1,real_type,MPI_MAX,0, &
                                                   MPI_COMM_WORLD,ierr)
        ! find the max div_hist_max between all processors
        call MPI_REDUCE(div_hist_max,div_hist_max_all,1,real_type,MPI_MAX,0, &
                                                   MPI_COMM_WORLD,ierr)

        ! print div to screen
        if (myid .eq. 0) then

            if ( isfinal .eq. 0) then

                write (*,104) div_avg_all,div_max_all
                104 format(' DIV-ERR INIT : DIV_AVG = ',es11.4,';  DIV_MAX = ',es11.4)

            elseif ( isfinal .eq. 1 .and. ishistout .eq. 0) then

                write (*,105) div_avg_all,div_max_all
                105 format(' DIV-ERR FINAL: DIV_AVG = ',es11.4,';  DIV_MAX = ',es11.4)

            elseif ( isfinal .eq. 1 .and. ishistout .eq. 1) then

                write (*,115) div_avg_all,div_max_all
                115 format(' DIV-ERR FINAL: DIV_AVG = ',es11.4,';  DIV_MAX = ',es11.4)
                write (*,106) div_hist_max_all
                106 format(' DIV-ERR HIST_MAX = ', es11.4)

                ! print ri if is_ri_var=1
                if (is_ri_var .eq. 1) then
                    write (*,258) got
                    258 format(' Ri_b = ', es12.5)
                endif

                ! print bulk velocity if dp_opt=2
                if (dp_opt .eq. 2) then
                    write (*,259) ubar+u_mrf
                    259 format(' Ubar = ', es12.5)
                endif

            else
                print *,'div out flag error!'
                stop
            endif

        endif


        return

    end subroutine



    ! output time-series
    subroutine output_time_series(u,v,w,p,t)

        use mpi

        implicit none

        real(mytype),dimension(:,:,:),intent(in)  :: u,v,w,t,p

        integer                           :: k,itmsr,jtmsr

        itmsr=ghst+1
        jtmsr=ghst+1

        open(99,file='results/tmsr_u.dat',status='old', &
            position='append',action='write')
        do k=kstr3,kend3,1

            if (k .eq. kstr3) then
                ! note that here we need to add the speed of moving reference frame
                write(99,'(2f14.8)',advance='no') time_sim,u(itmsr,jtmsr,k)+u_mrf
            elseif (k .eq. kend3) then
                write(99,'(1f14.8)') u(itmsr,jtmsr,k)+u_mrf
            else
                write(99,'(1f14.8)',advance='no') u(itmsr,jtmsr,k)+u_mrf
            endif

        enddo
        close(99)

        open(98,file='results/tmsr_v.dat',status='old', &
            position='append',action='write')
        do k=kstr3,kend3,1

            if (k .eq. kstr3) then
                write(98,'(2f14.8)',advance='no') time_sim,v(itmsr,jtmsr,k)
            elseif (k .eq. kend3) then
                write(98,'(1f14.8)') v(itmsr,jtmsr,k)
            else
                write(98,'(1f14.8)',advance='no') v(itmsr,jtmsr,k)
            endif

        enddo
        close(98)


        open(97,file='results/tmsr_w.dat',status='old', &
            position='append',action='write')
        do k=kstr3,kend3,1

            if (k .eq. kstr3) then
                write(97,'(2f14.8)',advance='no') time_sim,w(itmsr,jtmsr,k)
            elseif (k .eq. kend3) then
                write(97,'(1f14.8)') w(itmsr,jtmsr,k)
            else
                write(97,'(1f14.8)',advance='no') w(itmsr,jtmsr,k)
            endif

        enddo
        close(97)


        open(96,file='results/tmsr_t.dat',status='old', &
            position='append',action='write')
        do k=kstr3,kend3,1

            if (k .eq. kstr3) then
                write(96,'(2f14.8)',advance='no') time_sim,t(itmsr,jtmsr,k)
            elseif (k .eq. kend3) then
                write(96,'(1f14.8)') t(itmsr,jtmsr,k)
            else
                write(96,'(1f14.8)',advance='no') t(itmsr,jtmsr,k)
            endif

        enddo
        close(96)

        open(95,file='results/tmsr_p.dat',status='old', &
            position='append',action='write')
        do k=kstr3,kend3,1

            if (k .eq. kstr3) then
                write(95,'(2f14.8)',advance='no') time_sim,p(itmsr,jtmsr,k)
            elseif (k .eq. kend3) then
                write(95,'(1f14.8)') p(itmsr,jtmsr,k)
            else
                write(95,'(1f14.8)',advance='no') p(itmsr,jtmsr,k)
            endif

        enddo
        close(95)


        return

    end subroutine



    ! note that we store all the backup fields in staggered grids
    ! instead of cell center
    subroutine output_backup(u,v,w,p,t,ff,gg,hh,ss)

        use decomp_2d_io

        implicit none

        real(mytype),dimension(:,:,:),intent(in) :: u,v,w,p,t,ff,gg,hh,ss

        if (myid .eq. 0) then
            print *, 'Writing Backup Files...'
        endif

        ! note that here we need to add the speed of moving reference frame
        call decomp_2d_write_one(3,u(istr3:iend3,jstr3:jend3,kstr3:kend3)+u_mrf, &
            'results/bak_u.dat')

        call decomp_2d_write_one(3,v(istr3:iend3,jstr3:jend3,kstr3:kend3), &
            'results/bak_v.dat')

        call decomp_2d_write_one(3,w(istr3:iend3,jstr3:jend3,kstr3:kend3), &
            'results/bak_w.dat')

        call decomp_2d_write_one(3,p(istr3:iend3,jstr3:jend3,kstr3:kend3), &
            'results/bak_p.dat')

        call decomp_2d_write_one(3,t(istr3:iend3,jstr3:jend3,kstr3:kend3), &
            'results/bak_t.dat')

        ! for RK3, we don't need to write f,g,h,s
        ! note here the ff gg hh ss are based on the relative velocity (u-u_mrf)
        ! if running a restart simulation using AB2 or AB3, one should not change u_mrf
        if (dts .ne. 2) then

            call decomp_2d_write_one(3,ff(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                'results/bak_ff.dat')

            call decomp_2d_write_one(3,gg(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                'results/bak_gg.dat')

            call decomp_2d_write_one(3,hh(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                'results/bak_hh.dat')

            call decomp_2d_write_one(3,ss(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                'results/bak_ss.dat')

        endif

        return

    end subroutine



    ! note that we store all the instantaneous fields in staggered grids
    ! instead of cell center
    subroutine output_inst_fields(u,v,w,p,t,time_step)

        use decomp_2d_io

        implicit none

        real(mytype),dimension(:,:,:),intent(in) :: u, v, w, p, t
        integer, intent(in)               :: time_step

        character(len=8)                  :: fmt1
        character(len=6)                  :: xx1

        fmt1= '(i3.3)'
        write (xx1,fmt1) time_step/iinstfl

        if (myid .eq. 0) then
            print *, 'Writing Instantaneous Fields...'
        endif

        ! note that here we need to add the speed of moving reference frame
        call decomp_2d_write_one(3,u(istr3:iend3,jstr3:jend3,kstr3:kend3)+u_mrf, &
            'results/inst_u'//trim(xx1)//'.dat')

        call decomp_2d_write_one(3,v(istr3:iend3,jstr3:jend3,kstr3:kend3), &
            'results/inst_v'//trim(xx1)//'.dat')

        call decomp_2d_write_one(3,w(istr3:iend3,jstr3:jend3,kstr3:kend3), &
            'results/inst_w'//trim(xx1)//'.dat')

        call decomp_2d_write_one(3,p(istr3:iend3,jstr3:jend3,kstr3:kend3), &
            'results/inst_p'//trim(xx1)//'.dat')

        call decomp_2d_write_one(3,t(istr3:iend3,jstr3:jend3,kstr3:kend3), &
            'results/inst_t'//trim(xx1)//'.dat')

        return

    end subroutine



    ! write 2d slice
    subroutine output_2d_slices(u,v,w,p,t,time_step)

        use decomp_2d_io

        implicit none

        real(mytype),dimension(:,:,:),intent(in) :: u, v, w, p, t
        integer, intent(in)             :: time_step

        character(len=8)       :: fmt0
        character(len=8)       :: fmt1
        character(len=6)       :: ff1
        character(len=6)       :: ii1

        integer                :: size_xy=0,size_xz=0,size_yz=0,n,test_size,indx

        fmt0= '(i6.6)'
        fmt1= '(i4.4)'
        write (ff1,fmt0) time_step/intv_2d

        ! output xy slices
        if (isxy2d .eq. 1) then

            if (myid .eq. 0) then
                print *, 'Writing X-Y Slices...'
            endif

            ! we need to check how many slices to extract
            test_size=size(xy2d_id)
            indx=1
            do n=1,test_size,1
                if (xy2d_id(n) .eq. 0) then
                    size_xy=indx-1
                    exit
                endif
                indx=indx+1
            enddo

            do n=1,size_xy,1

                write (ii1,fmt1) xy2d_id(n)
                ! note that here we need to add the speed of moving reference frame
                call decomp_2d_write_plane(3,u(istr3:iend3,jstr3:jend3,kstr3:kend3)+u_mrf, &
                 3,xy2d_id(n),'results/slices/xy2d_u_k'//trim(ii1)//'_t'//trim(ff1) )
                call decomp_2d_write_plane(3,v(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                 3,xy2d_id(n),'results/slices/xy2d_v_k'//trim(ii1)//'_t'//trim(ff1) )
                call decomp_2d_write_plane(3,w(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                 3,xy2d_id(n),'results/slices/xy2d_w_k'//trim(ii1)//'_t'//trim(ff1) )
                call decomp_2d_write_plane(3,p(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                 3,xy2d_id(n),'results/slices/xy2d_p_k'//trim(ii1)//'_t'//trim(ff1) )
                call decomp_2d_write_plane(3,t(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                 3,xy2d_id(n),'results/slices/xy2d_t_k'//trim(ii1)//'_t'//trim(ff1) )

            enddo

        endif

        ! output xz slices
        if (isxz2d .eq. 1) then

            if (myid .eq. 0) then
                print *, 'Writing X-Z Slices...'
            endif

            ! we need to check how many slices to extract
            test_size=size(xz2d_id)
            indx=1
            do n=1,test_size,1
                if (xz2d_id(n) .eq. 0) then
                    size_xz=indx-1
                    exit
                endif
                indx=indx+1
            enddo

            do n=1,size_xz,1

                write (ii1,fmt1) xz2d_id(n)
                ! note that here we need to add the speed of moving reference frame
                call decomp_2d_write_plane(3,u(istr3:iend3,jstr3:jend3,kstr3:kend3)+u_mrf, &
                 2,xz2d_id(n),'results/slices/xz2d_u_j'//trim(ii1)//'_t'//trim(ff1) )
                call decomp_2d_write_plane(3,v(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                 2,xz2d_id(n),'results/slices/xz2d_v_j'//trim(ii1)//'_t'//trim(ff1) )
                call decomp_2d_write_plane(3,w(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                 2,xz2d_id(n),'results/slices/xz2d_w_j'//trim(ii1)//'_t'//trim(ff1) )
                call decomp_2d_write_plane(3,p(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                 2,xz2d_id(n),'results/slices/xz2d_p_j'//trim(ii1)//'_t'//trim(ff1) )
                call decomp_2d_write_plane(3,t(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                 2,xz2d_id(n),'results/slices/xz2d_t_j'//trim(ii1)//'_t'//trim(ff1) )

            enddo

        endif

        ! output yz slices
        if (isyz2d .eq. 1) then

            if (myid .eq. 0) then
                print *, 'Writing Y-Z Slices...'
            endif

            ! we need to check how many slices to extract
            test_size=size(yz2d_id)
            indx=1
            do n=1,test_size,1
                if (yz2d_id(n) .eq. 0) then
                    size_yz=indx-1
                    exit
                endif
                indx=indx+1
            enddo

            do n=1,size_yz,1

                write (ii1,fmt1) yz2d_id(n)
                ! note that here we need to add the speed of moving reference frame
                call decomp_2d_write_plane(3,u(istr3:iend3,jstr3:jend3,kstr3:kend3)+u_mrf, &
                 1,yz2d_id(n),'results/slices/yz2d_u_i'//trim(ii1)//'_t'//trim(ff1) )
                call decomp_2d_write_plane(3,v(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                 1,yz2d_id(n),'results/slices/yz2d_v_i'//trim(ii1)//'_t'//trim(ff1) )
                call decomp_2d_write_plane(3,w(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                 1,yz2d_id(n),'results/slices/yz2d_w_i'//trim(ii1)//'_t'//trim(ff1) )
                call decomp_2d_write_plane(3,p(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                 1,yz2d_id(n),'results/slices/yz2d_p_i'//trim(ii1)//'_t'//trim(ff1) )
                call decomp_2d_write_plane(3,t(istr3:iend3,jstr3:jend3,kstr3:kend3), &
                 1,yz2d_id(n),'results/slices/yz2d_t_i'//trim(ii1)//'_t'//trim(ff1) )

            enddo

        endif


        return

    end subroutine



    subroutine output_mean_fields(mean_u,mean_v,mean_w,mean_t,mean_uu, &
               mean_vv,mean_ww,mean_uw,mean_vw,mean_tt,mean_tw,time_step)

        use decomp_2d_io

        implicit none

        real(mytype),dimension(:,:,:),intent(in) :: &
        mean_u,mean_v,mean_w,mean_t,mean_uu,mean_vv,mean_ww,mean_uw,mean_vw, &
        mean_tt,mean_tw
        integer, intent(in)                      :: time_step

        if (myid .eq. 0) then
            print *, 'Writing Mean Fields...'
        endif

        call decomp_2d_write_one(3,mean_u(istr3:iend3,jstr3:jend3,kstr3:kend3)/time_step, &
         'results/mean_u.dat')

        call decomp_2d_write_one(3,mean_v(istr3:iend3,jstr3:jend3,kstr3:kend3)/time_step, &
        'results/mean_v.dat')

        call decomp_2d_write_one(3,mean_w(istr3:iend3,jstr3:jend3,kstr3:kend3)/time_step, &
        'results/mean_w.dat')

        call decomp_2d_write_one(3,mean_t(istr3:iend3,jstr3:jend3,kstr3:kend3)/time_step, &
        'results/mean_t.dat')

        call decomp_2d_write_one(3,mean_uu(istr3:iend3,jstr3:jend3,kstr3:kend3)/time_step, &
        'results/mean_uu.dat')

        call decomp_2d_write_one(3,mean_vv(istr3:iend3,jstr3:jend3,kstr3:kend3)/time_step, &
        'results/mean_vv.dat')

        call decomp_2d_write_one(3,mean_ww(istr3:iend3,jstr3:jend3,kstr3:kend3)/time_step, &
        'results/mean_ww.dat')

        call decomp_2d_write_one(3,mean_uw(istr3:iend3,jstr3:jend3,kstr3:kend3)/time_step, &
        'results/mean_uw.dat')

        call decomp_2d_write_one(3,mean_vw(istr3:iend3,jstr3:jend3,kstr3:kend3)/time_step, &
        'results/mean_vw.dat')

        call decomp_2d_write_one(3,mean_tt(istr3:iend3,jstr3:jend3,kstr3:kend3)/time_step, &
        'results/mean_tt.dat')

        call decomp_2d_write_one(3,mean_tw(istr3:iend3,jstr3:jend3,kstr3:kend3)/time_step, &
        'results/mean_tw.dat')


        return

    end subroutine



end module



