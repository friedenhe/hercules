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



module module_spectral

    use module_parameters
    implicit none



contains



    ! this function calculates spatial derivatives and the 
    ! right hand side terms for momentum equations
    subroutine get_momentum_cnvdiff_spec &
        (u,v,w,p,t,c_u,c_v,c_w,c_p,c_ff,c_gg,c_hh,c_qu,c_qv,c_qw)
                                           
        ! these are temporary arrays, their dimensions are: nx/p_row,ny/p_col,nz
        ! they are on z pencil
        use module_variables, only : c_src,c_var
        use module_tools
        
        implicit none

        real(mytype),dimension(:,:,:),intent(in)      :: u,v,w,p,t
        complex(mytype),dimension(:,:,:),intent(in)   :: c_u,c_v,c_w,c_p
        complex(mytype),dimension(:,:,:),intent(out)  :: c_ff,c_gg,c_hh
        complex(mytype),dimension(:,:,:),intent(out)  :: c_qu,c_qv,c_qw
        
        integer                :: i,j,k
        real(mytype)           :: tmp1,tmp2,a1,a2
        real(mytype)           :: coe_sk

        ! check if the skew-symmetric form is used
        if (issk .eq. 1) then
            coe_sk = 0.5d0
        else
            coe_sk = 1.d0
        endif 
        
        ! first calculate dpx_drive
        call update_dpx_drive

        !!!!!!!!!! U equation !!!!!!!!!!!
        !!!!!!!!!! now start building the right hand side !!!!!!!!!!
        ! qu=u^n+dt*(rk_a*ff^n+rk_b*ff^n-1) + &
        !               rk_c*0.5*dt*nu*(dudx2+dudy^2+dudz^2)^n-rk_c*dt*dpdx^n

        ! add u
        c_qu(:,:,:)=c_u(:,:,:)
        
        ! add dudx2
        call dev_fft_linear_cmplx(c_u,c_var,2,1) 
        c_qu(:,:,:)=c_qu(:,:,:)+rk_c*0.5d0*dt*nu*c_var(:,:,:)
        ! add dudy2
        call dev_fft_linear_cmplx(c_u,c_var,2,2) 
        c_qu(:,:,:)=c_qu(:,:,:)+rk_c*0.5d0*dt*nu*c_var(:,:,:)
        ! add dudz2, we store dudz2 to tmp1 in physical space
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            tmp1=  ( (u(i,j,k+1)-u(i,j,k))/dz_t(k) - &
                     (u(i,j,k)-u(i,j,k-1))/dz_b(k) ) / dz(k)
            c_var(i-ghst,j-ghst,k-ghst)=cmplx(tmp1, 0.d0, kind=mytype)

        enddo
        enddo
        enddo
        call mpi_2d_fft(c_var,-1)
        c_qu(:,:,:)=c_qu(:,:,:)+rk_c*0.5d0*dt*nu*c_var(:,:,:)
        
        ! add dpdx
        call dev_fft_linear_cmplx(c_p,c_var,1,1) 
        c_qu(:,:,:)=c_qu(:,:,:)-rk_c*dt*c_var(:,:,:)
        
        ! add ff^n-1
        c_qu(:,:,:)=c_qu(:,:,:)+rk_b*dt*c_ff(:,:,:)
        
        !!!!!!!!!!!!! calculate ff !!!!!!!!!!
        ! ff=-0.5*(duudx+duvdy+duwdz)^n - &
        !        0.5*(u*dudx+v*dudy+w*dudz)^n+dpx_drive-u_mrf*dudx^n
        c_ff(:,:,:)=cmplx( 0.d0,0.d0,kind=mytype )
        ! add u_mrf*dudx
        call dev_fft_linear_cmplx(c_u,c_var,1,1) 
        c_ff(:,:,:)=c_ff(:,:,:)-u_mrf*c_var(:,:,:)
        ! add duudx
        call dev_fft_nonlinear_cmplx(u,u,c_var,1,0)
        c_ff(:,:,:)=c_ff(:,:,:)-coe_sk*c_var(:,:,:)
        ! add duvdy
        call dev_fft_nonlinear_cmplx(u,v,c_var,2,0)
        c_ff(:,:,:)=c_ff(:,:,:)-coe_sk*c_var(:,:,:)
        ! add duwdz and dpx_drive, we store duwdz to tmp1 in physical space
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            a1=u(i,j,k+1)*l_t(k)+u(i,j,k)*(1.d0-l_t(k))
            a2=u(i,j,k-1)*l_b(k)+u(i,j,k)*(1.d0-l_b(k))
            tmp1   = ( a1 * w(i,j,k) - &
                       a2 * w(i,j,k-1) ) / dz(k)
                       
            if (dp_opt .eq. 3) then
                tmp2   = fc*( v(i,j,k) - vg )
            else
                tmp2   = dpx_drive
            endif
            
            ! note that here we use 1/coe_sk*tmp2
            ! since we will do a coe_sk*c_var later
            c_var(i-ghst,j-ghst,k-ghst)= cmplx(1.d0/coe_sk*tmp2 - tmp1, &
                                               0.d0, kind=mytype )

        enddo
        enddo
        enddo
        call mpi_2d_fft(c_var,-1)
        c_ff(:,:,:)=c_ff(:,:,:)+coe_sk*c_var(:,:,:)
        
        ! check if the skew-symmetric form is used, 
        ! if yes we need to add the other terms
        if (issk .eq. 1) then
        
            ! add u*dudx
            call dev_fft_linear_cmplx(c_u,c_src,1,1) ! c_dudx
            call mpi_2d_fft(c_src,1)
            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1
                
                ! multiply u and dudx in physical space
                c_var(i-ghst,j-ghst,k-ghst)= &
                cmplx( u(i,j,k) * real(c_src(i-ghst,j-ghst,k-ghst),kind=mytype) , &
                                                   0.d0, kind=mytype )
            
            enddo
            enddo
            enddo
            call mpi_2d_fft(c_var,-1)
            c_ff(:,:,:)=c_ff(:,:,:)-coe_sk*c_var(:,:,:)

            ! add v*dudy
            call dev_fft_linear_cmplx(c_u,c_src,1,2) ! c_dudy
            call mpi_2d_fft(c_src,1)
            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1
                
                ! multiply u and dudx in physical space
                c_var(i-ghst,j-ghst,k-ghst)= cmplx( v(i,j,k) * &
                           real(c_src(i-ghst,j-ghst,k-ghst),kind=mytype) , &
                                                   0.d0, kind=mytype )
            
            enddo
            enddo
            enddo
            call mpi_2d_fft(c_var,-1)
            c_ff(:,:,:)=c_ff(:,:,:)-coe_sk*c_var(:,:,:)

            ! add w*dudz
            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1
                
                a1=u(i,j,k+1)*l_t(k)+u(i,j,k)*(1.d0-l_t(k))
                a2=u(i,j,k-1)*l_b(k)+u(i,j,k)*(1.d0-l_b(k))
                tmp1   = ( a1  - a2 ) / dz(k)
                tmp2   = 0.5d0*( w(i,j,k) + w(i,j,k-1) )
                c_var(i-ghst,j-ghst,k-ghst)= cmplx( tmp1*tmp2,0.d0, kind=mytype )
            
            enddo
            enddo
            enddo
            call mpi_2d_fft(c_var,-1)
            c_ff(:,:,:)=c_ff(:,:,:)-coe_sk*c_var(:,:,:)
            
        endif
        !!!!!!!!!!!! finish ff !!!!!!!!!

        ! add ff to qu
        c_qu(:,:,:)=c_qu(:,:,:)+rk_a*dt*c_ff(:,:,:)
        !!!!!!!!!!!!!!! U equation finish !!!!!!

        
        
        !!!!!!!!!! V equation !!!!!!!!!!!
        !!!!!!!!!! now start building the right hand side !!!!!!!!!!
        ! qv=v^n+dt*(rk_a*gg^n+rk_b*gg^n-1) + &
        !                    rk_c*0.5*dt*nu*(dvdx2+dvdy^2+dvdz^2)^n-rk_c*dt*dpdy^n

        ! add v
        c_qv(:,:,:)=c_v(:,:,:)
        
        ! add dvdx2
        call dev_fft_linear_cmplx(c_v,c_var,2,1) 
        c_qv(:,:,:)=c_qv(:,:,:)+rk_c*0.5d0*dt*nu*c_var(:,:,:)
        ! add dvdy2
        call dev_fft_linear_cmplx(c_v,c_var,2,2) 
        c_qv(:,:,:)=c_qv(:,:,:)+rk_c*0.5d0*dt*nu*c_var(:,:,:)
        ! add dvdz2, we store dvdz2 to tmp1 in physical space
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            tmp1 = ( (v(i,j,k+1)-v(i,j,k))/dz_t(k) - &
                     (v(i,j,k)-v(i,j,k-1))/dz_b(k) ) / dz(k)
            c_var(i-ghst,j-ghst,k-ghst)=cmplx(tmp1, 0.d0, kind=mytype)

        enddo
        enddo
        enddo
        call mpi_2d_fft(c_var,-1)
        c_qv(:,:,:)=c_qv(:,:,:)+rk_c*0.5d0*dt*nu*c_var(:,:,:)
        
        ! add dpdy
        call dev_fft_linear_cmplx(c_p,c_var,1,2) 
        c_qv(:,:,:)=c_qv(:,:,:)-rk_c*dt*c_var(:,:,:)
        
        ! add gg0
        c_qv(:,:,:)=c_qv(:,:,:)+rk_b*dt*c_gg(:,:,:)
        
        !!!!!!!!!!!!! calculate gg !!!!!!!!!!
        ! gg=-0.5*(duvdx+dvvdy+dvwdz)^n - &
        !        0.5*(u*dvdx+v*dvdy+w*dvdz)^n-u_mrf*dvdx^n
        c_gg(:,:,:)=cmplx( 0.d0,0.d0,kind=mytype )
        ! add u_mrf*dvdx
        call dev_fft_linear_cmplx(c_v,c_var,1,1) 
        c_gg(:,:,:)=c_gg(:,:,:)-u_mrf*c_var(:,:,:)
        ! add duvdx
        call dev_fft_nonlinear_cmplx(u,v,c_var,1,0)
        c_gg(:,:,:)=c_gg(:,:,:)-coe_sk*c_var(:,:,:)
        ! add dvvdy
        call dev_fft_nonlinear_cmplx(v,v,c_var,2,0)
        c_gg(:,:,:)=c_gg(:,:,:)-coe_sk*c_var(:,:,:)
        ! add dvwdz, we store dvwdz to tmp1 in physical space
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1
            
            a1=v(i,j,k+1)*l_t(k)+v(i,j,k)*(1.d0-l_t(k))
            a2=v(i,j,k-1)*l_b(k)+v(i,j,k)*(1.d0-l_b(k))
            tmp1 = ( a1 * w(i,j,k) - &
                          a2 * w(i,j,k-1) ) / dz(k)
                          
            if (dp_opt .eq. 3) then
                tmp2   = fc*( ug- u(i,j,k) - u_mrf )
            else
                tmp2   = 0.d0
            endif
            
            c_var(i-ghst,j-ghst,k-ghst)= &
                  cmplx( 1.d0/coe_sk*tmp2-tmp1,0.d0, kind=mytype )

        enddo
        enddo
        enddo
        call mpi_2d_fft(c_var,-1)
        c_gg(:,:,:)=c_gg(:,:,:)+coe_sk*c_var(:,:,:)

        
        ! check if the skew-symmetric form is used, 
        ! if yes we need to add the other terms
        if (issk .eq. 1) then
        
            ! add u*dvdx
            call dev_fft_linear_cmplx(c_v,c_src,1,1) ! c_dvdx
            call mpi_2d_fft(c_src,1)
            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1
                
                ! multiply u and dvdx in physical space
                c_var(i-ghst,j-ghst,k-ghst)= cmplx( u(i,j,k) * &
                           real(c_src(i-ghst,j-ghst,k-ghst),kind=mytype) , &
                                                   0.d0, kind=mytype )
            
            enddo
            enddo
            enddo
            call mpi_2d_fft(c_var,-1)
            c_gg(:,:,:)=c_gg(:,:,:)-coe_sk*c_var(:,:,:)

            ! add v*dvdy
            call dev_fft_linear_cmplx(c_v,c_src,1,2) ! c_dvdy
            call mpi_2d_fft(c_src,1)
            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1
                
                ! multiply v and dvdy in physical space
                c_var(i-ghst,j-ghst,k-ghst)= cmplx( v(i,j,k) * &
                           real(c_src(i-ghst,j-ghst,k-ghst),kind=mytype) , &
                                                   0.d0, kind=mytype )
            
            enddo
            enddo
            enddo
            call mpi_2d_fft(c_var,-1)
            c_gg(:,:,:)=c_gg(:,:,:)-coe_sk*c_var(:,:,:)

            ! add w*dvdz
            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1
                
                a1=v(i,j,k+1)*l_t(k)+v(i,j,k)*(1.d0-l_t(k))
                a2=v(i,j,k-1)*l_b(k)+v(i,j,k)*(1.d0-l_b(k))
                tmp1   = ( a1  - a2 ) / dz(k)
                tmp2   = 0.5d0*( w(i,j,k) + w(i,j,k-1) )
                c_var(i-ghst,j-ghst,k-ghst)= cmplx( tmp1*tmp2,0.d0, kind=mytype )
            
            enddo
            enddo
            enddo
            call mpi_2d_fft(c_var,-1)
            c_gg(:,:,:)=c_gg(:,:,:)-coe_sk*c_var(:,:,:)
            
        endif
        !!!!!!!!!!!! finish gg !!!!!!!!!

        ! add gg to qv
        c_qv(:,:,:)=c_qv(:,:,:)+rk_a*dt*c_gg(:,:,:)
        !!!!!!!!!!!!!!! V equation finish !!!!!!

        
        
        !!!!!!!!!! W equation !!!!!!!!!!! ------ note that w is staggered
        !!!!!!!!!! now start building the right hand side !!!!!!!!!!
        ! qw=w^n+dt*(rk_a*hh^n+rk_b*hh^n-1) + &
        !                     rk_c*0.5*dt*nu*(dwdx2+dwdy^2+dwdz^2)^n-rk_c*dt*dpdz^n

        ! add w
        c_qw(:,:,:)=c_w(:,:,:)
        
        ! add dwdx2
        call dev_fft_linear_cmplx(c_w,c_var,2,1) 
        c_qw(:,:,:)=c_qw(:,:,:)+rk_c*0.5d0*dt*nu*c_var(:,:,:)
        ! add dwdy2
        call dev_fft_linear_cmplx(c_w,c_var,2,2) 
        c_qw(:,:,:)=c_qw(:,:,:)+rk_c*0.5d0*dt*nu*c_var(:,:,:)
        ! add dwdz2, we store dwdz2 to tmp1 in physical space
        do k=kstr3,kend3,1  !!!!! we might need a different range here
        do j=jstr3,jend3,1
        do i=istr3,iend3,1
        
            tmp1 =  ( (w(i,j,k+1)-w(i,j,k))/dz(k+1) - &
                      (w(i,j,k)-w(i,j,k-1))/dz(k) ) / dz_t(k)
            
            c_var(i-ghst,j-ghst,k-ghst)=cmplx(tmp1, 0.d0, kind=mytype)

        enddo
        enddo
        enddo
        call mpi_2d_fft(c_var,-1)
        c_qw(:,:,:)=c_qw(:,:,:)+rk_c*0.5d0*dt*nu*c_var(:,:,:)
        
        ! add dpdz
        do k=kstr3,kend3,1  !!!!! we might need a different range here
        do j=jstr3,jend3,1
        do i=istr3,iend3,1
        
            tmp1=  ( p(i,j,k+1)-p(i,j,k) ) / dz_t(k)
            
            c_var(i-ghst,j-ghst,k-ghst)=cmplx(tmp1, 0.d0, kind=mytype)

        enddo
        enddo
        enddo
        call mpi_2d_fft(c_var,-1)
        c_qw(:,:,:)=c_qw(:,:,:)-rk_c*dt*c_var(:,:,:)
        
        ! add hh0
        c_qw(:,:,:)=c_qw(:,:,:)+rk_b*dt*c_hh(:,:,:)
        
        !!!!!!!!!!!!! calculate hh !!!!!!!!!!
        ! hh=-0.5*(duwdx+dvwdy+dwwdz)^n - &
        !        0.5*(u*dwdx+v*dwdy+w*dwdz)^n-u_mrf*dwdx^n + g*(T-T0)/T0
        c_hh(:,:,:)=cmplx( 0.d0,0.d0,kind=mytype )
        ! add u_mrf*dwdx
        call dev_fft_linear_cmplx(c_w,c_var,1,1) 
        c_hh(:,:,:)=c_hh(:,:,:)-u_mrf*c_var(:,:,:)
        ! add duwdx
        call dev_fft_nonlinear_cmplx(u,w,c_var,1,1)
        c_hh(:,:,:)=c_hh(:,:,:)-coe_sk*c_var(:,:,:)
        ! add dvwdy
        call dev_fft_nonlinear_cmplx(v,w,c_var,2,1)
        c_hh(:,:,:)=c_hh(:,:,:)-coe_sk*c_var(:,:,:)
        ! add dwwdz and the buoyancy term, we store dwwdz to tmp1 in physical space
        do k=kstr3,kend3,1   ! range
        do j=jstr3,jend3,1
        do i=istr3,iend3,1
            
            !dwwdz
            a1 = ( 0.5d0*(w(i,j,k)+w(i,j,k+1)) )**2.d0
            a2 = ( 0.5d0*(w(i,j,k)+w(i,j,k-1)) )**2.d0
            tmp1= (a1-a2) / dz_t(k)

            if (isscalar .eq. 1) then
                ! buoyancy
                tmp2=got*( (1.d0-l_t(k))*t(i,j,k)+l_t(k)*t(i,j,k+1) -t_ref )
                ! note that here we use 1/coe_sk*buoyancy
                ! since we will do a coe_sk*c_var later
                c_var(i-ghst,j-ghst,k-ghst) = &
                         cmplx( 1.d0/coe_sk*tmp2 - tmp1,0.d0, kind=mytype )
            else
                c_var(i-ghst,j-ghst,k-ghst)= cmplx( -tmp1,0.d0, kind=mytype )
            endif

        enddo
        enddo
        enddo
        call mpi_2d_fft(c_var,-1)
        c_hh(:,:,:)=c_hh(:,:,:)+coe_sk*c_var(:,:,:)

        ! check if the skew-symmetric form is used, 
        ! if yes we need to add the other terms
        if (issk .eq. 1) then
        
            ! add u*dwdx
            call dev_fft_linear_cmplx(c_w,c_src,1,1) ! c_dwdx
            call mpi_2d_fft(c_src,1)
            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1
                
                ! multiply u and dwdx in physical space
                tmp1= l_t(k)*u(i,j,k+1)+(1.d0-l_t(k))*u(i,j,k)
                c_var(i-ghst,j-ghst,k-ghst)= cmplx( tmp1 * &
                           real(c_src(i-ghst,j-ghst,k-ghst),kind=mytype) , &
                                                   0.d0, kind=mytype )
            
            enddo
            enddo
            enddo
            call mpi_2d_fft(c_var,-1)
            c_hh(:,:,:)=c_hh(:,:,:)-coe_sk*c_var(:,:,:)

            ! add v*dwdy
            call dev_fft_linear_cmplx(c_w,c_src,1,2) ! c_dwdy
            call mpi_2d_fft(c_src,1)
            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1
                
                ! multiply v and dwdy in physical space
                tmp1= l_t(k)*v(i,j,k+1)+(1.d0-l_t(k))*v(i,j,k)
                c_var(i-ghst,j-ghst,k-ghst)= cmplx( tmp1 * &
                           real(c_src(i-ghst,j-ghst,k-ghst),kind=mytype) , &
                                                   0.d0, kind=mytype )
            
            enddo
            enddo
            enddo
            call mpi_2d_fft(c_var,-1)
            c_hh(:,:,:)=c_hh(:,:,:)-coe_sk*c_var(:,:,:)

            ! add w*dwdz
            do k=kstr3,kend3,1
            do j=jstr3,jend3,1
            do i=istr3,iend3,1
            
                a1 =  0.5d0*(w(i,j,k)+w(i,j,k+1)) 
                a2 =  0.5d0*(w(i,j,k)+w(i,j,k-1)) 
                tmp1= (a1-a2) / dz_t(k)
                c_var(i-ghst,j-ghst,k-ghst) = &
                           cmplx( w(i,j,k)*tmp1,0.d0, kind=mytype )
            
            enddo
            enddo
            enddo
            call mpi_2d_fft(c_var,-1)
            c_hh(:,:,:)=c_hh(:,:,:)-coe_sk*c_var(:,:,:)


        endif
        !!!!!!!!!!!! finish hh !!!!!!!!!

        ! add hh to qw
        c_qw(:,:,:)=c_qw(:,:,:)+rk_a*dt*c_hh(:,:,:)
        !!!!!!!!!!!!!!! W equation finish !!!!!!
        

        return
    end subroutine


    ! this function calculates spatial derivatives and the 
    ! right hand side terms for the T equation
    subroutine get_temperature_cnvdiff_spec &
        (u,v,w,t,c_t,c_ss,c_qt)
                                           
        ! these are temporary arrays, their dimensions are: nx/p_row,ny/p_col,nz
        ! they are on z pencil
        use module_variables, only : c_var
        use module_tools
        
        implicit none

        real(mytype),dimension(:,:,:),intent(in)      :: u,v,w,t
        complex(mytype),dimension(:,:,:),intent(in)   :: c_t
        complex(mytype),dimension(:,:,:),intent(out)  :: c_ss
        complex(mytype),dimension(:,:,:),intent(out)  :: c_qt
        
        integer                :: i,j,k
        real(mytype)           :: tmp1,a1,a2


        !!!!!!!!!! T equation !!!!!!!!!!!
        !!!!!!!!!! now start building the right hand side !!!!!!!!!!
        ! qt=T^n+dt*(rk_a*ss^n+rk_b*ss^n-1) + &
        !                       rk_c*0.5*dt*alpha*(dTdx2+dTdy^2+dTdz^2)^n

        ! add T
        c_qt(:,:,:)=c_t(:,:,:)
        
        ! add dTdx2
        call dev_fft_linear_cmplx(c_t,c_var,2,1) 
        c_qt(:,:,:)=c_qt(:,:,:)+rk_c*0.5d0*dt*alpha*c_var(:,:,:)
        ! add dTdy2
        call dev_fft_linear_cmplx(c_t,c_var,2,2) 
        c_qt(:,:,:)=c_qt(:,:,:)+rk_c*0.5d0*dt*alpha*c_var(:,:,:)
        ! add dTdz2, we store dTdz2 to tmp1 in physical space
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            tmp1=  ( (t(i,j,k+1)-t(i,j,k))/dz_t(k) - &
                     (t(i,j,k)-t(i,j,k-1))/dz_b(k) ) / dz(k)
            c_var(i-ghst,j-ghst,k-ghst)=cmplx(tmp1, 0.d0, kind=mytype)

        enddo
        enddo
        enddo
        call mpi_2d_fft(c_var,-1)
        c_qt(:,:,:)=c_qt(:,:,:)+rk_c*0.5d0*dt*alpha*c_var(:,:,:)
 
        ! add ss^n-1
        c_qt(:,:,:)=c_qt(:,:,:)+rk_b*dt*c_ss(:,:,:)
        
        !!!!!!!!!!!!! calculate ss !!!!!!!!!!
        ! ss= - (dTudx+dTvdy+dTwdz)^n - u_mrf*dTdx
        c_ss(:,:,:)=cmplx( 0.d0,0.d0,kind=mytype )
        ! add u_mrf*dTdx
        call dev_fft_linear_cmplx(c_t,c_var,1,1) 
        c_ss(:,:,:)=c_ss(:,:,:)-u_mrf*c_var(:,:,:)
        ! add dTudx
        call dev_fft_nonlinear_cmplx(u,t,c_var,1,0)
        c_ss(:,:,:)=c_ss(:,:,:)-c_var(:,:,:)
        ! add dTvdy
        call dev_fft_nonlinear_cmplx(t,v,c_var,2,0)
        c_ss(:,:,:)=c_ss(:,:,:)-c_var(:,:,:)
        ! add dTwdz, we store dTwdz to tmp1 in physical space
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            a1=t(i,j,k+1)*l_t(k)+t(i,j,k)*(1.d0-l_t(k))
            a2=t(i,j,k-1)*l_b(k)+t(i,j,k)*(1.d0-l_b(k))
            tmp1   = ( a1 * w(i,j,k) - &
                       a2 * w(i,j,k-1) ) / dz(k)
            c_var(i-ghst,j-ghst,k-ghst)= cmplx(tmp1,0.d0, kind=mytype )

        enddo
        enddo
        enddo
        call mpi_2d_fft(c_var,-1)
        c_ss(:,:,:)=c_ss(:,:,:)-c_var(:,:,:)

        ! add ss to qt
        c_qt(:,:,:)=c_qt(:,:,:)+rk_a*dt*c_ss(:,:,:)
        !!!!!!!!!!!!!!! U equation finish !!!!!!


        return
    end subroutine



    ! Poisson equation fft solver for spectral method,
    ! we call it Poisson equation solver,
    ! although we also use it to solve Helmholtz equation for U,V,W,T
    subroutine poisson_solver_fft_spec &
        (aa,bb,cc,c_src_in,c_var_out,zstagger,delnull,tbnd,ubnd)

        ! C standard, for using fftw
        use, intrinsic :: iso_c_binding
        use module_tools
        ! these are temporary arrays, their dimensions are: nx/p_row,ny/p_col,nz
        ! they are on z pencil
        use module_variables, only : c_src

        implicit none

        ! for using fftw
        include 'include/fftw3.f03'

        ! these are Thomas algorithm coefficients
        real(mytype),   dimension(:,:,:),intent(inout) :: aa,cc
        complex(mytype),dimension(:,:,:),intent(inout) :: bb
        ! 3d real array, input from source term
        complex(mytype),   dimension(:,:,:),   intent(in) :: c_src_in
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
        complex(mytype),dimension(:,:,:),     intent(out) :: c_var_out

        ! temporary variable used in Thomas algorithm
        complex(mytype)                                :: tmp
        ! temporary variable used in T and u eqn boundary condition
        real(mytype)                                   :: aa_t_zstr,cc_t_zend
        real(mytype)                                   :: aa_u_zstr,cc_u_zend
        integer   :: i,j,k,kmax3

        ! initiate variables
        c_var_out(:,:,:)=0.d0
        c_src(:,:,:)=c_src_in(:,:,:)

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

        ! Thomas algorithm to calculate c_var_out
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

            c_var_out(i,j,kmax3-ghst)=c_src(i,j,kmax3-ghst)/bb(i,j,kmax3)

        enddo
        enddo
        
        do k=kmax3-1,kstr3,-1
        do j=1,zsize(2),1
        do i=1,zsize(1),1

            c_var_out(i,j,k-ghst)=( c_src(i,j,k-ghst)- &
                    cc(i,j,k)*c_var_out(i,j,k-ghst+1) ) / bb(i,j,k)

        enddo
        enddo
        enddo


        return

    end subroutine



    
    ! this function calculates the divergence
    subroutine get_div_spec(c_u,c_v,w,c_div)
    
        ! these are temporary arrays, their dimensions are: nx/p_row,ny/p_col,nz
        ! they are on z pencil
        use module_variables, only : c_var
        use mpi
        use module_tools

        implicit none

        complex(mytype),dimension(:,:,:),intent(in)    :: c_u,c_v
        real(mytype),dimension(:,:,:),intent(in)       :: w
        complex(mytype),dimension(:,:,:),intent(out)   :: c_div

        integer              :: i,j,k
        real(mytype)         :: dwdz

        
        ! add dudx
        call dev_fft_linear_cmplx(c_u,c_var,1,1) 
        c_div(:,:,:)=c_var(:,:,:)

        ! add dvdy
        call dev_fft_linear_cmplx(c_v,c_var,1,2) 
        c_div(:,:,:)=c_div(:,:,:)+c_var(:,:,:)

        ! add dwdz
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            dwdz = ( w(i,j,k)-w(i,j,k-1) )/dz(k)
            c_var(i-ghst,j-ghst,k-ghst)= cmplx( dwdz,0.d0, kind=mytype )
            
        enddo
        enddo
        enddo
        call mpi_2d_fft(c_var,-1)
        c_div(:,:,:)=c_div(:,:,:)+c_var(:,:,:)
        

        return

    end subroutine
    
    
    
    
    ! calculate u^(n+1) using phi
    subroutine calculate_uvw_spec(phi,c_phi,c_u,c_v,c_w)
    
        ! temporary complex array to store derivatives
        use module_variables, only : c_var
        use module_tools
        
        implicit none
        
        real(mytype),dimension(:,:,:),intent(in)        :: phi
        complex(mytype),dimension(:,:,:),intent(in)     :: c_phi
        complex(mytype),dimension(:,:,:),intent(inout)  :: c_u,c_v,c_w

        integer                               :: i,j,k
        real(mytype)                          :: dphidz
        
        ! calculate new u
        call dev_fft_linear_cmplx(c_phi,c_var,1,1) 
        c_u(:,:,:)=c_u(:,:,:) - rk_c*dt*c_var(:,:,:)

        ! calculate new v
        call dev_fft_linear_cmplx(c_phi,c_var,1,2) 
        c_v(:,:,:)=c_v(:,:,:) - rk_c*dt*c_var(:,:,:)

        ! calculate new w
        do k=kstr3,kend3-1,1 ! note that we don't need to update w at nz
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            dphidz=  ( phi(i,j,k+1) - phi(i,j,k) ) / dz_t(k)
            c_var(i-ghst,j-ghst,k-ghst)= cmplx( dphidz,0.d0, kind=mytype )
            
        enddo
        enddo
        enddo
        call mpi_2d_fft(c_var,-1)

        ! we don't update w(nz), we will assign it to zeros when applying 
        ! boundary conidition
        do k=1,zsize(3)-1,1  ! notice this range!
        do j=1,zsize(2),1
        do i=1,zsize(1),1
            c_w(i,j,k)=c_w(i,j,k) - rk_c*dt*c_var(i,j,k)
        enddo
        enddo
        enddo

        return

    end subroutine
    
    
    
    ! calculate p^(n+1) using p^n and phi
    subroutine calculate_new_p_spec(phi,c_phi,c_p)
    
        ! temporary complex array to store derivatives
        use module_variables, only : c_var
        use module_tools
        
        implicit none
        
        real(mytype),dimension(:,:,:),intent(in)        :: phi
        complex(mytype),dimension(:,:,:),intent(in)     :: c_phi
        complex(mytype),dimension(:,:,:),intent(inout)  :: c_p

        integer                               :: i,j,k
        real(mytype)                          :: dphidz2
        
        ! add phi
        c_p(:,:,:)=c_p(:,:,:) + c_phi(:,:,:)
        ! add dphidx2
        call dev_fft_linear_cmplx(c_phi,c_var,2,1) 
        c_p(:,:,:)=c_p(:,:,:) - rk_c*0.5d0*dt*nu*c_var(:,:,:)
        ! add dphidy2
        call dev_fft_linear_cmplx(c_phi,c_var,2,2) 
        c_p(:,:,:)=c_p(:,:,:) - rk_c*0.5d0*dt*nu*c_var(:,:,:)
        ! add dphidz2
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            dphidz2=  ( (phi(i,j,k+1)-phi(i,j,k))/dz_t(k) -  &
                        (phi(i,j,k)-phi(i,j,k-1))/dz_b(k) ) / dz(k)
            c_var(i-ghst,j-ghst,k-ghst)= cmplx( dphidz2,0.d0, kind=mytype )
            
        enddo
        enddo
        enddo
        call mpi_2d_fft(c_var,-1)
        c_p(:,:,:)=c_p(:,:,:) - rk_c*0.5d0*dt*nu*c_var(:,:,:)

        return

    end subroutine
    
    

    ! print the div info to screen 
    subroutine screen_div_error_spec(u,v,w,c_u,c_v,isfinal,ishistout)

        use mpi
        use module_tools
        ! temporary complex array to store derivatives
        use module_variables, only :  c_src

        implicit none

        real(mytype),dimension(:,:,:),intent(in)     :: u,v,w
        complex(mytype),dimension(:,:,:),intent(in)  :: c_u,c_v
        
        integer,intent(in)                        :: isfinal,ishistout

        integer                           :: i,j,k
        real(mytype)                      :: div_avg, div_max
        real(mytype)                      :: div_avg_all,div_max_all
        real(mytype)                      :: div_hist_max_all,tmp1
        
        call get_div_spec(c_u,c_v,w,c_src)

        call mpi_2d_fft(c_src,1)

        ! calculate continuity error
        div_avg=0.d0
        div_max=0.d0
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            tmp1=real(c_src(i-ghst,j-ghst,k-ghst),kind=mytype)

            div_avg=div_avg+abs(tmp1)

            if ( abs(tmp1) .gt. div_max ) then
                div_max=abs(tmp1)
            endif

            if (    u(i,j,k) .ne. u(i,j,k) .or. &
                    v(i,j,k) .ne. v(i,j,k) .or. &
                    w(i,j,k) .ne. w(i,j,k) ) then
                print *, 'NaN error! ',i,' ',j,' ',k
                stop
            endif

            if (isfinal .eq. 1) then
                if ( abs(tmp1) .gt. div_hist_max ) then
                    div_hist_max=abs(tmp1)
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
   
    
    
    ! calculate mean fields
    subroutine calculate_mean_fields_spec &
    (u,v,w,t,mean_u,mean_v,mean_w,mean_t,mean_uu,mean_vv,mean_ww, &
     mean_uw,mean_vw,mean_tt,mean_tw)

        implicit none

        real(mytype),dimension(:,:,:),intent(in)    :: u,v,w,t
        real(mytype),dimension(:,:,:),intent(out)   :: mean_u,mean_v,mean_w,mean_t, &
         mean_uu,mean_vv,mean_ww,mean_uw, mean_vw, mean_tt,mean_tw

        integer                           :: i,j,k

        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            ! note that here we need to add the speed of moving reference frame
            mean_uu(i,j,k)=mean_uu(i,j,k)+(u(i,j,k)+u_mrf)**2.d0
            mean_vv(i,j,k)=mean_vv(i,j,k)+v(i,j,k)**2.d0
            mean_ww(i,j,k)=mean_ww(i,j,k)+w(i,j,k)**2.d0
            mean_uw(i,j,k)=mean_uw(i,j,k)+0.5d0*( w(i,j,k)+w(i,j,k-1) ) * &
                                                ( u(i,j,k)+u_mrf )
            mean_vw(i,j,k)=mean_vw(i,j,k)+0.5d0*( w(i,j,k)+w(i,j,k-1) ) * v(i,j,k) 
                                                 
            mean_tt(i,j,k)=mean_tt(i,j,k)+t(i,j,k)**2.d0
            mean_tw(i,j,k)=mean_tw(i,j,k)+t(i,j,k) * 0.5d0 * ( w(i,j,k)+w(i,j,k-1) )

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


    
    ! this function calculates derivatives for non-linear terms using spectral method
    ! dir=1: derivative in x direction
    ! dir=2: derivative in y direction
    ! weqn=1: for w equation. since w is staggered, we need special treatment
    ! if weqn=1, make sure to put w as r_in2
    subroutine dev_fft_nonlinear_cmplx(r_in1,r_in2,c_out3,dir,weqn)
    
        use module_tools

        implicit none

        real(mytype),dimension(:,:,:),intent(in)     :: r_in1,r_in2
        complex(mytype),dimension(:,:,:),intent(out) :: c_out3
        integer,intent(in)                        :: dir,weqn
        
        integer                      :: i,j,k, l, m
        real(mytype)     :: tmp1

        ! multiply these truncated variables
        do k=1,zsize(3),1
        do j=1,zsize(2),1
        do i=1,zsize(1),1

            if (weqn .eq. 1) then ! we need to interpolate u or v to w level

                tmp1=r_in1(i+ghst,j+ghst,k+1+ghst)*l_t(k+ghst) + &
                     r_in1(i+ghst,j+ghst,k+ghst)*(1.d0-l_t(k+ghst))

                c_out3(i,j,k)=cmplx( tmp1 * &
                                 r_in2(i+ghst,j+ghst,k+ghst),0.d0,kind=mytype )

            else

                c_out3(i,j,k)=cmplx( r_in1(i+ghst,j+ghst,k+ghst) * &
                                 r_in2(i+ghst,j+ghst,k+ghst),0.d0,kind=mytype )

            endif
            
            
 
        enddo
        enddo
        enddo

        ! mpi 2d fft for c_out3
        call mpi_2d_fft(c_out3,-1)
        
        ! calculate derivative of the non-linear term
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
            
            if (dir .eq. 1) then ! d?dx
            
                c_out3(i,j,k) =  c_out3(i,j,k)*iii*2.d0*pi/lx*l
                
            elseif (dir .eq. 2) then  ! d?dy
            
                c_out3(i,j,k) =  c_out3(i,j,k)*iii*2.d0*pi/ly*m
            
            else
            
                print *, 'flag error'
                stop
            
            endif

            ! remove Nyquist wavenumber
            if (l .eq. -nx/2 .or. m .eq. -ny/2) then
                c_out3(i,j,k)=cmplx(0.d0,0.d0,kind=mytype)
            endif

            
        enddo
        enddo
        enddo
        

        return

    end subroutine
    
    
    


    ! do fft and calculate Fourier coefficients for a variable
    subroutine real_2_cmplx(r_in,c_out)

        use module_tools

        implicit none

        real(mytype),dimension(:,:,:),intent(in)     :: r_in
        complex(mytype),dimension(:,:,:),intent(out) :: c_out
        
        integer                      :: i,j,k
        
        do k=1,zsize(3),1
        do j=1,zsize(2),1
        do i=1,zsize(1),1
            
            c_out(i,j,k)=cmplx(r_in(i+ghst,j+ghst,k+ghst),0.d0,kind=mytype)
            
        enddo
        enddo
        enddo

        call mpi_2d_fft(c_out,-1)

        return

    end subroutine




    ! do ifft and calculate a variable based on the Fourier coefficients
    subroutine cmplx_2_real(c_in,r_out)

        ! temporary complex array, their dimension is: nx/p_row,ny/p_col,nz
        use module_variables, only : c_var
        use module_tools
    
        implicit none

        real(mytype),dimension(:,:,:),intent(out)   :: r_out
        complex(mytype),dimension(:,:,:),intent(in) :: c_in
        
        integer                      :: i,j,k

        c_var(:,:,:)=c_in(:,:,:)

        call mpi_2d_fft(c_var,1)

        do k=1,zsize(3),1
        do j=1,zsize(2),1
        do i=1,zsize(1),1
            
            r_out(i+ghst,j+ghst,k+ghst)=real(c_var(i,j,k),kind=mytype)
            
        enddo
        enddo
        enddo

        return

    end subroutine



    ! truncate the Fourier coefficients using the 2/3 rule
    subroutine spectral_truncation(c_inout3)

        implicit none

        complex(mytype),dimension(:,:,:),intent(inout) :: c_inout3
        
        integer                      :: i,j,k,l,m

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
            
            if (l .ge. int(nx/3.d0) .or. l .le. -int(nx/3.d0) ) then
            
                c_inout3(i,j,k) =  cmplx(0.d0,0.d0,kind=mytype)

            endif
                
            
            if (m .ge. int(ny/3.d0) .or. m .le. -int(ny/3.d0) ) then
            
                c_inout3(i,j,k) =  cmplx(0.d0,0.d0,kind=mytype)

            endif

            
        enddo
        enddo
        enddo

        return

    end subroutine



    ! this function calculate derivatives for linear terms using spectral method
    ! ord: 1st or 2nd order derivatives; dir: direction? x or y?
    subroutine dev_fft_linear_cmplx(c_in3,c_out3,ord,dir)
    
        implicit none

        complex(mytype),dimension(:,:,:),intent(in)    :: c_in3
        complex(mytype),dimension(:,:,:),intent(out)   :: c_out3
        integer,intent(in)                             :: ord,dir
        
        integer                      :: i,j,k,l,m

        c_out3(:,:,:)=c_in3(:,:,:)

        ! calculate derivative
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
            
            if (ord .eq. 1 .and. dir .eq. 1) then ! d?dx
            
                c_out3(i,j,k) =  c_out3(i,j,k)*iii*2.d0*pi/lx*l
                
            elseif (ord .eq. 1 .and. dir .eq. 2) then  ! d?dy
            
                c_out3(i,j,k) =  c_out3(i,j,k)*iii*2.d0*pi/ly*m
            
            elseif (ord .eq. 2 .and. dir .eq. 1) then  ! d?dx2
            
                c_out3(i,j,k) = -c_out3(i,j,k)*(2.d0*pi/lx*l)**2.d0
            
            elseif (ord .eq. 2 .and. dir .eq. 2) then  ! d?dy2
            
                c_out3(i,j,k) = -c_out3(i,j,k)*(2.d0*pi/ly*m)**2.d0
            
            else
            
                print *, 'flag error'
                stop
            
            endif

            ! remove Nyquist wavenumber
            if (l .eq. -nx/2 .or. m .eq. -ny/2) then
                c_out3(i,j,k)=cmplx(0.d0,0.d0,kind=mytype)
            endif
            
        enddo
        enddo
        enddo
    

        return

    end subroutine
    
    
    
    ! Rayleigh damping layer at the top of the boundary layer
    subroutine rayleigh_damping_spec(c_varin)

        use mpi

        implicit none

        complex(mytype),dimension(:,:,:),intent(inout) :: c_varin

        real(mytype)                      :: zd,sigz
        complex(mytype)                   :: mean_local,mean_all
        integer                           :: i,j,k


        ! do the damping for the topmost nzdamp points
        do k=kend3-nzdamp,kend3,1
            
            ! first calculate the mean of the variable on the local processor
            mean_local=0.d0
            do j=jstr3,jend3,1
            do i=istr3,iend3,1
                mean_local=mean_local+c_varin(i-ghst,j-ghst,k-ghst)
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
                c_varin(i-ghst,j-ghst,k-ghst) = c_varin(i-ghst,j-ghst,k-ghst) - &
                                    sigz*( c_varin(i-ghst,j-ghst,k-ghst)-mean_all )
            enddo
            enddo

        enddo


        return

    end subroutine


end module
