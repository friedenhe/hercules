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



module module_cnvdiff

    use module_parameters

    implicit none

contains



    ! this function calculates the spatial derivatives and 
    ! the right hand side terms for the momentum equations
    ! i.e., f,g,h,qu,qv,qw etc
    subroutine get_momentum_cnvdiff &
    (u,v,w,p,t,uhx,uhy,uhz,vhx,vhy,vhz,whx,why,ff,ff0,gg,gg0,hh,hh0,qu,qv,qw)
                                           
        implicit none

        real(mytype),dimension(:,:,:),intent(in) ::u,v,w,p,t,ff0,gg0, &
                                           hh0,uhx,uhy,uhz,vhx,vhy,vhz,whx,why
        real(mytype),dimension(:,:,:),intent(out)::ff,gg,hh,qu,qv,qw

        integer              :: i,j,k
        real(mytype)         :: a1,a2
        real(mytype)         :: uux, uvy, uwz, uxx, uyy, uzz, dpdx, ux
        real(mytype)         :: vux, vvy, vwz, vxx, vyy, vzz, dpdy, vx
        real(mytype)         :: wux, wvy, wwz, wxx, wyy, wzz, dpdz, wx
        
        ! first calculate dpx_drive, it is value is dependent on dp_opt
        call update_dpx_drive

        
        !!!!!!!!!! U equation !!!!!!!!!!!
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1
        
            ! horizontal derivatives, check to use 2nd or 4th order schemes
            if ( cds .eq. 1) then
        
                uux= ( uhx(i+1,j,k)**2.d0 - &
                       uhx(i,j,k)**2.d0 ) / dx
        
                uvy= ( uhy(i,j,k)*vhx(i,j,k) - &
                       uhy(i,j-1,k)*vhx(i,j-1,k) ) / dy
        
                ux = ( uhx(i+1,j,k)-uhx(i,j,k) ) / dx
        
            elseif ( cds .eq. 2) then
        
                uux= ( -uhx(i+2,j,k)**2.d0 + &
                  27.d0*uhx(i+1,j,k)**2.d0 - &
                  27.d0*uhx(i,j,k)**2.d0   + &
                        uhx(i-1,j,k)**2.d0 ) / 24.d0 / dx
        
                uvy= ( -uhy(i,j+1,k)*vhx(i,j+1,k) + &
                  27.d0*uhy(i,j,k)  *vhx(i,j,k)   - &
                  27.d0*uhy(i,j-1,k)*vhx(i,j-1,k) + &
                        uhy(i,j-2,k)*vhx(i,j-2,k) ) / 24.d0 /dy
        
                ux = ( -uhx(i+2,j,k) + &
                  27.d0*uhx(i+1,j,k) - &
                  27.d0*uhx(i,j,k)   + &
                        uhx(i-1,j,k) ) / 24.d0 / dx
     
            endif
        
            ! we always use the 2nd order scheme for vertical derivatives
            uwz= ( uhz(i,j,k)  *whx(i,j,k) - &
                   uhz(i,j,k-1)*whx(i,j,k-1) ) / dz(k)

        
            ! horizontal derivatives, check to use 2nd or 4th schemes
            if ( cds .eq. 1) then
        
                uxx= ( u(i+1,j,k) - &
                  2.d0*u(i,j,k)   + &
                       u(i-1,j,k) ) / dx2
        
                uyy= ( u(i,j+1,k) - &
                  2.d0*u(i,j,k)   + &
                       u(i,j-1,k) ) / dy2
        
            elseif ( cds .eq. 2) then
        
                uxx= ( -u(i+2,j,k) + &
                  16.d0*u(i+1,j,k) - &
                  30.d0*u(i,j,k)   + &
                  16.d0*u(i-1,j,k) - &
                        u(i-2,j,k) ) / 12.d0 / dx2
        
                uyy= ( -u(i,j+2,k) + &
                  16.d0*u(i,j+1,k) - &
                  30.d0*u(i,j,k)   + &
                  16.d0*u(i,j-1,k) - &
                        u(i,j-2,k) ) / 12.d0 / dy2
                        
            endif
        
            ! we always use the 2nd order scheme for vertical derivatives
            uzz= ( (u(i,j,k+1)-u(i,j,k))/dz_t(k) - &
                   (u(i,j,k)-u(i,j,k-1))/dz_b(k) ) / dz(k)
        
            ! calculate the nonlinear terms for U
            ! check if it is Ekman layer, if yes, add the Coriolis term
            if (dp_opt .eq. 3) then

                ! constant geostrophic wind
                ff(i,j,k)  = -uux-uvy-uwz - u_mrf*ux + &
                    fc*( 0.5d0*(vhx(i,j,k)+vhx(i,j-1,k)) - vg )

            else

                ! + dpx_drive to drive the flow for dp_opt=1 or 2
                ff(i,j,k)  = -uux-uvy-uwz + dpx_drive - u_mrf*ux

            endif
        
            ! dpdx, check to use 2nd or 4th schemes
            if ( cds .eq. 1) then
        
                dpdx= ( p(i+1,j,k)-p(i,j,k) ) / dx
        
            elseif ( cds .eq. 2) then
        
                dpdx= ( -p(i+2,j,k) + &
                   27.d0*p(i+1,j,k) - &
                   27.d0*p(i,j,k)   + &
                         p(i-1,j,k) ) / 24.d0 / dx
                         
            endif

            ! the right hand side term for U eqn
            qu(i,j,k)= u(i,j,k) + dt*( rk_a*ff(i,j,k)+rk_b*ff0(i,j,k) ) + &
                               rk_c*0.5d0*dt*nu*(uxx+uyy+uzz) - rk_c*dt*dpdx
        
        
        enddo
        enddo
        enddo
        


        !!!!!!!!!! V equation !!!!!!!!!!!
        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1
        
            ! horizontal derivatives, check to use 2nd or 4th schemes
            if ( cds .eq. 1) then
        
                vux= ( vhx(i,j,k)  *uhy(i,j,k) - &
                       vhx(i-1,j,k)*uhy(i-1,j,k) ) / dx
        
                vvy= ( vhy(i,j+1,k)**2.d0 - &
                       vhy(i,j,k)**2.d0 ) / dy
        
                vx = ( vhx(i,j,k)-vhx(i-1,j,k) ) / dx
        
            elseif ( cds .eq. 2) then
        
                vux= ( -vhx(i+1,j,k)*uhy(i+1,j,k) + &
                  27.d0*vhx(i,j,k)  *uhy(i,j,k)   - &
                  27.d0*vhx(i-1,j,k)*uhy(i-1,j,k) +  &
                        vhx(i-2,j,k)*uhy(i-2,j,k) ) / 24.d0 / dx
        
                vvy= ( -vhy(i,j+2,k)**2.d0 + &
                  27.d0*vhy(i,j+1,k)**2.d0 - &
                  27.d0*vhy(i,j,k)**2.d0   + &
                        vhy(i,j-1,k)**2.d0 ) / 24.d0 / dy
        
                vx = ( -vhx(i+1,j,k) + &
                  27.d0*vhx(i,j,k)   - &
                  27.d0*vhx(i-1,j,k) +  &
                        vhx(i-2,j,k) ) / 24.d0 / dx
                        
            endif
        
            ! we always use the 2nd order scheme for vertical derivatives
            vwz= ( vhz(i,j,k)  *why(i,j,k) - &
                   vhz(i,j,k-1)*why(i,j,k-1) ) / dz(k)
        
            ! horizontal derivatives, check to use 2nd or 4th schemes
            if ( cds .eq. 1) then
        
                vxx= ( v(i+1,j,k) - &
                  2.d0*v(i,j,k)   + &
                       v(i-1,j,k) ) / dx2
        
                vyy= ( v(i,j+1,k) - &
                  2.d0*v(i,j,k)   + &
                       v(i,j-1,k) ) / dy2
        
            elseif ( cds .eq. 2) then
        
                vxx= ( -v(i+2,j,k) + &
                  16.d0*v(i+1,j,k) - &
                  30.d0*v(i,j,k)   + &
                  16.d0*v(i-1,j,k) - &
                        v(i-2,j,k) ) / 12.d0 / dx2
        
                vyy= ( -v(i,j+2,k) + &
                  16.d0*v(i,j+1,k) - &
                  30.d0*v(i,j,k)   + &
                  16.d0*v(i,j-1,k) - &
                        v(i,j-2,k) ) / 12.d0 / dy2
                        
            endif
        
            ! we always use the 2nd order scheme for vertical derivatives
            vzz= ( (v(i,j,k+1)-v(i,j,k))/dz_t(k) - &
                   (v(i,j,k)-v(i,j,k-1))/dz_b(k) ) / dz(k)
        
            ! calculate the nonlinear term for V
            ! check if it is Ekman layer
            if (dp_opt .eq. 3) then

                ! constant geostrophci wind
                gg(i,j,k) = -vux-vvy-vwz - u_mrf*vx + &
                 fc*( ug - 0.5d0*(uhy(i,j,k)+uhy(i-1,j,k)) - u_mrf )

            else

                gg(i,j,k) = -vux-vvy-vwz - u_mrf*vx

            endif
        
            ! dpdy, check to use 2nd or 4th schemes
            if ( cds .eq. 1) then
        
                dpdy = ( p(i,j+1,k)-p(i,j,k) ) / dy
        
            elseif ( cds .eq. 2) then
        
                dpdy= ( -p(i,j+2,k) + &
                   27.d0*p(i,j+1,k) - &
                   27.d0*p(i,j,k)   + &
                         p(i,j-1,k) ) / 24.d0 / dy

            endif
        
            ! the right hand side term for V eqn
            qv(i,j,k)=v(i,j,k) + dt*( rk_a*gg(i,j,k)+rk_b*gg0(i,j,k) ) +  &
                             rk_c*0.5d0*dt*nu*(vxx+vyy+vzz) - rk_c*dt*dpdy
        
        
        enddo
        enddo
        enddo

        
        !!!!!!!!!!!!!!!!!!!!!!!! W equation  !!!!!!!!!!!!!!!!!
        do k=kstr3,kend3-1,1 ! note that here the z range is different
        do j=jstr3,jend3,1
        do i=istr3,iend3,1
        
            ! horizontal derivatives, check to use 2nd or 4th schemes
            if ( cds .eq. 1) then
        
                wux= ( whx(i,j,k)  *uhz(i,j,k) - &
                       whx(i-1,j,k)*uhz(i-1,j,k) ) / dx
        
                wvy= ( why(i,j,k)  *vhz(i,j,k) - &
                       why(i,j-1,k)*vhz(i,j-1,k) ) / dy
        
                wx = ( whx(i,j,k)-whx(i-1,j,k) ) / dx
        
            elseif ( cds .eq. 2) then
        
                wux= ( -whx(i+1,j,k)*uhz(i+1,j,k) + &
                  27.d0*whx(i,j,k)  *uhz(i,j,k)   - &
                  27.d0*whx(i-1,j,k)*uhz(i-1,j,k) + &
                        whx(i-2,j,k)*uhz(i-2,j,k) ) / 24.d0 /dx
        
                wvy= ( -why(i,j+1,k)*vhz(i,j+1,k) + &
                  27.d0*why(i,j,k)  *vhz(i,j,k)   - &
                  27.d0*why(i,j-1,k)*vhz(i,j-1,k) + &
                        why(i,j-2,k)*vhz(i,j-2,k) ) / 24.d0 /dy
        
                wx = ( -whx(i+1,j,k) + &
                  27.d0*whx(i,j,k)   - &
                  27.d0*whx(i-1,j,k) + &
                        whx(i-2,j,k) ) / 24.d0 /dx
                        
            endif
        
            ! we always use the 2nd order scheme for vertical derivatives
            a1 = ( 0.5d0*(w(i,j,k)+w(i,j,k+1)) )**2.d0
            a2 = ( 0.5d0*(w(i,j,k)+w(i,j,k-1)) )**2.d0
            wwz= (a1-a2) / dz_t(k)
        
            ! horizontal derivatives, check to use 2nd or 4th schemes
            if ( cds .eq. 1) then
        
                wxx= ( w(i+1,j,k) - &
                  2.d0*w(i,j,k)   + &
                       w(i-1,j,k) ) / dx2
        
                wyy= ( w(i,j+1,k) - &
                  2.d0*w(i,j,k)   + &
                       w(i,j-1,k) ) / dy2
        
            elseif ( cds .eq. 2) then
        
                wxx= ( -w(i+2,j,k) + &
                  16.d0*w(i+1,j,k) - &
                  30.d0*w(i,j,k)   + &
                  16.d0*w(i-1,j,k) - &
                        w(i-2,j,k) ) / 12.d0 / dx2
        
                wyy= ( -w(i,j+2,k) + &
                  16.d0*w(i,j+1,k) - &
                  30.d0*w(i,j,k)   + &
                  16.d0*w(i,j-1,k) - &
                        w(i,j-2,k) ) / 12.d0 / dy2
                        
            endif
        
            ! we always use the 2nd order scheme for vertical derivatives
            wzz= ( (w(i,j,k+1)-w(i,j,k))/dz(k+1) - &
                   (w(i,j,k)-w(i,j,k-1))/dz(k) ) / dz_t(k)
        
            ! calculate the nonlinear term for W
            ! check whether to add the buoyancy term
            if (isscalar .eq. 1) then
        
                hh(i,j,k)  = -wux-wvy-wwz -u_mrf*wx + &
                 got*( (1.d0-l_t(k))*t(i,j,k)+l_t(k)*t(i,j,k+1) -t_ref )
        
            elseif (isscalar .eq. 0) then
        
                hh(i,j,k)  = -wux-wvy-wwz -u_mrf*wx
        
            else
                print *, 'isscalar error!'
                stop
            endif
        
            ! dpdz term
            dpdz = ( p(i,j,k+1)-p(i,j,k) ) / dz_t(k)
        
            ! the right hand side term for W eqn
            qw(i,j,k)=w(i,j,k) + dt*( rk_a*hh(i,j,k)+rk_b*hh0(i,j,k) ) + &
                                rk_c*0.5d0*dt*nu*(wxx+wyy+wzz) - rk_c*dt*dpdz
        
        
        enddo
        enddo
        enddo
        

        return
    end subroutine



    ! this function calculates the spatial derivatives and 
    ! the right hand side terms for the temperature equation
    ! i.e., ss, qt etc
    subroutine get_temperature_cnvdiff(u,v,w,t,thx,thy,thz,ss,ss0,qt)

        implicit none

        real(mytype),dimension(:,:,:),intent(in)  :: u,v,w,t,ss0,thx,thy,thz
        real(mytype),dimension(:,:,:),intent(out) :: ss,qt

        integer              :: i,j,k
        real(mytype)         :: tux, tvy, twz, txx, tyy, tzz, tx

        do k=kstr3,kend3,1
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            ! horizontal derivatives, check to use 2nd or 4th schemes
            if ( cds .eq. 1) then

                tux= ( thx(i,j,k)  *u(i,j,k) - &
                       thx(i-1,j,k)*u(i-1,j,k) ) / dx

                tvy= ( thy(i,j,k)  *v(i,j,k) - &
                       thy(i,j-1,k)*v(i,j-1,k) ) / dy

                tx = ( thx(i,j,k)-thx(i-1,j,k) ) / dx

            elseif ( cds .eq. 2) then

                tux= ( -thx(i+1,j,k)*u(i+1,j,k) + &
                  27.d0*thx(i,j,k)  *u(i,j,k)   - &
                  27.d0*thx(i-1,j,k)*u(i-1,j,k) + &
                        thx(i-2,j,k)*u(i-2,j,k) ) / 24.d0 / dx

                tvy= ( -thy(i,j+1,k)*v(i,j+1,k) + &
                  27.d0*thy(i,j,k)  *v(i,j,k)   - &
                  27.d0*thy(i,j-1,k)*v(i,j-1,k) + &
                        thy(i,j-2,k)*v(i,j-2,k) ) / 24.d0 / dy

                tx = ( -thx(i+1,j,k) + &
                  27.d0*thx(i,j,k)   - &
                  27.d0*thx(i-1,j,k) + &
                        thx(i-2,j,k) ) / 24.d0 / dx
                        
            endif

            ! we always use the 2nd order scheme for vertical derivatives
            twz= ( thz(i,j,k)*w(i,j,k) - thz(i,j,k-1)*w(i,j,k-1) ) / dz(k)

            ! horizontal derivative, check to use 2nd or 4th schemes
            if ( cds .eq. 1) then

                txx= ( t(i+1,j,k) - &
                  2.d0*t(i,j,k)   + &
                       t(i-1,j,k) ) / dx2

                tyy= ( t(i,j+1,k) - &
                  2.d0*t(i,j,k)   + &
                       t(i,j-1,k) ) / dy2

            elseif ( cds .eq. 2) then

                txx= ( -t(i+2,j,k) + &
                  16.d0*t(i+1,j,k) - &
                  30.d0*t(i,j,k)   + &
                  16.d0*t(i-1,j,k) - &
                        t(i-2,j,k) ) / 12.d0 / dx2

                tyy= ( -t(i,j+2,k) + &
                  16.d0*t(i,j+1,k) - &
                  30.d0*t(i,j,k)   + &
                  16.d0*t(i,j-1,k) - &
                        t(i,j-2,k) ) / 12.d0 / dy2
                        
            endif

            ! we always use the 2nd order scheme for vertical derivatives
            tzz= ( (t(i,j,k+1)-t(i,j,k))/dz_t(k) - &
                   (t(i,j,k)-t(i,j,k-1))/dz_b(k) ) / dz(k)

            ! nonlinear term for T
            ss(i,j,k)  = -tux-tvy-twz - u_mrf*tx

            ! the right hand side term for T
            qt(i,j,k)= t(i,j,k) + dt*( rk_a*ss(i,j,k)+rk_b*ss0(i,j,k) ) + &
                                  rk_c*0.5d0*dt * alpha * (txx+tyy+tzz)


        enddo
        enddo
        enddo

        return

    end subroutine




    


end module
