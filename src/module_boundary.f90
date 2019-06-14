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



module module_boundary

    use module_parameters

    implicit none

contains



    ! update boundary/ghost cell values for u,v,w in the z direction
    subroutine update_boundary_uvw(u,v,w)

        implicit none

        real(mytype),dimension(:,:,:),intent(inout) :: u, v, w

        integer                               :: i,j

        !!!!!!!!!!!!!!!! U Eqn !!!!!!!!!!!!!!!!!
        ! Top and bottom ghost cells for u
        ! check if it is open channel or closed channel
        if (ochannel .eq. 0) then

            do j=jstr3,jend3,1
            do i=istr3,iend3,1

                ! if moving reference frame is used, we need to 
                ! give the non-slip walls with a moving speed
                ! of -u_mrf in the x direction
                u(i,j,kend3+1)=-u(i,j,kend3)-2.d0*u_mrf ! non-slip
                u(i,j,kstr3-1)=-u(i,j,kstr3)-2.d0*u_mrf ! non-slip

            enddo
            enddo

        elseif (ochannel .eq. 1) then

            do j=jstr3,jend3,1
            do i=istr3,iend3,1

                ! if moving reference frame is used, we need to 
                ! give the non-slip walls with a moving speed
                ! of -u_mrf in the x direction
                u(i,j,kend3+1)= u(i,j,kend3)            !     slip
                u(i,j,kstr3-1)=-u(i,j,kstr3)-2.d0*u_mrf ! non-slip

            enddo
            enddo

        endif

        !!!!!!!!!!!!!!! V Eqn !!!!!!!!!!!!!!!!!
        ! Top and bottom ghost cells for v
        ! check if it is open channel or closed channel
        if (ochannel .eq. 0) then

            do j=jstr3,jend3,1
            do i=istr3,iend3,1

                v(i,j,kend3+1)=-v(i,j,kend3) ! non-slip
                v(i,j,kstr3-1)=-v(i,j,kstr3) ! non-slip

            enddo
            enddo

        elseif (ochannel .eq. 1) then

            do j=jstr3,jend3,1
            do i=istr3,iend3,1

                v(i,j,kend3+1)= v(i,j,kend3) !     slip
                v(i,j,kstr3-1)=-v(i,j,kstr3) ! non-slip

            enddo
            enddo

        endif

        !!!!!!!!!!!! W Equation !!!!!!!!!!!!!!
        ! Top and bottom boundary for w
        ! no matter it is open or closed channel, w=0 at the top and bottom
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            w(i,j,kend3)=0.d0
            w(i,j,kend3+1)=-w(i,j,kend3-1)
            w(i,j,kstr3-1)=0.d0

        enddo
        enddo

        ! communication between processors 
        ! for horizontal ghost cells
        ! we don't need this for spectral method
        if (cds .ne. 3) then
            
            call comm_bound(u)
            call comm_bound(v)
            call comm_bound(w)

        endif


        return

    end subroutine



    ! update boundary/ghost cell values for p
    subroutine update_boundary_p ( p )

        implicit none

        real(mytype),dimension(:,:,:),intent(inout) :: p

        integer                               :: i,j

        ! Top and bottom ghost cells for p, zero-gradient
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            p(i,j,kend3+1)=p(i,j,kend3)
            p(i,j,kstr3-1)=p(i,j,kstr3)

        enddo
        enddo
        
        ! communication between processors 
        ! for horizontal ghost cells
        ! we don't need this for spectral method
        if (cds .ne. 3) then

            call comm_bound(p)

        endif

        return

    end subroutine



    ! update boundary/ghost cell value for t
    subroutine update_boundary_t ( t )

        implicit none

        real(mytype),dimension(:,:,:),intent(inout) :: t

        integer                                     :: i,j

        ! Top and bottom ghost cells for t, fixed temperature
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            t(i,j,kend3+1)=2.d0*ttop-t(i,j,kend3)
            t(i,j,kstr3-1)=2.d0*tbot-t(i,j,kstr3)

        enddo
        enddo
        
        ! communication between processors 
        ! for horizontal ghost cells
        ! we don't need this for spectral method
        if (cds .ne. 3) then

            call comm_bound(t)

        endif

        return

    end subroutine



    ! update boundary/ghost cell values for interpolated variables uhx, vhx, etc.
    ! this is only needed for FD method
    subroutine update_boundary_uvwhxyz(u,v,uhx,uhy,uhz,vhx,vhy,vhz,whx,why)

        implicit none

        real(mytype),dimension(:,:,:),   intent(in):: u, v
        real(mytype),dimension(:,:,:),intent(inout):: uhx,uhy,uhz, &
                                                      vhx,vhy,vhz,whx,why

        integer                                    :: i,j

        ! check if it is open channel or closed channel
        if (ochannel .eq. 0) then

            ! Top boundary, non-slip
            do j=jstr3,jend3,1
            do i=istr3,iend3,1

                ! if moving reference frame is used, we need to 
                ! give the non-slip walls with a moving speed
                ! of -u_mrf in the x direction
                uhx(i,j,kend3+1)=-uhx(i,j,kend3)-2.d0*u_mrf
                uhy(i,j,kend3+1)=-uhy(i,j,kend3)-2.d0*u_mrf
                uhz(i,j,kend3)  =-u_mrf
                uhz(i,j,kend3+1)=-uhz(i,j,kend3-1)-2.d0*u_mrf

                vhx(i,j,kend3+1)=-vhx(i,j,kend3)
                vhy(i,j,kend3+1)=-vhy(i,j,kend3)
                vhz(i,j,kend3)  = 0.d0
                vhz(i,j,kend3+1)=-vhz(i,j,kend3-1)

                whx(i,j,kend3)  = 0.d0
                whx(i,j,kend3+1)=-whx(i,j,kend3-1)
                why(i,j,kend3)  = 0.d0
                why(i,j,kend3+1)=-why(i,j,kend3-1)

            enddo
            enddo

        elseif (ochannel .eq. 1) then

            ! open channel, slip for the top and non-slip for the bottom
            do j=jstr3,jend3,1
            do i=istr3,iend3,1

                uhx(i,j,kend3+1)= uhx(i,j,kend3)
                uhy(i,j,kend3+1)= uhy(i,j,kend3)
                uhz(i,j,kend3)  =   u(i,j,kend3)
                uhz(i,j,kend3+1)= uhz(i,j,kend3)

                vhx(i,j,kend3+1)= vhx(i,j,kend3)
                vhy(i,j,kend3+1)= vhy(i,j,kend3)
                vhz(i,j,kend3)  =   v(i,j,kend3)
                vhz(i,j,kend3+1)= vhz(i,j,kend3)

                whx(i,j,kend3)  = 0.d0
                whx(i,j,kend3+1)=-whx(i,j,kend3-1)
                why(i,j,kend3)  = 0.d0
                why(i,j,kend3+1)=-why(i,j,kend3-1)

            enddo
            enddo

        endif

        ! Bottom boundary, non-slip
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            ! if moving reference frame is used, we need to 
            ! give the non-slip walls with a moving speed
            ! of -u_mrf in the x direction
            uhx(i,j,kstr3-1)=-uhx(i,j,kstr3)-2.d0*u_mrf
            uhy(i,j,kstr3-1)=-uhy(i,j,kstr3)-2.d0*u_mrf
            uhz(i,j,kstr3-1)=-u_mrf

            vhx(i,j,kstr3-1)=-vhx(i,j,kstr3)
            vhy(i,j,kstr3-1)=-vhy(i,j,kstr3)
            vhz(i,j,kstr3-1)=0.d0

            whx(i,j,kstr3-1)=0.d0
            why(i,j,kstr3-1)=0.d0

        enddo
        enddo
        
        ! communication between processors 
        ! for horizontal ghost cells
        ! we don't need this for spectral method
        if (cds .ne. 3) then

            call comm_bound(uhx)
            call comm_bound(uhy)
            call comm_bound(uhz)
            call comm_bound(vhx)
            call comm_bound(vhy)
            call comm_bound(vhz)
            call comm_bound(whx)
            call comm_bound(why)

        endif

        return

    end subroutine



    ! update boundary/ghost cell values for interpolated variables thx, thx, etc.
    ! this is only needed for FD method
    subroutine update_boundary_thxyz (thx,thy,thz)

        implicit none

        real(mytype),dimension(:,:,:),intent(inout):: thx,thy,thz

        integer                               :: i,j

        ! Top ghost cells for t, fixed temperature
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            thx(i,j,kend3+1)=-thx(i,j,kend3)+2.d0*ttop
            thy(i,j,kend3+1)=-thy(i,j,kend3)+2.d0*ttop
            thz(i,j,kend3)  = ttop
            thz(i,j,kend3+1)=-thz(i,j,kend3-1)+2.d0*ttop

        enddo
        enddo

        ! Bottom ghost cells for t, fixed temperature
        do j=jstr3,jend3,1
        do i=istr3,iend3,1

            thx(i,j,kstr3-1)=-thx(i,j,kstr3)+2.d0*tbot
            thy(i,j,kstr3-1)=-thy(i,j,kstr3)+2.d0*tbot
            thz(i,j,kstr3-1)=tbot

        enddo
        enddo

        ! communication between processors 
        ! for horizontal ghost cells
        ! we don't need this for spectral method
        if (cds .ne. 3) then

            call comm_bound(thx)
            call comm_bound(thy)
            call comm_bound(thz)

        endif

        return

    end subroutine



    ! this function does the communication between processors for the
    ! East and West boundary ghost cells.
    subroutine comm_bound(var)

        use mpi

        implicit none

        real(mytype),dimension(:,:,:),intent(inout) :: var
        ! idsend: the processor id to send data to
        ! isrecv: the processor id to receive data from
        ! count_ns: how many cells to send/receive in the north and south directions
        ! count_ew: how many cells to send/receive in the east and west directions
        integer                           :: idsend,idrecv,count_ns,count_ew
        integer,dimension(MPI_STATUS_SIZE):: status1

        ! calculate how many grids to exchange, note that for north and south,
        ! we need to extend two ghost cells in the x direction
        count_ew=(ny/p_col)*(nz+2)
        count_ns=(nx/p_row+2*ghst)*(nz+2)

        ! DATA exchange
        ! first step: send data to EAST and receive data from WEST
        ! calculate send and receive id
        ! here myid_rowindx, pij is calculated in initiate_indices
        ! if it is the west boundary sub-domain, send data to
        ! the east boundary sub-domain
        ! otherwise, send data to the neighbor sub-domain to the west
        if (myid_rowindx .eq. p_row) then
            idsend=pij(1,myid_colindx)
        else
            idsend=myid+p_col
        endif
        if (myid_rowindx .eq. 1) then
            idrecv=pij(p_row,myid_colindx)
        else
            idrecv=myid-p_col
        endif
        ! do communication for the 1st ghost cells
        call MPI_SENDRECV(var(iend3,  jstr3:jend3,kstr3-1:kend3+1), &
                                 count_ew,real_type,idsend,99, &
                          var(istr3-1,jstr3:jend3,kstr3-1:kend3+1), &
                                 count_ew,real_type,idrecv,99, &
                                 MPI_COMM_WORLD,status1,ierr)
        ! do communication for the 2st ghost cells
        call MPI_SENDRECV(var(iend3-1,jstr3:jend3,kstr3-1:kend3+1), &
                                 count_ew,real_type,idsend,98, &
                          var(istr3-2,jstr3:jend3,kstr3-1:kend3+1), &
                                 count_ew,real_type,idrecv,98, &
                                 MPI_COMM_WORLD,status1,ierr)

        ! second step: send data to WEST and receive data from EAST
        ! calculate send and receive id
        if (myid_rowindx .eq. 1) then
            idsend=pij(p_row,myid_colindx)
        else
            idsend=myid-p_col
        endif
        if (myid_rowindx .eq. p_row) then
            idrecv=pij(1,myid_colindx)
        else
            idrecv=myid+p_col
        endif
        ! do communication for the 1st ghost cells
        call MPI_SENDRECV(var(istr3,  jstr3:jend3,kstr3-1:kend3+1), &
                                 count_ew,real_type,idsend,97, &
                          var(iend3+1,jstr3:jend3,kstr3-1:kend3+1), &
                                 count_ew,real_type,idrecv,97, &
                                 MPI_COMM_WORLD,status1,ierr)
        ! do communication for the 2st ghost cells
        call MPI_SENDRECV(var(istr3+1,jstr3:jend3,kstr3-1:kend3+1), &
                                 count_ew,real_type,idsend,96, &
                          var(iend3+2,jstr3:jend3,kstr3-1:kend3+1), &
                                 count_ew,real_type,idrecv,96, &
                                 MPI_COMM_WORLD,status1,ierr)

        ! third step: send data to NORTH and receive data from SOUTH
        ! calculate send and receive id
        if (myid_colindx .eq. p_col) then
            idsend=pij(myid_rowindx,1)
        else
            idsend=myid+1
        endif
        if (myid_colindx .eq. 1) then
            idrecv=pij(myid_rowindx,p_col)
        else
            idrecv=myid-1
        endif
        ! do communication for the 1st ghost cells
        call MPI_SENDRECV(var(istr3-ghst:iend3+ghst,jend3,  kstr3-1:kend3+1), &
                                 count_ns,real_type,idsend,95, &
                          var(istr3-ghst:iend3+ghst,jstr3-1,kstr3-1:kend3+1), &
                                 count_ns,real_type,idrecv,95, &
                                 MPI_COMM_WORLD,status1,ierr)
        ! do communication for the 2st ghost cells
        call MPI_SENDRECV(var(istr3-ghst:iend3+ghst,jend3-1,kstr3-1:kend3+1), &
                                 count_ns,real_type,idsend,94, &
                          var(istr3-ghst:iend3+ghst,jstr3-2,kstr3-1:kend3+1), &
                                 count_ns,real_type,idrecv,94, &
                                 MPI_COMM_WORLD,status1,ierr)

        ! fourth step: send data to SOUTH and receive data from NORTH
        ! calculate send and receive id
        if (myid_colindx .eq. 1) then
            idsend=pij(myid_rowindx,p_col)
        else
            idsend=myid-1
        endif
        if (myid_colindx .eq. p_col) then
            idrecv=pij(myid_rowindx,1)
        else
            idrecv=myid+1
        endif
        ! do communication for the 1st ghost cells
        call MPI_SENDRECV(var(istr3-ghst:iend3+ghst,jstr3,  kstr3-1:kend3+1), &
                                 count_ns,real_type,idsend,93, &
                          var(istr3-ghst:iend3+ghst,jend3+1,kstr3-1:kend3+1), &
                                 count_ns,real_type,idrecv,93, &
                                 MPI_COMM_WORLD,status1,ierr)
        ! do communication for the 2st ghost cells
        call MPI_SENDRECV(var(istr3-ghst:iend3+ghst,jstr3+1,kstr3-1:kend3+1), &
                                 count_ns,real_type,idsend,92, &
                          var(istr3-ghst:iend3+ghst,jend3+2,kstr3-1:kend3+1), &
                                 count_ns,real_type,idrecv,92, &
                                 MPI_COMM_WORLD,status1,ierr)

        return

    end subroutine



end module
