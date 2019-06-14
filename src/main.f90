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



program main

    use module_navier_stokes
    use module_parameters
    use mpi
    use decomp_2d

    implicit none

    ! read and initiate all the parameters
    call initiate_parameters

    ! initiate mpi
    call initiate_mpi
    
    ! create a folder to store the output results
    if (myid .eq. 0) then

        call system("mkdir -p ./results")

        if (isxy2d .eq. 1 .or. isxz2d .eq. 1 .or. isyz2d .eq. 1) then
            call system("mkdir -p ./results/slices")
        endif

    endif

    ! initiate the 2decomp lib
    call decomp_2d_init(nx,ny,nz,p_row,p_col)

    ! initiate some indices, e.g., istr3, iend3, i_offset, etc.
    ! we need to do this after decomp_2d_init
    call initiate_indices

    ! welcome info
    if (myid .eq. 0) then
        call welcome_info
    endif

    ! main subroutine for HERCULES
    call drive

    ! goodbye info
    if (myid .eq. 0) then
        call goodby_info
    endif

    ! finalize mpi
    call finalize_mpi
    call decomp_2d_finalize

    stop




contains


    ! this function prints some information on the screen
    subroutine welcome_info

        use mpi
        use module_parameters

        implicit none

        ! these are start date, time, etc
        character(len=8)     :: date1
        character(len=10)    :: time1
        character(len=5)     :: zone1
        integer,dimension(8) :: d1

        print *, '!**********************************************************!'
        print *, '!                      HERCULES V1.1                       !'
        print *, '!**********************************************************!'
        print *, ''

        ! channel configuration
        if (ochannel .eq. 0) then
            print *, 'Configurations: Closed Channel'
        elseif (ochannel .eq. 1) then
            print *, 'Configurations: Open Channel'
        else
            print *,'ochannel error'
            stop
        endif

        ! domain sizes
        write(*,'(A8,F8.3,A14,F8.3,A14,F8.3)'), &
        ' Lx    =',lx,' Ly    =',ly,' Lz    =',lz

        ! grid numbers
        write(*,'(A8,I8,A14,I8,A14,I8)'), &
        ' Nx    =',nx,' Ny    =',ny,' Nz    =',nz

        ! Reynolds number etc
        if (dp_opt .eq. 1) then ! the constant friction case

            write(*,'(A8,F8.1,A14,F8.1,A14,F8.2)'), &
            ' Re_t  =',1.d0/nu,' Ri_t  =',got,' Pr    =',nu/alpha

        elseif (dp_opt .eq. 2) then ! the constant mass flow rate case

            write(*,'(A8,F8.1,A14,F8.3,A14,F8.2)'), &
            ' Re_b  =',1.d0/nu,' Ri_b  =',got,' Pr    =',nu/alpha

        elseif (dp_opt .eq. 3) then ! the Ekman case

            write(*,'(A8,F8.5,A14,F8.2,A14,F8.2)'), &
            ' nu    =',nu,' g/t   =',got,' Pr    =',nu/alpha

            write(*,'(A8,F8.5,A14,F8.3,A14,F8.3)'), &
            ' f     =',fc,' Ug    =',ug,' Vg    =',vg

        endif

        ! temperature boundary conditions
        write(*,'(A8,F8.2,A14,F8.2,A14,F8.2)'), &
        ' Tb    =',tbot,' Tt    =',ttop,' T_ref =',t_ref

        ! simulation setup
        write(*,'(A8,I8,A14,F8.4)'), ' Its   =',imax,' dt    =',dt
        write(*,'(A11,I5)'), ' IsScalar: ',isscalar
        write(*,'(A11,I5)'), ' Isdamp  : ',isdamp
        write(*,'(A11,F5.2)'),' zstretch: ',zstretch
        write(*,'(A11,F5.2)'), ' u_mrf   : ',u_mrf
        write(*,'(A11,I5)'), ' Nprc    :  ',nprc

        ! spatial schemes
        select case (cds)
            case (1)
                write(*,'(A36)'), ' Spatial Derivatives : 2nd Order CDS'
            case (2)
                write(*,'(A36)'), ' Spatial Derivatives : 4th Order CDS'
            case (3)
                write(*,'(A36)'), ' Spatial Derivatives : Fourier      '
            case default
                print *, 'cds cheme error'
                stop
        end select

        ! print if the skew-symmetric form is used for nonlinear terms
        if (cds .eq. 3 .and. issk .eq. 1) then
            write(*,'(A25)'), ' Skew-Symmetric      : On'
        elseif (cds .eq. 3 .and. issk .ne. 1) then
            write(*,'(A26)'), ' Skew-Symmetric      : Off'
        endif

        ! temporal scheme
        select case (dts)
            case (1)
                write(*,'(A33)'), ' Time Scheme         : AB2 and CN'
            case (2)
                write(*,'(A33)'), ' Time Scheme         : RK3 and CN'
            case default
                print *, 'dts scheme error'
                stop
        end select

        ! how to drive the flow
        select case (dp_opt)
            case (1)
                print *, 'dP Option           : Constant Friction'
            case (2)
                print *, 'dP Option           : Constant Mass Flow Rate'
            case (3)
                print *, 'dP Option           : Constant Geostrophic Wind'
            case default
                print *, 'dp_opt error'
                stop
        end select

        if (is_ri_var .eq. 1) then
            print *,     'Cooling Surface     : On'
            write(*,'(A22,F8.4)'), &
               ' Cooling Rate        : ',(ri_end-ri_str)*1.d0/imax/dt
        elseif (is_ri_var .eq. 0) then
            print *,     'Cooling Surface     : Off'
        else
            print *, 'is_ri_var error'
            stop
        endif

        print *, ''

        ! print the start time
        call date_and_time(date1,time1,zone1,d1)
        write(*,'(A21,A8,A2,I2,A1,I2,A1,I2)'), &
              ' Simulation Started: ',date1, ', ',d1(5),':',d1(6),':',d1(7)
        print *, ''


        return

    end subroutine


    ! this function prints some information on the screen
    subroutine goodby_info

        use mpi
        use module_parameters

        implicit none

        ! these are start date, time, etc
        character(len=8)     :: date1
        character(len=10)    :: time1
        character(len=5)     :: zone1
        integer,dimension(8) :: d1

        ! print the end time
        call date_and_time(date1,time1,zone1,d1)
        write(*,'(A22,A8,A2,I2,A1,I2,A1,I2)'), ' Simulation finished: ',date1, &
            ', ',d1(5),':',d1(6),':',d1(7)

        print *, ''
        print *, '!***********************************************************!'
        print *, '!                    Simulation Finished!!                  !'
        print *, '!***********************************************************!'

        return

    end subroutine


end program
