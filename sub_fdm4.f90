module sub_fdm4
    implicit none
    real(8), parameter :: C20 = 1.0 
    real(8), parameter :: C40 = 9.0 / 8.0 
    real(8), parameter :: C41 = 1.0 / 24.0 
contains

    ! receive array Q and return x-derivative dxQ.
    ! differentiated quantity is located at voxcel center,
    ! derivative is located at voxcel boundary.
    subroutine FDIFF4_x_c2b(Q, dxQ, Nx, Nz, dx)
        implicit none
        integer, intent(in) :: Nx, Nz
        real(8), intent(in) :: dx
        real(8), intent(in) :: Q(Nx, Nz)
        real(8), intent(out) :: dxQ(Nx, Nz)

        integer :: i, j
        real(8) :: D20, D40, D41
        
        D20 = C20 / dx
        D40 = C40 / dx
        D41 = C41 / dx

        do j = 1, Nz
            dxQ(1, j) = ( Q(2, j) - Q(1, j) ) * D20
            do i = 2, Nx-2
                dxQ(i, j) = ( Q(i+1, j) - Q(i,   j) ) * D40 &
                        & - ( Q(i+2, j) - Q(i-1, j) ) * D41
            end do
            dxQ(Nx-1, j) = ( Q(Nx, j) - Q(Nx-1, j) ) * D20
            dxQ(Nx,   j) = ( 0.0      - Q(Nx  , j) ) * D20
            ! dxQ(Nx,   j) = ( Q(Nx, j) - Q(Nx-1, j) ) * D20
        end do

        return
        
    end subroutine FDIFF4_x_c2b


    ! receive array Q and return x-derivative dxQ.
    ! differentiated quantity is located at voxcel boundary,
    ! derivative is located at voxcel center.
    subroutine FDIFF4_x_b2c(Q, dxQ, Nx, Nz, dx)
        implicit none
        integer, intent(in) :: Nx, Nz
        real(8), intent(in) :: dx
        real(8), intent(in) :: Q(Nx, Nz)
        real(8), intent(out) :: dxQ(Nx, Nz)

        integer :: i, j
        real(8) :: D20, D40, D41
        
        D20 = C20 / dx
        D40 = C40 / dx
        D41 = C41 / dx

        do j = 1, Nz
            ! dxQ(1, j) = ( Q(2, j) - Q(1, j) ) * D20
            dxQ(1, j) = ( Q(1, j) - 0.0     ) * D20
            dxQ(2, j) = ( Q(2, j) - Q(1, j) ) * D20
            do i = 3, Nx-1
                dxQ(i, j) = ( Q(i,   j) - Q(i-1, j) ) * D40 &
                        & - ( Q(i+1, j) - Q(i-2, j) ) * D41
            end do
            dxQ(Nx, j) = ( Q(Nx, j) - Q(Nx-1, j) ) * D20
        end do

        return
        
    end subroutine FDIFF4_x_b2c



    ! symmetric derivative
    subroutine FDIFF4_x_c2b_sym(Q, dxQ, Nx, Nz, dx)
        implicit none
        integer, intent(in) :: Nx, Nz
        real(8), intent(in) :: dx
        real(8), intent(in) :: Q(Nx+1, Nz)
        real(8), intent(out) :: dxQ(Nx, Nz)

        integer :: i, j
        real(8) :: D20, D40, D41
        
        D20 = C20 / dx
        D40 = C40 / dx
        D41 = C41 / dx

        do j = 1, Nz
            dxQ(1, j) = ( Q(2, j) - Q(1, j) ) * D20
            do i = 2, Nx-1
                dxQ(i, j) = ( Q(i+1, j) - Q(i,   j) ) * D40 &
                        & - ( Q(i+2, j) - Q(i-1, j) ) * D41
            end do
            dxQ(Nx, j) = ( Q(Nx+1, j) - Q(Nx, j) ) * D20
        end do

        return
        
    end subroutine FDIFF4_x_c2b_sym

    subroutine FDIFF4_x_b2c_sym(Q, dxQ, Nx, Nz, dx)
        implicit none
        integer, intent(in) :: Nx, Nz
        real(8), intent(in) :: dx
        real(8), intent(in) :: Q(Nx, Nz)
        real(8), intent(out) :: dxQ(Nx+1, Nz)

        integer :: i, j
        real(8) :: D20, D40, D41
        
        D20 = C20 / dx
        D40 = C40 / dx
        D41 = C41 / dx

        do j = 1, Nz
            ! dxQ(1, j) = ( Q(2, j) - Q(1, j) ) * D20
            dxQ(1, j) = ( Q(1, j) - 0.0     ) * D20
            dxQ(2, j) = ( Q(2, j) - Q(1, j) ) * D20
            do i = 3, Nx-1
                dxQ(i, j) = ( Q(i,   j) - Q(i-1, j) ) * D40 &
                        & - ( Q(i+1, j) - Q(i-2, j) ) * D41
            end do
            dxQ(Nx  , j) = ( Q(Nx, j) - Q(Nx-1, j) ) * D20
            dxQ(Nx+1, j) = ( 0.0      - Q(Nx  , j) ) * D20
            ! dxQ(Nx+1, j) = ( Q(Nx, j) - Q(Nx-1, j) ) * D20
        end do

        return
        
    end subroutine FDIFF4_x_b2c_sym


    ! receive array Q and return z-derivative dzQ.
    ! differentiated quantity is located at voxcel center,
    ! derivative is located at voxcel boundary.
    subroutine FDIFF4_z_c2b(Q, dzQ, Nx, Nz, dz)
        implicit none
        integer, intent(in) :: Nx, Nz
        real(8), intent(in) :: dz
        real(8), intent(in) :: Q(Nx, Nz)
        real(8), intent(out) :: dzQ(Nx, Nz)

        integer :: i, j
        real(8) :: D20, D40, D41
        
        D20 = C20 / dz
        D40 = C40 / dz
        D41 = C41 / dz

        do i = 1, Nx
            dzQ(i, 1) = ( Q(i, 2) - Q(i, 1) ) * D20
        end do

        do j = 2, Nz-2
            do i = 1, Nx
                dzQ(i, j) = ( Q(i, j+1) - Q(i, j  ) ) * D40 &
                        & - ( Q(i, j+2) - Q(i, j-1) ) * D41
            end do
        end do

        do i = 1, Nx
            dzQ(i, Nz-1) = ( Q(i, Nz) - Q(i, Nz-1) ) * D20
            dzQ(i, Nz  ) = ( Q(i, Nz) - Q(i, Nz-1) ) * D20
        end do

        return
        
    end subroutine FDIFF4_z_c2b


    ! receive array Q and return z-derivative dzQ.
    ! differentiated quantity is located at voxcel boundary,
    ! derivative is located at voxcel center.
    subroutine FDIFF4_z_b2c(Q, dzQ, Nx, Nz, dz)
        implicit none
        integer, intent(in) :: Nx, Nz
        real(8), intent(in) :: dz
        real(8), intent(in) :: Q(Nx, Nz)
        real(8), intent(out) :: dzQ(Nx, Nz)

        integer :: i, j
        real(8) :: D20, D40, D41
        
        D20 = C20 / dz
        D40 = C40 / dz
        D41 = C41 / dz

        do i = 1, Nx
            dzQ(i, 1) = ( Q(i, 2) - Q(i, 1) ) * D20 ! This value is not refered during calculation
            dzQ(i, 2) = ( Q(i, 2) - Q(i, 1) ) * D20
        end do

        do j = 3, Nz-1
            do i = 1, Nx
                dzQ(i, j) = ( Q(i, j  ) - Q(i, j-1) ) * D40 &
                        & - ( Q(i, j+1) - Q(i, j-2) ) * D41
            end do
        end do

        do i = 1, Nx
            dzQ(i, Nz) = ( Q(i, Nz) - Q(i, Nz-1) ) * D20
        end do

        return
        
    end subroutine FDIFF4_z_b2c


    ! receive array Ux and return z-derivative dzUx.
    ! differentiated quantity Ux is located at voxcel center,
    ! derivative dzUx is located at voxcel boundary.
    subroutine FDIFF4_z_c2b_across(Q, dzQ, Nx, Nz, JSF, dz)
        implicit none
        integer, intent(in) :: Nx, Nz, JSF
        real(8), intent(in) :: dz
        real(8), intent(in) :: Q(Nx, Nz)
        real(8), intent(out) :: dzQ(Nx, Nz)

        integer :: i, j
        real(8) :: D20, D40, D41

        if ( JSF <= 4 .or. JSF >= Nz-4 ) then
            write(6, *) '[Error @ sub_fdm4] location of JSF is incorrect.'
            stop
        end if
        
        D20 = C20 / dz
        D40 = C40 / dz
        D41 = C41 / dz

        do i = 1, Nx
            dzQ(i, 1) = ( Q(i, 2) - Q(i, 1) ) * D20
        end do

        do j = 2, JSF-2
            do i = 1, Nx
                dzQ(i, j) = ( Q(i, j+1) - Q(i, j  ) ) * D40 &
                        & - ( Q(i, j+2) - Q(i, j-1) ) * D41
            end do
        end do

        do i = 1, Nx
            dzQ(i, JSF-1) = ( Q(i, JSF  ) - Q(i, JSF-1) ) * D20
            dzQ(i, JSF  ) = 0.0       ! ( Q(i, JSF+1) - Q(i, JSF) ) * D20
            dzQ(i, JSF+1) = ( Q(i, JSF+2) - Q(i, JSF+1) ) * D20
        end do

        do j = JSF+2, Nz-2
            do i = 1, Nx
                dzQ(i, j) = ( Q(i, j+1) - Q(i, j  ) ) * D40 &
                        & - ( Q(i, j+2) - Q(i, j-1) ) * D41
            end do
        end do

        do i = 1, Nx
            dzQ(i, Nz-1) = ( Q(i, Nz) - Q(i, Nz-1) ) * D20
            dzQ(i, Nz  ) = ( Q(i, Nz) - Q(i, Nz-1) ) * D20
        end do

        return
        
    end subroutine FDIFF4_z_c2b_across
    



    ! receive array Uz and return z-derivative dzUz.
    ! differentiated quantity Uz is located at voxcel boundary,
    ! derivative dzUz is located at voxcel center.
    subroutine FDIFF4_z_b2c_across(Uz, dzUz, Nx, Nz, JSF, dz)
        implicit none
        integer, intent(in) :: Nx, Nz, JSF
        real(8), intent(in) :: dz
        real(8), intent(in) :: Uz(1:Nx, 0:Nz)
        real(8), intent(out) :: dzUz(1:Nx, 1:Nz)

        integer :: i, j
        real(8) :: D20, D40, D41

        if ( JSF <= 5 .or. JSF >= Nz-5 ) then
            write(6, *) '[Error @ sub_fdm4] location of JSF is incorrect.'
            stop
        end if
        
        D20 = C20 / dz
        D40 = C40 / dz
        D41 = C41 / dz

        do i = 1, Nx
            dzUz(i, 1) = ( Uz(i, 1) - Uz(i, 0) ) * D20
        end do

        do j = 2, JSF-1
            do i = 1, Nx
                dzUz(i, j) = ( Uz(i, j  ) - Uz(i, j-1) ) * D40 &
                         & - ( Uz(i, j+1) - Uz(i, j-2) ) * D41
            end do
        end do

        do i = 1, Nx
            dzUz(i, JSF  ) = ( Uz(i, JSF  ) - Uz(i, JSF-1) ) * D20
            dzUz(i, JSF+1) = ( Uz(i, JSF+1) - Uz(i, JSF  ) ) * D20
        end do

        do j = JSF+2, Nz-1
            do i = 1, Nx
                dzUz(i, j) = ( Uz(i, j  ) - Uz(i, j-1) ) * D40 &
                         & - ( Uz(i, j+1) - Uz(i, j-2) ) * D41
            end do
        end do

        do i = 1, Nx
            dzUz(i, Nz) = ( Uz(i, Nz) - Uz(i, Nz-1) ) * D20
        end do

        return
        
    end subroutine FDIFF4_z_b2c_across

end module sub_fdm4