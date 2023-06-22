program MaedaFurumura
    use sub_fdm4 ! subroutines for finite diff.
    use sub_program ! other subroutines
    implicit none
    ! [Examples of parameters]
    ! When dx = dz = 200 m and dt = 0.01 s (These satisfy CFL conditon),
    ! Nx =  1000 & Nt =  50000 for short calculation ~ 500 s (e.g. just check the tsunami propagation),
    ! Nx = 10000 & Nt = 500000 for long calculation ~ 5000 s (e.g. investigate dispersion relations for long period).
    !---------- parameter (we use MKS system of units.) ----------!
    integer, parameter :: Nx = 10000, Nz = 200                   ! number of grid, must be even
    integer, parameter :: JSF = 20, Nzs = Nz-JSF                ! JSF = location of seafloor
    integer, parameter :: i0 = Nx/2+1, j0 = JSF+50                ! position of the epicenter (j0 is measured from sea surface)
    integer, parameter :: Nx_pml = 20, Nz_pml = 20              ! width of the perfectly matched layer
    integer, parameter :: Nt = 500000                            ! number of iteration
    integer, parameter :: NWRITE1 = 100, NWRITE2 = 600000        ! interval of record (1:sea surface & seafloor 2:whole)
    ! If you don't want to record whole waves, just set NWRITE2 > Nt.
    real(8), parameter :: dx = 200.0, dz = 200.0  ! space grid size [m]
    real(8), parameter :: dt = 0.01               ! time step [s], determined by CFL condition ( dt <= 0.3 * min(dx,dz) / Vmax )
    real(8), parameter :: T0 = 0.0, Trise = 2.0   ! when earthquake occurred / width of source time function (=rise time)
    real(8), parameter :: f0 = 2.0/Trise
    real(8), parameter :: M0 = 1.0d14             ! seismic moment (absolute)
    real(8), parameter :: mxx = -1.0/sqrt(2.0), mzz = 1.0/sqrt(2.0), mxz = 0.0 ! components of seismic moment tensor (relative)

    !---------- parameter (medium) ----------!
    ! a:air, w:water, s:solid
    real(8), parameter :: rhoa = 1.3
    real(8), parameter :: rhow_0 = 1.0d3, Cw = 1.5d3 , Qw = 1.0d6
    real(8), parameter :: Vp = 5.0d3, Vs = 3.0d3, Vps2 = Vp**2 - 2 * Vs**2
    real(8), parameter :: rhos_0 = 3.0d3, Qs = 200.0
    character*128, parameter :: title = 'MF_T5000s'


    !---------- constant ----------!
    real(8), parameter :: pi = acos(-1.0)
    real(8), parameter :: g = 9.81

    !---------- other variables ----------!
    integer :: ierror1, ierror2, ierror3, ierror4, ierror5
    integer :: i, j, ii, jj, ipml, jpml, il, ir, ir1, it
    real(8) :: Vmin, Vmax, Courant, lambda_min, PPW
    real(8) :: QQw, QQs, t, stime, sdrop, Hscale_w, Hscale_s
    real(8) :: dzUx_ij, dzUz_ij
    real(8) :: dxUx_l, dxUx_r, dxUz_l, dxUz_r, dzUx_l, dzUx_r, dzUz_l, dzUz_r
    real(8) :: accelx, accelx_l, accelx_r, accelz, accelz_l, accelz_r
    real(8) :: dxP_l, dxSxx_l, dxSxz_l, dzSxz_l, dzSzz_l
    real(8) :: dxP_r, dxSxx_r, dxSxz_r, dzSxz_r, dzSzz_r
    real(8) :: dzSxz_ij, dzSzz_ij

    !-------------------- arrays --------------------!
    !---------- Common Variables ----------!
    real(8) ::   Ux(1:Nx  , 1:Nz),   Vx(1:Nx  , 1:Nz)
    real(8) ::   Uz(1:Nx+1, 0:Nz),   Vz(1:Nx+1, 0:Nz)
    real(8) :: dxUx(1:Nx+1, 1:Nz), dzUx(1:Nx  , 1:Nz)
    real(8) :: dxUz(1:Nx  , 0:Nz), dzUz(1:Nx+1, 1:Nz)

    !---------- Variables for air & fluid ----------!
    real(8) :: P(1:Nx+1, 1:JSF), dxP(1:Nx, 1:JSF), dzP(1:Nx+1, 1:JSF)
    real(8) :: rhow_c(1:JSF), rhow_b(1:JSF)
    !---------- Variables for solid ----------!
    real(8) :: Sxx(1:Nx+1, JSF+1:Nz), Sxz(1:Nx, JSF:Nz), Szz(1:Nx+1, JSF+1:Nz)
    real(8) :: dxSxx(1:Nx, JSF+1:Nz), dxSxz(1:Nx+1, JSF:Nz), dzSxz(1:Nx, JSF:Nz), dzSzz(1:Nx+1, JSF+1:Nz)
    real(8) :: rhos_c(JSF+1:Nz), rhos_b(JSF+1:Nz)

    !---------- Variables for ADE CFS-PML ----------!
    real(8) ::  axUx_l(1:Nx_pml,     1:Nz ),  axUx_r(Nx-Nx_pml+2:Nx+1,     1:Nz )
    real(8) ::  axUz_l(1:Nx_pml,     0:Nz ),  axUz_r(Nx-Nx_pml+1:Nx  ,     0:Nz )
    real(8) ::   axP_l(1:Nx_pml,     1:JSF),   axP_r(Nx-Nx_pml+1:Nx  ,     1:JSF)
    real(8) :: axSxx_l(1:Nx_pml, JSF+1:Nz ), axSxx_r(Nx-Nx_pml+1:Nx  , JSF+1:Nz )
    real(8) :: axSxz_l(1:Nx_pml,   JSF:Nz ), axSxz_r(Nx-Nx_pml+2:Nx+1,   JSF:Nz )

    real(8) ::  azUx(1:Nx  , Nz-Nz_pml+1:Nz),  azUz(1:Nx+1, Nz-Nz_pml+1:Nz)
    real(8) :: azSzz(1:Nx+1, Nz-Nz_pml+1:Nz), azSxz(1:Nx  , Nz-Nz_pml+1:Nz)

    ! c : voxel's ceter / b : voxel's boundadry
    real(8) :: Gx_c(4, Nx_pml), Gx_b(4, Nx_pml)
    real(8) :: Gz_c(4, Nz_pml), Gz_b(4, Nz_pml)


    character*128 :: dir_name, fname_Ux_floor, fname_Uz_floor, fname_Uz_surface, fname_Uxall, fname_Uzall

    dir_name = './record_' // trim(title)
    call system('mkdir -p ' // trim(dir_name))

    fname_Uz_surface = trim(dir_name) // '/Uz_surface_' // trim(title) // '.dat'
    fname_Ux_floor = trim(dir_name) // '/Ux_floor_' // trim(title) // '.dat'
    fname_Uz_floor = trim(dir_name) // '/Uz_floor_' // trim(title) // '.dat'
    


    ! prepare file i/o
    open(11, file = fname_Uz_surface, form = 'unformatted', access='stream', &
    & status = 'replace', iostat = ierror1)
    if ( ierror1 /= 0 ) then
        write(6, *) 'failed to open file (11)'
        stop
    end if

    open(12, file = fname_Uz_floor, form = 'unformatted', access='stream', &
    & status = 'replace', iostat = ierror2)
    if ( ierror2 /= 0 ) then
        write(6, *) 'failed to open file (12)'
        stop
    end if

    open(13, file = fname_Ux_floor, form = 'unformatted', access='stream', &
    & status = 'replace', iostat = ierror3)
    if ( ierror3 /= 0 ) then
        write(6, *) 'failed to open file (13)'
        stop
    end if



    write(6, *) 'This is Maeda & Furumura (2013) formulation (symmetric grid system).'
    write(6, *) 'Width of the calculation area (x-direction):', Nx*dx/1e3, '[km]'
    write(6, *) 'Width of the calculation area (z-direction):', Nz*dz/1e3, '[km]'
    write(6, *) 'Depth of the ocean:', JSF*dz/1e3, '[km]'
    write(6, *) 'Quality factor in the solid:', Qs



    !---------- check CFL condition & wavelength condition ----------!
    Vmin = min(Cw, Vp, Vs)
    Vmax = max(Cw, Vp, Vs)
    Courant = Vmax * dt / min(dx, dz)
    lambda_min = Vmin / f0
    PPW = lambda_min / max(dx, dz)

    if ( 0.0 < Courant .and. Courant <= 1.0 ) then
        write(6, *) 'CFL condition is satisfied. ( Courant Number = ', Courant, ')'
    else
        write(6, *) 'CFL condition is NOT satisfied !!! ( Courant Number = ', Courant, ')'
        stop
    end if

    if ( PPW >= 7.0 ) then
        write(6, *) 'wavelength condition is satisfied. ( lambda / grid = ', PPW, ')'
    else
        write(6, *) 'wavelength condition is NOT satisfied !!! ( lambda / grid = ', PPW, ')'
        stop
    end if


    !-------------------- initialize --------------------!
    Ux(:, :) = 0.0 ; Uz(:, :) = 0.0
    Vx(:, :) = 0.0 ; Vz(:, :) = 0.0


    !-------------------- ADE CFS-PML ( Zhang and Shen 2010) --------------------!
     axUx_l(:, :) = 0.0 ;  axUx_r(:, :) = 0.0
     axUz_l(:, :) = 0.0 ;  axUz_r(:, :) = 0.0
      axP_l(:, :) = 0.0 ;   axP_r(:, :) = 0.0
    axSxx_l(:, :) = 0.0 ; axSxx_r(:, :) = 0.0
    axSxz_l(:, :) = 0.0 ; axSxz_r(:, :) = 0.0

     azUx(:, :) = 0.0 ;  azUz(:, :) = 0.0
    azSzz(:, :) = 0.0 ; azSxz(:, :) = 0.0

    ! i,j -> distance from the boundary between PML & physical domain / dx or dz
    do ipml = 1, Nx_pml 
        call adecfspml_profile( dble(ipml)-0.5, Nx_pml, dx, dt, Vp, f0, Gx_c(:, ipml) )
        call adecfspml_profile( dble(ipml)-1.0, Nx_pml, dx, dt, Vp, f0, Gx_b(:, ipml) )
    end do

    do jpml = 1, Nz_pml
        call adecfspml_profile( dble(jpml)-1.0, Nz_pml, dz, dt, Vp, f0, Gz_c(:, jpml) )
        call adecfspml_profile( dble(jpml)-0.5, Nz_pml, dz, dt, Vp, f0, Gz_b(:, jpml) )
    end do



    !--------------- ANELASTICITY CONDITION  Blanch et al. (1995) Geophys. ---------------#
    QQw = exp( -pi * f0 * dt / Qw )
    QQs = exp( -pi * f0 * dt / Qs )


    ! stratified density profile
    Hscale_w = Cw**2 / g
    Hscale_s = ( Vp**2 - Vs**2 * 4.0/3.0 ) / g
    do j = 1, JSF
        rhow_c(j) = rhow_0 * exp( ( j - 0.5 ) * dz / Hscale_w )
        rhow_b(j) = rhow_0 * exp( j * dz / Hscale_w )
    end do
    do j = JSF+1, Nz
        jj = j - JSF
        rhos_c(j) = rhos_0 * exp( ( jj - 0.5 ) * dz / Hscale_s )
        rhos_b(j) = rhos_0 * exp( jj * dz / Hscale_s )
    end do

    write(6, *) 'scale height of sea water :', Hscale_w/1e3, '[km]'
    write(6, *) 'water density at seafloor :', rhow_b(JSF), '[kg/m3]'
    write(6, *) 'scale height of solid :', Hscale_s/1e3, '[km]'
    write(6, *) 'solid density at the bottom :', rhos_b(Nz), '[kg/m3]'


    

    !---------------------------------------- TIME STEP START ----------------------------------------!
    do it = 0, Nt
        if ( mod(it, 500) == 0 ) then
            write(6, *) it
        end if

        t = it * dt

        !-------------------- record --------------------!
        if ( mod(it, NWRITE1) == 0 ) then
            write(11) -Uz(:, 0  )
            write(12) -Uz(:, JSF)
            write(13) 1.5 * Ux(:, JSF+1) - 0.5 * Ux(:, JSF+2)
        end if
        
        if ( it >= 1 .and. mod(it, NWRITE2) == 0 ) then

            write(fname_Uxall, '("/Ux_", i0, ".dat")') it
            write(fname_Uzall, '("/Uz_", i0, ".dat")') it

            fname_Uxall = trim(dir_name) // trim(fname_Uxall)
            fname_Uzall = trim(dir_name) // trim(fname_Uzall)

            open(14, file = fname_Uxall, form = 'unformatted', access='stream', &
            & status = 'replace', iostat = ierror4)
            open(15, file = fname_Uzall, form = 'unformatted', access='stream', &
            & status = 'replace', iostat = ierror5)

            if ( ierror4 /= 0 ) then
                write(6, *) 'failed to open file (14)'
                stop
            end if
            if ( ierror5 /= 0 ) then
                write(6, *) 'failed to open file (15)'
                stop
            end if

            write(14) Ux(:, :)
            write(15) Uz(:, :)
            
            close(14)
            close(15)
        end if


        !------------------------------ Ui -> diUj ------------------------------!
        call FDIFF4_x_b2c_sym(Ux, dxUx, Nx, Nz, dx) ! OK
        call FDIFF4_z_c2b_across(Ux, dzUx, Nx, Nz, JSF, dz) ! OK
        call FDIFF4_x_c2b_sym(Uz, dxUz, Nx, Nz+1, dx) ! OK
        call FDIFF4_z_b2c_across(Uz, dzUz, Nx+1, Nz, JSF, dz)
    



        !------------------------------ Constitutive Law ------------------------------!
        !-------------------- [water] strain -> pressure --------------------!        
        do j = 1, JSF
            do i = Nx_pml+1, Nx-Nx_pml+1
                P(i, j) = - rhow_c(j) * Cw**2 * ( dxUx(i, j) + dzUz(i, j) )
            end do
        end do

        !-------------------- [water] strain -> pressure (PML) --------------------!
        do j = 1, JSF
            do ipml = 1, Nx_pml
                il = Nx_pml+1-ipml ; ir = Nx-Nx_pml+ipml ; ir1 = ir+1

                dxUx_l = Gx_c(1, ipml) * dxUx(il , j) + Gx_c(2, ipml) * axUx_l(il , j) 
                dxUx_r = Gx_c(1, ipml) * dxUx(ir1, j) + Gx_c(2, ipml) * axUx_r(ir1, j)
                 
                P(il , j) = - rhow_c(j) * Cw**2 * ( dxUx_l + dzUz(il , j) )
                P(ir1, j) = - rhow_c(j) * Cw**2 * ( dxUx_r + dzUz(ir1, j) )
            end do
        end do


        !---------- boundary condition @ seafloor ----------!
        Sxz(:, JSF) = 0.0

        !-------------------- [solid] strain -> stress tensor --------------------!
        do j = JSF+1, Nz-Nz_pml
            do i = Nx_pml+1, Nx-Nx_pml
                Sxx(i, j) = ( Vp**2 * dxUx(i, j) + Vps2 * dzUz(i, j) ) * rhos_c(j) 
                Szz(i, j) = ( Vp**2 * dzUz(i, j) + Vps2 * dxUx(i, j) ) * rhos_c(j)
                Sxz(i, j) = Vs**2 * ( dxUz(i, j) + dzUx(i, j) ) * rhos_b(j)
            end do
            i = Nx-Nx_pml+1
            Sxx(i, j) = ( Vp**2 * dxUx(i, j) + Vps2 * dzUz(i, j) ) * rhos_c(j) 
            Szz(i, j) = ( Vp**2 * dzUz(i, j) + Vps2 * dxUx(i, j) ) * rhos_c(j)
        end do

        !-------------------- [solid] strain -> stress tensor (PML - left and right side) --------------------!
        do j = JSF+1, Nz-Nz_pml
            do ipml = 1, Nx_pml
                il = Nx_pml+1-ipml ; ir = Nx-Nx_pml+ipml ; ir1 = ir+1

                dxUx_l = Gx_c(1, ipml) * dxUx(il , j) + Gx_c(2, ipml) * axUx_l(il , j)
                dxUx_r = Gx_c(1, ipml) * dxUx(ir1, j) + Gx_c(2, ipml) * axUx_r(ir1, j)

                dxUz_l = Gx_b(1, ipml) * dxUz(il, j) + Gx_b(2, ipml) * axUz_l(il, j)
                dxUz_r = Gx_b(1, ipml) * dxUz(ir, j) + Gx_b(2, ipml) * axUz_r(ir, j)  
            
                Sxx(il, j) = ( Vp**2 * dxUx_l + Vps2 * dzUz(il, j) ) * rhos_c(j)
                Szz(il, j) = ( Vp**2 * dzUz(il, j) + Vps2 * dxUx_l ) * rhos_c(j)
                Sxz(il, j) = Vs**2 * ( dxUz_l + dzUx(il, j) ) * rhos_b(j)

                Sxx(ir1, j) = ( Vp**2 * dxUx_r + Vps2 * dzUz(ir1, j) ) * rhos_c(j) 
                Szz(ir1, j) = ( Vp**2 * dzUz(ir1, j) + Vps2 * dxUx_r ) * rhos_c(j) 
                Sxz(ir , j) = Vs**2 * ( dxUz_r + dzUx(ir, j) ) * rhos_b(j)
            end do
        end do

        !------------------------------ solid strain -> stress tensor (PML - bottom) ------------------------------!
        do jpml = 1, Nz_pml
            j = Nz-Nz_pml+jpml
            do i = Nx_pml+1, Nx-Nx_pml
                dzUx_ij = Gz_b(1, jpml) * dzUx(i, j) + Gz_b(2, jpml) * azUx(i, j) 
                dzUz_ij = Gz_c(1, jpml) * dzUz(i, j) + Gz_c(2, jpml) * azUz(i, j) 

                Sxx(i, j) = ( Vp**2 * dxUx(i, j) + Vps2 * dzUz_ij ) * rhos_c(j) 
                Szz(i, j) = ( Vp**2 * dzUz_ij + Vps2 * dxUx(i, j) ) * rhos_c(j)
                Sxz(i, j) = Vs**2 * ( dxUz(i, j) + dzUx_ij ) * rhos_b(j)
            end do
            i = Nx-Nx_pml+1
            dzUz_ij = Gz_c(1, jpml) * dzUz(i, j) + Gz_c(2, jpml) * azUz(i, j) 

            Sxx(i, j) = ( Vp**2 * dxUx(i, j) + Vps2 * dzUz_ij ) * rhos_c(j) 
            Szz(i, j) = ( Vp**2 * dzUz_ij + Vps2 * dxUx(i, j) ) * rhos_c(j)
        end do



        !------------------------------ solid strain -> stress tensor (PML - bottomleft and bottomright corner) ------------------------------!
        do jpml = 1, Nz_pml
            j = Nz-Nz_pml+jpml
            do ipml = 1, Nx_pml
                il = Nx_pml+1-ipml ; ir = Nx-Nx_pml+ipml ; ir1 = ir+1

                dxUx_l = Gx_c(1, ipml) * dxUx(il , j) + Gx_c(2, ipml) * axUx_l(il , j)
                dxUx_r = Gx_c(1, ipml) * dxUx(ir1, j) + Gx_c(2, ipml) * axUx_r(ir1, j)

                dxUz_l = Gx_b(1, ipml) * dxUz(il, j) + Gx_b(2, ipml) * axUz_l(il, j) 
                dxUz_r = Gx_b(1, ipml) * dxUz(ir, j) + Gx_b(2, ipml) * axUz_r(ir, j) 

                dzUx_l = Gz_b(1, jpml) * dzUx(il, j) + Gz_b(2, jpml) *   azUx(il, j) 
                dzUx_r = Gz_b(1, jpml) * dzUx(ir, j) + Gz_b(2, jpml) *   azUx(ir, j)

                dzUz_l = Gz_c(1, jpml) * dzUz(il , j) + Gz_c(2, jpml) *   azUz(il , j) 
                dzUz_r = Gz_c(1, jpml) * dzUz(ir1, j) + Gz_c(2, jpml) *   azUz(ir1, j)
        
                Sxx(il, j) = ( Vp**2 * dxUx_l + Vps2 * dzUz_l ) * rhos_c(j)
                Szz(il, j) = ( Vp**2 * dzUz_l + Vps2 * dxUx_l ) * rhos_c(j)
                Sxz(il, j) = Vs**2 * ( dxUz_l + dzUx_l ) * rhos_b(j)

                Sxx(ir1, j) = ( Vp**2 * dxUx_r + Vps2 * dzUz_r ) * rhos_c(j)
                Szz(ir1, j) = ( Vp**2 * dzUz_r + Vps2 * dxUx_r ) * rhos_c(j)
                Sxz(ir , j) = Vs**2 * ( dxUz_r + dzUx_r ) * rhos_b(j)
            end do
        end do






        !---------- inject source by stress drop ----------!
        if ( T0 <= t ) then
            stime = IKupper(t-T0, Trise)
            sdrop = M0 * stime / ( dx * dz )

            Sxx(i0, j0) = Sxx(i0, j0) - mxx * sdrop
            Szz(i0, j0) = Szz(i0, j0) - mzz * sdrop

            Sxz(i0-1, j0-1) = Sxz(i0-1, j0-1) - mxz * sdrop / 4.0
            Sxz(i0-1, j0  ) = Sxz(i0-1, j0  ) - mxz * sdrop / 4.0
            Sxz(i0  , j0-1) = Sxz(i0  , j0-1) - mxz * sdrop / 4.0
            Sxz(i0  , j0  ) = Sxz(i0  , j0  ) - mxz * sdrop / 4.0
        end if
        

        !------------------------------ P -> diP ------------------------------!
        call FDIFF4_x_c2b_sym(P, dxP, Nx, JSF, dx) ! OK
        call FDIFF4_z_c2b(P, dzP, Nx+1, JSF, dz) ! OK

        !------------------------------ Sij -> diSij ------------------------------!
        call FDIFF4_x_c2b_sym(Sxx, dxSxx, Nx, Nzs  , dx) ! OK
        call FDIFF4_x_b2c_sym(Sxz, dxSxz, Nx, Nzs+1, dx) ! OK
        call FDIFF4_z_b2c(Sxz, dzSxz, Nx, Nzs+1, dz) ! OK
        call FDIFF4_z_c2b(Szz, dzSzz, Nx+1, Nzs, dz) ! OK


        !------------------------------ update velocity & displacement ------------------------------!
        !---------------------------------- ( equation of motion ) ----------------------------------!        

        !------------------------------ sea surface (j = 0) ------------------------------!
        do i = Nx_pml+1, Nx-Nx_pml+1
            accelz = ( -2 * P(i, 1) / dz ) / ( rhow_0 + rhoa )
            Vz(i, 0) = ( Vz(i, 0) + dt * accelz ) * QQw
            Uz(i, 0) = ( Uz(i, 0) + dt * Vz(i, 0) ) * QQw
        end do
        

        !------------------------------ sea surface (j = 0 ; PML) ------------------------------!
        do ipml = 1, Nx_pml
            il = Nx_pml+1-ipml ; ir = Nx-Nx_pml+ipml ; ir1 = ir+1

            accelz_l = ( -2 * P(il , 1) / dz ) / ( rhow_0 + rhoa )
            accelz_r = ( -2 * P(ir1, 1) / dz ) / ( rhow_0 + rhoa )

            Vz(il , 0) = ( Vz(il , 0) + dt * accelz_l ) * QQw
            Vz(ir1, 0) = ( Vz(ir1, 0) + dt * accelz_r ) * QQw

            Uz(il , 0) = ( Uz(il , 0) + dt * Vz(il , 0) ) * QQw
            Uz(ir1, 0) = ( Uz(ir1, 0) + dt * Vz(ir1, 0) ) * QQw
        end do


        !------------------------------ water ------------------------------!
        do j = 1, JSF-1
            do i = Nx_pml+1, Nx-Nx_pml
                accelx = ( -dxP(i, j) + rhow_0 * g * dxUz(i, 0) ) / rhow_c(j)
                accelz = -dzP(i, j) / rhow_b(j)

                Vx(i, j) = ( Vx(i, j) + dt * accelx ) * QQw
                Vz(i, j) = ( Vz(i, j) + dt * accelz ) * QQw

                Ux(i, j) = ( Ux(i, j) + dt * Vx(i, j) ) * QQw
                Uz(i, j) = ( Uz(i, j) + dt * Vz(i, j) ) * QQw
            end do

            i = Nx-Nx_pml+1
            accelz = -dzP(i, j) / rhow_b(j)

            Vz(i, j) = ( Vz(i, j) + dt * accelz ) * QQw
            Uz(i, j) = ( Uz(i, j) + dt * Vz(i, j) ) * QQw
        end do


        !------------------------------ water (PML) ------------------------------!
        do j = 1, JSF-1
            do ipml = 1, Nx_pml
                il = Nx_pml+1-ipml ; ir = Nx-Nx_pml+ipml ; ir1 = ir+1

                dxUz_l = Gx_b(1, ipml) * dxUz(il, 0) + Gx_b(2, ipml) * axUz_l(il, 0) 
                dxUz_r = Gx_b(1, ipml) * dxUz(ir, 0) + Gx_b(2, ipml) * axUz_r(ir, 0) 

                dxP_l = Gx_b(1, ipml) * dxP(il, j) + Gx_b(2, ipml) * axP_l(il, j)
                dxP_r = Gx_b(1, ipml) * dxP(ir, j) + Gx_b(2, ipml) * axP_r(ir, j)

                accelx_l = ( -dxP_l + rhow_0 * g * dxUz_l ) / rhow_c(j)
                accelz_l = - dzP(il, j) / rhow_b(j)

                accelx_r = ( -dxP_r + rhow_0 * g * dxUz_r ) / rhow_c(j)
                accelz_r = - dzP(ir1, j) / rhow_b(j)

                Vx(il, j) = ( Vx(il, j) + dt * accelx_l ) * QQw
                Vx(ir, j) = ( Vx(ir, j) + dt * accelx_r ) * QQw

                Vz(il , j) = ( Vz(il , j) + dt * accelz_l ) * QQw
                Vz(ir1, j) = ( Vz(ir1, j) + dt * accelz_r ) * QQw

                Ux(il, j) = ( Ux(il, j) + dt * Vx(il, j) ) * QQw
                Ux(ir, j) = ( Ux(ir, j) + dt * Vx(ir, j) ) * QQw

                Uz(il , j) = ( Uz(il , j) + dt * Vz(il , j) ) * QQw
                Uz(ir1, j) = ( Uz(ir1, j) + dt * Vz(ir1, j) ) * QQw
            end do
        end do


        !------------------------------ seafloor (j = JSF) ------------------------------!
        do i = Nx_pml+1, Nx-Nx_pml
            accelx = ( -dxP(i, JSF) + rhow_0 * g * dxUz(i, 0) ) / rhow_c(JSF)
            accelz = ( 2 * ( Szz(i, JSF+1) + P(i, JSF) ) / dz ) / ( rhos_0 + rhow_b(JSF) )

            Vx(i, JSF) = ( Vx(i, JSF) + dt * accelx ) * QQw
            Vz(i, JSF) = ( Vz(i, JSF) + dt * accelz ) * QQs

            Ux(i, JSF) = ( Ux(i, JSF) + dt * Vx(i, JSF) ) * QQw
            Uz(i, JSF) = ( Uz(i, JSF) + dt * Vz(i, JSF) ) * QQs
        end do

        i = Nx-Nx_pml+1
        accelz = ( 2 * ( Szz(i, JSF+1) + P(i, JSF) ) / dz ) / ( rhos_0 + rhow_b(JSF) )

        Vz(i, JSF) = ( Vz(i, JSF) + dt * accelz ) * QQs
        Uz(i, JSF) = ( Uz(i, JSF) + dt * Vz(i, JSF) ) * QQs


        !------------------------------ seafloor (j = JSF ; PML) ------------------------------!
        do ipml = 1, Nx_pml
            il = Nx_pml+1-ipml ; ir = Nx-Nx_pml+ipml ; ir1 = ir+1

            dxUz_l = Gx_b(1, ipml) * dxUz(il, 0) + Gx_b(2, ipml) * axUz_l(il, 0) 
            dxUz_r = Gx_b(1, ipml) * dxUz(ir, 0) + Gx_b(2, ipml) * axUz_r(ir, 0) 

            dxP_l = Gx_b(1, ipml) * dxP(il, JSF) + Gx_b(2, ipml) * axP_l(il, JSF)
            dxP_r = Gx_b(1, ipml) * dxP(ir, JSF) + Gx_b(2, ipml) * axP_r(ir, JSF)

            accelx_l = ( -dxP_l + rhow_0 * g * dxUz_l ) / rhow_c(JSF) 
            accelx_r = ( -dxP_r + rhow_0 * g * dxUz_r ) / rhow_c(JSF)

            accelz_l = ( 2 * ( Szz(il , JSF+1) + P(il , JSF) ) / dz ) / ( rhos_0 + rhow_b(JSF) )
            accelz_r = ( 2 * ( Szz(ir1, JSF+1) + P(ir1, JSF) ) / dz ) / ( rhos_0 + rhow_b(JSF) )

            Vx(il, JSF) = ( Vx(il, JSF) + dt * accelx_l ) * QQw
            Vx(ir, JSF) = ( Vx(ir, JSF) + dt * accelx_r ) * QQw

            Vz(il , JSF) = ( Vz(il , JSF) + dt * accelz_l ) * QQs
            Vz(ir1, JSF) = ( Vz(ir1, JSF) + dt * accelz_r ) * QQs
            
            Ux(il, JSF) = ( Ux(il, JSF) + dt * Vx(il, JSF) ) * QQw
            Ux(ir, JSF) = ( Ux(ir, JSF) + dt * Vx(ir, JSF) ) * QQw

            Uz(il , JSF) = ( Uz(il , JSF) + dt * Vz(il , JSF) ) * QQs
            Uz(ir1, JSF) = ( Uz(ir1, JSF) + dt * Vz(ir1, JSF) ) * QQs
        end do



        !------------------------------ solid ------------------------------!
        do j = JSF+1, Nz-Nz_pml
            do i = Nx_pml+1, Nx-Nx_pml

                accelx = ( dxSxx(i, j) + dzSxz(i, j) + rhow_0 * g * dxUz(i, 0) ) / rhos_c(j)
                accelz = ( dxSxz(i, j) + dzSzz(i, j) ) / rhos_b(j)

                Vx(i, j) = ( Vx(i, j) + dt * accelx ) * QQs
                Vz(i, j) = ( Vz(i, j) + dt * accelz ) * QQs

                Ux(i, j) = ( Ux(i, j) + dt * Vx(i, j) ) * QQs
                Uz(i, j) = ( Uz(i, j) + dt * Vz(i, j) ) * QQs
            end do

            i = Nx-Nx_pml+1
            accelz = ( dxSxz(i, j) + dzSzz(i, j) ) / rhos_b(j)

            Vz(i, j) = ( Vz(i, j) + dt * accelz ) * QQs
            Uz(i, j) = ( Uz(i, j) + dt * Vz(i, j) ) * QQs
        end do


        !------------------------------ solid (PML - left and right side) ------------------------------!
        do j = JSF+1, Nz-Nz_pml
            do ipml = 1, Nx_pml
                il = Nx_pml+1-ipml ; ir = Nx-Nx_pml+ipml ; ir1 = ir+1

                dxUz_l = Gx_b(1, ipml) * dxUz(il, 0) + Gx_b(2, ipml) * axUz_l(il, 0) 
                dxUz_r = Gx_b(1, ipml) * dxUz(ir, 0) + Gx_b(2, ipml) * axUz_r(ir, 0)

                dxSxx_l = Gx_b(1, ipml) * dxSxx(il, j) + Gx_b(2, ipml) * axSxx_l(il, j) 
                dxSxx_r = Gx_b(1, ipml) * dxSxx(ir, j) + Gx_b(2, ipml) * axSxx_r(ir, j)

                dxSxz_l = Gx_c(1, ipml) * dxSxz(il , j) + Gx_c(2, ipml) * axSxz_l(il , j) 
                dxSxz_r = Gx_c(1, ipml) * dxSxz(ir1, j) + Gx_c(2, ipml) * axSxz_r(ir1, j) 

                accelx_l = ( dxSxx_l + dzSxz(il, j) + rhow_0 * g * dxUz_l ) / rhos_c(j)
                accelx_r = ( dxSxx_r + dzSxz(ir, j) + rhow_0 * g * dxUz_r ) / rhos_c(j)

                accelz_l = ( dxSxz_l + dzSzz(il , j) ) / rhos_b(j)
                accelz_r = ( dxSxz_r + dzSzz(ir1, j) ) / rhos_b(j)

                
                Vx(il, j) = ( Vx(il, j) + dt * accelx_l ) * QQs
                Vx(ir, j) = ( Vx(ir, j) + dt * accelx_r ) * QQs

                Vz(il , j) = ( Vz(il , j) + dt * accelz_l ) * QQs
                Vz(ir1, j) = ( Vz(ir1, j) + dt * accelz_r ) * QQs

                Ux(il, j) = ( Ux(il, j) + dt * Vx(il, j) ) * QQs
                Ux(ir, j) = ( Ux(ir, j) + dt * Vx(ir, j) ) * QQs

                Uz(il , j) = ( Uz(il , j) + dt * Vz(il , j) ) * QQs
                Uz(ir1, j) = ( Uz(ir1, j) + dt * Vz(ir1, j) ) * QQs
            end do
        end do


        !------------------------------ solid (PML - bottom) ------------------------------!
        do jpml = 1, Nz_pml
            j = Nz-Nz_pml+jpml
            do i = Nx_pml+1, Nx-Nx_pml

                dzSxz_ij = Gz_c(1, jpml) * dzSxz(i, j) + Gz_c(2, jpml) * azSxz(i, j) 
                dzSzz_ij = Gz_b(1, jpml) * dzSzz(i, j) + Gz_b(2, jpml) * azSzz(i, j)

                accelx = ( dxSxx(i, j) + dzSxz_ij + rhow_0 * g * dxUz(i, 0) ) / rhos_c(j)
                accelz = ( dxSxz(i, j) + dzSzz_ij ) / rhos_b(j)

                Vx(i, j) = ( Vx(i, j) + dt * accelx ) * QQs
                Vz(i, j) = ( Vz(i, j) + dt * accelz ) * QQs

                Ux(i, j) = ( Ux(i, j) + dt * Vx(i, j) ) * QQs
                Uz(i, j) = ( Uz(i, j) + dt * Vz(i, j) ) * QQs
            end do

            i = Nx-Nx_pml+1

            dzSzz_ij = Gz_b(1, jpml) * dzSzz(i, j) + Gz_b(2, jpml) * azSzz(i, j)
            accelz = ( dxSxz(i, j) + dzSzz_ij ) / rhos_b(j)

            Vz(i, j) = ( Vz(i, j) + dt * accelz ) * QQs
            Uz(i, j) = ( Uz(i, j) + dt * Vz(i, j) ) * QQs
        end do


        !------------------------------ solid (PML - bottomleft and bottomright corner) ------------------------------!
        do jpml = 1, Nz_pml
            j = Nz-Nz_pml+jpml
            do ipml = 1, Nx_pml
                il = Nx_pml+1-ipml ; ir = Nx-Nx_pml+ipml ; ir1 = ir+1

                dxUz_l = Gx_b(1, ipml) * dxUz(il, 0) + Gx_b(2, ipml) * axUz_l(il, 0) 
                dxUz_r = Gx_b(1, ipml) * dxUz(ir, 0) + Gx_b(2, ipml) * axUz_r(ir, 0)

                dxSxx_l = Gx_b(1, ipml) * dxSxx(il, j) + Gx_b(2, ipml) * axSxx_l(il, j) 
                dxSxx_r = Gx_b(1, ipml) * dxSxx(ir, j) + Gx_b(2, ipml) * axSxx_r(ir, j)

                dxSxz_l = Gx_c(1, ipml) * dxSxz(il , j) + Gx_c(2, ipml) * axSxz_l(il , j)
                dxSxz_r = Gx_c(1, ipml) * dxSxz(ir1, j) + Gx_c(2, ipml) * axSxz_r(ir1, j)

                dzSxz_l = Gz_c(1, jpml) * dzSxz(il, j) + Gz_c(2, jpml) *   azSxz(il, j)
                dzSxz_r = Gz_c(1, jpml) * dzSxz(ir, j) + Gz_c(2, jpml) *   azSxz(ir, j) 
                
                dzSzz_l = Gz_b(1, jpml) * dzSzz(il , j) + Gz_b(2, jpml) *   azSzz(il , j)
                dzSzz_r = Gz_b(1, jpml) * dzSzz(ir1, j) + Gz_b(2, jpml) *   azSzz(ir1, j) 

                accelx_l = ( dxSxx_l + dzSxz_l + rhow_0 * g * dxUz_l ) / rhos_c(j)
                accelx_r = ( dxSxx_r + dzSxz_r + rhow_0 * g * dxUz_r ) / rhos_c(j)

                accelz_l = ( dxSxz_l + dzSzz_l ) / rhos_b(j)
                accelz_r = ( dxSxz_r + dzSzz_r ) / rhos_b(j)
                

                Vx(il, j) = ( Vx(il, j) + dt * accelx_l ) * QQs
                Vx(ir, j) = ( Vx(ir, j) + dt * accelx_r ) * QQs

                Vz(il , j) = ( Vz(il , j) + dt * accelz_l ) * QQs
                Vz(ir1, j) = ( Vz(ir1, j) + dt * accelz_r ) * QQs

                Ux(il, j) = ( Ux(il, j) + dt * Vx(il, j) ) * QQs
                Ux(ir, j) = ( Ux(ir, j) + dt * Vx(ir, j) ) * QQs

                Uz(il , j) = ( Uz(il , j) + dt * Vz(il , j) ) * QQs
                Uz(ir1, j) = ( Uz(ir1, j) + dt * Vz(ir1, j) ) * QQs
            end do
        end do


        !------------------------------ update memory variables ------------------------------!
        !------------------------ ( auxiliary differential equation ) ------------------------!

        !-------------------- axP --------------------!
        do j = 1, JSF
            do ipml = 1, Nx_pml
                il = Nx_pml+1-ipml ; ir = Nx-Nx_pml+ipml ; ir1 = ir+1

                axP_l(il, j) = Gx_b(3, ipml) * axP_l(il, j) + Gx_b(4, ipml) * dxP(il, j)
                axP_r(ir, j) = Gx_b(3, ipml) * axP_r(ir, j) + Gx_b(4, ipml) * dxP(ir, j)
            end do
        end do

        !-------------------- axSij --------------------!
        do j = JSF+1, Nz
            do ipml = 1, Nx_pml 
                il = Nx_pml+1-ipml ; ir = Nx-Nx_pml+ipml ; ir1 = ir+1

                axSxx_l(il, j) = Gx_b(3, ipml) * axSxx_l(il, j) + Gx_b(4, ipml) * dxSxx(il, j)
                axSxx_r(ir, j) = Gx_b(3, ipml) * axSxx_r(ir, j) + Gx_b(4, ipml) * dxSxx(ir, j)
            end do
        end do

        do j = JSF, Nz 
            do ipml = 1, Nx_pml 
                il = Nx_pml+1-ipml ; ir = Nx-Nx_pml+ipml ; ir1 = ir+1

                axSxz_l(il , j) = Gx_c(3, ipml) * axSxz_l(il , j) + Gx_c(4, ipml) * dxSxz(il , j)
                axSxz_r(ir1, j) = Gx_c(3, ipml) * axSxz_r(ir1, j) + Gx_c(4, ipml) * dxSxz(ir1, j) 
            end do
        end do

        !-------------------- axUi --------------------!
        do j = 1, Nz
            do ipml = 1, Nx_pml 
                il = Nx_pml+1-ipml ; ir = Nx-Nx_pml+ipml ; ir1 = ir+1

                axUx_l(il , j) = Gx_c(3, ipml) * axUx_l(il , j) + Gx_c(4, ipml) * dxUx(il , j)
                axUx_r(ir1, j) = Gx_c(3, ipml) * axUx_r(ir1, j) + Gx_c(4, ipml) * dxUx(ir1, j)
            end do
        end do

        do j = 0, Nz 
            do ipml = 1, Nx_pml 
                il = Nx_pml+1-ipml ; ir = Nx-Nx_pml+ipml ; ir1 = ir+1

                axUz_l(il, j) = Gx_b(3, ipml) * axUz_l(il, j) + Gx_b(4, ipml) * dxUz(il, j)
                axUz_r(ir, j) = Gx_b(3, ipml) * axUz_r(ir, j) + Gx_b(4, ipml) * dxUz(ir, j)
            end do
        end do


        !-------------------- azSij, azUi --------------------!
        do jpml = 1, Nz_pml
            j = Nz-Nz_pml+jpml
            do i = 1, Nx 
                azSzz(i, j) = Gz_b(3, jpml) * azSzz(i, j) + Gz_b(4, jpml) * dzSzz(i, j)
                azSxz(i, j) = Gz_c(3, jpml) * azSxz(i, j) + Gz_c(4, jpml) * dzSxz(i, j)
                 azUx(i, j) = Gz_b(3, jpml) *  azUx(i, j) + Gz_b(4, jpml) *  dzUx(i, j)
                 azUz(i, j) = Gz_c(3, jpml) *  azUz(i, j) + Gz_c(4, jpml) *  dzUz(i, j)
            end do
            i = Nx+1
            azSzz(i, j) = Gz_b(3, jpml) * azSzz(i, j) + Gz_b(4, jpml) * dzSzz(i, j)
             azUz(i, j) = Gz_c(3, jpml) *  azUz(i, j) + Gz_c(4, jpml) *  dzUz(i, j)
        end do


    end do
    !---------------------------------------- TIME STEP END ----------------------------------------!

    close(11)
    close(12)
    close(13)

    write(6, *) 'Finished.'

    stop
    
end program MaedaFurumura