module sub_program
    implicit none
contains

    real(8) function Kupper(t, Trise)

        real(8), intent(in) :: t
        real(8), intent(in) :: Trise
        real(8), parameter :: pi = acos(-1.0)

        if ( 0.0 <= t .and. t <= Trise) then
            Kupper = 3*pi * ( sin(pi*t/Trise) )**3 / ( 4*Trise )
        else 
            Kupper = 0.0
        end if

        return
    
    end function Kupper


    real(8) function IKupper(t, Trise)

        real(8), intent(in) :: t, Trise
        real(8), parameter :: pi = acos(-1.0)

        if ( t < 0.0 ) then
            IKupper = 0.0
        else if ( 0.0 <= t .and. t <= Trise) then
            IKupper = ( cos(3*pi*t/Trise) - 9.0 * cos(pi*t/Trise) + 8.0 ) / 16.0
        else if ( t > Trise ) then
            IKupper = 1.0
        end if

        return

    end function IKupper



    !-------------------- dumping profile for ADE CFS-PML --------------------!
    subroutine adecfspml_profile( ipml, Npml, dx, dt, Vp, f0, G )
        integer, intent(in) :: Npml
        real(8), intent(in) :: ipml, dx, dt, Vp, f0
        real(8), intent(out) :: G(4)

        real(8) :: log10R0, R0, d0, alpha0, beta0, xpml
        real(8) :: d, alpha, beta, d_beta, alpha_d_beta
        integer, parameter :: pd = 1, palpha = 1, pbeta = 2
        real(8), parameter :: pi = acos(-1.0)

        ! distance from the boundary between PML & physical domain
        ! xpml = 0 for the boundary between PML & physical domain
        ! xpml = 1 for the edge of computatinal domain
        xpml = ipml / Npml 

        log10R0 = - ( log10(dble(Npml)) - 1.0 ) / log10(2.0)  - 3.0
        R0 = 10 ** log10R0
        d0 = - log(R0) * (pd + 1.0) * Vp / (2.0 * Npml * dx) 
        alpha0 = pi * f0
        beta0 = 8.0

        d = d0 * xpml ** pd
        alpha = alpha0 * ( 1.0 - xpml ** palpha )
        beta = 1.0 + (beta0 - 1.0) * xpml ** pbeta
        d_beta =  d / beta
        alpha_d_beta = alpha + d_beta

        G(1) = ( ( 2.0 + dt * alpha ) / beta ) / ( 2.0 + dt * alpha_d_beta )
        G(2) = ( -2.0 / beta                 ) / ( 2.0 + dt * alpha_d_beta )
        G(3) = ( 2.0 - dt * alpha_d_beta     ) / ( 2.0 + dt * alpha_d_beta )
        G(4) = ( 2.0 * dt * d_beta           ) / ( 2.0 + dt * alpha_d_beta )

        return

    end subroutine adecfspml_profile


end module sub_program