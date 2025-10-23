module classical_medium_mod

    use constants_mod

    implicit none

    type classic_medium

        integer         , allocatable :: indx(:) ! locates the material in the real grid
        double precision, allocatable :: PDx(:)
        
        character(len=20) :: medium_type
        character(len=10) :: metal
        integer           :: n_tot      ! total number of grid points occupide by the material
        double precision  :: w0
        double precision  :: Gamma
        double precision  :: eps_r
        double precision  :: A1
        double precision  :: A2
        double precision  :: C1
        double precision  :: C2
        double precision  :: C3
        double precision  :: C4

        double precision, allocatable :: AP_i(:)
        double precision, allocatable :: omegaP_i(:)
        double precision, allocatable :: GammaP_i(:)
        double precision, allocatable :: PLx(:,:)
        double precision, allocatable :: PLx_P(:,:)
        double precision, allocatable :: alpha_k(:)
        double precision, allocatable :: beta_k(:)
        double precision, allocatable :: gamma_k(:)
        double precision, allocatable :: B1_k(:)
        double precision, allocatable :: B2_k(:)
        double precision, allocatable :: Ex_P(:)
        double precision, allocatable :: tmpPL(:)

        
        integer          :: n_poles
        double precision :: omegaD
        double precision :: GammaD
        double precision :: omegaP

        contains
            procedure :: init_classical_medium, kill_classical_medium, get_polarization

    end type classic_medium

    contains

    subroutine init_classical_medium(this, z_coor, z_min, z_max, w0, Gamma, eps_r, dt, Nz, medium_type, metal)

        class(classic_medium) , intent(inout) :: this
        character(len=20)     , intent(in)    :: medium_type
        character(len=10)     , intent(in)    :: metal
        integer               , intent(in)    :: Nz
        double precision      , intent(in)    :: z_coor(Nz)
        double precision      , intent(in)    :: z_min
        double precision      , intent(in)    :: z_max
        double precision      , intent(in)    :: w0
        double precision      , intent(in)    :: Gamma
        double precision      , intent(in)    :: eps_r
        double precision      , intent(in)    :: dt

        integer :: kk, ii
        integer :: n_count


        this%medium_type = medium_type
        this%metal       = metal
        this%eps_r       = eps_r

        n_count = 0
        do kk=1, Nz
            if ((z_coor(kk)>=z_min) .and. (z_coor(kk)<=z_max)) then
                n_count = n_count + 1
            end if
        end do

        this%n_tot = n_count

        if (.not. allocated(this%indx)) allocate(this%PDx(n_count))
        if (.not. allocated(this%indx)) allocate(this%indx(n_count))

        this%PDx   = 0.0d0
        this%indx  = 0

        ii=1
        do kk=1, Nz
            if ((z_coor(kk)>=z_min) .and. (z_coor(kk)<=z_max)) then
                this%indx(ii) = kk
                ii = ii + 1
            end if
        end do

        select case(this%medium_type)

        case("dielectric")
        
        case("drude")

            this%w0    = w0
            this%Gamma = Gamma

            this%A1 = (2.0-Gamma*dt)/(2.0+Gamma*dt)
            this%A2 = eps0* w0**2 *dt/(2.0+Gamma*dt)
            this%C1 = (eps_r*eps0/dt-0.5*this%A2)/(eps_r*eps0/dt+0.5*this%A2)
            this%C3 = 1.0/(eps_r*eps0/dt+0.5*this%A2)
            this%C4 = 0.5*(this%A1+1.0)/(eps_r*eps0/dt+0.5*this%A2)

        case ("drude-lorentz")

            select case(this%metal)
            case("Ag")

                this%n_poles = 5

                if (.not. allocated(this%omegaP_i)) allocate(this%omegaP_i(this%n_poles))
                if (.not. allocated(this%GammaP_i)) allocate(this%GammaP_i(this%n_poles))
                if (.not. allocated(this%AP_i))     allocate(this%AP_i(this%n_poles))
                
                this%GammaD = 0.048 * ev_to_au ; this%omegaP = 9.01 * ev_to_au 
                this%omegaD = (9.01*0.919238815542512) * ev_to_au  !omega_p*sqrt(f0 = 0.845)
                this%omegaP_i(1) = 0.816 * ev_to_au ; this%omegaP_i(2) = 4.481 * ev_to_au 
                this%omegaP_i(3) = 8.185 * ev_to_au ; this%omegaP_i(4) = 9.083 * ev_to_au
                this%omegaP_i(5) = 20.29 * ev_to_au
                this%GammaP_i(1) = 3.886 * ev_to_au ; this%GammaP_i(2) = 0.452 * ev_to_au
                this%GammaP_i(3) = 0.065 * ev_to_au ; this%GammaP_i(4) = 0.916 * ev_to_au
                this%GammaP_i(5) = 2.419* ev_to_au
                this%AP_i(1)     = 0.065 * this%omegaP**2 ; this%AP_i(2)     = 0.124 * this%omegaP**2
                this%AP_i(3)     = 0.011 * this%omegaP**2 ; this%AP_i(4)     = 0.840 * this%omegaP**2
                this%AP_i(5)     = 5.646 * this%omegaP**2

            case("Al")
                this%n_poles = 4

                if (.not. allocated(this%omegaP_i)) allocate(this%omegaP_i(this%n_poles))
                if (.not. allocated(this%GammaP_i)) allocate(this%GammaP_i(this%n_poles))
                if (.not. allocated(this%AP_i))     allocate(this%AP_i(this%n_poles))
                
                this%GammaD = 0.047 * ev_to_au ; this%omegaP = 14.98 * ev_to_au 
                this%omegaD = (14.98 * 0.723187389270582) * ev_to_au  !omega_p*sqrt(f0 = 0.523)
                this%omegaP_i(1) = 0.162 * ev_to_au ; this%omegaP_i(2) = 1.544 * ev_to_au 
                this%omegaP_i(3) = 1.808 * ev_to_au ; this%omegaP_i(4) = 3.473 * ev_to_au
                this%GammaP_i(1) = 0.333 * ev_to_au ; this%GammaP_i(2) = 0.312 * ev_to_au
                this%GammaP_i(3) = 1.351 * ev_to_au ; this%GammaP_i(4) = 3.382 * ev_to_au
                this%AP_i(1)     = 0.227 * this%omegaP**2 ; this%AP_i(2)     = 0.050 * this%omegaP**2
                this%AP_i(3)     = 0.166 * this%omegaP**2 ; this%AP_i(4)     = 0.030 * this%omegaP**2
                
                
            case("Au")

                this%n_poles = 5

                if (.not. allocated(this%omegaP_i)) allocate(this%omegaP_i(this%n_poles))
                if (.not. allocated(this%GammaP_i)) allocate(this%GammaP_i(this%n_poles))
                if (.not. allocated(this%AP_i))     allocate(this%AP_i(this%n_poles))
                
                this%GammaD = 0.053 * ev_to_au ; this%omegaP = 9.03 * ev_to_au 
                this%omegaD = (9.03*0.871779788708135) * ev_to_au  !omega_p*sqrt(f0 = 0.760)
                this%omegaP_i(1) = 0.415 * ev_to_au ; this%omegaP_i(2) = 0.830 * ev_to_au 
                this%omegaP_i(3) = 2.969 * ev_to_au ; this%omegaP_i(4) = 4.304 * ev_to_au
                this%omegaP_i(5) = 13.32 * ev_to_au
                this%GammaP_i(1) = 0.241 * ev_to_au ; this%GammaP_i(2) = 0.345 * ev_to_au
                this%GammaP_i(3) = 0.870 * ev_to_au ; this%GammaP_i(4) = 2.494 * ev_to_au
                this%GammaP_i(5) = 2.214 * ev_to_au
                this%AP_i(1)     = 0.024 * this%omegaP**2 ; this%AP_i(2) = 0.010 * this%omegaP**2
                this%AP_i(3)     = 0.071 * this%omegaP**2 ; this%AP_i(4) = 0.601 * this%omegaP**2
                this%AP_i(5)     = 4.384 * this%omegaP**2

            end select
            
            if (.not. allocated(this%alpha_k))  allocate(this%alpha_k(this%n_poles))
            if (.not. allocated(this%beta_k))   allocate(this%beta_k(this%n_poles))
            if (.not. allocated(this%gamma_k))  allocate(this%gamma_k(this%n_poles))
            if (.not. allocated(this%B1_k))     allocate(this%B1_k(this%n_poles))
            if (.not. allocated(this%B2_k))     allocate(this%B2_k(this%n_poles))
            if (.not. allocated(this%PLx))      allocate(this%PLx(this%n_poles, n_count))
            if (.not. allocated(this%PLx_P))    allocate(this%PLx_P(this%n_poles, n_count))
            if (.not. allocated(this%Ex_P))     allocate(this%Ex_P(n_count))
            if (.not. allocated(this%tmpPL))    allocate(this%tmpPL(this%n_poles))

            this%PLx   = 0.0d0 
            this%PLx_P = 0.0d0
            this%Ex_P  = 0.0d0
            
            do ii =1, this%n_poles
                this%alpha_k(ii) = (2.0-this%omegaP_i(ii)**2*dt**2)/(1.0+0.5*this%GammaP_i(ii)*dt)
                this%beta_k(ii)  = (this%GammaP_i(ii)*dt-2.0)/(this%GammaP_i(ii)*dt+2.0)
                this%gamma_k(ii) = eps0*this%AP_i(ii)*dt/(this%GammaP_i(ii)*dt+2.0)
            end do

            this%A1 = (2.0-this%GammaD*dt)/(2.0+this%GammaD*dt)
            this%A2 = eps0*this%omegaD**2*dt/(2.0+this%GammaD*dt)

            this%C1=(eps_r*eps0/dt-0.5*this%A2)/(eps_r*eps0/dt+0.5*this%A2+0.5*SUM(this%gamma_k))
            this%C2=0.5*SUM(this%gamma_k)/(eps_r*eps0/dt+0.5*this%A2+0.5*SUM(this%gamma_k))
            this%C3=1.0/(eps_r*eps0/dt+0.5*this%A2+0.5*SUM(this%gamma_k))
            this%C4=0.5*(this%A1+1.0)/(eps_r*eps0/dt+0.5*this%A2+0.5*SUM(this%gamma_k))
            
            do ii = 1, this%n_poles
                this%B1_k(ii)=0.5*(1.0+this%alpha_k(ii))/(eps_r*eps0/dt+0.5*this%A2+0.5*SUM(this%gamma_k))
                this%B2_k(ii)=0.5*this%beta_k(ii)/(eps_r*eps0/dt+0.5*this%A2+0.5*SUM(this%gamma_k))
            enddo

        end select


    end subroutine init_classical_medium

    subroutine kill_classical_medium(this)

        class(classic_medium), intent(inout) :: this

        if (allocated(this%PDx))      deallocate(this%PDx)
        if (allocated(this%indx))     deallocate(this%indx)       
        if (allocated(this%omegaP_i)) deallocate(this%omegaP_i)
        if (allocated(this%GammaP_i)) deallocate(this%GammaP_i)
        if (allocated(this%AP_i))     deallocate(this%AP_i)
        if (allocated(this%alpha_k))  deallocate(this%alpha_k)
        if (allocated(this%beta_k))   deallocate(this%beta_k)
        if (allocated(this%gamma_k))  deallocate(this%gamma_k)
        if (allocated(this%B1_k))     deallocate(this%B1_k)
        if (allocated(this%B2_k))     deallocate(this%B2_k)
        if (allocated(this%PLx))      deallocate(this%PLx)
        if (allocated(this%PLx_P))    deallocate(this%PLx_P)
        if (allocated(this%Ex_P))     deallocate(this%Ex_P)
        if (allocated(this%tmpPL))    deallocate(this%tmpPL)

    end subroutine kill_classical_medium

    subroutine get_polarization(this, rotH, Ex, tmpE, kk)

        class(classic_medium), intent(inout) :: this
        double precision     , intent(in)    :: rotH
        double precision     , intent(in)    :: Ex
        double precision     , intent(out)   :: tmpE
        integer              , intent(in)    :: kk

        integer :: ii

        select case(this%medium_type)

        case ("dielectric")

        case ("drude")

           tmpE         = this%C1*Ex+this%C3*rotH-this%C4*this%PDx(kk)
           this%PDx(kk) = this%A1*this%PDx(kk)+this%A2*(tmpE+Ex) 

        case ("drude-lorentz")

            tmpE = this%C1 * Ex + this%C2 * this%Ex_P(kk) + this%C3 * rotH - this%C4 * this%PDx(kk)
            
            do ii = 1, this%n_poles
                tmpE = tmpE - this%B1_k(ii) * this%PLx(ii, kk) - this%B2_k(ii) * this%PLx_P(ii, kk)
            end do

            this%PDx(kk) = this%A1 * this%PDx(kk) + this%A2 * (tmpE+Ex)

            do ii = 1, this%n_poles
                this%tmpPL(ii) = this%alpha_k(ii) * this%PLx(ii,kk)      &
                               + this%beta_k(ii)  * this%PLx_P(ii, kk)   &
                               + this%gamma_k(ii) * (tmpE-this%Ex_P(kk))
            end do

            this%Ex_P(kk) = Ex

            do ii = 1, this%n_poles
                this%PLx_P(ii, kk) = this%PLx(ii, kk)
                this%PLx(ii, kk)   = this%tmpPL(ii)
            end do

        end select

    end subroutine get_polarization

end module classical_medium_mod