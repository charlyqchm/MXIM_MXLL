module classical_medium_mod

    use constants_mod

    implicit none

    type drude

        integer         , allocatable :: indx(:) ! locates the material in the real grid
        double precision, allocatable :: Px(:)
        
        integer           :: n_tot      ! total number of grid points occupide by the material
        double precision  :: w0
        double precision  :: Gamma
        double precision  :: eps_r
        double precision  :: A1
        double precision  :: A2
        double precision  :: C1
        double precision  :: C3
        double precision  :: C4

        contains
            procedure :: init_drude, kill_drude

    end type drude

    contains

    subroutine init_drude(this, z_coor, z_min, z_max, w0, Gamma, eps_r, dt, Nz)

        class(drude)    , intent(inout) :: this
        integer         , intent(in)    :: Nz
        double precision, intent(in)    :: z_coor(Nz)
        double precision, intent(in)    :: z_min
        double precision, intent(in)    :: z_max
        double precision, intent(in)    :: w0
        double precision, intent(in)    :: Gamma
        double precision, intent(in)    :: eps_r
        double precision, intent(in)    :: dt

        integer :: kk, ii
        integer :: n_count

        this%w0    = w0
        this%Gamma = Gamma
        this%eps_r  = eps_r

        this%A1 = (2.0-Gamma*dt)/(2.0+Gamma*dt)
        this%A2 = eps0* w0**2 *dt/(2.0+Gamma*dt)
        this%C1 = (eps_r*eps0/dt-0.5*this%A2)/(eps_r*eps0/dt+0.5*this%A2)
        this%C3 = 1.0/(eps_r*eps0/dt+0.5*this%A2)
        this%C4 = 0.5*(this%A1+1.0)/(eps_r*eps0/dt+0.5*this%A2)

        n_count = 0
        do kk=1, Nz
            if ((z_coor(kk)>z_min) .and. (z_coor(kk)<z_max)) then
                n_count = n_count + 1
            end if
        end do

        this%n_tot = n_count

        if (.not. allocated(this%indx)) allocate(this%Px(n_count))
        if (.not. allocated(this%indx)) allocate(this%indx(n_count))

        this%Px   = 0.0d0
        this%indx = 0

        ii=1
        do kk=1, Nz
            if ((z_coor(kk)>z_min) .and. (z_coor(kk)<z_max)) then
                this%indx(ii) = kk
                ii = ii + 1
            end if
        end do

    end subroutine init_drude

    subroutine kill_drude(this)

        class(drude), intent(inout) :: this

        if (allocated(this%Px)) deallocate(this%Px)
        if (allocated(this%indx)) deallocate(this%indx)       

    end subroutine kill_drude

end module classical_medium_mod