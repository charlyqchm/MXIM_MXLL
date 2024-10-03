module grid_mod

    type grid
        integer          :: Nz
        double precision, allocatable :: psi_Exz_1(:)
        double precision, allocatable :: psi_Exz_2(:)
        double precision, allocatable :: psi_Hyz_1(:)
        double precision, allocatable :: psi_Hyz_2(:)
        double precision, allocatable :: Ex(:)
        double precision, allocatable :: Ex_new(:)
        double precision, allocatable :: Hy(:)
        double precision, allocatable :: Jx(:)
        double precision, allocatable :: Jx_old(:)
        double precision, allocatable :: dJx(:)
        double precision, allocatable :: E_coor(:)
        double precision, allocatable :: H_coor(:)

        contains
            procedure :: init_grid, kill_grid

    end type grid

contains
    subroutine init_grid(this, Nz, npml, dz)
        class(grid)     , intent(inout) :: this
        integer         , intent(in)    :: Nz
        integer         , intent(in)    :: npml
        double precision, intent(in)    :: dz

        integer           :: ii
        double precision  :: z0

        this%Nz = Nz

        if (.not. allocated(this%Ex))        allocate(this%Ex(Nz))
        if (.not. allocated(this%Ex_new))    allocate(this%Ex_new(Nz))
        if (.not. allocated(this%Hy))        allocate(this%Hy(Nz-1))
        if (.not. allocated(this%Jx))        allocate(this%Jx(Nz))
        if (.not. allocated(this%Jx_old))    allocate(this%Jx_old(Nz))
        if (.not. allocated(this%dJx))       allocate(this%dJx(Nz))
        if (.not. allocated(this%E_coor))    allocate(this%E_coor(Nz))
        if (.not. allocated(this%H_coor))    allocate(this%H_coor(Nz-1))
        if (.not. allocated(this%psi_Exz_1)) allocate(this%psi_Exz_1(npml))
        if (.not. allocated(this%psi_Exz_2)) allocate(this%psi_Exz_2(npml))
        if (.not. allocated(this%psi_Hyz_1)) allocate(this%psi_Hyz_1(npml-1))
        if (.not. allocated(this%psi_Hyz_2)) allocate(this%psi_Hyz_2(npml-1))

        this%Ex        = 0.0d0
        this%Ex_new    = 0.0d0       
        this%Hy        = 0.0d0
        this%Jx        = 0.0d0
        this%Jx_old    = 0.0d0
        this%dJx       = 0.0d0
        this%psi_Exz_1 = 0.0d0
        this%psi_Exz_2 = 0.0d0
        this%psi_Hyz_1 = 0.0d0
        this%psi_Hyz_2 = 0.0d0
    
        z0 = -dz*Nz*0.5

        do ii=1,Nz
            this%E_coor(ii)=z0+dz*(ii-1)
        enddo
        do ii=1,Nz-1
            this%H_coor(ii)=z0+dz*0.5+dz*(ii-1)
        enddo

    end subroutine init_grid

    subroutine kill_grid(this)
        class(grid), intent(inout) :: this
        
        if (allocated(this%Ex))        deallocate(this%Ex)
        if (allocated(this%Ex_new))    deallocate(this%Ex_new)
        if (allocated(this%Hy))        deallocate(this%Hy)
        if (allocated(this%Jx))        deallocate(this%Jx)
        if (allocated(this%Jx_old))    deallocate(this%Jx_old)
        if (allocated(this%dJx))       deallocate(this%dJx)
        if (allocated(this%psi_Exz_1)) deallocate(this%psi_Exz_1)
        if (allocated(this%psi_Exz_2)) deallocate(this%psi_Exz_2)
        if (allocated(this%psi_Hyz_1)) deallocate(this%psi_Hyz_1)
        if (allocated(this%psi_Hyz_2)) deallocate(this%psi_Hyz_2)
        if (allocated(this%E_coor))    deallocate(this%E_coor)
        if (allocated(this%H_coor))    deallocate(this%H_coor)

    end subroutine kill_grid

end module grid_mod