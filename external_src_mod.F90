module external_src_mod

    use constants_mod

    implicit none

    type external_src

        integer                       :: n
        integer                       :: i0
        double precision              :: z0
        double precision, allocatable :: Ex_t(:)

        contains
            procedure :: init_external_src, kill_external_src

    end type external_src

contains

subroutine init_external_src(this, Nt, Nz, dz, z)

    class(external_src), intent(inout) :: this
    integer            , intent(in)    :: Nt 
    integer            , intent(in)    :: Nz
    double precision   , intent(in)    :: dz
    double precision   , intent(in)    :: z

    this%n  = Nt
    this%z0 = z
    this%i0 =  int(z/dz)  + int(Nz/2) + 1

    if (.not. allocated(this%Ex_t))  allocate(this%Ex_t(this%n))

    this%Ex_t = 0.0d0

end subroutine init_external_src

subroutine kill_external_src(this)
    
    class(external_src), intent(inout) :: this

    if (allocated(this%Ex_t))  deallocate(this%Ex_t)

end subroutine kill_external_src

subroutine read_and_init_ext_src(src, n_src, Nt, Nz, dz)

    type(external_src), allocatable, intent(inout) :: src(:)
    integer            , intent(in) :: n_src
    integer            , intent(in) :: Nt 
    integer            , intent(in) :: Nz
    double precision   , intent(in) :: dz

    logical           :: exists
    integer           :: nn, tt
    integer           :: io
    character(len=20) :: file_name
    character(len=20) :: file_number
    character(len=20) :: file_exten
    character(len=40) :: input_name
    double precision  :: z
    double precision  :: Ex

    if (n_src == 0) return

    if (.not. allocated(src)) allocate(src(n_src))

    file_exten = ".dat"

    do nn = 1, n_src
    
        write(file_number, '(I4.4)') nn
        file_name = 'external_src.'
        input_name = trim(file_name)//trim(file_number)//trim(file_exten) 

        inquire(file=input_name, exist=exists)

        if (.not. exists) then
            write(*,*) "Error. File   ", input_name, "does not exist."
            stop 
        end if

        open(newunit=io, file=input_name, status="old")

        read(io, *) z

        z = z * nm_to_au

        call src(nn)%init_external_src(Nt,Nz,dz,z)

        read(io, *)

        print*, Nt

        do tt=1, Nt
            read(io, *) Ex
            src(nn)%Ex_t(tt) = Ex
        end do 

        close(io)

    end do

end subroutine read_and_init_ext_src

end module external_src_mod