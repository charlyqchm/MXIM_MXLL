module IO_mod

    use constants_mod
    use grid_mod

    implicit none

    character(len = 2) :: dftb_atom_type(10)
    character(len = 2) :: dftb_max_ang_orb(10)
    integer            :: dftb_n_mol
    integer            :: dftb_n_atoms
    integer            :: dftb_euler_steps
    integer            :: dftb_n_types
    logical            :: dftb_periodic = .false.
    logical            :: dftb_scc      = .true.
    logical            :: dftb_ion_dyn  = .false.
    real(dp)           :: dftb_scc_tol  = 1.0e-10_dp
    real(dp)           :: dftb_td       = 0.1e0_dp
    
    integer   :: mxll_Nz       = 1000
    integer   :: mxll_Nt       = 10 !fs
    integer   :: mxll_n_media  = 0
    integer   :: mxll_n_pml
    real(dp)  :: mxll_dz       =  1.0 ! informed in nm
    real(dp)  :: mxll_dt       =  0.1 ! informed in au
    real(dp)  :: mxll_Ex_src   = 0.001
    real(dp)  :: mxll_density  = 1.0e-5_dp
    real(dp)  :: mxll_z_src
    real(dp)  :: mxll_z_detect
    real(dp)  :: mxll_w_src
    real(dp)  :: mxll_tau_src  !fs
    real(dp)  :: mxll_t_print_big   = 1.0 ! fs
    real(dp)  :: mxll_t_print_small = 0.1 ! fs
    real(dp)  :: mxll_w_drude     = 1.0 ! eV
    real(dp)  :: mxll_gamma_drude = 0.1 ! eV
    real(dp)  :: mxll_ep_drude    = 1.0e0_dp
    real(dp)  :: mxll_media_center(10)
    real(dp)  :: mxll_media_rad(10)

contains

subroutine read_input_variables()

    integer :: ierr, funit

    namelist /MXLL_DFTB/ dftb_atom_type, dftb_max_ang_orb, dftb_periodic, dftb_scc, &
    dftb_ion_dyn, dftb_scc_tol, dftb_td, dftb_n_mol, dftb_n_atoms, dftb_n_types,    &
    dftb_euler_steps, mxll_Nz, mxll_Nt, mxll_n_media, mxll_dz, mxll_dt,             &
    mxll_density, mxll_z_src, mxll_w_src, mxll_tau_src, mxll_t_print_big,           &
    mxll_w_drude, mxll_gamma_drude, mxll_ep_drude, mxll_n_pml, mxll_z_detect,       &
    mxll_Ex_src, mxll_media_center, mxll_media_rad, mxll_t_print_small

    
    ! Check whether file exists.
    inquire (file="inp", iostat=ierr)
    
    if (ierr /= 0) then
        write (*, '("Error: input file inp does not exist")')
        return
    end if
    
    ! Open and read Namelist file.
    open (action='read', file="inp", iostat=ierr, newunit=funit)
    read (nml=MXLL_DFTB, iostat=ierr, unit=funit)
    if (ierr /= 0) write (*, '("Error: invalid Namelist format")')
    
    close (funit)

end subroutine read_input_variables

subroutine write_grid(mxll_grid, t_step, time)

    type(grid)      , intent(in) :: mxll_grid
    integer         , intent(in) :: t_step
    double precision, intent(in) :: time

    character(len=20) :: file_name
    character(len=20) :: file_number
    character(len=20) :: file_exten
    character(len=40) :: output_name
    
    integer :: n_dim

    n_dim = mxll_grid%Nz

    write(file_number, '(I7.7)') t_step
    file_exten = ".dat"

    file_name = 'Ex_'
    output_name = trim(file_name)//trim(file_number)//trim(file_exten)

    call write_E_field(mxll_grid%Ex, mxll_grid%E_coor, n_dim, time, output_name)

    file_name = 'Hy_'
    output_name = trim(file_name)//trim(file_number)//trim(file_exten)

    call write_H_field(mxll_grid%Hy, mxll_grid%H_coor, n_dim-1, time, output_name)

end subroutine write_grid

subroutine write_E_field(vec, coor, n_dim, time, file_name)
    
    integer           , intent(in) :: n_dim
    double precision  , intent(in) :: time
    double precision  , intent(in) :: vec(n_dim)
    double precision  , intent(in) :: coor(n_dim)   
    character(len=40) , intent(in) :: file_name

    integer :: io
    integer :: ii

    open(newunit=io, file=file_name)

    write(io, *) "#  time:", time, "a.u."
    write(io, *) "#  z                    Ex  "

    do ii=1, n_dim
        write(io, *) vec(ii)
    end do

    close(io)

end subroutine write_E_field

subroutine write_H_field(vec, coor, n_dim, time, file_name)
    
    integer           , intent(in) :: n_dim
    double precision  , intent(in) :: time
    double precision  , intent(in) :: vec(n_dim)
    double precision  , intent(in) :: coor(n_dim)   
    character(len=40) , intent(in) :: file_name

    integer :: io
    integer :: ii

    open(newunit=io, file=file_name)

    write(io, *) "#  time:", time, "a.u."
    write(io, *) "#  z                    Hy  "

    do ii=1, n_dim
        write(io, *) vec(ii)
    end do

    close(io)

end subroutine write_H_field

end module IO_mod