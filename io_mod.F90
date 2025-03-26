module IO_mod

    use constants_mod
    use grid_mod

    implicit none

    !Atomic species
    character(len = 2) :: dftb_atom_type(10)             

    !Atomic orbital with the greatest angular momentum per atom.
    character(len = 2) :: dftb_max_ang_orb(10)            
    
    !Number of atomic species
    integer            :: dftb_n_types                    
    
    !Number of molecules (always placed in the center of the box).
    integer            :: dftb_n_mol                      
    
    !Number of atoms per molecule.
    integer            :: dftb_n_atoms                    
    
    !Steps frequency to take an Euler step in electron and Ehrenfest dynamics.
    integer            :: dftb_euler_steps                
    
    !Boundary conditions for the DFTB system (currently .false. option is only available)
    logical            :: dftb_periodic = .false.
    
    !Option to include self-consistency in DFTB systems.
    logical            :: dftb_scc      = .true.

    !Option to run Ehrenfest dynamics.
    logical            :: dftb_ion_dyn  = .false.

    !Option to run Born-Oppenheimer dynamics.
    logical            :: dftb_BO_dyn   = .false.

    !Option to include magnetic fields in BO dynamics.
    logical            :: dftb_B_field  = .false.

    !Tolerance of the energy change during the DFTB ground state calculation.
    real(dp)           :: dftb_scc_tol  = 1.0e-10_dp

    !Suggested time step for the DFTB system. This value is modified later
    ! to be a multiple of the Maxwell time step.
    real(dp)           :: dftb_td       = 0.1e0_dp
   
    !Number of grid points in Maxwell box.
    integer   :: mxll_Nz       = 1000

    !Simulation time in fs.
    integer   :: mxll_Nt       = 10

    !Number of Drude media
    integer   :: mxll_n_media  = 0

    !Number of PML points.
    integer   :: mxll_n_pml    = 20

    !Number of external sources that are included by 'external_src.xxxx.dat' files.
    integer   :: mxll_ext_src  = 0

    !Grid points separation in nm.
    real(dp)  :: mxll_dz       =  1.0

    !Maxwell time step in atomic units. The value can not be bigger than the corresponding
    ! to the Courant condition which is set by default.
    real(dp)  :: mxll_dt       =  1.0E10

    !Number of DFTB systems per atomic unit of volume.
    real(dp)  :: mxll_density  = 1.0e-5_dp
    
    !Amplitude of pulse source.
    real(dp)  :: mxll_Ex_src   = 0.001

    !Position of the pulse inside the box in nm.
    real(dp)  :: mxll_z_src

    !Frequency of the pulse in eV.
    real(dp)  :: mxll_w_src
    
    !Broadening of the pulse in fs.
    real(dp)  :: mxll_tau_src  !fs

    !Electric and magnetic fields in the position -mxll_z_detect(nm)  and mxll_z_detect(nm)
    ! are printed in 'detector.dat' file.
    real(dp)  :: mxll_z_detect

    !Plasmon frequency in eV of the Drude medium.
    real(dp)  :: mxll_w_drude     = 1.0

    !Damping factor in eV of the Drude medium.
    real(dp)  :: mxll_gamma_drude = 0.1 ! eV

    !Drude medium relative permitivity.
    real(dp)  :: mxll_ep_drude    = 1.0e0_dp

    !Position of the Drude medium center.
    real(dp)  :: mxll_media_center(10)

    !Radius of the Drude medium.
    real(dp)  :: mxll_media_rad(10)

    !Average energy density is collected between -mxll_energy_dz (nm) and mxll_energy_dz (nm).
    real(dp)  :: mxll_energy_dz = 0.0d0

    !Printing frequency in fs of the 'detector.dat' file.
    real(dp)  :: mxll_t_print_small = 0.1

    !Printing frequency in fs of the whole E- and B-fields
    real(dp)  :: mxll_t_print_big   = 1.0 ! fs

contains

subroutine read_input_variables()

    integer :: ierr, funit

    namelist /MXLL_DFTB/ dftb_atom_type, dftb_max_ang_orb, dftb_periodic, dftb_scc, &
    dftb_ion_dyn, dftb_scc_tol, dftb_td, dftb_n_mol, dftb_n_atoms, dftb_n_types,    &
    dftb_euler_steps, dftb_BO_dyn, dftb_B_field,                                    &
    mxll_Nz, mxll_Nt, mxll_n_media, mxll_dz, mxll_dt,                               &
    mxll_density, mxll_z_src, mxll_w_src, mxll_tau_src, mxll_t_print_big,           &
    mxll_w_drude, mxll_gamma_drude, mxll_ep_drude, mxll_n_pml, mxll_z_detect,       &
    mxll_Ex_src, mxll_media_center, mxll_media_rad, mxll_t_print_small,             &
    mxll_ext_src, mxll_energy_dz

    
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