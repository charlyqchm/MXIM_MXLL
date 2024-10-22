module q_medium_mod
    use constants_mod
    use, intrinsic :: iso_fortran_env, only : output_unit
    use dftbp_common_constants, only : AA__Bohr, V_m__au, eV__Hartree, fs__au, imag
    use dftbplus
    ! Only needed for the internal test system
    use testhelpers, only : writeAutotestTag

    implicit none
    
    type q_medium
        logical                            :: BO_dyn
        integer                            :: n_mol
        integer                            :: n_atoms
        double precision                   :: density
        character(len=2)    , allocatable  :: atom_names(:,:)
        integer             , allocatable  :: atom_type(:,:)
        integer             , allocatable  :: index(:)
        type(TDftbPlus)     , allocatable  :: dftbp(:)
        type(TDftbPlusInput), allocatable  :: input(:)
        double precision    , allocatable  :: dipole(:)
        double precision    , allocatable  :: dip_old(:)
        double precision    , allocatable  :: energy(:)
        real(dp)            , allocatable  :: coor(:,:,:)
        real(dp)            , allocatable  :: coor_old(:,:,:)
        real(dp)            , allocatable  :: coor_new(:,:,:)
        real(dp)            , allocatable  :: forces(:,:,:)
        real(dp)            , allocatable  :: velocities(:,:,:)
        real(dp)            , allocatable  :: at_masses(:,:)
        real(dp)            , allocatable  :: at_charges(:,:)

        type(fnode), pointer :: pRoot
        type(fnode), pointer :: pGeo 
        type(fnode), pointer :: pHam
        type(fnode), pointer :: pDftb
        type(fnode), pointer :: pMaxAng
        type(fnode), pointer :: pSlakos
        type(fnode), pointer :: pType2Files
        type(fnode), pointer :: pElecDyn
        type(fnode), pointer :: pElecStatic
        type(fnode), pointer :: pExternal
        type(fnode), pointer :: pPerturb
        type(fnode), pointer :: pLaser
        type(fnode), pointer :: pAnalysis

        contains
            procedure :: init_q_medium, kill_q_medium

    end type q_medium

contains

    subroutine init_q_medium(this, z_coor, len_mol, n_mol, n_atoms,n_type, density, &
                             td_step, n_steps, Nz, atom_type_list, max_ang_orb, scc,  &
                             scc_tol, periodic, ion_dyn, BO_dyn, euler_step)

        class(q_medium)    , intent(inout) :: this
        logical            , intent(in)    :: scc
        logical            , intent(in)    :: periodic
        logical            , intent(in)    :: ion_dyn
        logical            , intent(in)    :: BO_dyn
        character(len = 2) , intent(in)    :: atom_type_list(n_type)
        character(len = 2) , intent(in)    :: max_ang_orb(n_type)
        integer            , intent(in)    :: Nz
        integer            , intent(in)    :: n_steps
        integer            , intent(in)    :: n_mol
        integer            , intent(in)    :: n_type
        integer            , intent(in)    :: n_atoms
        integer            , intent(in)    :: euler_step
        double precision   , intent(in)    :: z_coor(Nz)
        double precision   , intent(in)    :: len_mol
        double precision   , intent(in)    :: td_step
        double precision   , intent(in)    :: density
        real(dp)           , intent(in)    :: scc_tol

        integer, parameter :: nAtom = 2

        real(dp) :: coor_aux(3,n_atoms)
        !Nitrogen atom types
        integer  :: species(n_atoms)
       
        character(len=20) :: file_name
        character(len=20) :: file_number
        character(len=20) :: file_exten
        character(len=20) :: dot
        character(len=40) :: input_name
        character(len=2)  :: at_name

        real(dp) :: merminEnergy

        integer :: kk, nn, ii
        integer :: n_atoms_aux
        integer :: io
        logical :: exists, atom_type_exists

        this%BO_dyn  = BO_dyn
        this%density = density
        this%n_atoms = n_atoms
        this%n_mol   = n_mol


        if (.not. allocated(this%index))       allocate(this%index(n_mol))
        if (.not. allocated(this%dftbp))       allocate(this%dftbp(n_mol))
        if (.not. allocated(this%input))       allocate(this%input(n_mol))
        if (.not. allocated(this%dipole))      allocate(this%dipole(n_mol))
        if (.not. allocated(this%energy))      allocate(this%energy(n_mol))
        if (.not. allocated(this%dip_old))     allocate(this%dip_old(n_mol))
        if (.not. allocated(this%coor))        allocate(this%coor(n_mol, n_atoms, 3))
        if (.not. allocated(this%atom_names))  allocate(this%atom_names(n_mol, n_atoms))
        if (.not. allocated(this%atom_type))   allocate(this%atom_type(n_mol, n_atoms))
        if (.not. allocated(this%at_charges))  allocate(this%at_charges(n_mol, n_atoms))

        if (this%BO_dyn) then
            if (.not. allocated(this%coor_old))   allocate(this%coor_old(n_mol, n_atoms, 3))
            if (.not. allocated(this%coor_new))   allocate(this%coor_new(n_mol, n_atoms, 3))
            if (.not. allocated(this%forces))     allocate(this%forces(n_mol, n_atoms, 3))
            if (.not. allocated(this%velocities)) allocate(this%velocities(n_mol, n_atoms, 3))
            if (.not. allocated(this%at_masses))  allocate(this%at_masses(n_mol, n_atoms))
            this%at_masses  = 0.0d0
            this%forces     = 0.0d0
            this%coor_old   = 0.0d0
            this%coor_new   = 0.0d0
            this%velocities = 0.0d0
        end if
        
        this%dipole     = 0.0d0
        this%dip_old    = 0.0d0
        this%energy     = 0.0d0
        this%at_charges = 0.0d0

        nn = 1
        do kk=1,Nz
            if((z_coor(kk)>(-len_mol/2.0) .and. nn<=n_mol)) then
                this%index(nn) = kk
                nn = nn+1
            end if
        enddo 

        !############## Reading xyz files ###############################################
        !we assume all the systems have the same composition

        dot = "."
        file_exten = ".xyz"
        do kk=1, n_mol
            
            write(file_number, '(I6.6)') kk
            file_name = 'molecule'
            input_name = trim(file_name)//trim(dot)//trim(file_number)//trim(file_exten)


            inquire(file=input_name, exist=exists)

            if (exists) then
            
                open(newunit=io, file=input_name, status="old")
                
                read(io, *) n_atoms_aux

                if (n_atoms_aux /= n_atoms) then
                    write(*,*) "Error. File   ", input_name, "does not have", n_atoms, "atoms."
                    stop   
                end if

                do nn=1, n_atoms
                    read(io, *) at_name, this%coor(kk,nn,1), this%coor(kk,nn,2), this%coor(kk,nn,3) 
                    this%atom_names(kk,nn) = at_name
                
                    atom_type_exists = .false.
                    do ii=1, n_type
                        if (at_name == atom_type_list(ii)) then
                            this%atom_type(kk,nn) = ii
                            atom_type_exists = .true.
                        end if
                    end do
                    
                    if (.not. atom_type_exists) then
                        write(*,*) "ERROR.", at_name, "is not an atom type."
                        stop
                    end if
                end do
            else
                write(*,*) "Error. File", file_name, "does not exist"
                stop 
            end if

            close(io)

        end do

        !################################################################################


        do nn=1, n_mol
            call TDftbPlus_init(this%dftbp(nn))
            call this%dftbp(nn)%getEmptyInput(this%input(nn))
            call this%input(nn)%getRootNode(this%pRoot)

            call setChild(this%pRoot, "Geometry", this%pGeo)
            call setChildValue(this%pGeo, "Periodic", periodic)
            call setChildValue(this%pGeo, "TypeNames", atom_type_list)
            coor_aux(:,:) = 0.0_dp
            species = this%atom_type(nn,:)
            call setChildValue(this%pGeo, "TypesAndCoordinates", reshape(species, [1, size(species)]), coor_aux)
            call setChild(this%pRoot, "Hamiltonian", this%pHam)
            call setChild(this%pHam, "Dftb", this%pDftb)
            call setChildValue(this%pDftb, "Scc", scc)
            call setChildValue(this%pDftb, "SccTolerance", scc_tol)

            call setChild(this%pDftb, "MaxAngularMomentum", this%pMaxAng)
            do kk=1, n_type
                call setChildValue(this%pMaxAng, atom_type_list(kk), max_ang_orb(kk))
            end do

            call setChild(this%pDftb, "SlaterKosterFiles", this%pSlakos)
            call setChild(this%pSlakos, "Type2FileNames", this%pType2Files)
            call setChildValue(this%pType2Files, "Prefix", "./")
            call setChildValue(this%pType2Files, "Separator", "-")
            call setChildValue(this%pType2Files, "Suffix", ".skf")

            if (this%BO_dyn) then
                call setChild(this%pRoot, "Analysis", this%pAnalysis)
                call setChildValue(this%pAnalysis, "PrintForces", .true.)

                call setChild(this%pDftb, "ElectricField", this%pElecStatic)
                call setChild(this%pElecStatic, "External", this%pExternal)
                call setChildValue(this%pExternal, "Strength", 0.0_dp)
                call setChildValue(this%pExternal, "Direction", [1.0_dp, 0.0_dp, 0.0_dp])
            else
                !  set up electron dynamics options
                call setChild(this%pRoot, "ElectronDynamics", this%pElecDyn)
                call setChildValue(this%pElecDyn, "Steps", n_steps)
                call setChildValue(this%pElecDyn, "TimeStep", td_step)
                call setChildValue(this%pElecDyn, "FieldStrength", 1.0_dp)
                call setChildValue(this%pElecDyn, "IonDynamics", ion_dyn)
                call setChildValue(this%pElecDyn, "InitialTemperature", 0.0_dp)
                call setChildValue(this%pElecDyn, "EulerFrequency", euler_step)

                call setChild(this%pElecDyn, "Perturbation", this%pPerturb)
                call setChild(this%pPerturb, "Laser", this%pLaser)
            ! these twovalues will be overriden
                call setChildValue(this%pLaser, "PolarisationDirection", [1.0_dp, 0.0_dp, 0.0_dp])
                call setChildValue(this%pLaser, "LaserEnergy", 1.0_dp)
            end if

        end do

        do nn=1, n_mol
            call this%dftbp(nn)%setupCalculator(this%input(nn))
            
            do kk=1, n_atoms
                coor_aux(:, kk)    = this%coor(nn,kk,:)*AA__Bohr
                this%coor(nn,kk,:) = this%coor(nn,kk,:)*AA__Bohr
            end do
            call this%dftbp(nn)%setGeometry(coor_aux)
            call this%dftbp(nn)%getEnergy(merminEnergy)
            if (this%BO_dyn) then
                call this%dftbp(nn)%getAtomicMasses(this%at_masses(nn,:))
            else
                call this%dftbp(nn)%initializeTimeProp(td_step, .true., .false.)
            end if
        end do

        this%coor_old = this%coor

    end subroutine init_q_medium

    subroutine kill_q_medium(this)
        
        class(q_medium), intent(inout) :: this
        integer :: ii

        do ii=1, this%n_mol
            call TDftbPlus_destruct(this%dftbp(ii))
        end do

        if (allocated(this%index))       deallocate(this%index)
        if (allocated(this%dftbp))       deallocate(this%dftbp)
        if (allocated(this%input))       deallocate(this%input)
        if (allocated(this%dipole))      deallocate(this%dipole)
        if (allocated(this%energy))      deallocate(this%energy)
        if (allocated(this%dip_old))     deallocate(this%dip_old)
        if (allocated(this%coor))        deallocate(this%coor)
        if (allocated(this%coor_old))    deallocate(this%coor_old)
        if (allocated(this%coor_new))    deallocate(this%coor_new)
        if (allocated(this%atom_names))  deallocate(this%atom_names)
        if (allocated(this%atom_type))   deallocate(this%atom_type)
        if (allocated(this%forces))      deallocate(this%forces)
        if (allocated(this%at_charges))  deallocate(this%at_charges)
        if (allocated(this%velocities))  deallocate(this%velocities)
        if (allocated(this%at_masses))   deallocate(this%at_masses)

    end subroutine kill_q_medium

end module q_medium_mod