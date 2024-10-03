program Maxwell_Maxim
    use constants_mod
    use grid_mod
    use q_medium_mod
    use td_propagator_mod
    use IO_mod
    use classical_medium_mod

    implicit none

    type(grid)       :: mxll_grid
    type(grid)       :: mxll_inc
    type(q_medium)   :: q_sys
    integer          :: Nz
    integer          :: Nt, Nt_fs
    integer          :: Nt_q
    integer          :: n_skip
    integer          :: tt
    integer          :: ii, kk
    integer          :: k_detect
    integer          :: k_source
    integer          :: out_unit = 111
    integer          :: t_print
    integer          :: t_print2
    integer          :: print_step
    integer          :: n_media
    integer          :: dip_unit
    double precision :: Ex_detect
    double precision :: Hy_detect
    double precision :: Ex_inc
    double precision :: Hy_inc
    double precision :: dz
    double precision :: dt
    double precision :: dt_q
    double precision :: dt_skip, dt_skipH
    double precision :: density
    double precision :: len_mol
    double precision :: z0
    double precision :: tau
    double precision :: omega = 17*ev_to_au! 17.44*ev_to_au H_2
    double precision :: E0
    double precision :: z_min, z_max
    double precision :: omegaD
    double precision :: GammaD
    double precision :: eps_r
    character(len=30):: output_file="output.dat"
    character(len=30):: dipole_file="dipole.dat"
    
    type(drude)     , allocatable :: media(:)
    type(drude)     , allocatable :: med_aux(:)
    double precision, allocatable :: z_coor(:)
    double precision, parameter, dimension (4) :: aBH=(/0.353222222d0,-0.488d0,0.145d0,-0.010222222d0/) !Blackman-Harris window for the pulse envelope


    call read_input_variables()

    Nz       = mxll_Nz
    dz       = mxll_dz * nm_to_au
    dt       = MIN(dz/(2.0*c0), mxll_dt)
    dt_q     = dftb_td
    len_mol  = dftb_n_mol * dz
    density  = mxll_density !matter per nm^3
    k_detect = int(mxll_Nz/2) + int(mxll_z_detect/mxll_dz)
    k_source = int(mxll_Nz/2) + int(mxll_z_src/mxll_dz)
    npml     = mxll_n_pml
    omega    = mxll_w_src * ev_to_au
    E0       = mxll_Ex_src
    tau      = mxll_tau_src * fs_to_au
    Nt_fs    = mxll_Nt
    Nt       = Nt_fs * int((1.0*fs_to_au)/(dt))
    t_print  = int(mxll_t_print_big*fs_to_au/dt) !1*int(1.0/(aut*dt))
    t_print2 = int(mxll_t_print_small*fs_to_au/dt)
    n_media  = mxll_n_media

    !Material data

    omegaD = mxll_w_drude*ev_to_au !7.039*ev_to_au
    GammaD = mxll_gamma_drude*ev_to_au !0.189*ev_to_au
    eps_r  = mxll_ep_drude    

    allocate(z_coor(Nz))
    
    if (n_media>0) allocate(media(n_media), med_aux(n_media))

    z0 = -dz*Nz*0.5

    do ii=1,Nz
        z_coor(ii)=z0+dz*(ii-1)
    enddo

    call mxll_grid%init_grid(Nz, npml, dz)
    call mxll_inc%init_grid(Nz, npml, dz)

    !Temporarly for two mirrors
    do ii=1, n_media
        z_min = (mxll_media_center(ii)-mxll_media_rad(ii))*nm_to_au
        z_max = (mxll_media_center(ii)+mxll_media_rad(ii))*nm_to_au
        call media(ii)%init_drude(z_coor, z_min,z_max, omegaD, GammaD, eps_r, dt, Nz)
        call med_aux(ii)%init_drude(z_coor, z_min,z_max, omegaD, GammaD, eps_r, dt, Nz)
    end do

    n_skip=1
    if(dt_q>dt)then
        n_skip=int(dt_q/dt)
    endif

    dt_skip=dble(n_skip)*dt
    dt_skipH=dt_skip/2.0d0

    Nt_q=1+(Nt/n_skip)

    call q_sys%init_q_medium(z_coor, len_mol, dftb_n_mol, dftb_n_atoms, dftb_n_types, density,  &
                             dt_skip, Nt_q, Nz, dftb_atom_type, dftb_max_ang_orb, dftb_scc,  &
                             dftb_scc_tol, dftb_periodic, dftb_ion_dyn, dftb_euler_steps)

    call init_td_propagator(Nz, dz, dt)



    open(unit=dip_unit, file=dipole_file, status='replace')

    write(dip_unit, '(a1,2x,a5,2x,100(a20,2x))') "#", "Time", "Dipoles"

    open(unit=out_unit, file=output_file, status='replace')

    write(out_unit, '(a1,2x,a5,2x,100(a20,2x))') "#", "Step", "Ex", "Hy", "Ex incident", "Hy incident"
    
    !#################################################################################
    !##### Here we init the propagation ##############################################
    !#################################################################################
    print_step = 0
    do tt=0, Nt
        
        time = (dble(tt)-0.5d0)*dt !average time (in s) during this iteration
        if (tt == 0) time = 0.0d0

        !temporary we define the pulse here#######

        if(time<=tau)then
            do kk = 2, Nz-1
                ! pulse(k)=E0*DEXP(-0.5*((t-0.5*tau)/(10*dt))**2)*DEXP(-0.5*(DBLE(k-k_source))**2)
                pulse(kk)=E0*cos(omega*time)*( &
                          aBH(1)+ &
                          aBH(2)*DCOS(2.0*pi0*time/tau)+ &
                          aBH(3)*DCOS(2.0*pi0*2.0*time/tau)+ &
                          aBH(4)*DCOS(2.0*pi0*3.0*time/tau)) * DEXP(-0.5*(DBLE(kk-k_source))**2)
            end do
        else
            pulse=0.0
        endif

        !#########################################

        if (n_media>0) then
            call td_propagate_field(mxll_grid, dt, dz, density, media, n_media)
            call td_propagate_field(mxll_inc, dt, dz, 0.0d0, med_aux, n_media)
        else
            call td_propagate_field(mxll_grid, dt, dz, density)
            call td_propagate_field(mxll_inc, dt, dz, 0.0d0)
        end if
        
        call td_propagate_q_medium(q_sys, mxll_grid, tt, n_skip, dt_skip, dip_unit)

        if (mod(tt, t_print) == 0) then
            call write_grid(mxll_grid, print_step, time)
            print_step = print_step + 1
        end if 

        Ex_detect = mxll_grid%Ex(k_detect)
        Hy_detect = 0.5*(mxll_grid%Hy(k_detect-1) + mxll_grid%Hy(k_detect))
        Ex_inc = mxll_inc%Ex(k_detect)
        Hy_inc = 0.5*(mxll_inc%Hy(k_detect-1) + mxll_inc%Hy(k_detect))

        if (mod(tt,t_print2) == 0) then
            write(out_unit, '(3x, 100(es20.12e3,2x))') time,  Ex_detect, Hy_detect, Ex_inc, Hy_inc
        end if
    end do

    write(*,*) "NUMBER OF MOLECULES:", q_sys%n_mol

    call mxll_grid%kill_grid()
    call mxll_inc%kill_grid()
    call q_sys%kill_q_medium()
    do ii =1, n_media
        call media(ii)%kill_drude()
        call med_aux(ii)%kill_drude()
    end do
    if (allocated(media))   deallocate(media)
    if (allocated(med_aux)) deallocate(med_aux)
    deallocate(z_coor)
    close(out_unit)
    close(dip_unit)

end program Maxwell_Maxim