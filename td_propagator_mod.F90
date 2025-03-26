module td_propagator_mod
    use constants_mod
    use grid_mod
    use q_medium_mod
    use classical_medium_mod

    implicit none

    ! integer, parameter :: npml=29,m=3,ma=1
    ! integer, parameter :: npml=19,m=3,ma=1                          !thickness of the PML is 19 spatial steps on both sides of the grid
    integer, parameter :: m=3,ma=1  
    double precision, parameter :: alphaCPML=0.05,kappaCPML=5.0
    double precision :: sigmaCPML
    double precision :: tskip,dtskip,dtskipsec,dtskipH,dtskipHsec, time   !New TDSE time Steps
    double precision :: dt_eps0,dt_mu0, dt_epsM
    integer          :: npml
    integer          :: nskip,nts, tq_step_old
    
    
    double precision, allocatable :: den_ez(:)
    double precision, allocatable :: den_hz(:)
    double precision, allocatable :: pulse(:)
    double precision, allocatable :: be_z(:),ce_z(:),alphae_z(:),sige_z(:),kappae_z(:)
    double precision, allocatable :: bh_z(:),ch_z(:),alphah_z(:),sigh_z(:),kappah_z(:)
    
    interface td_propagate_field
        module procedure td_prop_WI_media, td_prop_WO_media
    end interface td_propagate_field

    contains

    subroutine init_td_propagator(Nz, dz, dt)

        integer         , intent(in)    :: Nz
        double precision, intent(in)    :: dz
        double precision, intent(inout) :: dt

        integer k, kk

        if (.not. allocated(den_ez))   allocate(den_ez(Nz))
        if (.not. allocated(den_hz))   allocate(den_hz(Nz))
        if (.not. allocated(pulse))    allocate(pulse(Nz))
        if (.not. allocated(be_z))     allocate(be_z(npml))
        if (.not. allocated(ce_z))     allocate(ce_z(npml))
        if (.not. allocated(alphae_z)) allocate(alphae_z(npml))
        if (.not. allocated(sige_z))   allocate(sige_z(npml))
        if (.not. allocated(kappae_z)) allocate(kappae_z(npml))
        if (.not. allocated(bh_z))     allocate(bh_z(npml))
        if (.not. allocated(ch_z))     allocate(ch_z(npml))
        if (.not. allocated(alphah_z)) allocate(alphah_z(npml))
        if (.not. allocated(sigh_z))   allocate(sigh_z(npml))
        if (.not. allocated(kappah_z)) allocate(kappah_z(npml))

        pulse = 0.0d0

        time=0.0
        tskip=0.0

        tq_step_old=0 

        dt_eps0=dt/eps0; dt_mu0=dt/mu0
        dt_epsM=dt_eps0

        sigmaCPML=0.8*(m+1)/(dz*(mu0/eps0)**0.5)
        do k=1,npml
            sige_z(k)=sigmaCPML*((npml-k)/(npml-1.0))**m
            alphae_z(k)=alphaCPML*((k-1)/(npml-1.0))**ma
            kappae_z(k)=1.0+(kappaCPML-1.0)*((npml-k)/(npml-1.0))**m
            be_z(k)=exp(-(sige_z(k)/kappae_z(k)+alphae_z(k))*dt/eps0)
            if ((sige_z(k)==0.0).and. &
                (alphae_z(k)==0.0).and. &
                (k==npml))then
                ce_z(k)=0.0
            else
                ce_z(k)=sige_z(k)*(be_z(k)-1.0)/(sige_z(k)+kappae_z(k)*alphae_z(k))/kappae_z(k)
            endif
        enddo

        do k=1,npml-1
            sigh_z(k)=sigmaCPML*((npml-k-0.5)/(npml-1.0))**m
            alphah_z(k)=alphaCPML*((k-0.5)/(npml-1.0))**ma
            kappah_z(k)=1.0+(kappaCPML-1.0)*((npml-k-0.5)/(npml-1.0))**m
            bh_z(k)=exp(-(sigh_z(k)/kappah_z(k)+alphah_z(k))*dt/eps0)
            ch_z(k)=sigh_z(k)*(bh_z(k)-1.0)/(sigh_z(k)+kappah_z(k)*alphah_z(k))/kappah_z(k)
        enddo

        kk=npml
        do k=1,Nz-1
            if(k<=npml)then
                den_ez(k)=1.0/(kappae_z(k)*dz)
            elseif(k>=(Nz+1-npml))then
                den_ez(k)=1.0/(kappae_z(kk)*dz)
                kk=kk-1
            else
                den_ez(k)=1.0/dz
            endif
        enddo

        kk=npml-1
        do k=1,Nz-1
            if(k<=(npml-1))then
                den_hz(k)=1.0/(kappah_z(k)*dz)
            elseif(k>=(Nz+1-npml))then
                den_hz(k)=1.0/(kappah_z(kk)*dz)
                kk=kk-1
            else
                den_hz(k)=1.0/dz
            endif
        enddo

    end subroutine init_td_propagator

    subroutine td_prop_WI_media(mxll_grid, dt, dz, density, media, n_media)

        type(grid)       , intent(inout) :: mxll_grid
        type(drude)      , intent(inout) :: media(n_media)
        integer          , intent(in)    :: n_media
        double precision , intent(in)    :: dt
        double precision , intent(in)    :: dz
        double precision , intent(in)    :: density
        
        integer            :: ii, kk
        integer            :: id
        integer            :: Nz
        double precision   :: Px_av
        double precision   :: cx
        double precision   :: tmpE, A1, A2, C1, C3, C4

        cx = density*dt/eps0
        Nz = mxll_grid%Nz

        !~~~~~~ Hy ~~~~~~~~!
        do ii=1,Nz-1
            mxll_grid%Hy(ii) = mxll_grid%Hy(ii)+dt_mu0*(mxll_grid%Ex(ii)-mxll_grid%Ex(ii+1))*den_hz(ii)
        enddo
        !  PML for the left side Hy
        do ii=1,npml-1
            mxll_grid%psi_Hyz_1(ii)    = bh_z(ii)*mxll_grid%psi_Hyz_1(ii)+ch_z(ii)*(mxll_grid%Ex(ii)-mxll_grid%Ex(ii+1))/dz
            mxll_grid%Hy(ii) = mxll_grid%Hy(ii)+dt_mu0*mxll_grid%psi_Hyz_1(ii)
        enddo
        !  PML for the right side Hy
        kk=npml-1
        do ii=Nz+1-npml,Nz-1
            mxll_grid%psi_Hyz_2(kk)    = bh_z(kk)*mxll_grid%psi_Hyz_2(kk)+ch_z(kk)*(mxll_grid%Ex(ii)-mxll_grid%Ex(ii+1))/dz
            mxll_grid%Hy(ii) = mxll_grid%Hy(ii)+dt_mu0*mxll_grid%psi_Hyz_2(kk)
            kk=kk-1
        enddo

        do ii=2,Nz-1
            Px_av                = mxll_grid%Jx(ii)+(time-tskip)*mxll_grid%dJx(ii)
            mxll_grid%Ex_new(ii) = mxll_grid%Ex(ii)&
                                 +dt_epsM*(mxll_grid%Hy(ii-1)-mxll_grid%Hy(ii))*den_ez(ii) &
                                 - cx*Px_av + pulse(ii)
        enddo

        do ii=1, n_media
        do kk=1, media(ii)%n_tot
            id   = media(ii)%indx(kk)
            A1   = media(ii)%A1; A2   = media(ii)%A2; C1   = media(ii)%C1
            C3   = media(ii)%C3; C4   = media(ii)%C4

            tmpE = C1*mxll_grid%Ex(id)+C3*(mxll_grid%Hy(id-1)-mxll_grid%Hy(id))/dz-C4*media(ii)%Px(kk)
            media(ii)%Px(kk)=A1*media(ii)%Px(kk)+A2*(tmpE+mxll_grid%Ex(id))
            mxll_grid%Ex_new(id)=tmpE
        end do
        end do

        mxll_grid%Ex = mxll_grid%Ex_new

        !  PML for the left side Ex
        do ii=2,npml
            mxll_grid%psi_Exz_1(ii)=be_z(ii)*mxll_grid%psi_Exz_1(ii)+ce_z(ii)*(mxll_grid%Hy(ii-1)-mxll_grid%Hy(ii))/dz
            mxll_grid%Ex(ii)=mxll_grid%Ex(ii)+dt_eps0*mxll_grid%psi_Exz_1(ii)
        enddo
        !  PML for right side Ex
        kk=npml
        do ii=Nz+1-npml,Nz-1
            mxll_grid%psi_Exz_2(kk)=be_z(kk)*mxll_grid%psi_Exz_2(kk)+ce_z(kk)*(mxll_grid%Hy(ii-1)-mxll_grid%Hy(ii))/dz
            mxll_grid%Ex(ii)=mxll_grid%Ex(ii)+dt_eps0*mxll_grid%psi_Exz_2(kk)
            kk=kk-1
        enddo

    end subroutine td_prop_WI_media

    subroutine td_prop_WO_media(mxll_grid, dt, dz, density)

        type(grid)       , intent(inout) :: mxll_grid
        double precision , intent(in)    :: dt
        double precision , intent(in)    :: dz
        double precision , intent(in)    :: density
        
        integer            :: ii, kk
        integer            :: id
        integer            :: Nz
        double precision   :: Px_av
        double precision   :: cx
        double precision   :: tmpE, A1, A2, C1, C3, C4

        cx = density*dt/eps0
        Nz = mxll_grid%Nz

        !~~~~~~ Hy ~~~~~~~~!
        do ii=1,Nz-1
            mxll_grid%Hy(ii) = mxll_grid%Hy(ii)+dt_mu0*(mxll_grid%Ex(ii)-mxll_grid%Ex(ii+1))*den_hz(ii)
        enddo
        !  PML for the left side Hy
        do ii=1,npml-1
            mxll_grid%psi_Hyz_1(ii)    = bh_z(ii)*mxll_grid%psi_Hyz_1(ii)+ch_z(ii)*(mxll_grid%Ex(ii)-mxll_grid%Ex(ii+1))/dz
            mxll_grid%Hy(ii) = mxll_grid%Hy(ii)+dt_mu0*mxll_grid%psi_Hyz_1(ii)
        enddo
        !  PML for the right side Hy
        kk=npml-1
        do ii=Nz+1-npml,Nz-1
            mxll_grid%psi_Hyz_2(kk)    = bh_z(kk)*mxll_grid%psi_Hyz_2(kk)+ch_z(kk)*(mxll_grid%Ex(ii)-mxll_grid%Ex(ii+1))/dz
            mxll_grid%Hy(ii) = mxll_grid%Hy(ii)+dt_mu0*mxll_grid%psi_Hyz_2(kk)
            kk=kk-1
        enddo

        do ii=2,Nz-1
            Px_av                = mxll_grid%Jx(ii)+(time-tskip)*mxll_grid%dJx(ii)
            mxll_grid%Ex_new(ii) = mxll_grid%Ex(ii)&
                                 +dt_epsM*(mxll_grid%Hy(ii-1)-mxll_grid%Hy(ii))*den_ez(ii) &
                                 - cx*Px_av + pulse(ii)
        enddo

        mxll_grid%Ex = mxll_grid%Ex_new

        !  PML for the left side Ex
        do ii=2,npml
            mxll_grid%psi_Exz_1(ii)=be_z(ii)*mxll_grid%psi_Exz_1(ii)+ce_z(ii)*(mxll_grid%Hy(ii-1)-mxll_grid%Hy(ii))/dz
            mxll_grid%Ex(ii)=mxll_grid%Ex(ii)+dt_eps0*mxll_grid%psi_Exz_1(ii)
        enddo
        !  PML for right side Ex
        kk=npml
        do ii=Nz+1-npml,Nz-1
            mxll_grid%psi_Exz_2(kk)=be_z(kk)*mxll_grid%psi_Exz_2(kk)+ce_z(kk)*(mxll_grid%Hy(ii-1)-mxll_grid%Hy(ii))/dz
            mxll_grid%Ex(ii)=mxll_grid%Ex(ii)+dt_eps0*mxll_grid%psi_Exz_2(kk)
            kk=kk-1
        enddo

    end subroutine td_prop_WO_media    

    subroutine td_propagate_q_medium(q_sys, mxll_grid, t_step, n_skip, dt_skip, q_sys_unit, &
                                     coor_charge_unit)

        type(q_medium)   , intent(inout) :: q_sys
        type(grid)       , intent(inout) :: mxll_grid
        integer          , intent(in)    :: t_step
        integer          , intent(in)    :: n_skip
        integer          , intent(in)    :: q_sys_unit
        integer          , intent(in)    :: coor_charge_unit
        double precision , intent(in)    :: dt_skip

        integer               :: tq_step
        integer               :: n_mol
        integer               :: n_at
        integer               :: ii, jj, indx
        real(dp)              :: vv
        real(dp)              :: By
        real(dp)              :: E_field
        real(dp)              :: dipole(3, 1)
        real(dp)              :: aux_field(3)
        real(dp)              :: energy
        real(dp)              :: dip_new, dip_old
        real(dp), allocatable :: atomNetCharges(:, :)
        real(dp), allocatable :: coor_aux(:,:)
        real(dp), allocatable :: forces_aux(:,:)

        !temporal to calculate angles
        real(dp) :: mass
        real(dp) :: F_CM(3)
        ! real(dp) :: r_vec(3)
        ! real(dp) :: cos_teta2(3)

        aux_field    = 0.0d0
        
        tq_step = 1+(t_step/n_skip)

        if(tq_step <= tq_step_old) return

        tskip = time
        tq_step_old = tq_step
        
        mxll_grid%Jx_old = mxll_grid%Jx
        n_mol = q_sys%n_mol
        n_at  = q_sys%n_atoms
        
        allocate(atomNetCharges(n_at, 1), coor_aux(3, n_at), forces_aux(3, n_at))        
        
        if (q_sys%BO_dyn) then

            q_sys%Kin        = 0.0d0
            q_sys%dE_t       = -q_sys%E_gs
            q_sys%dip_tot    = 0.0d0
            q_sys%coor_av    = 0.0d0
            q_sys%charges_av = 0.0d0
            ! cos_teta2 = 0.0

            do ii=1, n_mol
                indx         = q_sys%index(ii)
                E_field      = mxll_grid%Ex(indx)
                By           = 0.5/mu0 * (mxll_grid%Hy(indx-1) + mxll_grid%Hy(indx))
                aux_field(1) = 1.0d0
                dipole       = 0.0d0

                q_sys%dip_old(ii) = q_sys%dipole(ii)
                
                do jj=1, n_at
                    coor_aux(:, jj) = q_sys%coor(ii,jj,:)
                end do
                call q_sys%dftbp(ii)%setGeometry(coor_aux)
                call q_sys%dftbp(ii)%setExternalEfield(E_field, aux_field)
                call q_sys%dftbp(ii)%getGradients(forces_aux)
                call q_sys%dftbp(ii)%getEnergy(energy)
                call q_sys%dftbp(ii)%getGrossCharges(atomNetCharges(:,1)) 

                q_sys%dE_t = q_sys%dE_t + energy

                q_sys%dipole(ii) = 0.0d0

                do jj=1, n_at
                    q_sys%at_charges(ii,jj) = atomNetCharges(jj,1)
                    q_sys%forces(ii,jj,:)   = -forces_aux(:,jj)  
                    if (q_sys%B_field) then
                        q_sys%forces(ii,jj,1) = q_sys%forces(ii,jj,1) - q_sys%at_charges(ii,jj)*q_sys%vel(ii,jj,3)*By
                        q_sys%forces(ii,jj,3) = q_sys%forces(ii,jj,3) + q_sys%at_charges(ii,jj)*q_sys%vel(ii,jj,1)*By
                    end if
                end do
                
                F_CM = 0.0d0
                mass = 0.0d0
                do jj=1, n_at
                    F_CM = F_CM + q_sys%at_masses(ii, jj)*q_sys%forces(ii,jj,:)
                    mass = mass + q_sys%at_masses(ii, jj)
                end do
                F_CM = F_CM/mass

                do jj=1, n_at
                   
                    if (tq_step == 1) then
                        q_sys%coor_new(ii,jj,:) =  q_sys%coor(ii,jj,:) + q_sys%vel(ii,jj,:)*dt_skip + &
                        0.5*(q_sys%forces(ii,jj,:)-F_CM(:))/q_sys%at_masses(ii,jj)*dt_skip**2
                    else
                        q_sys%coor_new(ii,jj,:) =  2.0*q_sys%coor(ii,jj,:) - q_sys%coor_old(ii,jj,:) + &
                        (q_sys%forces(ii,jj,:)-F_CM(:))/q_sys%at_masses(ii,jj)*dt_skip**2
                    end if
                   
                    q_sys%dipole(ii) = q_sys%dipole(ii)  + q_sys%at_charges(ii,jj) * q_sys%coor_new(ii,jj,1)
                
                    q_sys%vel(ii,jj,1) = (q_sys%coor_new(ii,jj,1)-q_sys%coor(ii,jj,1))/dt_skip
                    q_sys%vel(ii,jj,2) = (q_sys%coor_new(ii,jj,2)-q_sys%coor(ii,jj,2))/dt_skip
                    q_sys%vel(ii,jj,3) = (q_sys%coor_new(ii,jj,3)-q_sys%coor(ii,jj,3))/dt_skip

                    vv = q_sys%vel(ii,jj,1)**2 + q_sys%vel(ii,jj,2)**2 + q_sys%vel(ii,jj,3)**2

                    q_sys%Kin = q_sys%Kin  + 0.5 * vv * q_sys%at_masses(ii,jj)
                        
                    q_sys%charges_av(jj) = q_sys%charges_av(jj) + q_sys%at_charges(ii,jj)

                    q_sys%coor_av(jj, 1) = q_sys%coor_av(jj, 1) + q_sys%coor_new(ii,jj,1)
                    q_sys%coor_av(jj, 2) = q_sys%coor_av(jj, 2) + q_sys%coor_new(ii,jj,2)
                    q_sys%coor_av(jj, 3) = q_sys%coor_av(jj, 3) + q_sys%coor_new(ii,jj,3)

                    
                end do 
                
                !temporal to calculate theta
                ! r_vec(1) =  q_sys%coor_new(ii, 1, 1) - q_sys%coor_new(ii, 2, 1) 
                ! r_vec(2) =  q_sys%coor_new(ii, 1, 2) - q_sys%coor_new(ii, 2, 2)
                ! r_vec(3) =  q_sys%coor_new(ii, 1, 3) - q_sys%coor_new(ii, 2, 3)

                ! cos_teta2(1) = cos_teta2(1) + r_vec(1)**2/(r_vec(1)**2+ r_vec(2)**2+ r_vec(3)**2)
                ! cos_teta2(2) = cos_teta2(2) + r_vec(2)**2/(r_vec(1)**2+ r_vec(2)**2+ r_vec(3)**2)
                ! cos_teta2(3) = cos_teta2(3) + r_vec(3)**2/(r_vec(1)**2+ r_vec(2)**2+ r_vec(3)**2)

                q_sys%dip_tot = q_sys%dip_tot + q_sys%dipole(ii)
                
                q_sys%coor_old(ii,:,:) = q_sys%coor(ii,:,:)
                q_sys%coor(ii,:,:)     = q_sys%coor_new(ii,:,:)

                mxll_grid%Jx(indx)=(q_sys%dipole(ii)-q_sys%dip_old(ii))/dt_skip
            
                if(tq_step==1) mxll_grid%Jx(indx)= 0.0d0

                mxll_grid%dJx(indx)=(mxll_grid%Jx(indx)-mxll_grid%Jx_old(indx))/dt_skip

            end do

        else

            q_sys%dE_t       = -q_sys%E_gs
            q_sys%dip_tot    = 0.0d0
!TODO: The same information should be printed for non-BO dynamics
            do ii=1, n_mol
                indx         = q_sys%index(ii)
                E_field      = mxll_grid%Ex(indx)
                aux_field(1) = E_field
                dipole       = 0.0d0

                q_sys%dip_old(ii) = q_sys%dipole(ii)

                call q_sys%dftbp(ii)%setTdElectricField(aux_field)
                call q_sys%dftbp(ii)%doOneTdStep(tq_step, dipole=dipole, energy=energy, &
                                                atomNetCharges=atomNetCharges        , &
                                                coord=coor_aux, force= forces_aux)
                q_sys%dipole(ii) = dipole(1, 1)
                q_sys%energy(ii) = energy
                q_sys%dE_t       = q_sys%dE_t + energy
                q_sys%dip_tot    = q_sys%dip_tot + q_sys%dipole(ii)
                do jj=1, n_at
                    q_sys%at_charges(ii,jj) =  atomNetCharges(jj,1)
                    q_sys%coor(ii,jj,:) = coor_aux(:, jj)
                end do

                mxll_grid%Jx(indx)=(q_sys%dipole(ii)-q_sys%dip_old(ii))/dt_skip
            
                if(tq_step==1) mxll_grid%Jx(indx)= 0.0d0

                mxll_grid%dJx(indx)=(mxll_grid%Jx(indx)-mxll_grid%Jx_old(indx))/dt_skip

            end do
        
        end if


        if (q_sys%BO_dyn) then
            write(q_sys_unit,*) time,  q_sys%dE_t/q_sys%n_mol, q_sys%Kin/q_sys%n_mol, &
                                q_sys%dip_tot/q_sys%n_mol
        else
            write(q_sys_unit,*) time,  q_sys%dE_t/q_sys%n_mol, q_sys%dip_tot/q_sys%n_mol
        end if

        
        if (q_sys%BO_dyn) then
            write(coor_charge_unit,*) q_sys%n_atoms
            write(coor_charge_unit,*) "time :", time
            do ii=1, q_sys%n_atoms

                write(coor_charge_unit,*) q_sys%atom_names(1,ii), &
                q_sys%coor_av(ii, 1)/AA__Bohr/q_sys%n_mol,        &
                q_sys%coor_av(ii, 2)/AA__Bohr/q_sys%n_mol,        &
                q_sys%coor_av(ii, 3)/AA__Bohr/q_sys%n_mol,        &
                q_sys%charges_av(ii)/q_sys%n_mol

            end do

            ! write(777, *) time, cos_teta2/q_sys%n_mol

        end if

        deallocate(atomNetCharges)
        deallocate(coor_aux)
        deallocate(forces_aux)    

    end subroutine td_propagate_q_medium

end module td_propagator_mod