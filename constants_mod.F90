module constants_mod

    implicit none

    !costants

    double precision, parameter :: pi0   = 3.141592653589793d0
    double precision, parameter :: c0    = 137.0                 !299792458.0d0
    double precision, parameter :: sqrt2 = 1.414213562373095d0
    double precision, parameter :: sqrt3 = 1.732050807568877d0
    double precision, parameter :: sqrt6 = 2.449489742783178d0
    double precision, parameter :: eps0  = 1.0/(4.0*pi0)          !1.0d0/(c*c*mu0)
    double precision, parameter :: mu0   = 1.0/(eps0*c0**2)      !4.0d-7*pi
    double precision, parameter :: hbar  = 1.0                   !1.054571628d-34
    double complex  , parameter :: Z_I   = (0.0d0, 1.0d0)
    double complex  , parameter :: Z_0   = (0.0d0, 0.0d0)
    double complex  , parameter :: Z_ONE = (1.0d0, 0.0d0)

    !unit convertions

    double precision, parameter :: nm_to_au        = 18.897261259077823
    double precision, parameter :: ev_to_radsec    = 2.0*pi0*2.418d14
    double precision, parameter :: ev_to_au        = 1.0/27.2114
    double precision, parameter :: Hz_to_ev        = 4.1356691d-15
    double precision, parameter :: Debye_to_Cm     = 3.33564d-30
    double precision, parameter :: Debye_to_au     = 3.0/7.63
    double precision, parameter :: fs_to_au        = 1.0d0/2.4188843265864D-2  !0.0241888432650516d0
    double precision, parameter :: dipole_au_to_SI = 8.47835281d-30
    integer         , parameter :: dp = kind(1.0d0)

end module constants_mod