module types

implicit none
! by Nino Keuzenkamp

private

public :: calculation_preset, energy_type, atomic_orbital


! to store the user input & data about the calculation
! examples: user wants MP2. User wants UHF instead of RHF
type calculation_preset
    logical :: UHF = .false.        ! RHF by default
    logical :: MP2 = .false.
    logical :: CCS = .false.
    logical :: CCSD= .false.
    ! it may make sense to include n_cycles, out file path, tolerances in this type
    ! for tolerances, including them in this type would make it easier to get them from user input
end type

! collect all energy terms in the same type
type energy_type
    real(8)              :: HF         = 0  ! HF energy of current cycle
    real(8)              :: HF_old     = 0  ! HF energy of previous cycle
    real(8), allocatable :: all_SCF(:)      ! all ELECTRONIC energies from SCF cycle, will be printed to the output
    real(8)              :: nuc        = 0  ! nuclear repulsion energy
    real(8)              :: MP2        = 0  ! MP2 correction (correlation) energy
end type

! for a certain atom (defined by nuclear charge), assign all angular momentum to a exponent 
type atomic_orbital                                 ! example hydrogen with 3 s-functions:
    integer, allocatable :: angular(:)              ! angular   = (0,   0,   0)
    real(8), allocatable :: exponents(:)            ! exponents = (0.1, 1.0, 3.0)
end type

end module