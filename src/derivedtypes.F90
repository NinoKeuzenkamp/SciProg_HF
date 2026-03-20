module types

implicit none

private

public :: calculation_preset, energy_type


type calculation_preset
    logical :: UHF = .false.        ! RHF by default
    logical :: MP2 = .false.
    logical :: CCS = .false.
    logical :: CCSD= .false.
end type

! collect all energy terms in the same type
type energy_type
    real(8) :: HF       = 0
    real(8) :: HF_old   = 0         ! HF energy of previous cycle
    real(8) :: nuc      = 0
    real(8) :: MP2      = 0
end type

end module