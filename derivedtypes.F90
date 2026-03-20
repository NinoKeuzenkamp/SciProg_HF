module types

implicit none

private

public :: calculation_preset


type calculation_preset
    logical :: UHF = .false.        ! RHF by default
    logical :: MP2 = .false.
    logical :: CCS = .false.
    logical :: CCSD= .false.

end type

end module