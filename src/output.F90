module output
use molecular_structure
use types
implicit none

private

public :: write_to_file

contains

    subroutine write_to_file(molecule, energy, outfile, converged, cycles)
        type(molecular_structure_t), intent(in) :: molecule
        type(energy_type),           intent(in) :: energy

        character(32), intent(in) :: outfile
        logical,       intent(in) :: converged
        integer,       intent(in) :: cycles     ! the amount of cycles that were used

        integer , parameter :: unit = 32
        integer :: i

        open(unit, file=outfile, action = "write")

        write(unit, "(a, /)") "RESULTS OF CALCULATION"

        write(unit, *) "ATOM COORDINATES AND CHARGE"
        do i = 1, molecule%num_atoms
            write(unit, "(3f8.5, f5.2)") molecule%coord(:,i) * 0.529177D0, molecule%charge(i)
        enddo

        write(unit, "(/, a)") "HARTREE FOCK RESULTS (in Ha)"
        if (converged) then
            write(unit, "(a, i4, a)") "PROGRAM CONVERGED AFTER ", cycles, " CYCLES"
            write(unit, "(a, t30, f17.10)") "ELECTRONIC ENERGY: ",         energy%HF - energy%nuc
            write(unit, "(a, t30, f17.10)") "NUCLEAR REPULSION ENERGY: ",  energy%nuc
            write(unit, "(a, t30, f17.10)") "TOTAL HARTREE FOCK ENERGY: ",        energy%HF
        else
            write(unit, "(a, i4, a)") "PROGRAM DID NOT CONVERGE AFTER ", cycles, " CYCLES"
        endif




        close(unit)
        
    end subroutine write_to_file

end module