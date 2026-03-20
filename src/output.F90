module output
use molecular_structure
use types
implicit none

private

public :: write_to_file

contains

    subroutine write_to_file(molecule, energy, outfile, converged, cycles, eps, n_occ)
        type(molecular_structure_t), intent(in) :: molecule
        type(energy_type),           intent(in) :: energy

        character(32), intent(in) :: outfile
        logical,       intent(in) :: converged
        integer,       intent(in) :: cycles     ! the amount of cycles that were used
        real(8),       intent(in) :: eps(:)
        integer,       intent(in) :: n_occ

        integer , parameter :: unit = 32
        integer :: i

        open(unit, file=outfile, action = "write")

        write(unit, "(a, /)") "RESULTS OF CALCULATION"

        ! write all coords and their charge
        write(unit, "(a)") "ATOM COORDINATES AND CHARGE:"
        do i = 1, molecule%num_atoms
            write(unit, "(3f10.5, f5.2)") molecule%coord(:,i) * 0.529177D0, molecule%charge(i)
        enddo

        ! start of energy results
        write(unit, "(/, a)") "HARTREE FOCK RESULTS (in Ha)"

        ! write all energies from the SCF cycle 
        write(unit, "(a)") "ENERGIES FROM SCF CYCLES:"
        do i = 1, cycles
            write(unit, "(f17.10, a, i4)") energy%all_HF(i) + energy%nuc, " CYCLE: ", i
        enddo

        if (converged) then
            write(unit, "(/, a, i4, a)") "PROGRAM CONVERGED AFTER ", cycles, " CYCLES"

            write(unit, "(a, t30, f17.10)") "ELECTRONIC ENERGY: ",         energy%HF - energy%nuc
            write(unit, "(a, t30, f17.10)") "NUCLEAR REPULSION ENERGY: ",  energy%nuc
            write(unit, "(a, t30, f17.10)") "TOTAL HARTREE FOCK ENERGY: ", energy%HF

            write(unit, "(/, a)") "EIGENVALUES OF FILLED ORBITALS:"
            do i = 1, n_occ
                write(unit, "(a, i3, f17.10)") "Molecular orbital: ", i, eps(i)
            enddo

            write(unit, "(/, a, f17.10)") "HOMO-LUMO GAP:", abs(eps(n_occ) - eps(n_occ + 1))
        else
            write(unit, "(a, i4, a)") "PROGRAM DID NOT CONVERGE AFTER ", cycles, " CYCLES"
        endif




        close(unit)
        
    end subroutine write_to_file

end module