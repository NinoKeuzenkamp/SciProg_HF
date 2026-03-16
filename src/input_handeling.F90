module input_handeling

implicit none

private

public :: read_input

contains

    subroutine read_input(molecule, ao_basis, n_occ)
        use molecular_structure
        use ao_basis

        type(molecular_structure_t), intent(out) :: molecule
        type(basis_set_info_t),      intent(out) :: ao_basis
        integer, intent(out) :: n_occ
        
        character(32) :: inputfile
        integer, parameter :: unit = 21

        integer :: n_atoms, n_orbs, i
        real(8), allocatable :: coords(:, :), charge(:)

        ! for orbital adding loop
        integer :: atom_index, angular
        real(8) :: coefficient

        ! get user inputfile
        print *, "Please type the path to your input txt file. (max 32 characters)"
        read "(a32)", inputfile

        ! open file
        open(unit, file=inputfile)
        ! format file:  (contracted orbitals NOT implemented)
        ! n_atoms n_orbs
        ! x y z charge (for each atom)
        ! atom_index angular_momentum coefficient (for each orbital)

        ! read number of atoms and number of orbitals and allocate
        read(unit, *) n_atoms, n_orbs
        allocate(coords(3, n_atoms), charge(n_atoms))

        ! read atom data into arrays
        do i = 1, n_atoms
            read(unit, *) coords(1, i), coords(2, i), coords(3, i), charge(i)
        enddo

        ! turn arrays into molecule type
        call add_atoms_to_molecule(molecule, charge, coords)

        ! per orbital, add it to the ao_basis (contracted orbitals NOT implemented)
        do i = 1, n_orbs
            read(unit, *) atom_index, angular, coefficient

            ! add to ao_basis
            call add_shell_to_basis(ao_basis,                   &
                                    angular,                    &
                                    coords(:, atom_index),      &   ! array of x y z of atom numbered "atom_index"
                                    coefficient)                
        enddo

        close(unit)
        ! calculate n_occ
        n_occ = int(sum(charge)/2.D0)
    end subroutine

end module