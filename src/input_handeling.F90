module input_handeling

implicit none

private

public :: read_input_file, get_output_file

contains

    subroutine read_input_file(molecule, ao_basis, n_occ, n_cycles)
        use molecular_structure
        use ao_basis

        type(molecular_structure_t), intent(out) :: molecule
        type(basis_set_info_t),      intent(out) :: ao_basis
        integer, intent(out) :: n_occ, n_cycles
        
        character(32)      :: inputfile         ! user input file path
        integer, parameter :: unit = 21
        integer            :: stat

        integer :: n_atoms, total_charge, n_electron, i     ! total charge of SYSTEM
        real(8), allocatable :: coords(:, :), charge(:)     ! coords and charge of atoms to add to molecule-type

        ! used to read ao info from file
        integer :: atom_index, angular
        real(8) :: coefficient


        ! get user inputfile
        print "(/, a)", "Please type the path to your input txt file. (max 32 characters)"
        read "(a32)", inputfile


        open(unit, file=inputfile)
        ! FORMAT OF FILE:  (contracted orbitals NOT implemented)
        ! n_atoms total_charge_SYSTEM n_cycles
        ! x y z charge (for each atom in angstrom)
        ! atom_index angular_momentum coefficient (for each orbital)


        ! read number of atoms and number of orbitals and allocate
        read(unit, *) n_atoms, total_charge, n_cycles
        allocate(coords(3, n_atoms), charge(n_atoms))


        ! FOR MOLECULE
        ! read atom data into arrays
        do i = 1, n_atoms
            read(unit, *) coords(1, i), coords(2, i), coords(3, i), charge(i)
        enddo
        ! from angstrom to Bohr, 1 A =  0.529177 Bohr
        coords = coords * 1.889727D0    

        ! turn arrays into molecule type
        call add_atoms_to_molecule(molecule, charge, coords)


        ! FOR ATOMIC ORBITAL BASIS
        ! per orbital, add it to the ao_basis (contracted orbitals NOT implemented)
        do
            read(unit, *, iostat=stat) atom_index, angular, coefficient

            ! EOF reached
            if (stat /= 0) exit

            ! add ao to ao_basis
            call add_shell_to_basis(ao_basis,                   &
                                    angular,                    &
                                    coords(:, atom_index),      &   ! array of x y z of atom numbered "atom_index"
                                    coefficient)    
            
        enddo
        close(unit)


        ! calculate number of electrons
        n_electron = int(sum(charge)) - total_charge

        ! make sure n_electron is always a whole number
        if (mod(n_electron, 2) == 1) then   ! add implementation for unrestricted HF here
            stop "The number of electrons cannot be odd"
        else    ! -> restricted HF
            n_occ = n_electron/2
        endif
        
    end subroutine read_input_file

    ! ask where the user wants to place the results
    subroutine get_output_file(outfile)
        character(32), intent(out) :: outfile

        print "(/, a)", "Please input the path for your output txt file. (max 32 characters)"
        read "(a32)", outfile
    end subroutine get_output_file

end module