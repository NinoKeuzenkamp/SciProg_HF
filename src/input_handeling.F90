module input_handeling

implicit none
! by Nino Keuzenkamp

private

public :: read_input, get_output_file

contains

    ! read input from userfiles
    subroutine read_input(molecule, ao_basis, n_occ, n_cycles, preset)
        use molecular_structure
        use ao_basis
        use types

        type(molecular_structure_t), intent(out) :: molecule
        type(basis_set_info_t),      intent(out) :: ao_basis
        integer,                     intent(out) :: n_occ, n_cycles
        type(calculation_preset),    intent(out) :: preset

        integer :: n_electron
        character:: bool


        ! gets the following variables
        call read_molecule(molecule, n_cycles, n_electron)

        ! gets the ao_basis variable
        call read_basis(ao_basis, molecule)

        ! does the user want MP2?
        print "(/, a)", "Do you want to include MP2? (y/n)"
        read "(a1)", bool

        if (bool == "y") then
            preset%MP2 = .true.
        elseif (bool == "n") then
            ! false by default
            ! preset%MP2 = .false
        else 
            print "(a)", "The input was neither y nor n, the program will terminate"
            stop
        endif

        ! make sure n_electron is always a whole number
        if (mod(n_electron, 2) == 1) then   ! add implementation for unrestricted HF here
            stop "The number of electrons cannot be odd"
        else    ! -> restricted HF
            n_occ = n_electron/2
        endif
        
    end subroutine read_input


    ! FOR MOLECULE COORDINATES
    subroutine read_molecule(molecule, n_cycles, n_electron)
        use molecular_structure
        type(molecular_structure_t), intent(out) :: molecule
        integer,                     intent(out) :: n_cycles, n_electron

        integer :: n_atoms, total_charge_sys, i              
        real(8), allocatable :: coords(:, :), charge(:)     ! coords and charge of atoms to add to molecule-type

        character(32)      :: inputfile         ! user input file path, see README.txt to see format of file
        integer, parameter :: unit = 21


        ! get user inputfile for molecule
        print "(/, a)", "Please type the path to your input txt file. (max 32 characters)"
        read "(a32)", inputfile
        open(unit, file=inputfile)

        ! read number of atoms and allocate
        read(unit, *) n_atoms, total_charge_sys, n_cycles
        allocate(coords(3, n_atoms), charge(n_atoms))

        ! read atom data into arrays
        do i = 1, n_atoms
            read(unit, *) coords(1, i), coords(2, i), coords(3, i), charge(i)
        enddo

        ! from angstrom to Bohr, 1 A =  0.529177 Bohr
        coords = coords * 1.889727D0    

        ! turn arrays into molecule type
        call add_atoms_to_molecule(molecule, charge, coords)
        close(unit)

        ! calculate number of electrons
        n_electron = int(sum(charge)) - total_charge_sys

    end subroutine read_molecule


    ! FOR ATOMIC ORBITAL BASIS
    subroutine read_basis(ao_basis, molecule)
        use molecular_structure
        use ao_basis
        use types
        type(molecular_structure_t), intent(in) :: molecule
        type(basis_set_info_t),      intent(out):: ao_basis

        type(atomic_orbital),   allocatable :: AO(:)          ! store orbital exponents and angular momenta for each atom
        
        integer :: largest_charge    ! largest charge = largest atom defined in basis
        integer :: n_ao              ! number of AO's per atom

        integer :: i, j, charge
        character(32)      :: inputfile         ! user input file path, see README.txt to see format of file
        integer, parameter :: unit = 21


        ! get user inputfile
        print "(/, a)", "Please type the path to your basis txt file. (max 32 characters)"
        read "(a32)", inputfile
        open(unit, file=inputfile)

        ! create an AO-type for each atom defined in basis
        read(unit, *) largest_charge
        allocate(AO(largest_charge))

        ! per atom defined in basis
        do i = 1, largest_charge

            ! read how many AO's are present for this atom
            read(unit, *) n_ao
            allocate(AO(i)%angular(n_ao), AO(i)%exponents(n_ao))

            ! for each of those AO's
            do j = 1, n_ao
                
                ! place AO data in AO-type
                read(unit, *) AO(i)%angular(j), AO(i)%exponents(j)
            enddo
        enddo

        close(unit)


        ! add atomic basis functions to "ao_basis" for each atom in the molecule
        ! i.e. if there is a carbon atom in the molecule, then add the carbon AO's to that position
        do i = 1, molecule%num_atoms
            ! shorter notation, "charge of atom to which the AO's need to be assigned"
            charge = int(molecule%charge(i))

            ! add each AO to the ao_basis
            do j = 1, size(AO(charge)%angular)
                call add_shell_to_basis(ao_basis,               &
                                        AO(charge)%angular(j),  &
                                        molecule%coord(:, i),   &
                                        AO(charge)%exponents(j))
            enddo
        enddo

        
    end subroutine read_basis



    ! ask where the user wants to place the results
    subroutine get_output_file(outfile)
        character(32), intent(out) :: outfile

        print "(/, a)", "Please input the path for your output txt file. (max 32 characters)"
        print "(a)", "Note that the program cannot make new folders."
        read "(a32)", outfile
    end subroutine get_output_file

end module