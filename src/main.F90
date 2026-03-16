program HartreeFock

   ! Demonstration program that can be used as a starting point
   ! Lucas Visscher, March 2022

   use molecular_structure
   use ao_basis
   use compute_integrals
   use diagonalization
   use fock_implementation

     implicit none

     ! Variable containing the molecular structure
     type(molecular_structure_t) :: molecule
     ! Variable containing the atomic orbital basis
     type(basis_set_info_t) :: ao_basis

     ! Variable naming as in the description of the exercise
     integer  :: n_AO, n_occ
     integer  :: kappa, lambda, i
     real(8)  :: E_HF
     real(8), allocatable :: F(:,:),V(:,:),T(:,:),S(:,:), C(:,:), eps(:), D(:,:)

     integer :: n_cycles
     real(8) :: convergence
     real(8), allocatable :: hcore(:,:)

     ! The following large array can be eliminated when Fock matrix contruction is implemented
     real(8), allocatable :: ao_integrals (:,:,:,:)
   
     ! Definition of the molecule
     call define_molecule(molecule)

     ! Definition of the GTOs
     call define_basis(ao_basis)
     n_AO = ao_basis%nao
   
     ! Definition of the number of occupied orbitals
     n_occ = 3 ! hardwired for this demonstration program, should be set via input

     ! Compute the overlap matrix
     allocate (S(n_AO,n_AO))
     call   compute_1e_integrals ("OVL",ao_basis,ao_basis,S)

     ! Compute the kinetic matrix
     allocate (T(n_AO,n_AO))
     call   compute_1e_integrals ("KIN",ao_basis,ao_basis,T)

     ! Compute the potential matrix
     allocate (V(n_AO,n_AO))
     call   compute_1e_integrals ("POT",ao_basis,ao_basis,V,molecule)

     ! Compute the core Hamiltonian matrix (the potential is positive, we scale with -e = -1 to get to the potential energy matrix)
     allocate (F(n_AO,n_AO))
     allocate (hcore(n_AO,n_AO))
     hcore = T - V

     ! Diagonalize the Fock matrix
     allocate (C(n_AO,n_AO))
     allocate (eps(n_AO))
     call solve_genev (hcore,S,C,eps)
     print*, "Orbital energies for the core Hamiltonian:",eps

     ! Form the density matrix
     allocate (D(n_AO,n_AO))
     call density_matrix(D, C, n_ao, n_occ)

     
     ! Compute all 2-electron integrals
     allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))
     call generate_2int (ao_basis,ao_integrals)

    ! SCF loop
    n_cycles = 50
    convergence = 1.D-9
    do i = 1, n_cycles

      ! calculate Fock Matrix
      call fock_matrix(F, D, hcore, ao_integrals, n_ao)

      ! diagonalise fock matrix
      call solve_genev (F,S,C,eps)
      print*, "Orbital energies for the core Hamiltonian:",eps, "Cycle: ", i

      ! Form the density matrix
      call density_matrix(D, C, n_ao, n_occ)

      ! calculate HF energy for current cycle
      E_HF = sum((hcore + F) * D)
      print*, "The Hartree-Fock energy:    ", E_HF, "Cycle: ", i

      if (sqrt(sum((matmul(F,D) - matmul(D,F))**2)) < convergence) then
        print*, "The program has converged at cycle ", i
        stop
      endif

    enddo

    print *, "The program did not converged after ", n_cycles, " cycles"
  end


   subroutine define_molecule(molecule)
     ! This routine should be improved such that an arbitrary molecule can be given as input
     ! the coordinates below are for a be-he dimer oriented along the x-axis with a bond length of 2 au
     use molecular_structure
     type(molecular_structure_t), intent(inout) :: molecule
     real(8) :: charge(2),coord(3,2)
     charge(1)   = 4.D0
     charge(2)   = 2.D0
     coord       = 0.D0
     coord(1,2)  = 2.D0
     call add_atoms_to_molecule(molecule,charge,coord)
   end subroutine

   subroutine define_basis(ao_basis)
    ! This routine can be extended to use better basis sets 
    ! The coordinates of the shell centers are the nuclear coordinates
    ! Think of a refactoring of define_molecule and define_basis to ensure consistency 
     use ao_basis
     type(basis_set_info_t), intent(inout) :: ao_basis
     type(basis_func_info_t) :: gto
     ! Be:  2 uncontracted s-funs:    l      coord          exp      
     call add_shell_to_basis(ao_basis,0,(/0.D0,0.D0,0.D0/),4.D0)
     call add_shell_to_basis(ao_basis,0,(/0.D0,0.D0,0.D0/),1.D0)
     call add_shell_to_basis(ao_basis,0,(/0.D0,0.D0,0.D0/),2.D0)
     ! He:  1 uncontracted s-fun:     l      coord          exp      
     call add_shell_to_basis(ao_basis,0,(/2.D0,0.D0,0.D0/),1.D0)
     call add_shell_to_basis(ao_basis,0,(/2.D0,0.D0,0.D0/),2.D0)
   end subroutine

   
