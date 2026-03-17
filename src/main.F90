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
     real(8)  :: E_HF_old, E_HF_new
     real(8), allocatable :: F(:,:),V(:,:),T(:,:),S(:,:), C(:,:), eps(:)

     integer :: n_cycles
     real(8) :: tolerance_E, tolerance_D
     real(8), allocatable :: hcore(:,:), D_old(:,:), D_new(:,:)
     logical :: converged

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

     ! Form the density matrix
     allocate (D_new(n_AO,n_AO), D_old(n_AO,n_AO))
     call density_matrix(D_old, C, n_ao, n_occ)

     ! calculate energy of core only system
     E_HF_old = sum((hcore + F) * D_old)
     
     ! Compute all 2-electron integrals
     allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))
     call generate_2int (ao_basis,ao_integrals)


    ! set tolerances for convergence
    tolerance_E = 1.D-9
    tolerance_D = 1.D-3 
    ! SCF loop
    n_cycles = 50
    i = 1
    do while (i <= n_cycles)

      ! calculate Fock Matrix
      call fock_matrix(F, D_old, hcore, ao_integrals, n_ao)

      ! diagonalise fock matrix
      call solve_genev (F,S,C,eps)
      ! print*, "Orbital energies for the core Hamiltonian:",eps, "Cycle: ", i

      ! Form the density matrix
      call density_matrix(D_new, C, n_ao, n_occ)

      ! calculate HF energy for current cycle
      E_HF_new = sum((hcore + F) * D_new)
      ! print*, "The Hartree-Fock energy:    ", E_HF_new, "Cycle: ", i

      ! if converged, exit SCF loop
      converged = convergence_check(E_HF_new, E_HF_old, D_new, D_old, tolerance_E, tolerance_D)
      if (converged) exit

      E_HF_old = E_HF_new
      D_old    = D_new

      i = i + 1
    enddo

    if (converged) then
      print "(a, i3, a)", "The program has converged after ", i, " cycles"
      print "(a, f10.5, a)", "The energy has converged to ", E_HF_new, " Hartree"
    else
      print *, "The program did not converged after ", n_cycles, " cycles"
    endif
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

   
