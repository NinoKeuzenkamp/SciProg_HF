program HartreeFock

  ! Demonstration program that can be used as a starting point
  ! Lucas Visscher, March 2022

  use molecular_structure
  use ao_basis
  use compute_integrals
  use diagonalization
  use fock_implementation
  use input_handeling

  implicit none

  ! Variable containing the molecular structure
  type(molecular_structure_t) :: molecule
  ! Variable containing the atomic orbital basis
  type(basis_set_info_t) :: ao_basis

  integer  :: n_ao, n_occ
  integer  :: kappa, lambda, i          ! loop integers
  real(8)  :: E_HF_old, E_HF_new, E_nuc

  real(8), allocatable :: hcore(:,:), V(:,:), T(:,:), S(:,:), ao_integrals(:,:,:,:)
  real(8), allocatable :: D_old(:,:), D_new(:,:), F(:,:)
  real(8), allocatable :: C(:,:), eps(:)

  ! for SCF loop
  integer :: n_cycles
  real(8) :: tolerance_E, tolerance_D
  logical :: converged


  ! get molecule, ao_basis, and n_occ
  call read_input(molecule, ao_basis, n_occ, n_cycles)
  n_ao = ao_basis%nao

  ! allocate all arrays
  allocate(S(n_ao,n_ao))
  allocate(T(n_ao,n_ao))
  allocate(V(n_ao,n_ao))
  allocate(F(n_ao,n_ao))
  allocate(hcore(n_ao,n_ao))
  allocate(C(n_ao,n_ao))
  allocate(eps(n_ao))
  allocate(D_new(n_ao,n_ao))
  allocate(D_old(n_ao,n_ao))
  allocate(ao_integrals(n_ao,n_ao,n_ao,n_ao))

  ! Compute overlap matrix, "S", kinetic matrix, "T", and potential matrix, "V"
  call compute_1e_integrals ("OVL",ao_basis,ao_basis,S)
  call compute_1e_integrals ("KIN",ao_basis,ao_basis,T)
  call compute_1e_integrals ("POT",ao_basis,ao_basis,V,molecule)

  ! Compute the core Hamiltonian matrix (the potential is positive, so we scale with -e = -1 to get to the potential energy matrix)
  hcore = T - V

  ! Diagonalize the fock matrix
  call solve_genev (hcore,S,C,eps)

  ! Form the density matrix
  call density_matrix(D_old, C, n_ao, n_occ)

  ! calculate energy of core-only system
  E_HF_old = sum((hcore + F) * D_old)
  
  ! Compute all 2-electron integrals
  call generate_2int (ao_basis,ao_integrals)


  ! set tolerances for convergence
  tolerance_E = 1.D-9
  tolerance_D = 1.D-3 
  ! SCF loop
  i = 1
  do while (i <= n_cycles)

    ! calculate fock matrix
    call fock_matrix(F, D_old, hcore, ao_integrals, n_ao)

    ! diagonalise fock matrix
    call solve_genev (F,S,C,eps)

    ! form the density matrix
    call density_matrix(D_new, C, n_ao, n_occ)

    ! calculate HF energy for current cycle
    E_HF_new = sum((hcore + F) * D_new)
    ! print *, "The Hartree-Fock energy:    ", E_HF_new + E_nuc, " of cycle: ", i       ! for debugging energy

    ! if converged, exit SCF loop
    converged = convergence_check(E_HF_new, E_HF_old, D_new, D_old, tolerance_E, tolerance_D)
    if (converged) exit

    ! re-initialise "old" variables before entering the next cycle
    E_HF_old = E_HF_new
    D_old    = D_new

    i = i + 1
  enddo ! end of SCF loop

  E_nuc = nuclear_repulsion_energy(molecule)
  E_HF_new = E_HF_new + E_nuc
  E_HF_old = E_HF_old + E_nuc

  if (converged) then
    print "(a, i3, a)", "The program has converged after ", i, " cycles."
    print "(a, f10.5, a)", "The energy has converged to ", E_HF_new, " Hartree."
    print "(a, f10.5, a)", "The nuclear repulsion energy was ", E_nuc, " Hartree."
  else
    print "(a, i3, a)", "The program did not converged after ", n_cycles, " cycles."
    print "(a, f10.5, a)", "The energy of the last cycle was ", E_HF_new, " Hartree."
    print "(a, f10.5, a)", "The energy of the second last cycle was ", E_HF_old, " Hartree."
  endif

end program HartreeFock




   
