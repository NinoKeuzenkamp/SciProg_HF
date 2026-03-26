program HartreeFock

  ! Demonstration program that can be used as a starting point
  ! Lucas Visscher, March 2022

  ! modified by Nino Keuzenkamp

  use molecular_structure
  use ao_basis
  use compute_integrals
  use diagonalization
  use fock_implementation
  use input_handeling
  use output
  use types
  use MP2

  implicit none

  type(molecular_structure_t) :: molecule
  type(basis_set_info_t)      :: ao_basis

  integer           :: n_ao, n_occ
  integer           :: kappa, lambda, i          ! loop integers
  type(energy_type) :: energy                    ! collection of energy variables (results)
  type(calculation_preset) :: preset             ! see derivedtypes.F90,  options/settings might be a better name for this

  real(8), allocatable :: hcore(:,:), V(:,:), T(:,:), S(:,:), ao_integrals(:,:,:,:)
  real(8), allocatable :: D_old(:,:), D_new(:,:), F(:,:)
  real(8), allocatable :: C(:,:), eps(:)

  ! for SCF loop
  integer :: n_cycles
  real(8) :: tolerance_E, tolerance_D
  logical :: converged        ! true if the SCF cycle converged

  character(32) :: outfile
  real(8)       :: time_start, time_end ! to keep track of the time elapsed

  call cpu_time(time_start)

  ! get molecule, ao_basis, n_occ, n_cycles, outfile
  call read_input(molecule, ao_basis, n_occ, n_cycles, preset)
  n_ao = ao_basis%nao
  call get_output_file(outfile)

  ! allocate all arrays
  allocate(energy%all_SCF(n_cycles))
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
  ! Form the density matrix
  ! calculate energy of core-only system
  call solve_genev (hcore,S,C,eps)
  call density_matrix(D_old, C, n_ao, n_occ)
  energy%HF_old = sum(hcore * D_old)
  
  ! Compute all 2-electron integrals
  call generate_2int (ao_basis,ao_integrals)


  ! set tolerances for convergence
  tolerance_E = 1.D-6     ! in Hartree
  tolerance_D = 1.D-3 
  ! SCF loop
  i = 1
  do while (i <= n_cycles)

    ! calculate fock matrix
    ! diagonalise fock matrix
    ! form the density matrix
    call fock_matrix(F, D_old, hcore, ao_integrals, n_ao)
    call solve_genev (F,S,C,eps)
    call density_matrix(D_new, C, n_ao, n_occ)

    ! calculate HF energy for current cycle
    energy%HF = sum((hcore + F) * D_new)
    energy%all_SCF(i) = energy%HF
   
    ! convergence check
    converged = convergence_check(energy%HF, energy%HF_old, D_new, D_old, tolerance_E, tolerance_D)
    if (converged) then
      print "(/, a, i4, a)", "Program converged on cycle ", i, "."
      exit
    endif

    ! re-initialise "old" variables before entering the next cycle
    energy%HF_old = energy%HF
    D_old         = D_new

    i = i + 1
  enddo ! end of SCF loop

  ! also add nuclear repulsion energy to total energy
  energy%nuc    = nuclear_repulsion_energy(molecule)
  energy%HF     = energy%HF
  energy%HF_old = energy%HF_old

  call cpu_time(time_end)
  print "(a, f10.5, a)", "The program took ", time_end-time_start, " seconds to finish the SCF loop."
  call cpu_time(time_start)

  ! if MP2 is enabled
  if (preset%MP2) then
    print "(/, a)", "Now starting MP2 correlation energy calculation."
    energy%MP2 = MP2_energy(ao_integrals, C, eps, n_occ, n_ao)
    call cpu_time(time_end)
    print "(a, f10.5, a)", "The program took ", time_end-time_start, " seconds to finish the MP2 correlation energy calculation."
  endif

  ! i is the last cycle of the SCF loop!
  call write_to_file(molecule, energy, outfile, converged, i, eps, n_occ, preset)
  print "(/, a)", "Result file has been made."

end program HartreeFock




   
