module fock_implementation

implicit none
! by Nino Keuzenkamp

private

public :: fock_matrix, density_matrix, convergence_check, nuclear_repulsion_energy

contains

  ! make Fock matrix
  subroutine fock_matrix(F, D, hcore, ao_integrals, n_ao)
    real(8), intent(in)    :: D(:,:), hcore(:,:), ao_integrals(:,:,:,:)
    integer, intent(in)    :: n_ao

    real(8), intent(out) :: F(:,:)

    integer :: kappa, lambda

    F = hcore
    do lambda = 1, n_ao
      do kappa = 1, n_ao
        F(kappa,lambda) = F(kappa,lambda)              &
        + 2.D0 * sum(D*ao_integrals(kappa,lambda,:,:)) &
        - 1.D0 * sum(D*ao_integrals(kappa,:,:,lambda))
      end do
    end do
  end subroutine fock_matrix


  ! create density matrix
  subroutine density_matrix(D, C, n_ao, n_occ)
    real(8), intent(in)  :: C(:,:)
    integer, intent(in)  :: n_occ, n_ao

    real(8), intent(out) :: D(:,:)
    integer :: kappa, lambda 

    do lambda = 1, n_ao
      do kappa = 1, n_ao
          D(kappa,lambda) = sum(C(kappa,1:n_occ)*C(lambda,1:n_occ))
      end do
    end do

  end subroutine density_matrix


  ! check convergence
  ! if the energy difference and density matrix difference is within the tolerance, return true
  logical function convergence_check(E_new, E_old, D_new, D_old, tolerance_E, tolerance_D) result(bool)
  real(8), intent(in) :: D_new(:,:), D_old(:,:)
  real(8), intent(in) :: E_new, E_old
  real(8), intent(in) :: tolerance_E, tolerance_D

  bool = .false.

  if ((abs(E_new - E_old) < tolerance_E) .and. (sqrt(sum((D_new - D_old)**2)) < tolerance_D)) then
    bool = .true.
  endif
  
  end function convergence_check


  ! calculates the nuclear repulsion energy
  real(8) function nuclear_repulsion_energy(molecule) result(energy)
  use molecular_structure

  type(molecular_structure_t), intent(in) :: molecule

  real(8) :: distance
  integer :: i, j

  energy = 0
  ! for each unique atom pair
  do i = 1, molecule%num_atoms - 1
    do j = i + 1, molecule%num_atoms

      ! calculate the distance
      distance = sum((molecule%coord(:, i) - molecule%coord(:, j))**2)
      distance = sqrt(distance)

      ! add the energy to the total
      energy = energy + (molecule%charge(i) * molecule%charge(j))/(distance)
    enddo
  enddo

  end function nuclear_repulsion_energy

end module