module fock_implementation

implicit none

private

public :: fock_matrix, density_matrix

contains

  subroutine fock_matrix(F, D, hcore, ao_integrals, n_ao)
    real(8), intent(in)    :: D(:,:), hcore(:,:), ao_integrals(:,:,:,:)
    integer, intent(in)    :: n_ao

    real(8), intent(inout) :: F(:,:)

    integer :: kappa, lambda

    F = hcore
    do lambda = 1, n_ao
      do kappa = 1, n_ao
        F(kappa,lambda) = &
        F(kappa,lambda) &
        + 2.D0 * sum(D*ao_integrals(kappa,lambda,:,:)) &
        - 1.D0 * sum(D*ao_integrals(kappa,:,:,lambda))
      end do
    end do
  end subroutine fock_matrix

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
   
  !  subroutine temp()

  !  end subroutine temp

end module