module MP2

implicit none

private

public :: MP2_energy

contains

    real(8) function MP2_energy(ao_integrals, C, eps, n_occ, n_ao) result(energy)
        real(8), intent(in) :: ao_integrals(:,:,:,:), C(:,:), eps(:)
        integer, intent(in) :: n_occ, n_ao

        integer :: i, a, j, b
        energy = 0

        ! for all occupied orbitals
        do i = 1, n_occ
            do j = 1, n_occ

                ! with each virtual orbital
                do a = n_occ + 1, n_ao 
                    do  b = n_occ + 1, n_ao 

                        ! (ia||jb) / e_a + e_b - e_i - e_j

                        ! using physicist notation
                        ! [ <ij|ab> * [2<ij|ab> - <ij|ba>] ] / [e_i + e_j - e_a - e_b]
                        energy = energy + &
                                (AO_to_MO(ao_integrals, c, i, j, a, b) * &
                                (2 * AO_to_MO(ao_integrals, c, i, j, a, b) - AO_to_MO(ao_integrals, c, i, j, b, a))) / &
                                (eps(i) + eps(j) - eps(a) - eps(b))

                        ! using chemist notation
                        ! ! [ (ia|jb) * [2(ia|jb) - (ib|ja)] ] / [e_i + e_j - e_a - e_b]
                        ! energy = energy + &
                        !         (AO_to_MO(ao_integrals, c, i, a, j, b) * &
                        !         (2 * AO_to_MO(ao_integrals, c, i, a, j, b) - AO_to_MO(ao_integrals, c, j, i, b, a))) / &
                        !         (eps(i) + eps(j) - eps(a) - eps(b))
                                

                        ! print *, energy
                    enddo
                enddo

            enddo
        enddo


    end function MP2_energy

    real(8) function AO_to_MO(ao_integrals, C, p, q, r, s) result(integral)
        real(8), intent(in) :: ao_integrals(:,:,:,:), C(:,:)
        integer, intent(in) :: p, q, r, s

        integer :: kappa, lambda, mu, nu, n_ao

        n_ao = size(C, 1)
        integral = 0

        ! kappa: MO -> p
        do kappa = 1, n_ao
            ! lambda: MO -> q
            do lambda = 1, n_ao
                ! mu: MO -> r
                do mu = 1, n_ao
                    ! nu: MO -> s
                    do nu = 1, n_ao
                        ! ao_integrals uses chemist notation

                        ! ! keep chemist notation
                        ! integral = integral + ao_integrals(kappa, lambda, mu, nu) * &

                        ! indices swapped, chemist notation -> physicist notation
                        integral = integral + ao_integrals(kappa, mu, lambda, nu) * &

                        C(kappa, p) * C(lambda, q) * C(mu, r) * C(nu, s)
                    enddo
                enddo
            enddo
        enddo

    end function AO_to_MO

end module