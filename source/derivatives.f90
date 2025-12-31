module cea_derivatives
    use cea_param, only: dp, empty_dp, R=>gas_constant, &
                         snl=>species_name_len, &
                         Avgdr=>avogadro, &
                         Boltz=>boltzmann, &
                         pi
    use cea_equilibrium, only: EqSolver, EqSolution, gauss
    use cea_mixture, only: Mixture, MixtureThermo
    use fb_utils

    type :: EqTotals
        ! Inputs
        type(EqSolver) :: eq_solver
        type(EqSolution) :: eq_solution

        ! Final values
        real(dp), allocatable :: dnj_dstate1(:)      ! Total derivative of nj wrt state1 (h, s, t, or u)
        real(dp)              :: dn_dstate1          ! Total derivative of n wrt state1 (h, s, t, or u)
        real(dp)              :: dT_dstate1          ! Total derivative of T wrt state1 (h, s, t, or u)

        real(dp), allocatable :: dnj_dstate2(:)      ! Total derivative of nj wrt state2 (p or v)
        real(dp)              :: dn_dstate2          ! Total derivative of n wrt state2 (p or v)
        real(dp)              :: dT_dstate2          ! Total derivative of T wrt state2 (p or v)

        real(dp), allocatable :: dnj_dweights(:, :)  ! Total derivative of nj wrt reactant weights
        real(dp), allocatable :: dn_dweights(:)      ! Total derivative of n wrt reactant weights
        real(dp), allocatable :: dT_dweights(:)      ! Total derivative of T wrt reactant weights

        ! Solution workspace
        integer :: nx, nu   ! Sizes of x and u
        integer :: num_eqn  ! Number of equations in the implicit system (size of u_hat)
        real(dp), allocatable :: u_hat(:)            ! Solution of the converged implicit system
        real(dp), allocatable :: dR_du(:, :)         ! Partial derivative matrix of R wrt u
        real(dp), allocatable :: dR_dx(:, :)         ! Partial derivative matrix of R wrt x
        real(dp), allocatable :: du_dx(:, :)         ! Total derivative of u wrt x

        ! Slice indices
        integer :: m_eq(2)
        integer :: n_H0
        integer :: n_b0(2)
        integer :: n_P

        integer :: m_pi(2)
        integer :: m_dlnn
        integer :: m_dlnT
        integer :: m_lnnj(2)
        integer :: m_nj(2)
        integer :: m_n
        integer :: m_T
        integer :: m_b(2)
        integer :: m_Hj(2)
        integer :: m_Sj(2)
        integer :: m_Uj(2)
        integer :: m_Cpj(2)

        integer :: n_pi(2)
        integer :: n_dlnn
        integer :: n_dlnT
        integer :: n_lnnj(2)
        integer :: n_nj(2)
        integer :: n_n
        integer :: n_T
        integer :: n_b(2)
        integer :: n_Hj(2)
        integer :: n_Sj(2)
        integer :: n_Uj(2)
        integer :: n_Cpj(2)

        ! Sub-Jacobians of dR/dx
        real(dp), pointer :: dR_eq_dH0(:)
        real(dp), pointer :: dR_eq_db0(:,:)
        real(dp), pointer :: dR_eq_dP(:)
        real(dp), pointer :: dR_lnnj_dP(:)
        real(dp), pointer :: dR_n_dP
        real(dp), pointer :: dR_T_dP

        ! Sub-Jacobians of dR/du
        real(dp), pointer :: dR_eq_duhat(:,:)
        real(dp), pointer :: dR_eq_dnj(:,:)
        real(dp), pointer :: dR_eq_dn(:)
        real(dp), pointer :: dR_eq_db(:,:)
        real(dp), pointer :: dR_eq_dlnnj(:,:)
        real(dp), pointer :: dR_eq_dHj(:,:)
        real(dp), pointer :: dR_eq_dSj(:,:)
        real(dp), pointer :: dR_eq_dUj(:,:)
        real(dp), pointer :: dR_eq_dCpj(:,:)
        real(dp), pointer :: dR_eq_dCvj(:,:)
        real(dp), pointer :: dR_eq_dT(:)
        real(dp), pointer :: dR_nj_dlnnj(:,:)
        real(dp), pointer :: dR_nj_dnj(:,:)
        real(dp), pointer :: dR_lnnj_dpi(:,:)
        real(dp), pointer :: dR_lnnj_ddlnn(:)
        real(dp), pointer :: dR_lnnj_ddlnT(:)
        real(dp), pointer :: dR_lnnj_dn(:)
        real(dp), pointer :: dR_lnnj_dHj(:,:)
        real(dp), pointer :: dR_lnnj_dSj(:,:)
        real(dp), pointer :: dR_lnnj_dlnnj(:,:)
        real(dp), pointer :: dR_n_dpi(:)
        real(dp), pointer :: dR_n_ddlnn
        real(dp), pointer :: dR_n_ddlnT
        real(dp), pointer :: dR_n_dlnnj(:)
        real(dp), pointer :: dR_n_dn
        real(dp), pointer :: dR_n_dHj(:)
        real(dp), pointer :: dR_n_dSj(:)
        real(dp), pointer :: dR_T_dpi(:)
        real(dp), pointer :: dR_T_ddlnn
        real(dp), pointer :: dR_T_ddlnT
        real(dp), pointer :: dR_T_dlnnj(:)
        real(dp), pointer :: dR_T_dn
        real(dp), pointer :: dR_T_dT
        real(dp), pointer :: dR_T_dHj(:)
        real(dp), pointer :: dR_T_dSj(:)
        real(dp), pointer :: dR_b_dnj(:,:)
        real(dp), pointer :: dR_b_db(:,:)
        real(dp), pointer :: dR_Hj_dT(:)
        real(dp), pointer :: dR_Hj_dHj(:,:)
        real(dp), pointer :: dR_Sj_dT(:)
        real(dp), pointer :: dR_Sj_dSj(:,:)
        real(dp), pointer :: dR_Uj_dT(:)
        real(dp), pointer :: dR_Uj_dUj(:,:)
        real(dp), pointer :: dR_Cp_dT(:)
        real(dp), pointer :: dR_Cp_dCp(:,:)

    contains
        procedure :: compute_totals => EqTotals_compute_totals

        ! dR/dx functions
        procedure :: compute_dR_dx => EqTotals_compute_dR_dx
        procedure :: compute_dR_eq_dH0 => EqTotals_compute_dR_eq_dH0
        procedure :: compute_dR_eq_db0 => EqTotals_compute_dR_eq_db0
        procedure :: compute_dR_eq_dP => EqTotals_compute_dR_eq_dP
        procedure :: compute_dR_lnnj_dP => EqTotals_compute_dR_lnnj_dP
        procedure :: compute_dR_n_dP => EqTotals_compute_dR_n_dP
        procedure :: compute_dR_T_dP => EqTotals_compute_dR_T_dP

        ! dR/du functions
        procedure :: compute_dR_du => EqTotals_compute_dR_du
        procedure :: compute_dR_eq_dlnnj => EqTotals_compute_dR_eq_dlnnj
        procedure :: compute_dR_eq_dnj => EqTotals_compute_dR_eq_dnj
        procedure :: compute_dR_eq_dn => EqTotals_compute_dR_eq_dn
        procedure :: compute_dR_eq_dT => EqTotals_compute_dR_eq_dT
        procedure :: compute_dR_eq_db => EqTotals_compute_dR_eq_db
        procedure :: compute_dR_eq_dHj => EqTotals_compute_dR_eq_dHj
        procedure :: compute_dR_eq_dSj => EqTotals_compute_dR_eq_dSj
        procedure :: compute_dR_eq_dUj => EqTotals_compute_dR_eq_dUj
        procedure :: compute_dR_eq_dCpj => EqTotals_compute_dR_eq_dCpj
        procedure :: compute_dR_eq_dCvj => EqTotals_compute_dR_eq_dCvj

        procedure :: compute_dR_lnnj_dpi => EqTotals_compute_dR_lnnj_dpi
        procedure :: compute_dR_lnnj_ddlnn => EqTotals_compute_dR_lnnj_ddlnn
        procedure :: compute_dR_lnnj_ddlnT => EqTotals_compute_dR_lnnj_ddlnT
        procedure :: compute_dR_lnnj_dlnnj => EqTotals_compute_dR_lnnj_dlnnj
        procedure :: compute_dR_lnnj_dn => EqTotals_compute_dR_lnnj_dn
        procedure :: compute_dR_lnnj_dHj => EqTotals_compute_dR_lnnj_dHj
        procedure :: compute_dR_lnnj_dSj => EqTotals_compute_dR_lnnj_dSj

        procedure :: compute_dR_nj_dlnnj => EqTotals_compute_dR_nj_dlnnj

        procedure :: compute_dR_n_dpi => EqTotals_compute_dR_n_dpi
        procedure :: compute_dR_n_ddlnn => EqTotals_compute_dR_n_ddlnn
        procedure :: compute_dR_n_ddlnT => EqTotals_compute_dR_n_ddlnT
        procedure :: compute_dR_n_dlnnj => EqTotals_compute_dR_n_dlnnj
        procedure :: compute_dR_n_dHj => EqTotals_compute_dR_n_dHj
        procedure :: compute_dR_n_dSj => EqTotals_compute_dR_n_dSj

        procedure :: compute_dR_T_dpi => EqTotals_compute_dR_T_dpi
        procedure :: compute_dR_T_ddlnn => EqTotals_compute_dR_T_ddlnn
        procedure :: compute_dR_T_ddlnT => EqTotals_compute_dR_T_ddlnT
        procedure :: compute_dR_T_dlnnj => EqTotals_compute_dR_T_dlnnj
        procedure :: compute_dR_T_dn => EqTotals_compute_dR_T_dn
        procedure :: compute_dR_T_dHj => EqTotals_compute_dR_T_dHj
        procedure :: compute_dR_T_dSj => EqTotals_compute_dR_T_dSj

        procedure :: compute_dR_b_dnj => EqTotals_compute_dR_b_dnj

        procedure :: compute_dR_Hj_dT => EqTotals_compute_dR_Hj_dT

        procedure :: compute_dR_Sj_dT => EqTotals_compute_dR_Sj_dT

        procedure :: compute_dR_Uj_dT => EqTotals_compute_dR_Uj_dT

        procedure :: compute_dR_Cp_dT => EqTotals_compute_dR_Cp_dT

        ! Helper functions to compute smaller partials
        procedure :: dG_dP => EqTotals_dG_dP
        procedure :: dG_dnj => EqTotals_dG_dnj
        procedure :: dG_dn => EqTotals_dG_dn
        procedure :: dG_dlnnj => EqTotals_dG_dlnnj
        procedure :: dG_dHj => EqTotals_dG_dHj
        procedure :: dG_dSj => EqTotals_dG_dSj
        procedure :: dG_dUj => EqTotals_dG_dUj

        ! procedure :: dlambda_dP => EqTotals_dlambda_dP
        ! procedure :: dlambda_dn => EqTotals_dlambda_dn
        ! procedure :: dlambda_dlnnj => EqTotals_dlambda_dlnnj
        ! procedure :: dlambda_ddlnnj => EqTotals_dlambda_ddlnnj
        ! procedure :: dlambda_ddlnn => EqTotals_dlambda_ddlnn
        ! procedure :: dlambda_ddlnT => EqTotals_dlambda_ddlnT

        ! procedure :: ddlnnj_dpi => EqTotals_ddlnnj_dpi
        ! procedure :: ddlnnj_ddlnn => EqTotals_ddlnnj_ddlnn
        ! procedure :: ddlnnj_ddlnT => EqTotals_ddlnnj_ddlnT
        ! procedure :: ddlnnj_dlnnj => EqTotals_ddlnnj_dlnnj
        ! procedure :: ddlnnj_dn => EqTotals_ddlnnj_dn
        ! procedure :: ddlnnj_dP => EqTotals_ddlnnj_dP
        ! procedure :: ddlnnj_dHj => EqTotals_ddlnnj_dHj
        ! procedure :: ddlnnj_dSj => EqTotals_ddlnnj_dSj

        ! Finite difference checks
        procedure :: finite_difference_totals => EqTotals_finite_difference_totals
        procedure :: fd_dR_eq_dH0 => EqTotals_fd_dR_eq_dH0
        procedure :: fd_dR_eq_db0 => EqTotals_fd_dR_eq_db0
        procedure :: fd_dR_eq_dP => EqTotals_fd_dR_eq_dP
        procedure :: fd_dR_eq_dnj => EqTotals_fd_dR_eq_dnj
        procedure :: fd_dR_eq_dn => EqTotals_fd_dR_eq_dn
        procedure :: fd_dR_eq_dT => EqTotals_fd_dR_eq_dT
        procedure :: fd_dR_eq_db => EqTotals_fd_dR_eq_db
        procedure :: fd_dR_eq_dlnnj => EqTotals_fd_dR_eq_dlnnj
        procedure :: fd_dR_eq_dHj => EqTotals_fd_dR_eq_dHj
        procedure :: fd_dR_eq_dSj => EqTotals_fd_dR_eq_dSj
        procedure :: fd_dR_eq_dUj => EqTotals_fd_dR_eq_dUj
        procedure :: fd_dR_eq_dCpj => EqTotals_fd_dR_eq_dCpj

        procedure :: fd_dR_lnnj_dP => EqTotals_fd_dR_lnnj_dP
        procedure :: fd_dR_lnnj_dpi => EqTotals_fd_dR_lnnj_dpi
        procedure :: fd_dR_lnnj_ddlnn => EqTotals_fd_dR_lnnj_ddlnn
        procedure :: fd_dR_lnnj_ddlnT => EqTotals_fd_dR_lnnj_ddlnT
        procedure :: fd_dR_lnnj_dlnnj => EqTotals_fd_dR_lnnj_dlnnj
        procedure :: fd_dR_lnnj_dn => EqTotals_fd_dR_lnnj_dn
        procedure :: fd_dR_lnnj_dHj => EqTotals_fd_dR_lnnj_dHj
        procedure :: fd_dR_lnnj_dSj => EqTotals_fd_dR_lnnj_dSj

        procedure :: fd_dR_nj_dlnnj => EqTotals_fd_dR_nj_dlnnj

        procedure :: fd_dR_n_dP => EqTotals_fd_dR_n_dP
        procedure :: fd_dR_n_dpi => EqTotals_fd_dR_n_dpi
        procedure :: fd_dR_n_ddlnn => EqTotals_fd_dR_n_ddlnn
        procedure :: fd_dR_n_ddlnT => EqTotals_fd_dR_n_ddlnT
        procedure :: fd_dR_n_dlnnj => EqTotals_fd_dR_n_dlnnj
        procedure :: fd_dR_n_dn => EqTotals_fd_dR_n_dn
        procedure :: fd_dR_n_dHj => EqTotals_fd_dR_n_dHj
        procedure :: fd_dR_n_dSj => EqTotals_fd_dR_n_dSj

        procedure :: fd_dR_T_dP => EqTotals_fd_dR_T_dP
        procedure :: fd_dR_T_dpi => EqTotals_fd_dR_T_dpi
        procedure :: fd_dR_T_ddlnn => EqTotals_fd_dR_T_ddlnn
        procedure :: fd_dR_T_ddlnT => EqTotals_fd_dR_T_ddlnT
        procedure :: fd_dR_T_dlnnj => EqTotals_fd_dR_T_dlnnj
        procedure :: fd_dR_T_dn => EqTotals_fd_dR_T_dn
        procedure :: fd_dR_T_dHj => EqTotals_fd_dR_T_dHj
        procedure :: fd_dR_T_dSj => EqTotals_fd_dR_T_dSj

        procedure :: fd_dR_b_dnj => EqTotals_fd_dR_b_dnj

        procedure :: fd_dlambda_dP => EqTotals_fd_dlambda_dP
        procedure :: fd_dlambda_ddlnn => EqTotals_fd_dlambda_ddlnn
        procedure :: fd_dlambda_ddlnT => EqTotals_fd_dlambda_ddlnT
        procedure :: fd_dlambda_dpi => EqTotals_fd_dlambda_dpi
        procedure :: fd_dlambda_dn => EqTotals_fd_dlambda_dn
        procedure :: fd_dlambda_dlnnj => EqTotals_fd_dlambda_dlnnj
        procedure :: fd_dlambda_ddlnnj => EqTotals_fd_dlambda_ddlnnj
        procedure :: fd_dlambda_dHj => EqTotals_fd_dlambda_dHj
        procedure :: fd_dlambda_dSj => EqTotals_fd_dlambda_dSj

        procedure :: check_dR_dx => EqTotals_check_dR_dx
        procedure :: check_dR_du => EqTotals_check_dR_du

    end type
    interface EqTotals
        module procedure :: EqTotals_init
    end interface

contains

    ! ----------------------------------------------------------------------
    ! EqTotals
    ! ----------------------------------------------------------------------
    function EqTotals_init(eq_solver, eq_solution) result(eq_totals)
        type(EqTotals), target :: eq_totals
        type(EqSolver), intent(in) :: eq_solver
        type(EqSolution), intent(in) :: eq_solution

        ! Locals
        integer :: nx, nu   ! Size of x and u
        integer :: num_eqn  ! Number of equations in the implicit system (size of u_hat)
        integer :: ne       ! Number of elements
        integer :: ng       ! Number of gas species
        integer :: np       ! Number of products
        integer :: nr       ! Number of reactants

        ! Define shorthand
        nx = 2 + eq_solver%num_elements
        nu = 2*eq_solver%num_elements + 6*eq_solver%num_gas + 4
        num_eqn = eq_solution%num_equations(eq_solver)
        ne = eq_solver%num_elements
        ng = eq_solver%num_gas
        np = eq_solver%num_products
        nr = eq_solver%num_reactants

        eq_totals%eq_solver = eq_solver
        eq_totals%eq_solution = eq_solution
        eq_totals%num_eqn = num_eqn
        eq_totals%nx = nx
        eq_totals%nu = nu

        allocate(eq_totals%dnj_dstate1(np), &
                 eq_totals%dnj_dstate2(np), &
                 eq_totals%dnj_dweights(np, nr), &
                 eq_totals%dn_dweights(nr), &
                 eq_totals%dT_dweights(nr), &
                 eq_totals%du_dx(nu, nx), &
                 eq_totals%dR_du(nu, nu), &
                 eq_totals%dR_dx(nu, nx), &
                 eq_totals%u_hat(num_eqn))

        ! Set the slice indices
        eq_totals%n_H0 = 1
        eq_totals%n_b0 = [2, ne+1]
        eq_totals%n_P  = ne + 2

        eq_totals%m_eq   = [1, num_eqn]
        eq_totals%m_lnnj = [num_eqn+1, num_eqn+ng]
        eq_totals%m_nj   = [num_eqn+ng+1, num_eqn+2*ng]
        eq_totals%m_n    = num_eqn+2*ng+1
        eq_totals%m_T    = num_eqn+2*ng+2
        eq_totals%m_b    = [num_eqn+2*ng+3,num_eqn+2*ng+ne+2]
        eq_totals%m_Hj   = [num_eqn+2*ng+ne+3, num_eqn+3*ng+ne+2]
        eq_totals%m_Sj   = [num_eqn+3*ng+ne+3, num_eqn+4*ng+ne+2]
        eq_totals%m_Uj   = [num_eqn+4*ng+ne+3, num_eqn+5*ng+ne+2]
        eq_totals%m_Cpj  = [num_eqn+5*ng+ne+3, num_eqn+6*ng+ne+2]

        eq_totals%n_pi   = [1, ne]
        eq_totals%n_dlnn = ne + 1
        eq_totals%n_dlnT = ne + 2
        eq_totals%n_lnnj = [ne+3, ne+ng+2]
        eq_totals%n_nj   = [ne+ng+3, ne+2*ng+2]
        eq_totals%n_n    = ne+2*ng+3
        eq_totals%n_T    = ne+2*ng+4
        eq_totals%n_b    = [ne+2*ng+5,2*ne+2*ng+4]
        eq_totals%n_Hj   = [2*ne+2*ng+5, 2*ne+3*ng+4]
        eq_totals%n_Sj   = [2*ne+3*ng+5, 2*ne+4*ng+4]
        eq_totals%n_Uj   = [2*ne+4*ng+5, 2*ne+5*ng+4]
        eq_totals%n_Cpj  = [2*ne+5*ng+5, 2*ne+6*ng+4]

        ! Associate subarray pointers for dR/dx
        eq_totals%dR_eq_dH0 => eq_totals%dR_dx(eq_totals%m_eq(1):eq_totals%m_eq(2), eq_totals%n_H0)
        eq_totals%dR_eq_db0 => eq_totals%dR_dx(eq_totals%m_eq(1):eq_totals%m_eq(2), &
                                               eq_totals%n_b0(1):eq_totals%n_b0(2))
        eq_totals%dR_eq_dP  => eq_totals%dR_dx(eq_totals%m_eq(1):eq_totals%m_eq(2), eq_totals%n_P)

        eq_totals%dR_lnnj_dP => eq_totals%dR_dx(eq_totals%m_lnnj(1):eq_totals%m_lnnj(2), eq_totals%n_P)

        eq_totals%dR_n_dP => eq_totals%dR_dx(eq_totals%m_n, eq_totals%n_P)

        eq_totals%dR_T_dP => eq_totals%dR_dx(eq_totals%m_T, eq_totals%n_P)

        ! Associate subarray pointers for dR/du
        eq_totals%dR_eq_duhat => eq_totals%dR_du(eq_totals%m_eq(1):eq_totals%m_eq(2), &
                                                 eq_totals%m_eq(1):eq_totals%m_eq(2))
        eq_totals%dR_eq_dnj   => eq_totals%dR_du(eq_totals%m_eq(1):eq_totals%m_eq(2), &
                                                 eq_totals%m_nj(1):eq_totals%m_nj(2))
        eq_totals%dR_eq_dn    => eq_totals%dR_du(eq_totals%m_eq(1):eq_totals%m_eq(2), &
                                                 eq_totals%m_n)
        eq_totals%dR_eq_dT    => eq_totals%dR_du(eq_totals%m_eq(1):eq_totals%m_eq(2), &
                                                 eq_totals%m_T)
        eq_totals%dR_eq_db    => eq_totals%dR_du(eq_totals%m_eq(1):eq_totals%m_eq(2), &
                                                 eq_totals%m_b(1):eq_totals%m_b(2))
        eq_totals%dR_eq_dlnnj => eq_totals%dR_du(eq_totals%m_eq(1):eq_totals%m_eq(2), &
                                                 eq_totals%m_lnnj(1):eq_totals%m_lnnj(2))
        eq_totals%dR_eq_dHj   => eq_totals%dR_du(eq_totals%m_eq(1):eq_totals%m_eq(2), &
                                                 eq_totals%m_Hj(1):eq_totals%m_Hj(2))
        eq_totals%dR_eq_dSj   => eq_totals%dR_du(eq_totals%m_eq(1):eq_totals%m_eq(2), &
                                                 eq_totals%m_Sj(1):eq_totals%m_Sj(2))
        eq_totals%dR_eq_dUj   => eq_totals%dR_du(eq_totals%m_eq(1):eq_totals%m_eq(2), &
                                                 eq_totals%m_Uj(1):eq_totals%m_Uj(2))
        eq_totals%dR_eq_dCpj  => eq_totals%dR_du(eq_totals%m_eq(1):eq_totals%m_eq(2), &
                                                 eq_totals%m_Cpj(1):eq_totals%m_Cpj(2))

        eq_totals%dR_nj_dlnnj => eq_totals%dR_du(eq_totals%m_nj(1):eq_totals%m_nj(2), &
                                                 eq_totals%n_lnnj(1):eq_totals%n_lnnj(2))
        eq_totals%dR_nj_dnj   => eq_totals%dR_du(eq_totals%m_nj(1):eq_totals%m_nj(2), &
                                                 eq_totals%n_nj(1):eq_totals%n_nj(2))

        eq_totals%dR_lnnj_dpi   => eq_totals%dR_du(eq_totals%m_lnnj(1):eq_totals%m_lnnj(2), &
                                                   eq_totals%n_pi(1):eq_totals%n_pi(2))
        eq_totals%dR_lnnj_ddlnn => eq_totals%dR_du(eq_totals%m_lnnj(1):eq_totals%m_lnnj(2), &
                                                   eq_totals%n_dlnn)
        eq_totals%dR_lnnj_ddlnT => eq_totals%dR_du(eq_totals%m_lnnj(1):eq_totals%m_lnnj(2), &
                                                   eq_totals%n_dlnT)
        eq_totals%dR_lnnj_dn    => eq_totals%dR_du(eq_totals%m_lnnj(1):eq_totals%m_lnnj(2), &
                                                   eq_totals%n_n)
        eq_totals%dR_lnnj_dHj   => eq_totals%dR_du(eq_totals%m_lnnj(1):eq_totals%m_lnnj(2), &
                                                   eq_totals%n_Hj(1):eq_totals%n_Hj(2))
        eq_totals%dR_lnnj_dSj   => eq_totals%dR_du(eq_totals%m_lnnj(1):eq_totals%m_lnnj(2), &
                                                   eq_totals%n_Sj(1):eq_totals%n_Sj(2))
        eq_totals%dR_lnnj_dlnnj => eq_totals%dR_du(eq_totals%m_lnnj(1):eq_totals%m_lnnj(2), &
                                                   eq_totals%n_lnnj(1):eq_totals%n_lnnj(2))

        eq_totals%dR_n_dpi   => eq_totals%dR_du(eq_totals%m_n, eq_totals%n_pi(1):eq_totals%n_pi(2))
        eq_totals%dR_n_ddlnn => eq_totals%dR_du(eq_totals%m_n, eq_totals%n_dlnn)
        eq_totals%dR_n_ddlnT => eq_totals%dR_du(eq_totals%m_n, eq_totals%n_dlnT)
        eq_totals%dR_n_dlnnj => eq_totals%dR_du(eq_totals%m_n, eq_totals%n_lnnj(1):eq_totals%n_lnnj(2))
        eq_totals%dR_n_dn    => eq_totals%dR_du(eq_totals%m_n, eq_totals%n_n)
        eq_totals%dR_n_dHj   => eq_totals%dR_du(eq_totals%m_n, eq_totals%n_Hj(1):eq_totals%n_Hj(2))
        eq_totals%dR_n_dSj   => eq_totals%dR_du(eq_totals%m_n, eq_totals%n_Sj(1):eq_totals%n_Sj(2))

        eq_totals%dR_T_dpi   => eq_totals%dR_du(eq_totals%m_T, eq_totals%n_pi(1):eq_totals%n_pi(2))
        eq_totals%dR_T_ddlnn => eq_totals%dR_du(eq_totals%m_T, eq_totals%n_dlnn)
        eq_totals%dR_T_ddlnT => eq_totals%dR_du(eq_totals%m_T, eq_totals%n_dlnT)
        eq_totals%dR_T_dlnnj => eq_totals%dR_du(eq_totals%m_T, eq_totals%n_lnnj(1):eq_totals%n_lnnj(2))
        eq_totals%dR_T_dn    => eq_totals%dR_du(eq_totals%m_T, eq_totals%n_n)
        eq_totals%dR_T_dT    => eq_totals%dR_du(eq_totals%m_T, eq_totals%n_T)
        eq_totals%dR_T_dHj   => eq_totals%dR_du(eq_totals%m_T, eq_totals%n_Hj(1):eq_totals%n_Hj(2))
        eq_totals%dR_T_dSj   => eq_totals%dR_du(eq_totals%m_T, eq_totals%n_Sj(1):eq_totals%n_Sj(2))

        eq_totals%dR_b_dnj => eq_totals%dR_du(eq_totals%m_b(1):eq_totals%m_b(2), &
                                              eq_totals%n_nj(1):eq_totals%n_nj(2))
        eq_totals%dR_b_db  => eq_totals%dR_du(eq_totals%m_b(1):eq_totals%m_b(2), &
                                              eq_totals%n_b(1):eq_totals%n_b(2))

        eq_totals%dR_Hj_dT  => eq_totals%dR_du(eq_totals%m_Hj(1):eq_totals%m_Hj(2), &
                                               eq_totals%n_T)
        eq_totals%dR_Hj_dHj => eq_totals%dR_du(eq_totals%m_Hj(1):eq_totals%m_Hj(2), &
                                               eq_totals%n_Hj(1):eq_totals%n_Hj(2))

        eq_totals%dR_Sj_dT  => eq_totals%dR_du(eq_totals%m_Sj(1):eq_totals%m_Sj(2), &
                                               eq_totals%n_T)
        eq_totals%dR_Sj_dSj => eq_totals%dR_du(eq_totals%m_Sj(1):eq_totals%m_Sj(2), &
                                               eq_totals%n_Sj(1):eq_totals%n_Sj(2))

        eq_totals%dR_Uj_dT  => eq_totals%dR_du(eq_totals%m_Uj(1):eq_totals%m_Uj(2), &
                                               eq_totals%n_T)
        eq_totals%dR_Uj_dUj => eq_totals%dR_du(eq_totals%m_Uj(1):eq_totals%m_Uj(2), &
                                               eq_totals%n_Uj(1):eq_totals%n_Uj(2))

        eq_totals%dR_Cp_dT  => eq_totals%dR_du(eq_totals%m_Cpj(1):eq_totals%m_Cpj(2), &
                                               eq_totals%n_T)
        eq_totals%dR_Cp_dCp => eq_totals%dR_du(eq_totals%m_Cpj(1):eq_totals%m_Cpj(2), &
                                               eq_totals%n_Cpj(1):eq_totals%n_Cpj(2))

    end function EqTotals_init

    subroutine EqTotals_compute_totals(self, reactant_weights)
        ! Compute the total derivatives of the equilibrium solution
        ! Computes derivatives of: nj, n, and T with respect to state1, state2, and reactant weights
        !
        ! We create a full residual vector R, including both the implicit and explicit equations of
        ! the equilibrium solve. Importantly, we add *many* variables to the state vector u to compute
        ! the total derivatives. R_eq is the implicit system of equations that we use to solve equilibrium.
        ! The remaining R terms are explicit equations for quantities that are computed after the R_eq system
        ! solved, but that R_eq depend on as well. We have:
        ! R = [R_eq, R_nj, R_lnnj, R_n, R_T, R_b, R_Hj, R_Sj, R_Uj, R_Cpj]
        ! x = [H0/S0/U0, b0, P]
        ! u = [pi, ln(n), ln(T), ln(nj), nj, n, T, b, Hj, Sj, Uj, Cp_j]
        !
        ! R_eq = G*u_hat - f = 0; where u_hat = [pi, ln(n), ln(T)]
        ! R_nj = nj - exp(ln(nj) + lambda*dln(nj)) = 0
        ! R_lnnj = dln(nj) + mu_j - dln(n) - Aij*pi - dln(T)*H_j = 0
        ! R_n = n - exp(ln(n) + lambda*dln(n)) = 0
        ! R_T = T - exp(ln(T) + lambda*dln(T)) = 0
        ! R_b = b - Aij*nj = 0
        ! R_Hj = Hj - f(T) = 0
        ! R_Sj = Sj - f(T) = 0
        ! R_Uj = Uj - f(T) = 0
        ! R_Cpj = Cp_j - f(T) = 0
        !
        ! The required total derivatives found by:
        ! du/dx = -dR/du^-1 * dR/dx

        ! Arguments
        class(EqTotals) :: self
        real(dp), intent(in) :: reactant_weights(:)

        ! Locals
        logical :: const_t, const_p ! Flag that is true if problem is constant temperature
        real(dp) :: tmp(self%nu, self%nu+1)  ! Temporary matrix for solving the total derivative matrix equation
        real(dp), allocatable :: dnj_db0(:,:), dn_db0(:), dT_db0(:)  ! Total derivatives of nj, n, and T wrt reactant weights
        integer :: ne, ng  ! Number of elements and gas species
        integer :: i, j, ierr

        ! Define shorthand
        ne = self%eq_solver%num_elements
        ng = self%eq_solver%num_gas
        const_t = self%eq_solution%constraints%is_constant_temperature()
        const_p = self%eq_solution%constraints%is_constant_pressure()

        ! Save the solution of the converged implicit system
        self%u_hat(:) = self%eq_solution%G(:self%num_eqn, self%num_eqn+1)

        ! Compute the partial derivative matrix of R wrt x
        call self%compute_dR_dx()

        ! Compute the partial derivative matrix of R wrt u
        call self%compute_dR_du()

        ! Solve the total derivative matrix equation du/dx = -dR/du^-1 * dR/dx
        do i = 1, self%nx
            tmp(1:self%nu, 1:self%nu) = -self%dR_du
            tmp(1:self%nu, self%nu+1) = self%dR_dx(:, i)
            call gauss(tmp, ierr)
            self%du_dx(:, i) = tmp(1:self%nu, self%nu+1)
        end do

        ! Assign the total derivatives
        self%dnj_dstate1 = self%du_dx(ne+2:ne+1+2*ng, 1)
        self%dn_dstate1 = self%du_dx(ne+3+2*ng, 1)
        self%dT_dstate1 = self%du_dx(ne+4+2*ng, 1)

        self%dnj_dstate2 = self%du_dx(ne+2:ne+1+2*ng, self%nx)
        self%dn_dstate2 = self%du_dx(ne+3+2*ng, self%nx)
        self%dT_dstate2 = self%du_dx(ne+4+2*ng, self%nx)

        dnj_db0 = self%du_dx(ne+2:ne+1+2*ng, 1:self%nx-1)
        dn_db0 = self%du_dx(ne+3+2*ng, 1:self%nx-1)
        dT_db0 = self%du_dx(ne+4+2*ng, 1:self%nx-1)

        ! write(*,*) "du_dx = "
        ! do i = 1, self%nu
        !     write(*,*) (self%du_dx(i, j), j = 1,self%nx)
        ! end do

        write(*,*) "dnj_state1 = ", self%dnj_dstate1
        write(*,*) "dn_state1 = ", self%dn_dstate1
        write(*,*) "dT_state1 = ", self%dT_dstate1

        write(*,*) "dnj_state2 = ", self%dnj_dstate2
        write(*,*) "dn_state2 = ", self%dn_dstate2
        write(*,*) "dT_state2 = ", self%dT_dstate2

        ! write(*,*) "dnj_db0 = "
        ! do i = 1, ne
        !     write(*,*) dnj_db0(:, i)
        ! end do
        ! write(*,*) "dn_db0 = ", dn_db0
        ! write(*,*) "dT_db0 = ", dT_db0

        ! TODO: Get the desired total derivatives

        ! self%dnj_dstate1 = 0.0d0
        ! if (const_p) self%dn_dstate1 = 0.0d0
        ! if (.not. const_t) then
        !     self%dT_dstate1 = 0.0d0
        ! else
        !     self%dT_dstate1 = 1.0d0   ! Derivative of T wrt T is 1
        ! end if

        ! self%dnj_dstate2 = 0.0d0
        ! if (const_p) self%dn_dstate2 = 0.0d0
        ! if (.not. const_t) self%dT_dstate2 = 0.0d0

    end subroutine

    subroutine EqTotals_compute_dR_dx(self)
        ! Compute the partial derivative matrix of R wrt x

        ! Arguments
        class(EqTotals), target :: self

        ! Build the partial derivative matrix of R wrt x
        self%dR_dx = 0.0d0
        call self%compute_dR_eq_dH0()
        call self%compute_dR_eq_db0()
        call self%compute_dR_eq_dP()
        call self%compute_dR_lnnj_dP()
        call self%compute_dR_n_dP()
        call self%compute_dR_T_dP()

        call self%check_dR_dx()

    end subroutine

    subroutine EqTotals_compute_dR_eq_dH0(self)
        ! Compute the partial derivative matrix of R_eq wrt H0

        ! Arguments
        class(EqTotals), target :: self

        ! Compute the partial derivative matrix of R_eq wrt H0
        self%dR_eq_dH0(self%num_eqn) = -1.0d0/self%eq_solution%T

    end subroutine

    subroutine EqTotals_compute_dR_eq_db0(self)
        ! Compute the partial derivative matrix of R_eq wrt b0

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ne       ! Number of elements
        integer :: i        ! Loop index

        ! Define shorthand
        ne = self%eq_solver%num_elements

        ! Compute the partial derivative matrix of R_eq wrt b0
        do i = 1, ne
            self%dR_eq_db0(i, i) = -1.0d0
        end do

    end subroutine

    subroutine EqTotals_compute_dR_eq_dP(self)
        ! Compute the partial derivative matrix of R_eq wrt P

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: dG_dP(self%num_eqn, self%num_eqn+1) ! Partial derivative of G wrt P
        real(dp) :: dA_dP(self%num_eqn, self%num_eqn)   ! Partial derivative of A wrt P
        real(dp) :: df_dP(self%num_eqn)                 ! Partial derivative of f wrt P

        call self%dG_dP(dG_dP)
        dA_dP = dG_dP(:, :self%num_eqn)
        df_dP = dG_dP(:, self%num_eqn+1)
        self%dR_eq_dP(:) = matmul(dA_dP, self%u_hat) - df_dP

    end subroutine

    subroutine EqTotals_compute_dR_lnnj_dP(self)
        ! Compute the partial derivative matrix of R_lnnj wrt P

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: lambda  ! Damped update factor
        real(dp) :: P
        real(dp), pointer :: dln_nj(:)
        real(dp) :: dlambda_dP ! Partial derivative of lambda wrt P
        real(dp) :: ddlnnj_dP(self%eq_solver%num_gas)  ! Partial derivative of dln(nj) wrt P

        ! Define shorthand
        dln_nj => self%eq_solution%dln_nj
        lambda = self%eq_solution%lambda
        P = self%eq_solution%calc_pressure()

        ! dR_lnnj_dP = -(dlambda_dP*dlnnj + lambda*ddlnnj_dP)
        ! dlambda_dP = dlambda_ddlnnj*ddlnnj_dP
        call self%fd_dlambda_dP(dlambda_dP)  ! ** TODO **
        ddlnnj_dP = -1.0d0/P

        self%dR_lnnj_dP = -(dlambda_dP*dln_nj + lambda*ddlnnj_dP)

    end subroutine

    subroutine EqTotals_compute_dR_n_dP(self)
        ! Compute the partial derivative of R_n wrt P

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: lambda  ! Damped update factor
        real(dp) :: dln_n, n
        real(dp) :: dlambda_dP ! Partial derivative of lambda wrt P

        ! Define shorthand
        lambda = self%eq_solution%lambda
        dln_n = self%eq_solution%dln_n
        n = self%eq_solution%n

        ! dR_n_dP = -(dlambda_dP*dln_n)*n
        ! dlambda_dP = dlambda_ddlnn*ddlnn_dP
        call self%fd_dlambda_dP(dlambda_dP)  ! ** TODO **
        self%dR_n_dP = -dlambda_dP*dln_n*n

    end subroutine

    subroutine EqTotals_compute_dR_T_dP(self)
        ! Compute the partial derivative of R_n wrt P

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: lambda  ! Damped update factor
        real(dp) :: dln_T, T
        real(dp) :: dlambda_dP ! Partial derivative of lambda wrt P

        ! Define shorthand
        lambda = self%eq_solution%lambda
        dln_T = self%eq_solution%dln_T
        T = self%eq_solution%T

        ! dR_T_dP = -(dlambda_dP*dln_T)*T
        ! dlambda_dP = dlambda_ddlnn*ddlnn_dP
        call self%fd_dlambda_dP(dlambda_dP)  ! ** TODO **
        self%dR_T_dP = -dlambda_dP*dln_T*T

    end subroutine

    subroutine EqTotals_compute_dR_du(self)
        ! Compute the partial derivative matrix of R wrt u

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: num_eqn
        integer :: ng
        integer :: ne
        integer :: i

        ! Define shorthand
        num_eqn = self%num_eqn
        ng = self%eq_solver%num_gas
        ne = self%eq_solver%num_elements

        ! Build the partial derivative matrix of R wrt u
        self%dR_du = 0.0d0
        call self%eq_solver%assemble_matrix(self%eq_solution)
        self%dR_eq_duhat = self%eq_solution%G(:num_eqn, :num_eqn)
        call self%compute_dR_eq_dnj()
        call self%compute_dR_eq_dn()
        call self%compute_dR_eq_dT()
        call self%compute_dR_eq_db()
        call self%compute_dR_eq_dlnnj()
        call self%compute_dR_eq_dHj()
        call self%compute_dR_eq_dSj()
        call self%compute_dR_eq_dUj()
        call self%compute_dR_eq_dCpj()

        call self%compute_dR_nj_dlnnj()
        ! Add diagonal terms
        do i = 1, ng
            self%dR_nj_dnj(i, i) = 1.0d0
        end do

        call self%compute_dR_lnnj_dpi()
        call self%compute_dR_lnnj_ddlnn()
        call self%compute_dR_lnnj_ddlnT()
        call self%compute_dR_lnnj_dn()
        call self%compute_dR_lnnj_dHj()
        call self%compute_dR_lnnj_dSj()
        call self%compute_dR_lnnj_dlnnj() ! Add diagonal terms
        ! do i = 1, ng
        !     dR_lnnj_dlnnj(i, i) = 1.0d0
        ! end do

        call self%compute_dR_n_dpi()
        call self%compute_dR_n_ddlnn()
        call self%compute_dR_n_ddlnT()
        call self%compute_dR_n_dlnnj()
        call self%compute_dR_n_dHj()
        call self%compute_dR_n_dSj()
        self%dR_n_dn = 1.0d0  ! Add diagonal term

        call self%compute_dR_T_dpi()
        call self%compute_dR_T_ddlnn()
        call self%compute_dR_T_ddlnT()
        call self%compute_dR_T_dlnnj()
        call self%compute_dR_T_dHj()
        call self%compute_dR_T_dSj()
        self%dR_T_dT = 1.0d0  ! Add diagonal term

        call self%compute_dR_b_dnj()
        ! Add diagonal terms
        do i = 1, ne
            self%dR_b_db(i, i) = 1.0d0
        end do

        call self%compute_dR_Hj_dT()
        ! Add diagonal terms
        do i = 1, ng
            self%dR_Hj_dHj(i, i) = 1.0d0
        end do

        call self%compute_dR_Sj_dT()
        ! Add diagonal terms
        do i = 1, ng
            self%dR_Sj_dSj(i, i) = 1.0d0
        end do

        call self%compute_dR_Uj_dT()
        ! Add diagonal terms
        do i = 1, ng
            self%dR_Uj_dUj(i, i) = 1.0d0
        end do

        call self%compute_dR_Cp_dT()
        ! Add diagonal terms
        do i = 1, ng
            self%dR_Cp_dCp(i, i) = 1.0d0
        end do

        ! FINITE DIFFERENCE CHECK
        call self%check_dR_du()

    end subroutine

    subroutine EqTotals_compute_dR_eq_dnj(self)
        ! Compute the partial derivative matrix of R_eq wrt nj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng                                    ! Number of gas species
        integer :: i                                     ! Index variables
        real(dp) :: dG_dnj(self%num_eqn, self%num_eqn+1) ! Partial derivative of G wrt nj
        real(dp) :: dA_dnj(self%num_eqn, self%num_eqn)   ! Partial derivative of A wrt nj
        real(dp) :: df_dnj(self%num_eqn)                 ! Partial derivative of f wrt nj

        ! Define shorthand
        ng = self%eq_solver%num_gas

        do i = 1, ng
            call self%dG_dnj(dG_dnj, i)
            dA_dnj = dG_dnj(:, :self%num_eqn)
            df_dnj = dG_dnj(:, self%num_eqn+1)
            self%dR_eq_dnj(:, i) = matmul(dA_dnj, self%u_hat) - df_dnj
        end do

    end subroutine

    subroutine EqTotals_compute_dR_eq_dn(self)
        ! Compute the partial derivative matrix of R_eq wrt n

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: dG_dn(self%num_eqn, self%num_eqn+1) ! Partial derivative of G wrt nj
        real(dp) :: dA_dn(self%num_eqn, self%num_eqn)   ! Partial derivative of A wrt nj
        real(dp) :: df_dn(self%num_eqn)                 ! Partial derivative of f wrt nj

        call self%dG_dn(dG_dn)
        dA_dn = dG_dn(:, :self%num_eqn)
        df_dn = dG_dn(:, self%num_eqn+1)
        self%dR_eq_dn = matmul(dA_dn, self%u_hat) - df_dn

    end subroutine

    subroutine EqTotals_compute_dR_eq_dT(self)
        ! Compute the partial derivative matrix of R_eq wrt T

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: state1, d_hsu_delta_dT
        logical :: const_s, const_h, const_u

        ! Define shorthand
        state1 = self%eq_solution%constraints%state1
        const_s = self%eq_solution%constraints%is_constant_entropy()
        const_h = self%eq_solution%constraints%is_constant_enthalpy()
        const_u = self%eq_solution%constraints%is_constant_energy()

        if (const_s) then
            d_hsu_delta_dT = 0.0d0
        else if (const_h) then
            d_hsu_delta_dT = -state1/(self%eq_solution%T*self%eq_solution%T)
        else if (const_u) then
            d_hsu_delta_dT = -state1/(self%eq_solution%T*self%eq_solution%T)
        end if

        self%dR_eq_dT(self%num_eqn) = -d_hsu_delta_dT

    end subroutine

    subroutine EqTotals_compute_dR_eq_db(self)
        ! Compute the partial derivative matrix of R_eq wrt b

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ne       ! Number of elements
        integer :: i

        ! Define shorthand
        ne = self%eq_solver%num_elements

        ! Compute the partial derivative matrix of R_eq wrt b
        do i = 1, ne
            self%dR_eq_db(i, i) = 1.0d0
        end do

    end subroutine

    subroutine EqTotals_compute_dR_eq_dlnnj(self)
        ! Compute the partial derivative matrix of R_eq wrt nj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng                                      ! Number of gas species
        integer :: i                                       ! Index variables
        real(dp) :: dG_dlnnj(self%num_eqn, self%num_eqn+1) ! Partial derivative of G wrt nj
        real(dp) :: dA_dlnnj(self%num_eqn, self%num_eqn)   ! Partial derivative of A wrt nj
        real(dp) :: df_dlnnj(self%num_eqn)                 ! Partial derivative of f wrt nj

        ! Define shorthand
        ng = self%eq_solver%num_gas

        do i = 1, ng
            call self%dG_dlnnj(dG_dlnnj, i)
            dA_dlnnj = dG_dlnnj(:, :self%num_eqn)
            df_dlnnj = dG_dlnnj(:, self%num_eqn+1)
            self%dR_eq_dlnnj(:, i) = matmul(dA_dlnnj, self%u_hat) - df_dlnnj
        end do

    end subroutine

    subroutine EqTotals_compute_dR_eq_dHj(self)
        ! Compute the partial derivative matrix of R_eq wrt Hj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng                                      ! Number of gas species
        integer :: i                                       ! Index variables
        real(dp) :: dG_dHj(self%num_eqn, self%num_eqn+1)   ! Partial derivative of G wrt hj
        real(dp) :: dA_dHj(self%num_eqn, self%num_eqn)     ! Partial derivative of A wrt hj
        real(dp) :: df_dHj(self%num_eqn)                   ! Partial derivative of f wrt hj

        ! Define shorthand
        ng = self%eq_solver%num_gas

        do i = 1, ng
            call self%dG_dHj(dG_dHj, i)
            dA_dHj = dG_dHj(:, :self%num_eqn)
            df_dHj = dG_dHj(:, self%num_eqn+1)
            self%dR_eq_dHj(:, i) = matmul(dA_dHj, self%u_hat) - df_dHj
        end do

    end subroutine

    subroutine EqTotals_compute_dR_eq_dSj(self)
        ! Compute the partial derivative matrix of R_eq wrt Sj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng                                      ! Number of gas species
        integer :: i                                       ! Index variables
        real(dp) :: dG_dSj(self%num_eqn, self%num_eqn+1)   ! Partial derivative of G wrt Sj
        real(dp) :: dA_dSj(self%num_eqn, self%num_eqn)     ! Partial derivative of A wrt Sj
        real(dp) :: df_dSj(self%num_eqn)                   ! Partial derivative of f wrt Sj

        ! Define shorthand
        ng = self%eq_solver%num_gas

        do i = 1, ng
            call self%dG_dSj(dG_dSj, i)
            dA_dSj = dG_dSj(:, :self%num_eqn)
            df_dSj = dG_dSj(:, self%num_eqn+1)
            self%dR_eq_dSj(:, i) = matmul(dA_dSj, self%u_hat) - df_dSj
        end do

    end subroutine

    subroutine EqTotals_compute_dR_eq_dUj(self)
        ! Compute the partial derivative matrix of R_eq wrt Uj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng                                      ! Number of gas species
        integer :: i                                       ! Index variables
        real(dp) :: dG_dUj(self%num_eqn, self%num_eqn+1)   ! Partial derivative of G wrt hj
        real(dp) :: dA_dUj(self%num_eqn, self%num_eqn)     ! Partial derivative of A wrt hj
        real(dp) :: df_dUj(self%num_eqn)                   ! Partial derivative of f wrt hj

        ! Define shorthand
        ng = self%eq_solver%num_gas

        do i = 1, ng
            call self%dG_dUj(dG_dUj, i)
            dA_dUj = dG_dUj(:, :self%num_eqn)
            df_dUj = dG_dUj(:, self%num_eqn+1)
            self%dR_eq_dUj(:, i) = matmul(dA_dUj, self%u_hat) - df_dUj
        end do

    end subroutine

    subroutine EqTotals_compute_dR_eq_dCpj(self)
        ! Compute the partial derivative matrix of R_eq wrt Cpj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species
        logical :: const_p  ! Isobaric flag
        real(dp) :: dA_dCpj(self%num_eqn, self%num_eqn)  ! Partial derivative of A wrt Sj
        integer :: i, j, k

        ! Define shorthand
        ng = self%eq_solver%num_gas
        const_p = self%eq_solution%constraints%is_constant_pressure()

        ! Define shorthand
        ng = self%eq_solver%num_gas

        ! Compute the partial derivative matrix of R_eq wrt Cpj
        if (const_p) then
            do i = 1, ng
                dA_dCpj = 0.0d0
                dA_dCpj(self%num_eqn, self%num_eqn) = self%eq_solution%nj(i)
                self%dR_eq_dCpj(:, i) = matmul(dA_dCpj, self%u_hat)
            end do
        end if

    end subroutine

    subroutine EqTotals_compute_dR_eq_dCvj(self)
        ! Compute the partial derivative matrix of R_eq wrt Cpj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species
        logical :: const_p  ! Isobaric flag
        integer :: i, j

        ! Define shorthand
        ng = self%eq_solver%num_gas
        const_p = self%eq_solution%constraints%is_constant_pressure()

        ! Compute the partial derivative matrix of R_eq wrt Cpj
        if (.not. const_p) then
            do i = 1, ng
                self%dR_eq_dCvj(self%num_eqn, i) = self%eq_solution%nj(i)*self%u_hat(self%num_eqn)
            end do
        end if

    end subroutine

    subroutine EqTotals_compute_dR_nj_dlnnj(self)
        ! Compute the partial derivative matrix of R_nj wrt ln(nj)

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species
        integer :: i

        ! Define shorthand
        ng = self%eq_solver%num_gas

        ! Compute the partial derivative matrix of R_nj wrt ln(nj)
        do i = 1, ng
            self%dR_nj_dlnnj(i, i) = -exp(self%eq_solution%ln_nj(i))
        end do

    end subroutine

    subroutine EqTotals_compute_dR_lnnj_dpi(self)
        ! Compute the partial derivative matrix of R_lnnj wrt pi

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species

        ! Define shorthand
        ng = self%eq_solver%num_gas

        ! Compute the partial derivative matrix of R_lnnj wrt pi
        self%dR_lnnj_dpi = self%eq_solver%products%stoich_matrix(:ng,:)

    end subroutine

    subroutine EqTotals_compute_dR_lnnj_ddlnn(self)
        ! Compute the partial derivative matrix of R_lnnj wrt ln(n)

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species
        integer :: ne       ! Number of elements
        integer :: i
        real(dp) :: dlambda_ddlnn  ! Partial of the damped update factor wrt dln(n)
        real(dp) :: ddlnnj_ddlnn(self%eq_solver%num_gas)  !

        ! Define shorthand
        ng = self%eq_solver%num_gas
        ne = self%eq_solver%num_elements

        dlambda_ddlnn = 1.0d0  ! TODO

        ! Compute the partial derivative matrix of R_lnnj wrt ln(n)
        ddlnnj_ddlnn = -1.0d0

        self%dR_lnnj_ddlnn = dlambda_ddlnn*self%eq_solution%dln_n + self%eq_solution%lambda*ddlnnj_ddlnn

    end subroutine

    subroutine EqTotals_compute_dR_lnnj_ddlnT(self)
        ! Compute the partial derivative matrix of R_lnnj wrt ln(T)

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species
        integer :: i

        ! Define shorthand
        ng = self%eq_solver%num_gas

        ! Compute the partial derivative matrix of R_lnnj wrt ln(T)
        do i = 1, ng
            self%dR_lnnj_ddlnT(i) = -self%eq_solution%thermo%enthalpy(i)
        end do

    end subroutine

    subroutine EqTotals_compute_dR_lnnj_dlnnj(self)
        ! Compute the partial derivative matrix of R_lnnj wrt lnnj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species
        integer :: i, j
        real(dp) :: lambda
        real(dp) :: n
        real(dp), pointer :: dln_nj(:)
        real(dp) :: dlambda_dlnnj(self%eq_solver%num_gas)

        ! Define shorthand
        ng = self%eq_solver%num_gas
        lambda = self%eq_solution%lambda
        n = self%eq_solution%n
        dln_nj => self%eq_solution%dln_nj

        ! Compute the partial derivative of R_lnnj wrt n:
        ! dR_lnnj/dlnnj = 1.0 - [d(lambda)/dln(nj) * dln_nj - lambda]
        call self%fd_dlambda_dlnnj(dlambda_dlnnj)
        do i = 1, ng
            self%dR_lnnj_dlnnj(i, i) = 1.0d0 - (dlambda_dlnnj(i)*dln_nj(i) - lambda)
        end do

    end subroutine

    subroutine EqTotals_compute_dR_lnnj_dn(self)
        ! Compute the partial derivative matrix of R_lnnj wrt n

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: lambda
        real(dp) :: n
        real(dp), pointer :: dln_nj(:)
        real(dp) :: dlambda_dn

        ! Define shorthand
        lambda = self%eq_solution%lambda
        n = self%eq_solution%n
        dln_nj => self%eq_solution%dln_nj

        ! Compute the partial derivative of R_lnnj wrt n:
        ! dR_lnnj/dn = -[d(lambda)/dn * dln_nj + lambda/n]
        call self%fd_dlambda_dn(dlambda_dn)
        self%dR_lnnj_dn(:) = -(dlambda_dn*dln_nj + lambda/n)

    end subroutine

    subroutine EqTotals_compute_dR_lnnj_dHj(self)
        ! Compute the partial derivative matrix of R_lnnj wrt Hj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species
        integer :: i, j

        ! Define shorthand
        ng = self%eq_solver%num_gas

        ! Compute the partial derivative matrix of R_lnnj wrt Hj
        do i = 1, ng
            self%dR_lnnj_dHj(i, i) = self%eq_solution%lambda*(1.0d0 - self%eq_solution%dln_T)
        end do

    end subroutine

    subroutine EqTotals_compute_dR_lnnj_dSj(self)
        ! Compute the partial derivative matrix of R_lnnj wrt Sj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species
        integer :: i, j

        ! Define shorthand
        ng = self%eq_solver%num_gas

        ! Compute the partial derivative matrix of R_lnnj wrt Sj
        do i = 1, ng
            self%dR_lnnj_dSj(i, i) = -1.0d0*self%eq_solution%lambda
        end do

    end subroutine

    subroutine EqTotals_compute_dR_n_dpi(self)
        ! Compute the partial derivative matrix of R_n wrt pi

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: lambda
        real(dp) :: n
        real(dp) :: dln_n
        real(dp) :: dlambda_dpi(self%eq_solver%num_elements)
        real(dp) :: ddln_n_dpi(self%eq_solver%num_elements)

        lambda = self%eq_solution%lambda
        n = self%eq_solution%n
        dln_n = self%eq_solution%dln_n

        ! dR_n/dpi = -(dlambda/dpi*dlnn + lambda*dlnn/dpi)*n
        call self%eq_solver%assemble_matrix(self%eq_solution)
        call self%fd_dlambda_dpi(dlambda_dpi)  ! TODO
        ddln_n_dpi = self%eq_solution%G(:self%eq_solver%num_elements, self%eq_solver%num_elements+1)
        self%dR_n_dpi = -(dlambda_dpi*dln_n + lambda*ddln_n_dpi)*n

    end subroutine

    subroutine EqTotals_compute_dR_n_ddlnn(self)
        ! Compute the partial derivative matrix of R_n wrt dln(n)

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: lambda
        real(dp) :: n
        real(dp) :: dln_n
        real(dp) :: dlambda_ddlnn

        lambda = self%eq_solution%lambda
        n = self%eq_solution%n
        dln_n = self%eq_solution%dln_n

        ! dR_n/ddlnn = -(dlambda/ddlnn*dlnn + lambda)*n
        call self%fd_dlambda_ddlnn(dlambda_ddlnn)  ! TODO
        self%dR_n_ddlnn = -(dlambda_ddlnn*dln_n + lambda)*n

    end subroutine

    subroutine EqTotals_compute_dR_n_ddlnT(self)
        ! Compute the partial derivative matrix of R_n wrt dlnT

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: lambda
        real(dp) :: n
        real(dp) :: dln_n
        real(dp) :: dlambda_ddlnT
        real(dp) :: ddln_n_ddlnT

        lambda = self%eq_solution%lambda
        n = self%eq_solution%n
        dln_n = self%eq_solution%dln_n

        ! dR_n/ddlnT = -(dlambda/ddlnT*dlnn + lambda*dlnn/ddlnT)*n
        call self%eq_solver%assemble_matrix(self%eq_solution)
        call self%fd_dlambda_ddlnT(dlambda_ddlnT)  ! TODO
        ddln_n_ddlnT = self%eq_solution%G(self%eq_solver%num_elements+2, self%eq_solver%num_elements+1)
        self%dR_n_ddlnT = -(dlambda_ddlnT*dln_n + lambda*ddln_n_ddlnT)*n

    end subroutine

    subroutine EqTotals_compute_dR_n_dlnnj(self)
        ! Compute the partial derivative matrix of R_n wrt lnnj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: n
        real(dp) :: dln_n
        real(dp) :: dlambda_dlnnj(self%eq_solver%num_gas)

        n = self%eq_solution%n
        dln_n = self%eq_solution%dln_n

        ! dR_n/dlnnj = -(dlambda/dlnnj*dlnn)*n
        call self%fd_dlambda_dlnnj(dlambda_dlnnj)  ! TODO
        self%dR_n_dlnnj = -(dlambda_dlnnj*dln_n)*n

    end subroutine

    subroutine EqTotals_compute_dR_n_dHj(self)
        ! Compute the partial derivative matrix of R_n wrt Hj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: n
        real(dp) :: dln_n
        real(dp) :: dlambda_dHj(self%eq_solver%num_gas)

        n = self%eq_solution%n
        dln_n = self%eq_solution%dln_n

        ! dR_n/dHj = -(dlambda/dHj*dlnn)*n
        call self%fd_dlambda_dHj(dlambda_dHj)  ! TODO
        self%dR_n_dHj = -(dlambda_dHj*dln_n)*n

    end subroutine

    subroutine EqTotals_compute_dR_n_dSj(self)
        ! Compute the partial derivative matrix of R_n wrt Sj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: n
        real(dp) :: dln_n
        real(dp) :: dlambda_dSj(self%eq_solver%num_gas)

        n = self%eq_solution%n
        dln_n = self%eq_solution%dln_n

        ! dR_n/dSj = -(dlambda/dSj*dlnn)*n
        call self%fd_dlambda_dSj(dlambda_dSj)  ! TODO
        self%dR_n_dSj = -(dlambda_dSj*dln_n)*n

    end subroutine

    subroutine EqTotals_compute_dR_T_dpi(self)
        ! Compute the partial derivative matrix of R_T wrt pi

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: lambda
        real(dp) :: T
        real(dp) :: dln_T
        real(dp) :: dlambda_dpi(self%eq_solver%num_elements)
        real(dp) :: ddln_T_dpi(self%eq_solver%num_elements)

        lambda = self%eq_solution%lambda
        T = self%eq_solution%T
        dln_T = self%eq_solution%dln_T

        ! dR_T/dpi = -(dlambda/dpi*dlnT + lambda*dlnT/dpi)*T
        call self%eq_solver%assemble_matrix(self%eq_solution)
        call self%fd_dlambda_dpi(dlambda_dpi)  ! TODO
        ddln_T_dpi = self%eq_solution%G(:self%eq_solver%num_elements, self%eq_solver%num_elements+2)
        self%dR_T_dpi = -(dlambda_dpi*dln_T + lambda*ddln_T_dpi)*T

    end subroutine

    subroutine EqTotals_compute_dR_T_ddlnn(self)
        ! Compute the partial derivative matrix of R_T wrt dln(n)

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: lambda
        real(dp) :: T
        real(dp) :: dln_T
        real(dp) :: dlambda_ddlnn
        real(dp) :: ddlnT_ddlnn

        lambda = self%eq_solution%lambda
        T = self%eq_solution%T
        dln_T = self%eq_solution%dln_T

        ! dR_T/ddlnn = -(dlambda/ddlnn*dlnT + lambda*ddlnnT/ddlnn)*T
        call self%eq_solver%assemble_matrix(self%eq_solution)
        call self%fd_dlambda_ddlnn(dlambda_ddlnn)  ! TODO
        ddlnT_ddlnn = self%eq_solution%G(self%eq_solver%num_elements+2, self%eq_solver%num_elements+1)
        self%dR_T_ddlnn = -(dlambda_ddlnn*dln_T + lambda*ddlnT_ddlnn)*T

    end subroutine

    subroutine EqTotals_compute_dR_T_ddlnT(self)
        ! Compute the partial derivative matrix of R_T wrt dlnT

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: lambda
        real(dp) :: T
        real(dp) :: dln_T
        real(dp) :: dlambda_ddlnT

        lambda = self%eq_solution%lambda
        T = self%eq_solution%T
        dln_T = self%eq_solution%dln_T

        ! dR_T/ddlnT = -(dlambda/ddlnT*dlnT + lambda)*T
        call self%fd_dlambda_ddlnT(dlambda_ddlnT)  ! TODO
        self%dR_T_ddlnT = -(dlambda_ddlnT*dln_T + lambda)*T

    end subroutine

    subroutine EqTotals_compute_dR_T_dlnnj(self)
        ! Compute the partial derivative matrix of R_T wrt lnnj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: T
        real(dp) :: dln_T
        real(dp) :: dlambda_dlnnj(self%eq_solver%num_gas)

        T = self%eq_solution%n
        dln_T = self%eq_solution%dln_n

        ! dR_T/dlnnj = -(dlambda/dlnnj*dlnT)*T
        call self%fd_dlambda_dlnnj(dlambda_dlnnj)  ! TODO
        self%dR_T_dlnnj = -(dlambda_dlnnj*dln_T)*T

    end subroutine

    subroutine EqTotals_compute_dR_T_dn(self)
        ! Compute the partial derivative matrix of R_T wrtdn

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: T
        real(dp) :: dln_T
        real(dp) :: dlambda_dn

        T = self%eq_solution%n
        dln_T = self%eq_solution%dln_n

        ! dR_T/dlnnj = -(dlambda/dn*dlnT)*T
        call self%fd_dlambda_dn(dlambda_dn)  ! TODO
        self%dR_T_dn = -(dlambda_dn*dln_T)*T

    end subroutine

    subroutine EqTotals_compute_dR_T_dHj(self)
        ! Compute the partial derivative matrix of R_T wrt Hj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: T
        real(dp) :: dln_T
        real(dp) :: dlambda_dHj(self%eq_solver%num_gas)

        T = self%eq_solution%n
        dln_T = self%eq_solution%dln_T

        ! dR_T/dHj = -(dlambda/dHj*dlnT)*T
        call self%fd_dlambda_dHj(dlambda_dHj)  ! TODO
        self%dR_T_dHj = -(dlambda_dHj*dln_T)*T

    end subroutine

    subroutine EqTotals_compute_dR_T_dSj(self)
        ! Compute the partial derivative matrix of R_T wrt Sj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        real(dp) :: T
        real(dp) :: dln_T
        real(dp) :: dlambda_dSj(self%eq_solver%num_gas)

        T = self%eq_solution%n
        dln_T = self%eq_solution%dln_T

        ! dR_T/dSj = -(dlambda/dSj*dlnT)*T
        call self%fd_dlambda_dSj(dlambda_dSj)  ! TODO
        self%dR_T_dSj = -(dlambda_dSj*dln_T)*T

    end subroutine

    subroutine EqTotals_compute_dR_b_dnj(self)
        ! Compute the partial derivative matrix of R_b wrt nj

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species

        ! Define shorthand
        ng = self%eq_solver%num_gas

        ! Compute the partial derivative matrix of R_b wrt nj
        self%dR_b_dnj = -transpose(self%eq_solver%products%stoich_matrix(:ng,:))

    end subroutine

    subroutine EqTotals_compute_dR_Hj_dT(self)
        ! Compute the partial derivative matrix of R_Hj wrt T

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species
        integer :: i

        ! Define shorthand
        ng = self%eq_solver%num_gas

        ! Compute the partial derivative matrix of R_Hj wrt T
        do i = 1, ng
            self%dR_Hj_dT(i) = -self%eq_solver%products%species(i)%calc_denthalpy_dT(self%eq_solution%T)
        end do

    end subroutine

    subroutine EqTotals_compute_dR_Sj_dT(self)
        ! Compute the partial derivative matrix of R_Sj wrt T

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species
        integer :: i

        ! Define shorthand
        ng = self%eq_solver%num_gas

        ! Compute the partial derivative matrix of R_Sj wrt T
        do i = 1, ng
            self%dR_Sj_dT(i) = -self%eq_solver%products%species(i)%calc_dentropy_dT(self%eq_solution%T)
        end do

    end subroutine

    subroutine EqTotals_compute_dR_Uj_dT(self)
        ! Compute the partial derivative matrix of R_Uj wrt T

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species
        integer :: i

        ! Define shorthand
        ng = self%eq_solver%num_gas

        ! Compute the partial derivative matrix of R_Uj wrt T
        do i = 1, ng
            self%dR_Uj_dT(i) = -self%eq_solver%products%species(i)%calc_denergy_dT(self%eq_solution%T)
        end do

    end subroutine

    subroutine EqTotals_compute_dR_Cp_dT(self)
        ! Compute the partial derivative matrix of R_Cp wrt T

        ! Arguments
        class(EqTotals), target :: self

        ! Locals
        integer :: ng       ! Number of gas species
        integer :: i

        ! Define shorthand
        ng = self%eq_solver%num_gas

        ! Compute the partial derivative matrix of R_Cp wrt T
        do i = 1, ng
            self%dR_Cp_dT(i) = -self%eq_solver%products%species(i)%calc_dcp_dT(self%eq_solution%T)
        end do

    end subroutine

    subroutine EqTotals_compute_db0_dw(self, reactant_weights, db0_dw)
        ! Compute the derivative of b0 wrt the reactant weights
        ! This is used to compute the total derivatives of the equilibrium solution
        ! Computes derivatives of: nj, n, and T with respect to state1, state2, and reactant weights

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(in) :: reactant_weights(:)
        real(dp), intent(inout), allocatable :: db0_dw(:,:)

        ! Locals
        integer  :: nr        ! Number of reactants
        integer  :: ne        ! Number of elements
        real(dp) :: mw        ! Molecular weight
        integer :: i, j       ! Indices
        real(dp) :: dmoles_dw ! Derivative of moles wrt weights

        ! Define shorthand
        nr = self%eq_solver%num_reactants
        ne = self%eq_solver%num_elements

        ! Compute derivative of b0 wrt weights, d(b0)/dw
        db0_dw = 0.0d0
        do i = 1, ne
            do j = 1, nr
                mw = self%eq_solver%reactants%species(j)%molecular_weight
                dmoles_dw = (sum(reactant_weights) - reactant_weights(i)) / (mw*sum(reactant_weights)**2)
                db0_dw(i, j) = (self%eq_solver%reactants%stoich_matrix(j, i)**2-sum(self%eq_solver%reactants%stoich_matrix(:, i)))*dmoles_dw
            end do
        end do

    end subroutine

    ! ---------------------------------------------------------------------------
    ! Helper functions: partials of assemble_matrix
    ! ---------------------------------------------------------------------------

    subroutine EqTotals_dG_dP(self, dG_dP)

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout), target :: dG_dP(:,:)

        ! Locals
        integer  :: ng                          ! Number of gas species
        integer  :: nc                          ! Number of condensed species
        integer  :: na                          ! Number of active condensed species
        integer  :: ne                          ! Number of elements
        integer  :: num_eqn                     ! Active number of equations
        real(dp) :: tmp(self%eq_solver%num_gas)   ! Common sub-expression storage
        real(dp) :: d_tmp_dP(self%eq_solver%num_gas)   ! Common sub-expression storage
        real(dp) :: mu_g(self%eq_solver%num_gas)  ! Gas phase chemical potentials [unitless]
        real(dp) :: d_mu_g_dP  !
        real(dp) :: d_hsu_delta_dP             ! Residual for enthalpy / entropy constraint
        real(dp) :: n                           ! Total moles of mixture
        real(dp) :: P                           ! Pressure of mixture (bar)
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)           ! Log of gas species concentrations [kmol-per-kg]
        real(dp), pointer :: cp(:), cv(:)       ! Species heat capacities [unitless]
        real(dp), pointer :: h_g(:), h_c(:)     ! Gas/condensed enthalpies [unitless]
        real(dp), pointer :: s_g(:), s_c(:)     ! Gas/condensed entropies [unitless]
        real(dp), pointer :: u_g(:), u_c(:)     ! Gas/condensed energies [unitless]
        real(dp), pointer :: A_g(:,:), A_c(:,:) ! Gas/condensed stoichiometric matrices
        real(dp), pointer :: G(:,:)             ! Augmented Newton iteration matrix
        integer :: r, c                         ! Iteration matrix row/column indices
        integer :: i, j                         ! Loop counters
        logical :: const_p, const_t, const_s, const_h, const_u  ! Flags enabling/disabling matrix equations
        logical :: rhs_converged                ! Flag for RHS convergence

        ! Define shorthand
        ng = self%eq_solver%num_gas
        nc = self%eq_solver%num_condensed
        ne = self%eq_solver%num_elements
        na = count(self%eq_solution%is_active)
        num_eqn = self%num_eqn
        const_p = self%eq_solution%constraints%is_constant_pressure()
        const_t = self%eq_solution%constraints%is_constant_temperature()
        const_s = self%eq_solution%constraints%is_constant_entropy()
        const_h = self%eq_solution%constraints%is_constant_enthalpy()
        const_u = self%eq_solution%constraints%is_constant_energy()

        ! Associate subarray pointers
        G   => dG_dP
        A_g => self%eq_solver%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A_c => self%eq_solver%products%stoich_matrix(ng+1:,:)
        n = self%eq_solution%n
        nj  => self%eq_solution%nj
        nj_g => self%eq_solution%nj(:ng)
        ln_nj => self%eq_solution%ln_nj
        cp  => self%eq_solution%thermo%cp
        cv  => self%eq_solution%thermo%cv
        h_g => self%eq_solution%thermo%enthalpy(:ng)
        h_c => self%eq_solution%thermo%enthalpy(ng+1:)
        s_g => self%eq_solution%thermo%entropy(:ng)
        s_c => self%eq_solution%thermo%entropy(ng+1:)
        u_g => self%eq_solution%thermo%energy(:ng)
        u_c => self%eq_solution%thermo%energy(ng+1:)

        ! Get the mixture pressure
        P = self%eq_solution%calc_pressure()

        ! Compute gas phase chemical potentials
        mu_g = h_g - s_g + ln_nj + log(P/n)
        d_mu_g_dP = 1.0d0/P

        if (const_s) then
            d_hsu_delta_dP = (-1.0d0/P)*sum(nj_g)
        else if (const_h) then
            d_hsu_delta_dP = 0.0d0
        else if (const_u) then
            d_hsu_delta_dP = 0.0d0
        end if

        ! Initialize the iteration matrix
        G = 0.0d0
        r = 0
        c = 0

        ! Set a flag if the entropy and temperature residuals are already converged
        rhs_converged = .false.

        !-------------------------------------------------------
        ! Equation (2.24/2.45): Element constraints
        !-------------------------------------------------------
        do i = 1,ne
            tmp = nj_g*A_g(:,i)
            r = r+1
            c = 0

            ! Pi derivatives
            do j = 1,ne
                c = c+1
            end do

            ! Condensed derivatives
            ! Handled in (2.25) below (symmetric)
            c = c+na

            ! Delta ln(n) derivative
            ! Symmetric with (2.26) pi derivative
            if (const_p) then
                c = c+1
            end if

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
            end if

            ! Right hand side
            if (.not. rhs_converged) then
                G(r,c+1) = d_mu_g_dP*sum(tmp)
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.25/2.46): Condensed phase constraints
        !-------------------------------------------------------
        do i = 1,nc
            if (.not. self%eq_solution%is_active(i)) cycle
            r = r+1
            c = 0

            ! Pi derivatives
            ! Symmetric with (2.24) condensed derivatives
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            ! Handled via zero initialization
            if (const_p) c = c+1

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.26)
        !-------------------------------------------------------
        if (const_p) then
            r = r+1
            c = 0

            ! Pi derivatives
            ! Handled in Eq (2.24) above.
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            c = c+1

            ! Delta ln(T) derviative
            if (.not. const_t) then
                c = c+1
            end if

            ! Right-hand-side
            if (.not. rhs_converged) then
                G(r,c+1) = sum(nj_g)*d_mu_g_dP
            end if

        end if

        !---------------------------------------------------------
        ! Equation (2.27)/(2.28)/(2.47)/(2.48): Energy constraints
        !---------------------------------------------------------
        if (.not. const_t) then
            r = r+1
            c = 0

            ! Select entropy/enthalpy constraint
            if (const_s) then
                tmp = nj_g*(s_g-ln_nj-log(P/n))
                d_tmp_dP = nj_g*(-1.0d0/P)
            else if (const_h) then
                tmp = nj_g*h_g
                d_tmp_dP = 0.0d0
            else if (const_u) then
                tmp = nj_g*u_g
                d_tmp_dP = 0.0d0
            end if

            ! Pi derivatives
            if (.not. rhs_converged) then
                do j = 1,ne
                    c = c+1
                    G(r,c) = dot_product(d_tmp_dP, A_g(:,j))
                end do
            end if

            ! Condensed derivatives
            if (.not. rhs_converged) then
                do j = 1,nc
                    if (.not. self%eq_solution%is_active(j)) cycle
                    c = c+1
                end do
            end if

            ! Delta ln(n) derivative
            if (.not. rhs_converged) then
                if (const_p) then
                    c = c+1
                    G(r,c) = sum(d_tmp_dP)
                end if
            end if

            ! Delta ln(T) derivative
            c = c+1
            if (const_p) then
                if (.not. rhs_converged) then
                    G(r,c) = dot_product(d_tmp_dP, h_g)
                else
                    c = c + ne + na + 1
                end if
            else
                G(r,c) = dot_product(d_tmp_dP, u_g)
            end if

            ! Right-hand-side
            if (.not. rhs_converged) then
                G(r,c+1) = d_hsu_delta_dP + dot_product(d_tmp_dP, mu_g) + d_mu_g_dP*sum(tmp)
                if (const_s) then
                    if (.not. const_p) then
                        G(r,c+1) = G(r,c+1) - d_mu_g_dP*sum(nj_g)
                    end if
                end if
            end if

        end if

    end subroutine

    subroutine EqTotals_dG_dnj(self, dG_dnj, k)

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout), target :: dG_dnj(:,:)
        integer, intent(in) :: k  ! Index of nj to differentiate with respect to

        ! Locals
        integer  :: ng                          ! Number of gas species
        integer  :: nc                          ! Number of condensed species
        integer  :: na                          ! Number of active condensed species
        integer  :: ne                          ! Number of elements
        integer  :: num_eqn                     ! Active number of equations
        real(dp) :: tmp                         ! Common sub-expression storage
        real(dp) :: mu_g(self%eq_solver%num_gas)  ! Gas phase chemical potentials [unitless]
        real(dp) :: d_n_delta_dnj               ! Residual for total moles / pressure constraint
        real(dp) :: d_hsu_delta_dnj             ! Residual for enthalpy / entropy constraint
        real(dp) :: n                           ! Total moles of mixture
        real(dp) :: P                           ! Pressure of mixture (bar)
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)           ! Log of gas species concentrations [kmol-per-kg]
        real(dp), pointer :: cp(:), cv(:)       ! Species heat capacities [unitless]
        real(dp), pointer :: h_g(:), h_c(:)     ! Gas/condensed enthalpies [unitless]
        real(dp), pointer :: s_g(:), s_c(:)     ! Gas/condensed entropies [unitless]
        real(dp), pointer :: u_g(:), u_c(:)     ! Gas/condensed energies [unitless]
        real(dp), pointer :: A_g(:,:), A_c(:,:) ! Gas/condensed stoichiometric matrices
        real(dp), pointer :: G(:,:)             ! Augmented Newton iteration matrix
        real(dp), pointer :: h_or_s_or_u(:)     ! For evaluating Eq 2.27/2.28
        integer :: r, c                         ! Iteration matrix row/column indices
        integer :: i, j                         ! Loop counters
        logical :: const_p, const_t, const_s, const_h, const_u  ! Flags enabling/disabling matrix equations
        logical :: rhs_converged                ! Flag for RHS convergence

        ! Define shorthand
        ng = self%eq_solver%num_gas
        nc = self%eq_solver%num_condensed
        ne = self%eq_solver%num_elements
        na = count(self%eq_solution%is_active)
        num_eqn = self%num_eqn
        const_p = self%eq_solution%constraints%is_constant_pressure()
        const_t = self%eq_solution%constraints%is_constant_temperature()
        const_s = self%eq_solution%constraints%is_constant_entropy()
        const_h = self%eq_solution%constraints%is_constant_enthalpy()
        const_u = self%eq_solution%constraints%is_constant_energy()

        ! Associate subarray pointers
        G   => dG_dnj
        A_g => self%eq_solver%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A_c => self%eq_solver%products%stoich_matrix(ng+1:,:)
        n = self%eq_solution%n
        nj  => self%eq_solution%nj
        nj_g => self%eq_solution%nj(:ng)
        ln_nj => self%eq_solution%ln_nj
        cp  => self%eq_solution%thermo%cp
        cv  => self%eq_solution%thermo%cv
        h_g => self%eq_solution%thermo%enthalpy(:ng)
        h_c => self%eq_solution%thermo%enthalpy(ng+1:)
        s_g => self%eq_solution%thermo%entropy(:ng)
        s_c => self%eq_solution%thermo%entropy(ng+1:)
        u_g => self%eq_solution%thermo%energy(:ng)
        u_c => self%eq_solution%thermo%energy(ng+1:)

        ! Get the mixture pressure
        P = self%eq_solution%calc_pressure()

        ! Compute gas phase chemical potentials
        mu_g = h_g - s_g + ln_nj + log(P/n)

        if (const_s) then
            d_hsu_delta_dnj = 0.d0 !TODO: (-self%calc_derivative_entropy_sum_wrt_nj(k))
        else if (const_h) then
            d_hsu_delta_dnj = (-self%eq_solution%thermo%enthalpy(k))
        else if (const_u) then
            d_hsu_delta_dnj = (-self%eq_solution%thermo%energy(k))
        end if

        d_n_delta_dnj = -1.0d0

        ! Initialize the iteration matrix
        G = 0.0d0
        r = 0
        c = 0

        ! Set a flag if the entropy and temperature residuals are already converged
        rhs_converged = .false.

        !-------------------------------------------------------
        ! Equation (2.24/2.45): Element constraints
        !-------------------------------------------------------
        do i = 1,ne
            r = r+1
            c = 0

            ! Pi derivatives
            do j = 1,ne
                c = c+1
                G(r,c) = A_g(k,i)*A_g(k,j)
            end do

            ! Condensed derivatives
            ! Handled in (2.25) below (symmetric)
            c = c+na

            ! Delta ln(n) derivative
            ! Symmetric with (2.26) pi derivative
            if (const_p) then
                c = c+1
                G(r,c) = A_g(k,i)
                G(c,r) = G(r,c)
            end if

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
                if (const_p) then
                    G(r,c) = A_g(k,i)*h_g(k)
                else
                    G(r,c) = A_g(k,i)*u_g(k)
                end if
            end if

            ! Right hand side
            if (.not. rhs_converged) then
                G(r,c+1) = A_g(k,i)*mu_g(k)
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.25/2.46): Condensed phase constraints
        !-------------------------------------------------------
        do i = 1,nc
            if (.not. self%eq_solution%is_active(i)) cycle
            r = r+1
            c = 0

            ! Pi derivatives
            ! Symmetric with (2.24) condensed derivatives
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            ! Handled via zero initialization
            if (const_p) c = c+1

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.26)
        !-------------------------------------------------------
        if (const_p) then
            r = r+1
            c = 0

            ! Pi derivatives
            ! Handled in Eq (2.24) above.
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            c = c+1
            G(r,c) = -d_n_delta_dnj

            ! Delta ln(T) derviative
            if (.not. const_t) then
                c = c+1
                G(r,c) = h_g(k)
            end if

            ! Right-hand-side
            if (.not. rhs_converged) then
                G(r,c+1) = d_n_delta_dnj + mu_g(k)
            end if

        end if

        !---------------------------------------------------------
        ! Equation (2.27)/(2.28)/(2.47)/(2.48): Energy constraints
        !---------------------------------------------------------
        if (.not. const_t) then
            r = r+1
            c = 0

            ! Select entropy/enthalpy constraint
            if (const_s) then
                tmp = (s_g(k)-ln_nj(k)-log(P/n))
                h_or_s_or_u => self%eq_solution%thermo%entropy(ng+1:)
            else if (const_h) then
                tmp = h_g(k)
                h_or_s_or_u => self%eq_solution%thermo%enthalpy(ng+1:)
            else if (const_u) then
                tmp = u_g(k)
                h_or_s_or_u => self%eq_solution%thermo%energy(ng+1:)
            end if

            ! Pi derivatives
            if (.not. rhs_converged) then
                do j = 1,ne
                    c = c+1
                    G(r,c) = tmp*A_g(k,j)
                    if (.not. const_p .and. const_s) then
                        G(r,c) = G(r,c) - A_g(k, j)
                    end if
                end do
            end if

            ! Condensed derivatives
            if (.not. rhs_converged) then
                do j = 1,nc
                    if (.not. self%eq_solution%is_active(j)) cycle
                    c = c+1
                    G(r,c) = h_or_s_or_u(j)
                end do
            end if

            ! Delta ln(n) derivative
            if (.not. rhs_converged) then
                if (const_p) then
                    c = c+1
                    G(r,c) = tmp
                end if
            end if

            ! Delta ln(T) derivative
            c = c+1
            if (const_p) then
                if (.not. rhs_converged) then
                    G(r,c) = cp(k) + tmp*h_g(k)
                else
                    c = c + ne + na + 1
                    G(r,c) = cp(k) + h_g(k)*h_g(k)
                end if
            else
                G(r,c) = cv(k) + tmp*u_g(k)
                if (const_s) then
                    G(r,c) = G(r,c) - u_g(k)
                end if
            end if

            ! Right-hand-side
            if (.not. rhs_converged) then
                G(r,c+1) = d_hsu_delta_dnj + tmp*mu_g(k)
                if (const_s) then
                    if (const_p) then
                        G(r,c+1) = G(r,c+1) + d_n_delta_dnj
                    else
                        G(r,c+1) = G(r,c+1) - mu_g(k)
                    end if
                end if
            end if

        end if

    end subroutine

    subroutine EqTotals_dG_dn(self, dG_dn)

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout), target :: dG_dn(:,:)

        ! Locals
        integer  :: ng                          ! Number of gas species
        integer  :: nc                          ! Number of condensed species
        integer  :: na                          ! Number of active condensed species
        integer  :: ne                          ! Number of elements
        integer  :: num_eqn                     ! Active number of equations
        real(dp) :: tmp(self%eq_solver%num_gas)      ! Common sub-expression storage
        real(dp) :: d_tmp_dn(self%eq_solver%num_gas) !
        real(dp) :: mu_g(self%eq_solver%num_gas)  ! Gas phase chemical potentials [unitless]
        real(dp) :: d_mu_g_dn                   !
        real(dp) :: d_n_delta_dn                !
        real(dp) :: d_hsu_delta_dn              !
        real(dp) :: n                           ! Total moles of mixture
        real(dp) :: P                           ! Pressure of mixture (bar)
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)           ! Log of gas species concentrations [kmol-per-kg]
        real(dp), pointer :: cp(:), cv(:)       ! Species heat capacities [unitless]
        real(dp), pointer :: h_g(:), h_c(:)     ! Gas/condensed enthalpies [unitless]
        real(dp), pointer :: s_g(:), s_c(:)     ! Gas/condensed entropies [unitless]
        real(dp), pointer :: u_g(:), u_c(:)     ! Gas/condensed energies [unitless]
        real(dp), pointer :: A_g(:,:), A_c(:,:) ! Gas/condensed stoichiometric matrices
        real(dp), pointer :: G(:,:)             ! Augmented Newton iteration matrix
        integer :: r, c                         ! Iteration matrix row/column indices
        integer :: i, j                         ! Loop counters
        logical :: const_p, const_t, const_s, const_h, const_u  ! Flags enabling/disabling matrix equations
        logical :: rhs_converged                ! Flag for RHS convergence

        ! Define shorthand
        ng = self%eq_solver%num_gas
        nc = self%eq_solver%num_condensed
        ne = self%eq_solver%num_elements
        na = count(self%eq_solution%is_active)
        num_eqn = self%num_eqn
        const_p = self%eq_solution%constraints%is_constant_pressure()
        const_t = self%eq_solution%constraints%is_constant_temperature()
        const_s = self%eq_solution%constraints%is_constant_entropy()
        const_h = self%eq_solution%constraints%is_constant_enthalpy()
        const_u = self%eq_solution%constraints%is_constant_energy()

        ! Associate subarray pointers
        G   => dG_dn
        A_g => self%eq_solver%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A_c => self%eq_solver%products%stoich_matrix(ng+1:,:)
        n = self%eq_solution%n
        nj  => self%eq_solution%nj
        nj_g => self%eq_solution%nj(:ng)
        ln_nj => self%eq_solution%ln_nj
        cp  => self%eq_solution%thermo%cp
        cv  => self%eq_solution%thermo%cv
        h_g => self%eq_solution%thermo%enthalpy(:ng)
        h_c => self%eq_solution%thermo%enthalpy(ng+1:)
        s_g => self%eq_solution%thermo%entropy(:ng)
        s_c => self%eq_solution%thermo%entropy(ng+1:)
        u_g => self%eq_solution%thermo%energy(:ng)
        u_c => self%eq_solution%thermo%energy(ng+1:)

        ! Get the mixture pressure
        P = self%eq_solution%calc_pressure()

        ! Compute gas phase chemical potentials
        mu_g = h_g - s_g + ln_nj + log(P/n)
        d_mu_g_dn = -1.d0/n

        d_hsu_delta_dn = 0.0d0
        if (const_s) then
            d_hsu_delta_dn = -sum(nj_g)/n
        end if

        d_n_delta_dn = 1.0d0

        ! Initialize the iteration matrix
        G = 0.0d0
        r = 0
        c = 0

        ! Set a flag if the entropy and temperature residuals are already converged
        rhs_converged = .false.

        !-------------------------------------------------------
        ! Equation (2.24/2.45): Element constraints
        !-------------------------------------------------------
        do i = 1,ne
            tmp = nj_g*A_g(:,i)
            r = r+1
            c = 0

            ! Pi derivatives
            do j = 1,ne
                c = c+1
            end do

            ! Condensed derivatives
            ! Handled in (2.25) below (symmetric)
            c = c+na

            ! Delta ln(n) derivative
            ! Symmetric with (2.26) pi derivative
            if (const_p) then
                c = c+1
            end if

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
            end if

            ! Right hand side
            if (.not. rhs_converged) then
                G(r,c+1) = d_mu_g_dn*sum(tmp)
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.25/2.46): Condensed phase constraints
        !-------------------------------------------------------
        do i = 1,nc
            if (.not. self%eq_solution%is_active(i)) cycle
            r = r+1
            c = 0

            ! Pi derivatives
            ! Symmetric with (2.24) condensed derivatives
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            ! Handled via zero initialization
            if (const_p) c = c+1

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.26)
        !-------------------------------------------------------
        if (const_p) then
            r = r+1
            c = 0

            ! Pi derivatives
            ! Handled in Eq (2.24) above.
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            c = c+1
            G(r,c) = -d_n_delta_dn

            ! Delta ln(T) derviative
            if (.not. const_t) then
                c = c+1
            end if

            ! Right-hand-side
            if (.not. rhs_converged) then
                G(r,c+1) = d_n_delta_dn + d_mu_g_dn*sum(nj_g)
            end if

        end if

        !---------------------------------------------------------
        ! Equation (2.27)/(2.28)/(2.47)/(2.48): Energy constraints
        !---------------------------------------------------------
        if (.not. const_t) then
            r = r+1
            c = 0

            ! Select entropy/enthalpy constraint
            if (const_s) then
                tmp = nj_g*(s_g-ln_nj-log(P/n))
                d_tmp_dn = nj_g/n
            else if (const_h) then
                tmp = nj_g*h_g
                d_tmp_dn = 0.0d0
            else if (const_u) then
                tmp = nj_g*u_g
                d_tmp_dn = 0.0d0
            end if

            ! Pi derivatives
            if (.not. rhs_converged) then
                do j = 1,ne
                    c = c+1
                    G(r,c) = dot_product(d_tmp_dn, A_g(:,j))
                end do
            end if

            ! Condensed derivatives
            if (.not. rhs_converged) then
                do j = 1,nc
                    if (.not. self%eq_solution%is_active(j)) cycle
                    c = c+1
                end do
            end if

            ! Delta ln(n) derivative
            if (.not. rhs_converged) then
                if (const_p) then
                    c = c+1
                    G(r,c) = sum(d_tmp_dn)
                end if
            end if

            ! Delta ln(T) derivative
            c = c+1
            if (const_p) then
                if (.not. rhs_converged) then
                    G(r,c) = dot_product(d_tmp_dn, h_g)
                else
                    c = c + ne + na + 1
                end if
            else
                G(r,c) = dot_product(d_tmp_dn, u_g)
            end if

            ! Right-hand-side
            if (.not. rhs_converged) then
                G(r,c+1) = d_hsu_delta_dn + dot_product(d_tmp_dn, mu_g) + d_mu_g_dn*sum(tmp)
                if (const_s) then
                    if (const_p) then
                        G(r,c+1) = G(r,c+1) + d_n_delta_dn
                    else
                        G(r,c+1) = G(r,c+1) - d_mu_g_dn*sum(nj_g)
                    end if
                end if
            end if

        end if

    end subroutine

    subroutine EqTotals_dG_dlnnj(self, dG_dlnnj, k)

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout), target :: dG_dlnnj(:,:)
        integer, intent(in) :: k  ! Index of ln(nj) to differentiate with respect to

        ! Locals
        integer  :: ng                          ! Number of gas species
        integer  :: nc                          ! Number of condensed species
        integer  :: na                          ! Number of active condensed species
        integer  :: ne                          ! Number of elements
        integer  :: num_eqn                     ! Active number of equations
        real(dp) :: tmp(self%eq_solver%num_gas) ! Common sub-expression storage
        real(dp) :: d_tmp_dlnnj                 !
        real(dp) :: mu_g(self%eq_solver%num_gas)  ! Gas phase chemical potentials [unitless]
        real(dp) :: d_mu_g_dlnnj                !
        real(dp) :: d_hsu_delta_dlnnj           !
        real(dp) :: n                           ! Total moles of mixture
        real(dp) :: P                           ! Pressure of mixture (bar)
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)           ! Log of gas species concentrations [kmol-per-kg]
        real(dp), pointer :: cp(:), cv(:)       ! Species heat capacities [unitless]
        real(dp), pointer :: h_g(:), h_c(:)     ! Gas/condensed enthalpies [unitless]
        real(dp), pointer :: s_g(:), s_c(:)     ! Gas/condensed entropies [unitless]
        real(dp), pointer :: u_g(:), u_c(:)     ! Gas/condensed energies [unitless]
        real(dp), pointer :: A_g(:,:), A_c(:,:) ! Gas/condensed stoichiometric matrices
        real(dp), pointer :: G(:,:)             ! Augmented Newton iteration matrix
        integer :: r, c                         ! Iteration matrix row/column indices
        integer :: i, j                         ! Loop counters
        logical :: const_p, const_t, const_s, const_h, const_u  ! Flags enabling/disabling matrix equations
        logical :: rhs_converged                ! Flag for RHS convergence

        ! Define shorthand
        ng = self%eq_solver%num_gas
        nc = self%eq_solver%num_condensed
        ne = self%eq_solver%num_elements
        na = count(self%eq_solution%is_active)
        num_eqn = self%num_eqn
        const_p = self%eq_solution%constraints%is_constant_pressure()
        const_t = self%eq_solution%constraints%is_constant_temperature()
        const_s = self%eq_solution%constraints%is_constant_entropy()
        const_h = self%eq_solution%constraints%is_constant_enthalpy()
        const_u = self%eq_solution%constraints%is_constant_energy()

        ! Associate subarray pointers
        G   => dG_dlnnj
        A_g => self%eq_solver%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A_c => self%eq_solver%products%stoich_matrix(ng+1:,:)
        n = self%eq_solution%n
        nj  => self%eq_solution%nj
        nj_g => self%eq_solution%nj(:ng)
        ln_nj => self%eq_solution%ln_nj
        cp  => self%eq_solution%thermo%cp
        cv  => self%eq_solution%thermo%cv
        h_g => self%eq_solution%thermo%enthalpy(:ng)
        h_c => self%eq_solution%thermo%enthalpy(ng+1:)
        s_g => self%eq_solution%thermo%entropy(:ng)
        s_c => self%eq_solution%thermo%entropy(ng+1:)
        u_g => self%eq_solution%thermo%energy(:ng)
        u_c => self%eq_solution%thermo%energy(ng+1:)

        ! Get the mixture pressure
        P = self%eq_solution%calc_pressure()

        ! Compute gas phase chemical potentials
        mu_g = h_g - s_g + ln_nj + log(P/n)
        d_mu_g_dlnnj = 1.0d0

        d_hsu_delta_dlnnj = 0.0d0
        if (const_s) then
            d_hsu_delta_dlnnj = 1.0d0 !(cons%state1 - soln%calc_entropy_sum(self))
        end if

        ! Initialize the iteration matrix
        G = 0.0d0
        r = 0
        c = 0

        ! Set a flag if the entropy and temperature residuals are already converged
        rhs_converged = .false.

        !-------------------------------------------------------
        ! Equation (2.24/2.45): Element constraints
        !-------------------------------------------------------
        do i = 1,ne
            tmp = nj_g*A_g(:,i)
            r = r+1
            c = 0

            ! Pi derivatives
            do j = 1,ne
                c = c+1
            end do

            ! Condensed derivatives
            ! Handled in (2.25) below (symmetric)
            c = c+na

            ! Delta ln(n) derivative
            ! Symmetric with (2.26) pi derivative
            if (const_p) then
                c = c+1
            end if

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
            end if

            ! Right hand side
            if (.not. rhs_converged) then
                G(r,c+1) = d_mu_g_dlnnj*tmp(k)
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.25/2.46): Condensed phase constraints
        !-------------------------------------------------------
        do i = 1,nc
            if (.not. self%eq_solution%is_active(i)) cycle
            r = r+1
            c = 0

            ! Pi derivatives
            ! Symmetric with (2.24) condensed derivatives
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            ! Handled via zero initialization
            if (const_p) c = c+1

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.26)
        !-------------------------------------------------------
        if (const_p) then
            r = r+1
            c = 0

            ! Pi derivatives
            ! Handled in Eq (2.24) above.
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            c = c+1

            ! Delta ln(T) derviative
            if (.not. const_t) then
                c = c+1
            end if

            ! Right-hand-side
            if (.not. rhs_converged) then
                G(r,c+1) = d_mu_g_dlnnj*nj_g(k)
            end if

        end if

        !---------------------------------------------------------
        ! Equation (2.27)/(2.28)/(2.47)/(2.48): Energy constraints
        !---------------------------------------------------------
        if (.not. const_t) then
            r = r+1
            c = 0

            ! Select entropy/enthalpy constraint
            if (const_s) then
                tmp = nj_g*(s_g-ln_nj-log(P/n))
                d_tmp_dlnnj = -nj_g(k)
            else if (const_h) then
                tmp = nj_g*h_g
                d_tmp_dlnnj = 0.0d0
            else if (const_u) then
                tmp = nj_g*u_g
                d_tmp_dlnnj = 0.0d0
            end if

            ! Pi derivatives
            if (.not. rhs_converged) then
                do j = 1,ne
                    c = c+1
                    G(r,c) = d_tmp_dlnnj*A_g(k,j)
                end do
            end if

            ! Condensed derivatives
            if (.not. rhs_converged) then
                do j = 1,nc
                    if (.not. self%eq_solution%is_active(j)) cycle
                    c = c+1
                end do
            end if

            ! Delta ln(n) derivative
            if (.not. rhs_converged) then
                if (const_p) then
                    c = c+1
                    G(r,c) = d_tmp_dlnnj
                end if
            end if

            ! Delta ln(T) derivative
            c = c+1
            if (const_p) then
                if (.not. rhs_converged) then
                    G(r,c) = d_tmp_dlnnj*h_g(k)
                else
                    c = c + ne + na + 1
                end if
            else
                G(r,c) = d_tmp_dlnnj*u_g(k)
            end if

            ! Right-hand-side
            if (.not. rhs_converged) then
                G(r,c+1) = d_hsu_delta_dlnnj + d_tmp_dlnnj*mu_g(k) + d_mu_g_dlnnj*tmp(k)
                if (const_s) then
                    if (.not. const_p) then
                        G(r,c+1) = G(r,c+1) - d_mu_g_dlnnj*sum(nj_g)
                    end if
                end if
            end if

        end if

    end subroutine

    subroutine EqTotals_dG_dHj(self, dG_dHj, k)

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout), target :: dG_dHj(:,:)
        integer, intent(in) :: k  ! Index of Hj to differentiate with respect to

        ! Locals
        integer  :: ng                          ! Number of gas species
        integer  :: nc                          ! Number of condensed species
        integer  :: na                          ! Number of active condensed species
        integer  :: ne                          ! Number of elements
        integer  :: num_eqn                     ! Active number of equations
        real(dp) :: tmp(self%eq_solver%num_gas) ! Common sub-expression storage
        real(dp) :: d_tmp_dhj                 !
        real(dp) :: mu_g(self%eq_solver%num_gas)  ! Gas phase chemical potentials [unitless]
        real(dp) :: d_mu_g_dhj               !
        real(dp) :: d_hsu_delta_dhj          !
        real(dp) :: n                           ! Total moles of mixture
        real(dp) :: P                           ! Pressure of mixture (bar)
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)           ! Log of gas species concentrations [kmol-per-kg]
        real(dp), pointer :: cp(:), cv(:)       ! Species heat capacities [unitless]
        real(dp), pointer :: h_g(:), h_c(:)     ! Gas/condensed enthalpies [unitless]
        real(dp), pointer :: s_g(:), s_c(:)     ! Gas/condensed entropies [unitless]
        real(dp), pointer :: u_g(:), u_c(:)     ! Gas/condensed energies [unitless]
        real(dp), pointer :: A_g(:,:), A_c(:,:) ! Gas/condensed stoichiometric matrices
        real(dp), pointer :: G(:,:)             ! Augmented Newton iteration matrix
        integer :: r, c                         ! Iteration matrix row/column indices
        integer :: i, j                         ! Loop counters
        logical :: const_p, const_t, const_s, const_h, const_u  ! Flags enabling/disabling matrix equations
        logical :: rhs_converged                ! Flag for RHS convergence

        ! Define shorthand
        ng = self%eq_solver%num_gas
        nc = self%eq_solver%num_condensed
        ne = self%eq_solver%num_elements
        na = count(self%eq_solution%is_active)
        num_eqn = self%num_eqn
        const_p = self%eq_solution%constraints%is_constant_pressure()
        const_t = self%eq_solution%constraints%is_constant_temperature()
        const_s = self%eq_solution%constraints%is_constant_entropy()
        const_h = self%eq_solution%constraints%is_constant_enthalpy()
        const_u = self%eq_solution%constraints%is_constant_energy()

        ! Associate subarray pointers
        G   => dG_dHj
        A_g => self%eq_solver%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A_c => self%eq_solver%products%stoich_matrix(ng+1:,:)
        n = self%eq_solution%n
        nj  => self%eq_solution%nj
        nj_g => self%eq_solution%nj(:ng)
        ln_nj => self%eq_solution%ln_nj
        cp  => self%eq_solution%thermo%cp
        cv  => self%eq_solution%thermo%cv
        h_g => self%eq_solution%thermo%enthalpy(:ng)
        h_c => self%eq_solution%thermo%enthalpy(ng+1:)
        s_g => self%eq_solution%thermo%entropy(:ng)
        s_c => self%eq_solution%thermo%entropy(ng+1:)
        u_g => self%eq_solution%thermo%energy(:ng)
        u_c => self%eq_solution%thermo%energy(ng+1:)

        ! Get the mixture pressure
        P = self%eq_solution%calc_pressure()

        ! Compute gas phase chemical potentials
        mu_g = h_g - s_g + ln_nj + log(P/n)
        d_mu_g_dhj = 1.0d0

        if (const_s) then
            d_hsu_delta_dhj = -1.0d0
        else if (const_h) then
            d_hsu_delta_dhj = -nj(k)
        else if (const_u) then
            d_hsu_delta_dhj = 0.0d0
        end if

        ! Initialize the iteration matrix
        G = 0.0d0
        r = 0
        c = 0

        ! Set a flag if the entropy and temperature residuals are already converged
        rhs_converged = .false.
        !if (const_s .and. soln%entropy_converged .and. soln%temperature_converged) rhs_converged = .true.

        !-------------------------------------------------------
        ! Equation (2.24/2.45): Element constraints
        !-------------------------------------------------------
        do i = 1,ne
            tmp = nj_g*A_g(:,i)
            r = r+1
            c = 0

            ! Pi derivatives
            do j = 1,ne
                c = c+1
            end do

            ! Condensed derivatives
            ! Handled in (2.25) below (symmetric)
            c = c+na

            ! Delta ln(n) derivative
            ! Symmetric with (2.26) pi derivative
            if (const_p) then
                c = c+1
            end if

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
                if (const_p) then
                    G(r,c) = tmp(k)
                end if
            end if

            ! Right hand side
            if (.not. rhs_converged) then
                G(r,c+1) = d_mu_g_dhj*tmp(k)
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.25/2.46): Condensed phase constraints
        !-------------------------------------------------------
        do i = 1,nc
            if (.not. self%eq_solution%is_active(i)) cycle
            r = r+1
            c = 0

            ! Pi derivatives
            ! Symmetric with (2.24) condensed derivatives
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            ! Handled via zero initialization
            if (const_p) c = c+1

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.26)
        !-------------------------------------------------------
        if (const_p) then
            r = r+1
            c = 0

            ! Pi derivatives
            ! Handled in Eq (2.24) above.
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            c = c+1

            ! Delta ln(T) derviative
            if (.not. const_t) then
                c = c+1
                G(r,c) = nj_g(k)
            end if

            ! Right-hand-side
            if (.not. rhs_converged) then
                G(r,c+1) = nj_g(k)*d_mu_g_dhj
            end if

        end if

        !---------------------------------------------------------
        ! Equation (2.27)/(2.28)/(2.47)/(2.48): Energy constraints
        !---------------------------------------------------------
        if (.not. const_t) then
            r = r+1
            c = 0

            ! Select entropy/enthalpy constraint
            if (const_s) then
                tmp = nj_g*(s_g-ln_nj-log(P/n))
                d_tmp_dhj = 0.0d0
            else if (const_h) then
                tmp = nj_g*h_g
                d_tmp_dhj = nj_g(k)
            else if (const_u) then
                tmp = nj_g*u_g
                d_tmp_dhj = 0.0d0
            end if

            ! Pi derivatives
            if (.not. rhs_converged) then
                do j = 1,ne
                    c = c+1
                    G(r,c) = d_tmp_dhj*A_g(k,j)
                end do
            end if

            ! Condensed derivatives
            if (.not. rhs_converged) then
                do j = 1,nc
                    if (.not. self%eq_solution%is_active(j)) cycle
                    c = c+1
                end do
            end if

            ! Delta ln(n) derivative
            if (.not. rhs_converged) then
                if (const_p) then
                    c = c+1
                    G(r,c) = d_tmp_dhj
                end if
            end if

            ! Delta ln(T) derivative
            c = c+1
            if (const_p) then
                if (.not. rhs_converged) then
                    G(r,c) = d_tmp_dhj*h_g(k) + tmp(k)
                else
                    c = c + ne + na + 1
                    G(r,c) = 2.0d0*nj_g(k)*h_g(k)
                end if
            else
                G(r,c) = d_tmp_dhj*u_g(k)
            end if

            ! Right-hand-side
            if (.not. rhs_converged) then
                G(r,c+1) = d_hsu_delta_dhj + tmp(k)*d_mu_g_dhj + d_tmp_dhj*mu_g(k)
                if (const_s) then
                    if (.not. const_p) then
                        G(r,c+1) = G(r,c+1) - nj_g(k)*d_mu_g_dhj
                    end if
                end if
            end if

        end if

    end subroutine

    subroutine EqTotals_dG_dSj(self, dG_dSj, k)

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout), target :: dG_dSj(:,:)
        integer, intent(in) :: k  ! Index of Hj to differentiate with respect to

        ! Locals
        integer  :: ng                          ! Number of gas species
        integer  :: nc                          ! Number of condensed species
        integer  :: na                          ! Number of active condensed species
        integer  :: ne                          ! Number of elements
        integer  :: num_eqn                     ! Active number of equations
        real(dp) :: tmp(self%eq_solver%num_gas) ! Common sub-expression storage
        real(dp) :: d_tmp_dsj                 !
        real(dp) :: mu_g(self%eq_solver%num_gas)  ! Gas phase chemical potentials [unitless]
        real(dp) :: d_mu_g_dsj               !
        real(dp) :: d_hsu_delta_dsj          !
        real(dp) :: n                           ! Total moles of mixture
        real(dp) :: P                           ! Pressure of mixture (bar)
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)           ! Log of gas species concentrations [kmol-per-kg]
        real(dp), pointer :: cp(:), cv(:)       ! Species heat capacities [unitless]
        real(dp), pointer :: h_g(:), h_c(:)     ! Gas/condensed enthalpies [unitless]
        real(dp), pointer :: s_g(:), s_c(:)     ! Gas/condensed entropies [unitless]
        real(dp), pointer :: u_g(:), u_c(:)     ! Gas/condensed energies [unitless]
        real(dp), pointer :: A_g(:,:), A_c(:,:) ! Gas/condensed stoichiometric matrices
        real(dp), pointer :: G(:,:)             ! Augmented Newton iteration matrix
        integer :: r, c                         ! Iteration matrix row/column indices
        integer :: i, j                         ! Loop counters
        logical :: const_p, const_t, const_s, const_h, const_u  ! Flags enabling/disabling matrix equations
        logical :: rhs_converged                ! Flag for RHS convergence

        ! Define shorthand
        ng = self%eq_solver%num_gas
        nc = self%eq_solver%num_condensed
        ne = self%eq_solver%num_elements
        na = count(self%eq_solution%is_active)
        num_eqn = self%num_eqn
        const_p = self%eq_solution%constraints%is_constant_pressure()
        const_t = self%eq_solution%constraints%is_constant_temperature()
        const_s = self%eq_solution%constraints%is_constant_entropy()
        const_h = self%eq_solution%constraints%is_constant_enthalpy()
        const_u = self%eq_solution%constraints%is_constant_energy()

        ! Associate subarray pointers
        G   => dG_dSj
        A_g => self%eq_solver%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A_c => self%eq_solver%products%stoich_matrix(ng+1:,:)
        n = self%eq_solution%n
        nj  => self%eq_solution%nj
        nj_g => self%eq_solution%nj(:ng)
        ln_nj => self%eq_solution%ln_nj
        cp  => self%eq_solution%thermo%cp
        cv  => self%eq_solution%thermo%cv
        h_g => self%eq_solution%thermo%enthalpy(:ng)
        h_c => self%eq_solution%thermo%enthalpy(ng+1:)
        s_g => self%eq_solution%thermo%entropy(:ng)
        s_c => self%eq_solution%thermo%entropy(ng+1:)
        u_g => self%eq_solution%thermo%energy(:ng)
        u_c => self%eq_solution%thermo%energy(ng+1:)

        ! Get the mixture pressure
        P = self%eq_solution%calc_pressure()

        ! Compute gas phase chemical potentials
        mu_g = h_g - s_g + ln_nj + log(P/n)
        d_mu_g_dsj = -1.0d0

        d_hsu_delta_dsj = 0.0d0
        if (const_s) then
            d_hsu_delta_dsj = -nj(k)
        end if

        ! Initialize the iteration matrix
        G = 0.0d0
        r = 0
        c = 0

        ! Set a flag if the entropy and temperature residuals are already converged
        rhs_converged = .false.

        !-------------------------------------------------------
        ! Equation (2.24/2.45): Element constraints
        !-------------------------------------------------------
        do i = 1,ne
            tmp = nj_g*A_g(:,i)
            r = r+1
            c = 0

            ! Pi derivatives
            do j = 1,ne
                c = c+1
            end do

            ! Condensed derivatives
            ! Handled in (2.25) below (symmetric)
            c = c+na

            ! Delta ln(n) derivative
            ! Symmetric with (2.26) pi derivative
            if (const_p) then
                c = c+1
            end if

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
            end if

            ! Right hand side
            if (.not. rhs_converged) then
                G(r,c+1) = tmp(k)*d_mu_g_dsj
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.25/2.46): Condensed phase constraints
        !-------------------------------------------------------
        do i = 1,nc
            if (.not. self%eq_solution%is_active(i)) cycle
            r = r+1
            c = 0

            ! Pi derivatives
            ! Symmetric with (2.24) condensed derivatives
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            ! Handled via zero initialization
            if (const_p) c = c+1

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.26)
        !-------------------------------------------------------
        if (const_p) then
            r = r+1
            c = 0

            ! Pi derivatives
            ! Handled in Eq (2.24) above.
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            c = c+1

            ! Delta ln(T) derviative
            if (.not. const_t) then
                c = c+1
            end if

            ! Right-hand-side
            if (.not. rhs_converged) then
                G(r,c+1) = nj_g(k)*d_mu_g_dsj
            end if

        end if

        !---------------------------------------------------------
        ! Equation (2.27)/(2.28)/(2.47)/(2.48): Energy constraints
        !---------------------------------------------------------
        if (.not. const_t) then
            r = r+1
            c = 0

            ! Select entropy/enthalpy constraint
            if (const_s) then
                tmp = nj_g*(s_g-ln_nj-log(P/n))
                d_tmp_dsj = nj_g(k)
            else if (const_h) then
                tmp = nj_g*h_g
                d_tmp_dsj = 0.0d0
            else if (const_u) then
                tmp = nj_g*u_g
                d_tmp_dsj = 0.0d0
            end if

            ! Pi derivatives
            if (.not. rhs_converged) then
                do j = 1,ne
                    c = c+1
                    G(r,c) = d_tmp_dsj*A_g(k,j)
                end do
            end if

            ! Condensed derivatives
            if (.not. rhs_converged) then
                do j = 1,nc
                    if (.not. self%eq_solution%is_active(j)) cycle
                    c = c+1
                end do
            end if

            ! Delta ln(n) derivative
            if (.not. rhs_converged) then
                if (const_p) then
                    c = c+1
                    G(r,c) = d_tmp_dsj
                end if
            end if

            ! Delta ln(T) derivative
            c = c+1
            if (const_p) then
                if (.not. rhs_converged) then
                    G(r,c) = d_tmp_dsj*h_g(k)
                else
                    c = c + ne + na + 1
                end if
            else
                G(r,c) = d_tmp_dsj*u_g(k)
            end if

            ! Right-hand-side
            if (.not. rhs_converged) then
                G(r,c+1) = d_hsu_delta_dsj + d_tmp_dsj*mu_g(k) + tmp(k)*d_mu_g_dsj
                if (const_s) then
                    if (.not. const_p) then
                        G(r,c+1) = G(r,c+1) - nj_g(k)*d_mu_g_dsj
                    end if
                end if
            end if

        end if

    end subroutine

    subroutine EqTotals_dG_dUj(self, dG_dUj, k)

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout), target :: dG_dUj(:,:)
        integer, intent(in) :: k  ! Index of Uj to differentiate with respect to

        ! Locals
        integer  :: ng                          ! Number of gas species
        integer  :: nc                          ! Number of condensed species
        integer  :: na                          ! Number of active condensed species
        integer  :: ne                          ! Number of elements
        integer  :: num_eqn                     ! Active number of equations
        real(dp) :: tmp(self%eq_solver%num_gas) ! Common sub-expression storage
        real(dp) :: d_tmp_duj                 !
        real(dp) :: mu_g(self%eq_solver%num_gas)  ! Gas phase chemical potentials [unitless]
        real(dp) :: d_mu_g_duj               !
        real(dp) :: d_hsu_delta_duj          !
        real(dp) :: n                           ! Total moles of mixture
        real(dp) :: P                           ! Pressure of mixture (bar)
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)           ! Log of gas species concentrations [kmol-per-kg]
        real(dp), pointer :: cp(:), cv(:)       ! Species heat capacities [unitless]
        real(dp), pointer :: h_g(:), h_c(:)     ! Gas/condensed enthalpies [unitless]
        real(dp), pointer :: s_g(:), s_c(:)     ! Gas/condensed entropies [unitless]
        real(dp), pointer :: u_g(:), u_c(:)     ! Gas/condensed energies [unitless]
        real(dp), pointer :: A_g(:,:), A_c(:,:) ! Gas/condensed stoichiometric matrices
        real(dp), pointer :: G(:,:)             ! Augmented Newton iteration matrix
        integer :: r, c                         ! Iteration matrix row/column indices
        integer :: i, j                         ! Loop counters
        logical :: const_p, const_t, const_s, const_h, const_u  ! Flags enabling/disabling matrix equations
        logical :: rhs_converged                ! Flag for RHS convergence

        ! Define shorthand
        ng = self%eq_solver%num_gas
        nc = self%eq_solver%num_condensed
        ne = self%eq_solver%num_elements
        na = count(self%eq_solution%is_active)
        num_eqn = self%num_eqn
        const_p = self%eq_solution%constraints%is_constant_pressure()
        const_t = self%eq_solution%constraints%is_constant_temperature()
        const_s = self%eq_solution%constraints%is_constant_entropy()
        const_h = self%eq_solution%constraints%is_constant_enthalpy()
        const_u = self%eq_solution%constraints%is_constant_energy()

        ! Associate subarray pointers
        G   => dG_dUj
        A_g => self%eq_solver%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A_c => self%eq_solver%products%stoich_matrix(ng+1:,:)
        n = self%eq_solution%n
        nj  => self%eq_solution%nj
        nj_g => self%eq_solution%nj(:ng)
        ln_nj => self%eq_solution%ln_nj
        cp  => self%eq_solution%thermo%cp
        cv  => self%eq_solution%thermo%cv
        h_g => self%eq_solution%thermo%enthalpy(:ng)
        h_c => self%eq_solution%thermo%enthalpy(ng+1:)
        s_g => self%eq_solution%thermo%entropy(:ng)
        s_c => self%eq_solution%thermo%entropy(ng+1:)
        u_g => self%eq_solution%thermo%energy(:ng)
        u_c => self%eq_solution%thermo%energy(ng+1:)

        ! Get the mixture pressure
        P = self%eq_solution%calc_pressure()

        ! Compute gas phase chemical potentials
        mu_g = h_g - s_g + ln_nj + log(P/n)

        if (const_u) then
            d_hsu_delta_duj = -nj(k)
        else
            d_hsu_delta_duj = 0.0d0
        end if

        ! Initialize the iteration matrix
        G = 0.0d0
        r = 0
        c = 0

        ! Set a flag if the entropy and temperature residuals are already converged
        rhs_converged = .false.

        !-------------------------------------------------------
        ! Equation (2.24/2.45): Element constraints
        !-------------------------------------------------------
        do i = 1,ne
            tmp = nj_g*A_g(:,i)
            r = r+1
            c = 0

            ! Pi derivatives
            do j = 1,ne
                c = c+1
            end do

            ! Condensed derivatives
            ! Handled in (2.25) below (symmetric)
            c = c+na

            ! Delta ln(n) derivative
            ! Symmetric with (2.26) pi derivative
            if (const_p) then
                c = c+1
            end if

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
                if (.not. const_p) then
                    G(r,c) = tmp(k)
                end if
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.25/2.46): Condensed phase constraints
        !-------------------------------------------------------
        do i = 1,nc
            if (.not. self%eq_solution%is_active(i)) cycle
            r = r+1
            c = 0

            ! Pi derivatives
            ! Symmetric with (2.24) condensed derivatives
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            ! Handled via zero initialization
            if (const_p) c = c+1

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
            end if

        end do

        !-------------------------------------------------------
        ! Equation (2.26)
        !-------------------------------------------------------
        if (const_p) then
            r = r+1
            c = 0

            ! Pi derivatives
            ! Handled in Eq (2.24) above.
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            c = c+1

            ! Delta ln(T) derviative
            if (.not. const_t) then
                c = c+1
            end if

        end if

        !---------------------------------------------------------
        ! Equation (2.27)/(2.28)/(2.47)/(2.48): Energy constraints
        !---------------------------------------------------------
        if (.not. const_t) then
            r = r+1
            c = 0

            ! Select entropy/enthalpy constraint
            if (const_s) then
                tmp = nj_g*(s_g-ln_nj-log(P/n))
                d_tmp_duj = 0.0d0
            else if (const_h) then
                tmp = nj_g*h_g
                d_tmp_duj = 0.0d0
            else if (const_u) then
                tmp = nj_g*u_g
                d_tmp_duj = nj_g(k)
            end if

            ! Pi derivatives
            if (.not. rhs_converged) then
                do j = 1,ne
                    c = c+1
                    G(r,c) = d_tmp_duj*A_g(k,j)
                end do
            end if

            ! Condensed derivatives
            if (.not. rhs_converged) then
                do j = 1,nc
                    if (.not. self%eq_solution%is_active(j)) cycle
                    c = c+1
                end do
            end if

            ! Delta ln(n) derivative
            if (.not. rhs_converged) then
                if (const_p) then
                    c = c+1
                    G(r,c) = d_tmp_duj
                end if
            end if

            ! Delta ln(T) derivative
            c = c+1
            if (.not. const_p) then
                G(r,c) = tmp(k) + d_tmp_duj*u_g(k)
                if (const_s) then
                    G(r,c) = G(r,c) - nj_g(k)
                end if
            end if

            ! Right-hand-side
            if (.not. rhs_converged) then
                G(r,c+1) = d_hsu_delta_duj + d_tmp_duj
            end if

        end if

    end subroutine

    ! -----------------------------------------------------------------------
    ! Finite difference routines
    ! -----------------------------------------------------------------------

    subroutine EqTotals_finite_difference_totals(self, reactant_weights, h)
        ! Compute the total derivatives of the equilibrium solution using finite differences
        ! Computes derivatives of: nj, n, and T with respect to state1, state2, and reactant weights
        ! This uses a forward finite difference scheme with a step size of h

        ! TODO: implement backward and central difference options
        ! TODO: get rid of reactant_weights as an argument here

        ! Arguments
        class(EqTotals) :: self
        real(dp), intent(inout) :: reactant_weights(:)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_   ! Step size
        logical :: const_t, const_p ! Flag that is true if problem is constant temperature
        integer :: i, j, k  ! Loop index

        const_t = self%eq_solution%constraints%is_constant_temperature()
        const_p = self%eq_solution%constraints%is_constant_pressure()

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ! *** Assume the problem has been solved once and the solution is stored in self%eq_solution ***

        ! Perturb the variables and solve the problem again
        pert_soln = self%eq_solution
        pert_soln%constraints%state1 = self%eq_solution%constraints%state1 + h_
        call self%eq_solver%solve(pert_soln, self%eq_solution%constraints%type, &
                                  self%eq_solution%constraints%state1, &
                                  self%eq_solution%constraints%state2, &
                                  reactant_weights)
        self%dnj_dstate1 = (pert_soln%nj - self%eq_solution%nj) / h_
        if (const_p) self%dn_dstate1 = (pert_soln%n - self%eq_solution%n) / h_
        if (.not. const_t) then
            self%dT_dstate1 = (pert_soln%T - self%eq_solution%T) / h_
        else
            self%dT_dstate1 = 1.0d0   ! Derivative of T wrt T is 1
        end if
        pert_soln%constraints%state1 = self%eq_solution%constraints%state1 - h_

        pert_soln = self%eq_solution
        pert_soln%constraints%state1 = self%eq_solution%constraints%state2 + h_
        call self%eq_solver%solve(pert_soln, self%eq_solution%constraints%type, &
                                  self%eq_solution%constraints%state1, &
                                  self%eq_solution%constraints%state2, &
                                  reactant_weights)
        self%dnj_dstate2 = (pert_soln%nj - self%eq_solution%nj) / h_
        if (const_p) self%dn_dstate2 = (pert_soln%n - self%eq_solution%n) / h_
        if (.not. const_t) self%dT_dstate2 = (pert_soln%T - self%eq_solution%T) / h_
        pert_soln%constraints%state1 = self%eq_solution%constraints%state2 - h_

        pert_soln = self%eq_solution
        do i = 1, self%eq_solver%num_reactants
            reactant_weights(i) = reactant_weights(i) + h_
            call self%eq_solver%solve(pert_soln, self%eq_solution%constraints%type, &
                                    self%eq_solution%constraints%state1, &
                                    self%eq_solution%constraints%state2, &
                                    reactant_weights)
            self%dnj_dweights(:, i) = (pert_soln%nj - self%eq_solution%nj) / h_
            if (const_p) self%dn_dweights(i) = (pert_soln%n - self%eq_solution%n) / h_
            if (.not. const_t) self%dT_dweights(i) = (pert_soln%T - self%eq_solution%T) / h_
            reactant_weights(i) = reactant_weights(i) - h_

        end do

    end subroutine

    subroutine EqTotals_fd_dR_eq_dH0(self, dR_eq_dH0, h)
        ! Finite difference the partial derivative matrix of R_eq wrt state1 (H0, S0, T0, or U0)

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_eq_dH0(:)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        logical :: const_h
        real(dp) :: h_  ! Step size
        real(dp) :: R_eq(self%num_eqn), R_eq_pert(self%num_eqn)  ! R_eq residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        const_h = self%eq_solution%constraints%is_constant_enthalpy()

        if (const_h) then
            call self%eq_solver%assemble_matrix(self%eq_solution)
            R_eq = -self%eq_solution%G(:, self%num_eqn+1)
            pert_soln = self%eq_solution
            pert_soln%constraints%state1 = pert_soln%constraints%state1 + h_
            call self%eq_solver%assemble_matrix(pert_soln)
            R_eq_pert = -pert_soln%G(:, self%num_eqn+1)

            dR_eq_dH0 = (R_eq_pert - R_eq) / h_
        else
            dR_eq_dH0 = 0.0d0
        end if

    end subroutine

    subroutine EqTotals_fd_dR_eq_db0(self, dR_eq_db0, h)
        ! Finite difference the partial derivative matrix of R_eq wrt b0

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_eq_db0(:, :)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_eq(self%num_eqn), R_eq_pert(self%num_eqn)  ! R_eq residuals
        integer :: i

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        do i = 1, self%eq_solver%num_elements
            call self%eq_solver%assemble_matrix(self%eq_solution)
            R_eq = -self%eq_solution%G(:, self%num_eqn+1)
            pert_soln = self%eq_solution
            pert_soln%constraints%b0(i) = pert_soln%constraints%b0(i) + h_
            call self%eq_solver%assemble_matrix(pert_soln)
            R_eq_pert = -pert_soln%G(:, self%num_eqn+1)

            dR_eq_db0(:, i) = (R_eq_pert - R_eq) / h_

            pert_soln%constraints%b0(i) = pert_soln%constraints%b0(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_eq_dP(self, dR_eq_dP, h)
        ! Finite difference the partial derivative matrix of R_eq wrt state1 (H0, S0, T0, or U0)

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_eq_dP(:)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        logical :: const_p
        real(dp) :: h_  ! Step size
        real(dp) :: R_eq(self%num_eqn), R_eq_pert(self%num_eqn)  ! R_eq residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        const_p = self%eq_solution%constraints%is_constant_pressure()

        if (const_p) then
            call self%eq_solver%assemble_matrix(self%eq_solution)
            R_eq = matmul(self%eq_solution%G(:, :self%num_eqn), self%u_hat) - self%eq_solution%G(:, self%num_eqn+1)
            pert_soln = self%eq_solution
            pert_soln%constraints%state2 = pert_soln%constraints%state2 + h_
            call self%eq_solver%assemble_matrix(pert_soln)
            R_eq_pert = matmul(pert_soln%G(:, :self%num_eqn), self%u_hat) - pert_soln%G(:, self%num_eqn+1)

            dR_eq_dP = (R_eq_pert - R_eq) / h_
        else
            ! TODO: implement this for volume problems
            dR_eq_dP = 0.0d0
        end if

    end subroutine

    subroutine EqTotals_fd_dR_eq_dnj(self, dR_eq_dnj, h)
        ! Finite difference the partial derivative matrix of R_eq wrt nj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_eq_dnj(:, :)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_eq(self%num_eqn), R_eq_pert(self%num_eqn)  ! R_eq residuals
        integer :: i, j, k

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        do i = 1, self%eq_solver%num_gas
            call self%eq_solver%assemble_matrix(self%eq_solution)
            R_eq = matmul(self%eq_solution%G(:, :self%num_eqn), self%u_hat) - self%eq_solution%G(:, self%num_eqn+1)

            pert_soln = self%eq_solution
            pert_soln%nj(i) = pert_soln%nj(i) + h_
            call self%eq_solver%assemble_matrix(pert_soln)
            R_eq_pert = matmul(pert_soln%G(:, :self%num_eqn), self%u_hat) - pert_soln%G(:, self%num_eqn+1)

            dR_eq_dnj(:, i) = (R_eq_pert - R_eq) / h_

            pert_soln%nj(i) = pert_soln%nj(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_eq_dn(self, dR_eq_dn, h)
        ! Finite difference the partial derivative matrix of R_eq wrt n

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_eq_dn(:)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_eq(self%num_eqn), R_eq_pert(self%num_eqn)  ! R_eq residuals
        integer :: i

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        call self%eq_solver%assemble_matrix(self%eq_solution)
        R_eq = matmul(self%eq_solution%G(:, :self%num_eqn), self%u_hat) - self%eq_solution%G(:, self%num_eqn+1)

        pert_soln = self%eq_solution
        pert_soln%n = pert_soln%n + h_
        call self%eq_solver%assemble_matrix(pert_soln)
        R_eq_pert = matmul(pert_soln%G(:, :self%num_eqn), self%u_hat) - pert_soln%G(:, self%num_eqn+1)

        dR_eq_dn = (R_eq_pert - R_eq) / h_

    end subroutine

    subroutine EqTotals_fd_dR_eq_dT(self, dR_eq_dT, h)
        ! Finite difference the partial derivative matrix of R_eq wrt n

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_eq_dT(:)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_eq(self%num_eqn), R_eq_pert(self%num_eqn)  ! R_eq residuals
        integer :: i

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        call self%eq_solver%assemble_matrix(self%eq_solution)
        R_eq = matmul(self%eq_solution%G(:, :self%num_eqn), self%u_hat) - self%eq_solution%G(:, self%num_eqn+1)

        pert_soln = self%eq_solution
        pert_soln%T = pert_soln%T + h_
        call self%eq_solver%assemble_matrix(pert_soln)
        R_eq_pert = matmul(pert_soln%G(:, :self%num_eqn), self%u_hat) - pert_soln%G(:, self%num_eqn+1)

        dR_eq_dT = (R_eq_pert - R_eq) / h_

    end subroutine

    subroutine EqTotals_fd_dR_eq_db(self, dR_eq_db, h)
        ! Finite difference the partial derivative matrix of R_eq wrt b

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_eq_db(:, :)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_eq(self%num_eqn), R_eq_pert(self%num_eqn)  ! R_eq residuals
        integer :: i

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        do i = 1, self%eq_solver%num_elements
            call self%eq_solver%assemble_matrix(self%eq_solution)
            R_eq = -self%eq_solution%G(:, self%num_eqn+1)
            pert_soln = self%eq_solution
            pert_soln%b(i) = pert_soln%b(i) + h_
            call self%eq_solver%assemble_matrix(pert_soln)
            R_eq_pert = -pert_soln%G(:, self%num_eqn+1)

            dR_eq_db(:, i) = (R_eq_pert - R_eq) / h_

            pert_soln%b(i) = pert_soln%b(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_eq_dlnnj(self, dR_eq_dlnnj, h)
        ! Finite difference the partial derivative matrix of R_eq wrt nj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_eq_dlnnj(:, :)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_eq(self%num_eqn), R_eq_pert(self%num_eqn)  ! R_eq residuals
        integer :: i

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        do i = 1, self%eq_solver%num_gas
            call self%eq_solver%assemble_matrix(self%eq_solution)
            R_eq = matmul(self%eq_solution%G(:, :self%num_eqn), self%u_hat) - self%eq_solution%G(:, self%num_eqn+1)

            pert_soln = self%eq_solution
            pert_soln%ln_nj(i) = pert_soln%ln_nj(i) + h_
            call self%eq_solver%assemble_matrix(pert_soln)
            R_eq_pert = matmul(pert_soln%G(:, :self%num_eqn), self%u_hat) - pert_soln%G(:, self%num_eqn+1)

            dR_eq_dlnnj(:, i) = (R_eq_pert - R_eq) / h_

            pert_soln%ln_nj(i) = pert_soln%ln_nj(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_eq_dHj(self, dR_eq_dHj, h)
        ! Finite difference the partial derivative matrix of R_eq wrt Hj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_eq_dHj(:, :)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_eq(self%num_eqn), R_eq_pert(self%num_eqn)  ! R_eq residuals
        integer :: i, j, k

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        do i = 1, self%eq_solver%num_gas
            call self%eq_solver%assemble_matrix(self%eq_solution)
            R_eq = matmul(self%eq_solution%G(:, :self%num_eqn), self%u_hat) - self%eq_solution%G(:, self%num_eqn+1)

            pert_soln = self%eq_solution
            pert_soln%thermo%enthalpy(i) = pert_soln%thermo%enthalpy(i) + h_
            call self%eq_solver%assemble_matrix(pert_soln)
            R_eq_pert = matmul(pert_soln%G(:, :self%num_eqn), self%u_hat) - pert_soln%G(:, self%num_eqn+1)

            dR_eq_dHj(:, i) = (R_eq_pert - R_eq) / h_

            pert_soln%thermo%enthalpy(i) = pert_soln%thermo%enthalpy(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_eq_dSj(self, dR_eq_dSj, h)
        ! Finite difference the partial derivative matrix of R_eq wrt Sj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_eq_dSj(:, :)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_eq(self%num_eqn), R_eq_pert(self%num_eqn)  ! R_eq residuals
        integer :: i, j, k

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        do i = 1, self%eq_solver%num_gas
            call self%eq_solver%assemble_matrix(self%eq_solution)
            R_eq = matmul(self%eq_solution%G(:, :self%num_eqn), self%u_hat) - self%eq_solution%G(:, self%num_eqn+1)

            pert_soln = self%eq_solution
            pert_soln%thermo%entropy(i) = pert_soln%thermo%entropy(i) + h_
            call self%eq_solver%assemble_matrix(pert_soln)
            R_eq_pert = matmul(pert_soln%G(:, :self%num_eqn), self%u_hat) - pert_soln%G(:, self%num_eqn+1)

            dR_eq_dSj(:, i) = (R_eq_pert - R_eq) / h_

            pert_soln%thermo%entropy(i) = pert_soln%thermo%entropy(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_eq_dUj(self, dR_eq_dUj, h)
        ! Finite difference the partial derivative matrix of R_eq wrt Uj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_eq_dUj(:, :)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_eq(self%num_eqn), R_eq_pert(self%num_eqn)  ! R_eq residuals
        integer :: i, j, k

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        do i = 1, self%eq_solver%num_gas
            call self%eq_solver%assemble_matrix(self%eq_solution)
            R_eq = matmul(self%eq_solution%G(:, :self%num_eqn), self%u_hat) - self%eq_solution%G(:, self%num_eqn+1)

            pert_soln = self%eq_solution
            pert_soln%thermo%energy(i) = pert_soln%thermo%energy(i) + h_
            call self%eq_solver%assemble_matrix(pert_soln)
            R_eq_pert = matmul(pert_soln%G(:, :self%num_eqn), self%u_hat) - pert_soln%G(:, self%num_eqn+1)

            dR_eq_dUj(:, i) = (R_eq_pert - R_eq) / h_

            pert_soln%thermo%energy(i) = pert_soln%thermo%energy(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_eq_dCpj(self, dR_eq_dCpj, h)
        ! Finite difference the partial derivative matrix of R_eq wrt Sj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_eq_dCpj(:, :)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_eq(self%num_eqn), R_eq_pert(self%num_eqn)  ! R_eq residuals
        integer :: i, j, k

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        do i = 1, self%eq_solver%num_gas
            call self%eq_solver%assemble_matrix(self%eq_solution)
            R_eq = matmul(self%eq_solution%G(:, :self%num_eqn), self%u_hat) - self%eq_solution%G(:, self%num_eqn+1)

            pert_soln = self%eq_solution
            pert_soln%thermo%cp(i) = pert_soln%thermo%cp(i) + h_
            call self%eq_solver%assemble_matrix(pert_soln)
            R_eq_pert = matmul(pert_soln%G(:, :self%num_eqn), self%u_hat) - pert_soln%G(:, self%num_eqn+1)

            ! write(*,*) "R_eq = ", R_eq
            ! write(*,*) "R_eq_pert = ", R_eq_pert
            ! write(*,*) "dR_eq = ", R_eq_pert - R_eq

            ! write(*,*) "cp = ", self%eq_solution%thermo%cp(i)
            ! write(*,*) "cp_pert = ", pert_soln%thermo%cp(i)
            ! write(*,*) "dcp = ", pert_soln%thermo%cp(i) - self%eq_solution%thermo%cp(i)

            dR_eq_dCpj(:, i) = (R_eq_pert - R_eq) / h_

            ! write(*,*) "FD dG_dCpj(",i,") = "
            ! do j = 1, self%num_eqn
            !     write(*,*) (pert_soln%G(j, k)-self%eq_solution%G(j, k), k = 1, self%num_eqn+1)
            ! end do

            pert_soln%thermo%cp(i) = pert_soln%thermo%cp(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_lnnj_dP(self, dR_lnnj_dP, h)
        ! Finite difference the partial derivative matrix of R_lnnj wrt P

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_lnnj_dP(:)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        logical :: const_p
        real(dp) :: h_  ! Step size
        real(dp) :: R_lnnj(self%eq_solver%num_gas), R_lnnj_pert(self%eq_solver%num_gas)  ! R_eq residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        const_p = self%eq_solution%constraints%is_constant_pressure()

        if (const_p) then
            R_lnnj = 0.0d0

            write(*,*) "d_lnnj = ", self%eq_solution%dln_nj
            pert_soln = self%eq_solution
            pert_soln%constraints%state2 = pert_soln%constraints%state2 + h_

            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_lnnj_pert = -pert_soln%lambda*pert_soln%dln_nj

            dR_lnnj_dP = (R_lnnj_pert - R_lnnj) / h_
        else
            ! TODO: implement this for volume problems
            dR_lnnj_dP = 0.0d0
        end if

    end subroutine

    subroutine EqTotals_fd_dR_lnnj_dpi(self, dR_lnnj_dpi, h)
        ! Finite difference the partial derivative matrix of R_lnnj wrt pi

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_lnnj_dpi(:,:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        integer :: ne
        real(dp) :: lambda
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_lnnj(self%eq_solver%num_gas), R_lnnj_pert(self%eq_solver%num_gas)

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ne = self%eq_solver%num_elements
        lambda = self%eq_solution%lambda

        R_lnnj = 0.0d0 ! -self%eq_solution%lambda*self%eq_solution%dln_nj
        do i = 1, ne
            pert_soln = self%eq_solution
            pert_soln%pi(i) = pert_soln%pi(i) + h_

            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_lnnj_pert = -pert_soln%lambda*pert_soln%dln_nj

            dR_lnnj_dpi(:, i) = (R_lnnj_pert - R_lnnj) / h_

            pert_soln%pi(i) = pert_soln%pi(i) - h_
            pert_soln%G(i,self%num_eqn+1) = pert_soln%G(i,self%num_eqn+1) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_lnnj_ddlnn(self, dR_lnnj_ddlnn, h)
        ! Finite difference the partial derivative matrix of R_lnnj wrt dln(n)

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_lnnj_ddlnn(:)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_lnnj(self%eq_solver%num_gas), R_lnnj_pert(self%eq_solver%num_gas)

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        R_lnnj = -self%eq_solution%lambda*self%eq_solution%dln_nj

        pert_soln = self%eq_solution
        pert_soln%dln_n = pert_soln%dln_n + h_

        call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
        R_lnnj_pert = -pert_soln%lambda*pert_soln%dln_nj

        dR_lnnj_ddlnn = (R_lnnj_pert - R_lnnj) / h_

    end subroutine

    subroutine EqTotals_fd_dR_lnnj_ddlnT(self, dR_lnnj_ddlnT, h)
        ! Finite difference the partial derivative matrix of R_lnnj wrt dln(T)

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_lnnj_ddlnT(:)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_lnnj(self%eq_solver%num_gas), R_lnnj_pert(self%eq_solver%num_gas)

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        R_lnnj = 0.0d0! -self%eq_solution%lambda*self%eq_solution%dln_nj

        pert_soln = self%eq_solution
        pert_soln%dln_T = pert_soln%dln_T + h_

        call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
        R_lnnj_pert = -pert_soln%lambda*pert_soln%dln_nj

        dR_lnnj_ddlnT = (R_lnnj_pert - R_lnnj) / h_

    end subroutine

    subroutine EqTotals_fd_dR_lnnj_dlnnj(self, dR_lnnj_dlnnj, h)
        ! Finite difference the partial derivative matrix of R_lnnj wrt lnnj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_lnnj_dlnnj(:,:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_lnnj(self%eq_solver%num_gas), R_lnnj_pert(self%eq_solver%num_gas)

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas

        R_lnnj = 0.0d0!-self%eq_solution%lambda*self%eq_solution%dln_nj
        do i = 1, ng
            pert_soln = self%eq_solution
            pert_soln%ln_nj(i) = pert_soln%ln_nj(i) + h_

            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_lnnj_pert = -pert_soln%lambda*pert_soln%dln_nj

            dR_lnnj_dlnnj(:, i) = (R_lnnj_pert - R_lnnj) / h_

            pert_soln%ln_nj(i) = pert_soln%ln_nj(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_lnnj_dn(self, dR_lnnj_dn, h)
        ! Finite difference the partial derivative matrix of R_lnnj wrt n

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_lnnj_dn(:)
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_lnnj(self%eq_solver%num_gas), R_lnnj_pert(self%eq_solver%num_gas)

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        R_lnnj = 0.0d0!-self%eq_solution%lambda*self%eq_solution%dln_nj

        pert_soln = self%eq_solution
        pert_soln%n = pert_soln%n + h_

        call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
        R_lnnj_pert = -pert_soln%lambda*pert_soln%dln_nj

        dR_lnnj_dn = (R_lnnj_pert - R_lnnj) / h_

    end subroutine

    subroutine EqTotals_fd_dR_lnnj_dHj(self, dR_lnnj_dHj, h)
        ! Finite difference the partial derivative matrix of R_lnnj wrt Hj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_lnnj_dHj(:,:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_lnnj(self%eq_solver%num_gas), R_lnnj_pert(self%eq_solver%num_gas)

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas

        R_lnnj = 0.0d0!-self%eq_solution%lambda*self%eq_solution%dln_nj
        do i = 1, ng
            pert_soln = self%eq_solution
            pert_soln%thermo%enthalpy(i) = pert_soln%thermo%enthalpy(i) + h_

            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_lnnj_pert = -pert_soln%lambda*pert_soln%dln_nj

            dR_lnnj_dHj(:, i) = (R_lnnj_pert - R_lnnj) / h_

            pert_soln%thermo%enthalpy(i) = pert_soln%thermo%enthalpy(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_lnnj_dSj(self, dR_lnnj_dSj, h)
        ! Finite difference the partial derivative matrix of R_lnnj wrt Sj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_lnnj_dSj(:,:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_lnnj(self%eq_solver%num_gas), R_lnnj_pert(self%eq_solver%num_gas)

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas

        R_lnnj = 0.0d0!-self%eq_solution%lambda*self%eq_solution%dln_nj
        do i = 1, ng
            pert_soln = self%eq_solution
            pert_soln%thermo%entropy(i) = pert_soln%thermo%entropy(i) + h_

            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_lnnj_pert = -pert_soln%lambda*pert_soln%dln_nj

            dR_lnnj_dSj(:, i) = (R_lnnj_pert - R_lnnj) / h_

            pert_soln%thermo%entropy(i) = pert_soln%thermo%entropy(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_nj_dlnnj(self, dR_nj_dlnnj, h)
        ! Finite difference the partial derivative matrix of R_nj wrt lnnj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_nj_dlnnj(:,:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_nj(self%eq_solver%num_gas), R_nj_pert(self%eq_solver%num_gas)

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas

        R_nj = 0.0d0!self%eq_solution%nj - exp(self%eq_solution%ln_nj)
        do i = 1, ng
            pert_soln = self%eq_solution
            pert_soln%ln_nj(i) = pert_soln%ln_nj(i) + h_

            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_nj_pert = self%eq_solution%nj - exp(pert_soln%ln_nj)

            dR_nj_dlnnj(:, i) = (R_nj_pert - R_nj) / h_

            pert_soln%ln_nj(i) = pert_soln%ln_nj(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_n_dP(self, dR_n_dP, h)
        ! Finite difference the partial derivative matrix of R_n wrt P

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_n_dP
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        logical :: const_p
        real(dp) :: h_  ! Step size
        real(dp) :: R_n, R_n_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        const_p = self%eq_solution%constraints%is_constant_pressure()

        if (const_p) then
            R_n = -self%eq_solution%lambda*self%eq_solution%dln_n

            pert_soln = self%eq_solution
            pert_soln%n = pert_soln%n + h_
            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_n_pert = -pert_soln%lambda*pert_soln%dln_n

            dR_n_dP = (R_n_pert - R_n) / h_
        else
            ! TODO: implement this for volume problems
            dR_n_dP = 0.0d0
        end if

    end subroutine

    subroutine EqTotals_fd_dR_n_dpi(self, dR_n_dpi, h)
        ! Finite difference the partial derivative matrix of R_n wrt pi

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_n_dpi(:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_n, R_n_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        R_n = -self%eq_solution%lambda*self%eq_solution%dln_n

        do i = 1, self%eq_solver%num_elements
            pert_soln = self%eq_solution
            pert_soln%pi(i) = pert_soln%pi(i) + h_
            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_n_pert = -pert_soln%lambda*pert_soln%dln_n

            dR_n_dpi(i) = (R_n_pert - R_n) / h_

            pert_soln%pi(i) = pert_soln%pi(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_n_ddlnn(self, dR_n_ddlnn, h)
        ! Finite difference the partial derivative matrix of R_n wrt dlnn

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_n_ddlnn
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_n, R_n_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        R_n = -self%eq_solution%lambda*self%eq_solution%dln_n

        pert_soln = self%eq_solution
        pert_soln%dln_n = pert_soln%dln_n + h_
        call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
        R_n_pert = -pert_soln%lambda*pert_soln%dln_n

        dR_n_ddlnn = (R_n_pert - R_n) / h_

    end subroutine

    subroutine EqTotals_fd_dR_n_ddlnT(self, dR_n_ddlnT, h)
        ! Finite difference the partial derivative matrix of R_n wrt dlnT

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_n_ddlnT
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_n, R_n_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        R_n = -self%eq_solution%lambda*self%eq_solution%dln_n

        pert_soln = self%eq_solution
        pert_soln%dln_T = pert_soln%dln_T + h_
        call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
        R_n_pert = -pert_soln%lambda*pert_soln%dln_n

        dR_n_ddlnT = (R_n_pert - R_n) / h_

    end subroutine

    subroutine EqTotals_fd_dR_n_dn(self, dR_n_dn, h)
        ! Finite difference the partial derivative matrix of R_n wrt n

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_n_dn
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_n, R_n_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        R_n = -self%eq_solution%lambda*self%eq_solution%dln_n

        pert_soln = self%eq_solution
        pert_soln%n = pert_soln%n + h_
        call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
        R_n_pert = -pert_soln%lambda*pert_soln%dln_n

        dR_n_dn = (R_n_pert - R_n) / h_

    end subroutine

    subroutine EqTotals_fd_dR_n_dlnnj(self, dR_n_dlnnj, h)
        ! Finite difference the partial derivative matrix of R_n wrt lnnj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_n_dlnnj(:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_n, R_n_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas

        R_n = -self%eq_solution%lambda*self%eq_solution%dln_n

        do i = 1, ng
            pert_soln = self%eq_solution
            pert_soln%ln_nj(i) = pert_soln%ln_nj(i) + h_
            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_n_pert = -pert_soln%lambda*pert_soln%dln_n

            dR_n_dlnnj(i) = (R_n_pert - R_n) / h_

            pert_soln%ln_nj(i) = pert_soln%ln_nj(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_n_dHj(self, dR_n_dHj, h)
        ! Finite difference the partial derivative matrix of R_n wrt Hj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_n_dHj(:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_n, R_n_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas

        R_n = -self%eq_solution%lambda*self%eq_solution%dln_n

        do i = 1, ng
            pert_soln = self%eq_solution
            pert_soln%thermo%enthalpy(i) = pert_soln%thermo%enthalpy(i) + h_
            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_n_pert = -pert_soln%lambda*pert_soln%dln_n

            dR_n_dHj(i) = (R_n_pert - R_n) / h_

            pert_soln%thermo%enthalpy(i) = pert_soln%thermo%enthalpy(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_n_dSj(self, dR_n_dSj, h)
        ! Finite difference the partial derivative matrix of R_n wrt Sj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_n_dSj(:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_n, R_n_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas

        R_n = -self%eq_solution%lambda*self%eq_solution%dln_n

        do i = 1, ng
            pert_soln = self%eq_solution
            pert_soln%thermo%entropy(i) = pert_soln%thermo%entropy(i) + h_
            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_n_pert = -pert_soln%lambda*pert_soln%dln_n

            dR_n_dSj(i) = (R_n_pert - R_n) / h_

            pert_soln%thermo%entropy(i) = pert_soln%thermo%entropy(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_T_dP(self, dR_T_dP, h)
        ! Finite difference the partial derivative matrix of R_T wrt P

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_T_dP
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        logical :: const_p
        real(dp) :: h_  ! Step size
        real(dp) :: R_T, R_T_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        const_p = self%eq_solution%constraints%is_constant_pressure()

        if (const_p) then
            R_T = -self%eq_solution%lambda*self%eq_solution%dln_T

            pert_soln = self%eq_solution
            pert_soln%T = pert_soln%T + h_
            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_T_pert = -pert_soln%lambda*pert_soln%dln_T

            dR_T_dP = (R_T_pert - R_T) / h_
        else
            ! TODO: implement this for volume problems
            dR_T_dP = 0.0d0
        end if

    end subroutine

    subroutine EqTotals_fd_dR_T_dpi(self, dR_T_dpi, h)
        ! Finite difference the partial derivative matrix of R_T wrt pi

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_T_dpi(:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_T, R_T_pert  ! R_T residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        R_T = -self%eq_solution%lambda*self%eq_solution%dln_T

        do i = 1, self%eq_solver%num_elements
            pert_soln = self%eq_solution
            pert_soln%pi(i) = pert_soln%pi(i) + h_
            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_T_pert = -pert_soln%lambda*pert_soln%dln_T

            dR_T_dpi(i) = (R_T_pert - R_T) / h_

            pert_soln%pi(i) = pert_soln%pi(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_T_ddlnn(self, dR_T_ddlnn, h)
        ! Finite difference the partial derivative matrix of R_T wrt dlnn

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_T_ddlnn
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_T, R_T_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        R_T = -self%eq_solution%lambda*self%eq_solution%dln_n

        pert_soln = self%eq_solution
        pert_soln%dln_n = pert_soln%dln_n + h_
        call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
        R_T_pert = -pert_soln%lambda*pert_soln%dln_T

        dR_T_ddlnn = (R_T_pert - R_T) / h_

    end subroutine

    subroutine EqTotals_fd_dR_T_ddlnT(self, dR_T_ddlnT, h)
        ! Finite difference the partial derivative matrix of R_T wrt dlnT

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_T_ddlnT
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_T, R_T_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        R_T = -self%eq_solution%lambda*self%eq_solution%dln_T

        pert_soln = self%eq_solution
        pert_soln%dln_T = pert_soln%dln_T + h_
        call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
        R_T_pert = -pert_soln%lambda*pert_soln%dln_T

        dR_T_ddlnT = (R_T_pert - R_T) / h_

    end subroutine

    subroutine EqTotals_fd_dR_T_dn(self, dR_T_dn, h)
        ! Finite difference the partial derivative matrix of R_T wrt n

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_T_dn
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_T, R_T_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        R_T = -self%eq_solution%lambda*self%eq_solution%dln_T

        pert_soln = self%eq_solution
        pert_soln%n = pert_soln%n + h_
        call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
        R_T_pert = -pert_soln%lambda*pert_soln%dln_T

        dR_T_dn = (R_T_pert - R_T) / h_

    end subroutine

    subroutine EqTotals_fd_dR_T_dlnnj(self, dR_T_dlnnj, h)
        ! Finite difference the partial derivative matrix of R_T wrt lnnj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_T_dlnnj(:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_T, R_T_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas

        R_T = -self%eq_solution%lambda*self%eq_solution%dln_T

        do i = 1, ng
            pert_soln = self%eq_solution
            pert_soln%ln_nj(i) = pert_soln%ln_nj(i) + h_
            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_T_pert = -pert_soln%lambda*pert_soln%dln_T

            dR_T_dlnnj(i) = (R_T_pert - R_T) / h_

            pert_soln%ln_nj(i) = pert_soln%ln_nj(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_T_dHj(self, dR_T_dHj, h)
        ! Finite difference the partial derivative matrix of R_T wrt Hj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_T_dHj(:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_T, R_T_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas

        R_T = -self%eq_solution%lambda*self%eq_solution%dln_T

        do i = 1, ng
            pert_soln = self%eq_solution
            pert_soln%thermo%enthalpy(i) = pert_soln%thermo%enthalpy(i) + h_
            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_T_pert = -pert_soln%lambda*pert_soln%dln_T

            dR_T_dHj(i) = (R_T_pert - R_T) / h_

            pert_soln%thermo%enthalpy(i) = pert_soln%thermo%enthalpy(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_T_dSj(self, dR_T_dSj, h)
        ! Finite difference the partial derivative matrix of R_T wrt Sj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_T_dSj(:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_T, R_T_pert  ! R_n residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas

        R_T = -self%eq_solution%lambda*self%eq_solution%dln_T

        do i = 1, ng
            pert_soln = self%eq_solution
            pert_soln%thermo%entropy(i) = pert_soln%thermo%entropy(i) + h_
            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            R_T_pert = -pert_soln%lambda*pert_soln%dln_T

            dR_T_dSj(i) = (R_T_pert - R_T) / h_

            pert_soln%thermo%entropy(i) = pert_soln%thermo%entropy(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dR_b_dnj(self, dR_b_dnj, h)
        ! Finite difference the partial derivative matrix of R_b wrt nj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dR_b_dnj(:,:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: R_b(self%eq_solver%num_elements), R_b_pert(self%eq_solver%num_elements)  ! R_b residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas

        R_b = 0.0d0

        do i = 1, ng
            pert_soln = self%eq_solution
            pert_soln%nj(i) = pert_soln%nj(i) + h_
            R_b_pert = self%eq_solution%b - self%eq_solver%products%elements_from_species(pert_soln%nj)

            dR_b_dnj(:, i) = (R_b_pert - R_b) / h_

            pert_soln%nj(i) = pert_soln%nj(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dlambda_dP(self, dlambda_dP, h)
        ! Finite difference the partial derivative matrix of lambda wrt P

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dlambda_dP
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        logical :: const_p
        real(dp) :: lambda, lambda_pert  ! R_eq residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        lambda = self%eq_solution%lambda
        const_p = self%eq_solution%constraints%is_constant_pressure()
        pert_soln = self%eq_solution

        if (const_p) then
            pert_soln%constraints%state2 = pert_soln%constraints%state2 + h_
            lambda_pert = self%eq_solver%compute_damped_update_factor(pert_soln)

            dlambda_dP = (lambda_pert - lambda) / h_
        else
            ! TODO: implement this for volume problems
            dlambda_dP = 0.0d0
        end if


    end subroutine

    subroutine EqTotals_fd_dlambda_dn(self, dlambda_dn, h)
        ! Finite difference the partial derivative matrix of lambda wrt n

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dlambda_dn
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: lambda, lambda_pert  ! R_eq residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        lambda = self%eq_solution%lambda

        pert_soln = self%eq_solution
        pert_soln%n = pert_soln%n + h_
        lambda_pert = self%eq_solver%compute_damped_update_factor(pert_soln)

        dlambda_dn = (lambda_pert - lambda) / h_

    end subroutine

    subroutine EqTotals_fd_dlambda_dlnnj(self, dlambda_dlnnj, h)
        ! Finite difference the partial derivative matrix of lambda wrt ln_nj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dlambda_dlnnj(:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: lambda, lambda_pert  ! R_eq residuals
        integer :: i

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas
        lambda = self%eq_solution%lambda

        pert_soln = self%eq_solution
        do i = 1, ng
            pert_soln%ln_nj(i) = pert_soln%ln_nj(i) + h_
            lambda_pert = self%eq_solver%compute_damped_update_factor(pert_soln)

            dlambda_dlnnj(i) = (lambda_pert - lambda) / h_
            pert_soln%ln_nj(i) = pert_soln%ln_nj(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dlambda_ddlnnj(self, dlambda_ddlnnj, h)
        ! Finite difference the partial derivative matrix of lambda wrt dln_nj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dlambda_ddlnnj(:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: ng  ! Number of gas species
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: lambda, lambda_pert  ! R_eq residuals
        integer :: i

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas
        lambda = self%eq_solution%lambda

        pert_soln = self%eq_solution
        do i = 1, ng
            pert_soln%dln_nj(i) = pert_soln%dln_nj(i) + h_
            lambda_pert = self%eq_solver%compute_damped_update_factor(pert_soln)

            dlambda_ddlnnj(i) = (lambda_pert - lambda) / h_
            pert_soln%dln_nj(i) = pert_soln%dln_nj(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dlambda_ddlnn(self, dlambda_ddlnn, h)
        ! Finite difference the partial derivative matrix of lambda wrt dln(n)

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dlambda_ddlnn
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: lambda, lambda_pert  ! R_eq residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        lambda = self%eq_solution%lambda

        pert_soln = self%eq_solution
        pert_soln%dln_n = pert_soln%dln_n + h_
        lambda_pert = self%eq_solver%compute_damped_update_factor(pert_soln)

        dlambda_ddlnn = (lambda_pert - lambda) / h_

    end subroutine

    subroutine EqTotals_fd_dlambda_ddlnT(self, dlambda_ddlnT, h)
        ! Finite difference the partial derivative matrix of lambda wrt dln(T)

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dlambda_ddlnT
        real(dp), intent(in), optional :: h

        ! Locals
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: lambda, lambda_pert  ! R_eq residuals

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        lambda = self%eq_solution%lambda

        pert_soln = self%eq_solution
        pert_soln%dln_T = pert_soln%dln_T + h_
        lambda_pert = self%eq_solver%compute_damped_update_factor(pert_soln)

        dlambda_ddlnT = (lambda_pert - lambda) / h_

    end subroutine

    subroutine EqTotals_fd_dlambda_dpi(self, dlambda_dpi, h)
        ! Finite difference the partial derivative matrix of lambda wrt pi

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dlambda_dpi(:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: ne
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: lambda, lambda_pert
        integer :: i

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ne = self%eq_solver%num_elements
        lambda = self%eq_solution%lambda

        pert_soln = self%eq_solution
        do i = 1, ne
            pert_soln%pi(i) = pert_soln%pi(i) + h_
            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            lambda_pert = pert_soln%lambda!self%eq_solver%compute_damped_update_factor(pert_soln)

            dlambda_dpi(i) = (lambda_pert - lambda) / h_
            pert_soln%pi(i) = pert_soln%pi(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dlambda_dHj(self, dlambda_dHj, h)
        ! Finite difference the partial derivative matrix of lambda wrt Hj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dlambda_dHj(:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: lambda, lambda_pert
        integer :: i

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas
        lambda = self%eq_solution%lambda

        do i = 1, ng
            pert_soln = self%eq_solution
            pert_soln%thermo%enthalpy(i) = pert_soln%thermo%enthalpy(i) + h_
            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            lambda_pert = pert_soln%lambda!self%eq_solver%compute_damped_update_factor(pert_soln)

            dlambda_dHj(i) = (lambda_pert - lambda) / h_
            pert_soln%thermo%enthalpy(i) = pert_soln%thermo%enthalpy(i) - h_
        end do

    end subroutine

    subroutine EqTotals_fd_dlambda_dSj(self, dlambda_dSj, h)
        ! Finite difference the partial derivative matrix of lambda wrt Sj

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(inout) :: dlambda_dSj(:)
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: ng
        type(EqSolution) :: pert_soln  ! Perturbed solution
        real(dp) :: h_  ! Step size
        real(dp) :: lambda, lambda_pert
        integer :: i

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ng = self%eq_solver%num_gas
        lambda = self%eq_solution%lambda

        do i = 1, ng
            pert_soln = self%eq_solution
            pert_soln%thermo%entropy(i) = pert_soln%thermo%entropy(i) + h_
            call self%eq_solver%update_solution(pert_soln, fd_mode=.true.)
            lambda_pert = pert_soln%lambda!self%eq_solver%compute_damped_update_factor(pert_soln)

            dlambda_dSj(i) = (lambda_pert - lambda) / h_
            pert_soln%thermo%entropy(i) = pert_soln%thermo%entropy(i) - h_
        end do

    end subroutine

    subroutine EqTotals_check_dR_dx(self, h)
        ! Finite difference each sub-Jacobian in dR/dx, compare the analytic derivatives and print the result

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(in), optional :: h

        ! Locals
        real(dp) :: h_  ! Step size
        integer :: i, j
        integer :: ne

        ! More locals: sub-Jacobians of dR/dx
        real(dp), pointer :: dR_eq_dH0(:)
        real(dp), pointer :: dR_eq_db0(:,:)
        real(dp), pointer :: dR_eq_dP(:)
        real(dp), pointer :: dR_lnnj_dP(:)
        real(dp), pointer :: dR_n_dP
        real(dp), pointer :: dR_T_dP

        ! More locals: finite difference sub-Jacobians of dR/dx
        real(dp) :: dR_eq_dH0_fd(self%num_eqn)
        real(dp) :: dR_eq_db0_fd(self%num_eqn, self%eq_solver%num_elements)
        real(dp) :: dR_eq_dP_fd(self%num_eqn)
        real(dp) :: dR_lnnj_dP_fd(self%eq_solver%num_gas)
        real(dp) :: dR_n_dP_fd
        real(dp) :: dR_T_dP_fd

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ! Define shorthand
        ne = self%eq_solver%num_elements

        ! Get the analytic partials (** assumes compute_dR_dx has been called **)
        dR_eq_dH0 => self%dR_eq_dH0
        dR_eq_db0 => self%dR_eq_db0
        dR_eq_dP => self%dR_eq_dP
        dR_lnnj_dP => self%dR_lnnj_dP
        dR_n_dP => self%dR_n_dP
        dR_T_dP => self%dR_T_dP

        ! FINITE DIFFERENCE CHECK
        call self%fd_dR_eq_dH0(dR_eq_dH0_fd)
        call self%fd_dR_eq_db0(dR_eq_db0_fd)
        call self%fd_dR_eq_dP(dR_eq_dP_fd)
        call self%fd_dR_lnnj_dP(dR_lnnj_dP_fd)
        call self%fd_dR_n_dP(dR_n_dP_fd)
        call self%fd_dR_T_dP(dR_T_dP_fd)

        write(*,*) "----------------------------------------------"
        write(*,*) "Finite difference check of dR/du"
        write(*,*) "----------------------------------------------"

        write(*,*)
        write(*,*) "dR_eq/dH0"
        write(*, 2200) "analytic:", "check:"
        do i = 1, self%num_eqn
            write(*,2100) dR_eq_dH0(i), dR_eq_dH0_fd(i)
        end do

        ! TODO: why aren't these printing??
        write(*,*) ""
        write(*,*) "dR_eq/db0"
        write(*,*) "analytic:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_db0(i, j), j = 1,ne)
        end do

        write(*,*) "check:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_db0_fd(i, j), j = 1,ne)
        end do

        write(*,*)
        write(*,*) "dR_eq/dP"
        write(*, 2200) "analytic:", "check:"
        do i = 1, self%num_eqn
            write(*,2100) dR_eq_dP(i), dR_eq_dP_fd(i)
        end do

        write(*,*)
        write(*,*) "dR_lnnj/dP"
        write(*, 2200) "analytic:", "check:"
        do i = 1, self%num_eqn
            write(*,2100) dR_lnnj_dP(i), dR_lnnj_dP_fd(i)
        end do

        write(*,*)
        write(*,*) "dR_n/dP"
        write(*, 2200) "analytic:", "check:"
        write(*,2100) dR_n_dP, dR_n_dP_fd

        write(*,*)
        write(*,*) "dR_T/dP"
        write(*, 2200) "analytic:", "check:"
        write(*,2100) dR_T_dP, dR_T_dP_fd

        10100 format(10F20.8)
        2100 format(2F20.8)
        2200 format(2A15)

    end subroutine

    subroutine EqTotals_check_dR_du(self, h)
        ! Finite difference each sub-Jacobian in dR/du, compare the analytic derivatives and print the result

        ! Arguments
        class(EqTotals), target :: self
        real(dp), intent(in), optional :: h

        ! Locals
        integer :: i, j
        integer :: ne, ng
        real(dp) :: h_  ! Step size

        ! More locals: sub-Jacobians of dR/du
        real(dp), pointer :: dR_eq_duhat(:,:)
        real(dp), pointer :: dR_eq_dnj(:,:)
        real(dp), pointer :: dR_eq_dn(:)
        real(dp), pointer :: dR_eq_db(:,:)
        real(dp), pointer :: dR_eq_dlnnj(:,:)
        real(dp), pointer :: dR_eq_dHj(:,:)
        real(dp), pointer :: dR_eq_dSj(:,:)
        real(dp), pointer :: dR_eq_dUj(:,:)
        real(dp), pointer :: dR_eq_dCpj(:,:)
        real(dp), pointer :: dR_eq_dT(:)
        real(dp), pointer :: dR_nj_dlnnj(:,:)
        real(dp), pointer :: dR_nj_dnj(:,:)
        real(dp), pointer :: dR_lnnj_dpi(:,:)
        real(dp), pointer :: dR_lnnj_ddlnn(:)
        real(dp), pointer :: dR_lnnj_ddlnT(:)
        real(dp), pointer :: dR_lnnj_dn(:)
        real(dp), pointer :: dR_lnnj_dHj(:,:)
        real(dp), pointer :: dR_lnnj_dSj(:,:)
        real(dp), pointer :: dR_lnnj_dlnnj(:,:)
        real(dp), pointer :: dR_n_dpi(:)
        real(dp), pointer :: dR_n_ddlnn
        real(dp), pointer :: dR_n_ddlnT
        real(dp), pointer :: dR_n_dlnnj(:)
        real(dp), pointer :: dR_n_dn
        real(dp), pointer :: dR_n_dHj(:)
        real(dp), pointer :: dR_n_dSj(:)
        real(dp), pointer :: dR_T_dpi(:)
        real(dp), pointer :: dR_T_ddlnn
        real(dp), pointer :: dR_T_ddlnT
        real(dp), pointer :: dR_T_dlnnj(:)
        real(dp), pointer :: dR_T_dn
        real(dp), pointer :: dR_T_dT
        real(dp), pointer :: dR_T_dHj(:)
        real(dp), pointer :: dR_T_dSj(:)
        real(dp), pointer :: dR_b_dnj(:,:)
        real(dp), pointer :: dR_b_db(:,:)
        real(dp), pointer :: dR_Hj_dT(:)
        real(dp), pointer :: dR_Hj_dHj(:,:)
        real(dp), pointer :: dR_Sj_dT(:)
        real(dp), pointer :: dR_Sj_dSj(:,:)
        real(dp), pointer :: dR_Uj_dT(:)
        real(dp), pointer :: dR_Uj_dUj(:,:)
        real(dp), pointer :: dR_Cp_dT(:)
        real(dp), pointer :: dR_Cp_dCp(:,:)

        ! More locals: finite difference values
        real(dp) :: dR_eq_dlnnj_fd(self%num_eqn, self%eq_solver%num_gas)
        real(dp) :: dR_eq_dnj_fd(self%num_eqn, self%eq_solver%num_gas)
        real(dp) :: dR_eq_dn_fd(self%num_eqn)
        real(dp) :: dR_eq_dT_fd(self%num_eqn)
        real(dp) :: dR_eq_db_fd(self%num_eqn, self%eq_solver%num_elements)
        real(dp) :: dR_eq_dHj_fd(self%num_eqn, self%eq_solver%num_gas)
        real(dp) :: dR_eq_dSj_fd(self%num_eqn, self%eq_solver%num_gas)
        real(dp) :: dR_eq_dUj_fd(self%num_eqn, self%eq_solver%num_gas)
        real(dp) :: dR_eq_dCpj_fd(self%num_eqn, self%eq_solver%num_gas)
        !real(dp) :: dR_eq_dCvj_fd(self%num_eqn, self%eq_solver%num_gas)

        real(dp) :: dR_lnnj_dpi_fd(self%eq_solver%num_gas, self%eq_solver%num_elements)
        real(dp) :: dR_lnnj_ddlnn_fd(self%eq_solver%num_gas)
        real(dp) :: dR_lnnj_ddlnT_fd(self%eq_solver%num_gas)
        real(dp) :: dR_lnnj_dlnnj_fd(self%eq_solver%num_gas, self%eq_solver%num_gas)
        real(dp) :: dR_lnnj_dn_fd(self%eq_solver%num_gas)
        real(dp) :: dR_lnnj_dHj_fd(self%eq_solver%num_gas, self%eq_solver%num_gas)
        real(dp) :: dR_lnnj_dSj_fd(self%eq_solver%num_gas, self%eq_solver%num_gas)

        real(dp) :: dR_nj_dlnnj_fd(self%eq_solver%num_gas, self%eq_solver%num_gas)

        real(dp) :: dR_n_dpi_fd(self%eq_solver%num_elements)
        real(dp) :: dR_n_ddlnn_fd
        real(dp) :: dR_n_ddlnT_fd
        real(dp) :: dR_n_dn_fd
        real(dp) :: dR_n_dlnnj_fd(self%eq_solver%num_gas)
        real(dp) :: dR_n_dHj_fd(self%eq_solver%num_gas)
        real(dp) :: dR_n_dSj_fd(self%eq_solver%num_gas)

        real(dp) :: dR_T_dpi_fd(self%eq_solver%num_elements)
        real(dp) :: dR_T_ddlnn_fd
        real(dp) :: dR_T_ddlnT_fd
        real(dp) :: dR_T_dlnnj_fd(self%eq_solver%num_gas)
        real(dp) :: dR_T_dn_fd
        real(dp) :: dR_T_dHj_fd(self%eq_solver%num_gas)
        real(dp) :: dR_T_dSj_fd(self%eq_solver%num_gas)

        real(dp) :: dR_b_dnj_fd(self%eq_solver%num_elements, self%eq_solver%num_gas)

        ! Optional argument handling
        h_ = 1.0d-8
        if (present(h)) h_ = h

        ! Define shorthand
        ne = self%eq_solver%num_elements
        ng = self%eq_solver%num_gas

        ! Get the analytic partials (** assumes compute_dR_dx has been called **)
        dR_eq_duhat        => self%dR_eq_duhat
        dR_eq_dnj          => self%dR_eq_dnj
        dR_eq_dn           => self%dR_eq_dn
        dR_eq_db           => self%dR_eq_db
        dR_eq_dlnnj        => self%dR_eq_dlnnj
        dR_eq_dHj          => self%dR_eq_dHj
        dR_eq_dSj          => self%dR_eq_dSj
        dR_eq_dUj          => self%dR_eq_dUj
        dR_eq_dCpj         => self%dR_eq_dCpj
        dR_eq_dT           => self%dR_eq_dT
        dR_nj_dlnnj        => self%dR_nj_dlnnj
        dR_nj_dnj          => self%dR_nj_dnj
        dR_lnnj_dpi        => self%dR_lnnj_dpi
        dR_lnnj_ddlnn      => self%dR_lnnj_ddlnn
        dR_lnnj_ddlnT      => self%dR_lnnj_ddlnT
        dR_lnnj_dn         => self%dR_lnnj_dn
        dR_lnnj_dHj        => self%dR_lnnj_dHj
        dR_lnnj_dSj        => self%dR_lnnj_dSj
        dR_lnnj_dlnnj      => self%dR_lnnj_dlnnj
        dR_n_dpi           => self%dR_n_dpi
        dR_n_ddlnn         => self%dR_n_ddlnn
        dR_n_ddlnT         => self%dR_n_ddlnT
        dR_n_dlnnj         => self%dR_n_dlnnj
        dR_n_dn            => self%dR_n_dn
        dR_n_dHj           => self%dR_n_dHj
        dR_n_dSj           => self%dR_n_dSj
        dR_T_dpi           => self%dR_T_dpi
        dR_T_ddlnn         => self%dR_T_ddlnn
        dR_T_ddlnT         => self%dR_T_ddlnT
        dR_T_dlnnj         => self%dR_T_dlnnj
        dR_T_dn            => self%dR_T_dn
        dR_T_dT            => self%dR_T_dT
        dR_T_dHj           => self%dR_T_dHj
        dR_T_dSj           => self%dR_T_dSj
        dR_b_dnj           => self%dR_b_dnj
        dR_b_db            => self%dR_b_db
        dR_Hj_dT           => self%dR_Hj_dT
        dR_Hj_dHj          => self%dR_Hj_dHj
        dR_Sj_dT           => self%dR_Sj_dT
        dR_Sj_dSj          => self%dR_Sj_dSj
        dR_Uj_dT           => self%dR_Uj_dT
        dR_Uj_dUj          => self%dR_Uj_dUj
        dR_Cp_dT           => self%dR_Cp_dT
        dR_Cp_dCp          => self%dR_Cp_dCp

        ! FINITE DIFFERENCE CHECK
        call self%fd_dR_eq_dnj(dR_eq_dnj_fd)
        call self%fd_dR_eq_dn(dR_eq_dn_fd)
        call self%fd_dR_eq_db(dR_eq_db_fd)
        call self%fd_dR_eq_dlnnj(dR_eq_dlnnj_fd)
        call self%fd_dR_eq_dHj(dR_eq_dHj_fd)
        call self%fd_dR_eq_dSj(dR_eq_dSj_fd)
        call self%fd_dR_eq_dUj(dR_eq_dUj_fd)
        call self%fd_dR_eq_dCpj(dR_eq_dCpj_fd)
        call self%fd_dR_eq_dT(dR_eq_dT_fd)
        call self%fd_dR_lnnj_dpi(dR_lnnj_dpi_fd)
        call self%fd_dR_lnnj_ddlnn(dR_lnnj_ddlnn_fd)
        call self%fd_dR_lnnj_ddlnT(dR_lnnj_ddlnT_fd)
        call self%fd_dR_lnnj_dn(dR_lnnj_dn_fd)
        call self%fd_dR_lnnj_dlnnj(dR_lnnj_dlnnj_fd)
        call self%fd_dR_lnnj_dHj(dR_lnnj_dHj_fd)
        call self%fd_dR_lnnj_dSj(dR_lnnj_dSj_fd)
        call self%fd_dR_nj_dlnnj(dR_nj_dlnnj_fd)
        call self%fd_dR_n_dpi(dR_n_dpi_fd)
        call self%fd_dR_n_ddlnn(dR_n_ddlnn_fd)
        call self%fd_dR_n_ddlnT(dR_n_ddlnT_fd)
        call self%fd_dR_n_dlnnj(dR_n_dlnnj_fd)
        call self%fd_dR_n_dn(dR_n_dn_fd)
        call self%fd_dR_n_dHj(dR_n_dHj_fd)
        call self%fd_dR_n_dSj(dR_n_dSj_fd)
        call self%fd_dR_T_dpi(dR_T_dpi_fd)
        call self%fd_dR_T_ddlnn(dR_T_ddlnn_fd)
        call self%fd_dR_T_ddlnT(dR_T_ddlnT_fd)
        call self%fd_dR_T_dlnnj(dR_T_dlnnj_fd)
        call self%fd_dR_T_dn(dR_T_dn_fd)
        call self%fd_dR_T_dHj(dR_T_dHj_fd)
        call self%fd_dR_T_dSj(dR_T_dSj_fd)
        call self%fd_dR_b_dnj(dR_b_dnj_fd)

        write(*,*) "----------------------------------------------"
        write(*,*) "Finite difference check of dR/du"
        write(*,*) "----------------------------------------------"

        write(*,*) ""
        write(*,*) "dR_eq/dnj"
        write(*,*) "analytic:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_dnj(i, j), j = 1,ng)
        end do

        write(*,*) "check:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_dnj_fd(i, j), j = 1,ng)
        end do

        write(*,*)
        write(*,*) "dR_eq/dn"
        write(*, 2200) "analytic:", "check:"
        do i = 1, self%num_eqn
            write(*,2100) dR_eq_dn(i), dR_eq_dn_fd(i)
        end do

        write(*,*)
        write(*,*) "dR_eq/dT"
        write(*, 2200) "analytic:", "check:"
        do i = 1, self%num_eqn
            write(*,2100) dR_eq_dT(i), dR_eq_dT_fd(i)
        end do

        write(*,*) ""
        write(*,*) "dR_eq/db"
        write(*,*) "analytic:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_db(i, j), j = 1,ne)
        end do

        write(*,*) "check:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_db_fd(i, j), j = 1,ne)
        end do

        write(*,*) ""
        write(*,*) "dR_eq/dlnnj"
        write(*,*) "analytic:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_dlnnj(i, j), j = 1,ng)
        end do

        write(*,*) "check:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_dlnnj_fd(i, j), j = 1,ng)
        end do

        write(*,*) ""
        write(*,*) "dR_eq/dHj"
        write(*,*) "analytic:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_dHj(i, j), j = 1,ng)
        end do

        write(*,*) "check:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_dHj_fd(i, j), j = 1,ng)
        end do

        write(*,*) ""
        write(*,*) "dR_eq/dSj"
        write(*,*) "analytic:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_dSj(i, j), j = 1,ng)
        end do

        write(*,*) "check:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_dSj_fd(i, j), j = 1,ng)
        end do

        write(*,*) ""
        write(*,*) "dR_eq/dUj"
        write(*,*) "analytic:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_dUj(i, j), j = 1,ng)
        end do

        write(*,*) "check:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_dUj_fd(i, j), j = 1,ng)
        end do

        write(*,*) ""
        write(*,*) "dR_eq/dCpj"
        write(*,*) "analytic:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_dCpj(i, j), j = 1,ng)
        end do

        write(*,*) "check:"
        do i = 1, self%num_eqn
            write(*,10100) (dR_eq_dCpj_fd(i, j), j = 1,ng)
        end do

        write(*,*) ""
        write(*,*) "dR_lnnj/dpi"
        write(*,*) "analytic:"
        do i = 1, ng
            write(*,10100) (dR_lnnj_dpi(i, j), j = 1,ne)
        end do

        write(*,*) "check:"
        do i = 1, ng
            write(*,10100) (dR_lnnj_dpi_fd(i, j), j = 1,ne)
        end do

        write(*,*)
        write(*,*) "dR_lnnj/ddlnn"
        write(*, 2200) "analytic:", "check:"
        do i = 1, ng
            write(*,2100) dR_lnnj_ddlnn(i), dR_lnnj_ddlnn_fd(i)
        end do

        write(*,*)
        write(*,*) "dR_lnnj/ddlnT"
        write(*, 2200) "analytic:", "check:"
        do i = 1, ng
            write(*,2100) dR_lnnj_ddlnT(i), dR_lnnj_ddlnT_fd(i)
        end do

!         write(*,*) "FD dR_lnnj_dlnnj = "
!         do i = 1, ng
!             write(*,*) (dR_lnnj_dlnnj_fd(i,j), j = 1, ng)
!         end do

        write(*,*)
        write(*,*) "dR_lnnj/dn"
        write(*, 2200) "analytic:", "check:"
        do i = 1, ng
            write(*,2100) dR_lnnj_dn(i), dR_lnnj_dn_fd(i)
        end do

        write(*,*) ""
        write(*,*) "dR_lnnj/dHj"
        write(*,*) "analytic:"
        do i = 1, ng
            write(*,10100) (dR_lnnj_dHj(i, j), j = 1,ng)
        end do

        write(*,*) "check:"
        do i = 1, ng
            write(*,10100) (dR_lnnj_dHj_fd(i, j), j = 1,ng)
        end do

        write(*,*) ""
        write(*,*) "dR_lnnj/dSj"
        write(*,*) "analytic:"
        do i = 1, ng
            write(*,10100) (dR_lnnj_dSj(i, j), j = 1,ng)
        end do

        write(*,*) "check:"
        do i = 1, ng
            write(*,10100) (dR_lnnj_dSj_fd(i, j), j = 1,ng)
        end do

        write(*,*) ""
        write(*,*) "dR_nj/dlnnj"
        write(*,*) "analytic:"
        do i = 1, ng
            write(*,10100) (dR_nj_dlnnj(i, j), j = 1,ng)
        end do

        write(*,*) "check:"
        do i = 1, ng
            write(*,10100) (dR_nj_dlnnj_fd(i, j), j = 1,ng)
        end do

        write(*,*)
        write(*,*) "dR_n/dpi"
        write(*, 2200) "analytic:", "check:"
        do i = 1, ne
            write(*,2100) dR_n_dpi(i), dR_n_dpi_fd(i)
        end do

        write(*,*)
        write(*,*) "dR_n/ddlnn"
        write(*, 2200) "analytic:", "check:"
        write(*,2100) dR_n_ddlnn, dR_n_ddlnn_fd

        write(*,*)
        write(*,*) "dR_n/ddlnT"
        write(*, 2200) "analytic:", "check:"
        write(*,2100) dR_n_ddlnT, dR_n_ddlnT_fd

        write(*,*)
        write(*,*) "dR_n/dlnnj"
        write(*, 2200) "analytic:", "check:"
        do i = 1, ng
            write(*,2100) dR_n_dlnnj(i), dR_n_dlnnj_fd(i)
        end do

        write(*,*)
        write(*,*) "dR_n/dHj"
        write(*, 2200) "analytic:", "check:"
        do i = 1, ng
            write(*,2100) dR_n_dHj(i), dR_n_dHj_fd(i)
        end do

        write(*,*)
        write(*,*) "dR_n/dSj"
        write(*, 2200) "analytic:", "check:"
        do i = 1, ng
            write(*,2100) dR_n_dSj(i), dR_n_dSj_fd(i)
        end do

        write(*,*)
        write(*,*) "dR_T/dpi"
        write(*, 2200) "analytic:", "check:"
        do i = 1, ne
            write(*,2100) dR_T_dpi(i), dR_T_dpi_fd(i)
        end do

        write(*,*)
        write(*,*) "dR_T/ddlnn"
        write(*, 2200) "analytic:", "check:"
        write(*,2100) dR_T_ddlnn, dR_T_ddlnn_fd

        write(*,*)
        write(*,*) "dR_T/ddlnT"
        write(*, 2200) "analytic:", "check:"
        write(*,2100) dR_T_ddlnT, dR_T_ddlnT_fd

        write(*,*)
        write(*,*) "dR_T/dlnnj"
        write(*, 2200) "analytic:", "check:"
        do i = 1, ng
            write(*,2100) dR_T_dlnnj(i), dR_T_dlnnj_fd(i)
        end do

        write(*,*)
        write(*,*) "dR_T/dn"
        write(*, 2200) "analytic:", "check:"
        write(*,2100) dR_T_dn, dR_T_dn_fd

        write(*,*)
        write(*,*) "dR_T/dHj"
        write(*, 2200) "analytic:", "check:"
        do i = 1, ng
            write(*,2100) dR_T_dHj(i), dR_T_dHj_fd(i)
        end do

        write(*,*)
        write(*,*) "dR_T/dSj"
        write(*, 2200) "analytic:", "check:"
        do i = 1, ng
            write(*,2100) dR_T_dSj(i), dR_T_dSj_fd(i)
        end do

        write(*,*) ""
        write(*,*) "dR_b/dnj"
        write(*,*) "analytic:"
        do i = 1, ne
            write(*,10100) (dR_b_dnj(i, j), j = 1,ng)
        end do

        write(*,*) "check:"
        do i = 1, ne
            write(*,10100) (dR_b_dnj_fd(i, j), j = 1,ng)
        end do

        10100 format(10F20.8)
        2100 format(2F20.8)
        2200 format(2A15)

    end subroutine

end module