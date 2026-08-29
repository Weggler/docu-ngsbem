from math import log

import pytest

from convergence_timing.Interpolation_Convergence import (
    run_pi1_pi0_convergence,
    run_pid_convergence,
)


def convergence_rate(results, error, ndof):
    """For fixed order on a two dimensional surface mesh, ndof ~ h**(-2).
    Thus error ~ h**p is equivalent to error ~ ndof**(-p/2).
    """

    coarse, fine = results[0], results[-1]
    return log(coarse[error] / fine[error]) / log(fine[ndof] / coarse[ndof])


def assert_convergence(results, error, ndof, expected_rate):
    errors = [result[error] for result in results]

    assert errors[0] > errors[1] > errors[2]
    assert convergence_rate(results, error, ndof) >= 0.75 * expected_rate


@pytest.mark.parametrize(
    ("order", "expected_pi1_rate", "expected_pi0_rate"),
    [
        (1, 1.0, 0.5),
        (3, 2.0, 1.5),
    ],
)
def test_pi1_pi0_convergence(order, expected_pi1_rate, expected_pi0_rate):
    results = run_pi1_pi0_convergence(
        orders=(order,),
        refinement_levels=range(4),
    )[1:]

    assert_convergence(results, "err_int_1", "ndof_Pi1", expected_pi1_rate)
    assert_convergence(results, "err_int_0", "ndof_Pi0", expected_pi0_rate)


@pytest.mark.parametrize(
    ("order", "expected_rate"),
    [
        (1, 0.5),
        (3, 1.5),
    ],
)
def test_pid_convergence(order, expected_rate):
    results = run_pid_convergence(
        orders=(order,),
        refinement_levels=range(4),
    )[1:]

    assert_convergence(results, "err_int_d", "ndof", expected_rate)
