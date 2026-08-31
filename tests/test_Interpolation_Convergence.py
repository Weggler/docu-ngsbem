import pytest

from convergence_timing.Interpolation_Convergence import (
    run_pi1_pi0_convergence,
    run_pid_convergence,
)


def assert_convergence(results, error, expected_order):
    errors = [result[error] for result in results]
    alpha = f"alpha_{error.removeprefix('err_')}"
    observed_orders = [result[alpha] for result in results]

    assert errors[0] > errors[1] > errors[2]
    assert min(observed_orders) >= 0.75 * expected_order


@pytest.mark.parametrize(
    ("order", "expected_pi1_rate", "expected_pi0_rate"),
    [
        (1, 2.0, 1.0),
        (3, 4.0, 3.0),
    ],
)
def test_pi1_pi0_convergence(order, expected_pi1_rate, expected_pi0_rate):
    results = run_pi1_pi0_convergence(
        orders=(order,),
        refinement_levels=range(4),
    )[1:]

    assert_convergence(results, "err_int_1", expected_pi1_rate)
    assert_convergence(results, "err_int_0", expected_pi0_rate)


@pytest.mark.parametrize(
    ("order", "expected_rate"),
    [
        (1, 1.0),
        (3, 3.0),
    ],
)
def test_pid_convergence(order, expected_rate):
    results = run_pid_convergence(
        orders=(order,),
        refinement_levels=range(4),
    )[1:]

    assert_convergence(results, "err_int_d", expected_rate)
