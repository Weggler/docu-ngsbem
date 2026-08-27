from math import log

import pytest

from convergence_timing.Interpolation_Convergence import (
    run_pi0,
    run_pi1,
    run_pid,
)


def convergence_rate(results):
    """For fixed order on a two dimensional surface mesh, ndof ~ h**(-2).
    Thus error ~ h**p is equivalent to error ~ ndof**(-p/2).
    """

    coarse, fine = results[0], results[-1]
    return log(coarse["err"] / fine["err"]) / log(fine["ndof"] / coarse["ndof"])


@pytest.mark.parametrize(
    ("runner", "order", "expected_rate"),
    [
        pytest.param(run_pi1, 1, 1.0, id="pi1"),
        pytest.param(run_pi0, 1, 0.5, id="pi0"),
        pytest.param(run_pid, 1, 0.5, id="pid"),
        pytest.param(run_pi1, 3, 2.0, id="pi1-high-order"),
        pytest.param(run_pi0, 3, 1.5, id="pi0-high-order"),
        pytest.param(run_pid, 3, 1.5, id="pid-high-order"),
    ],
)
def test_convergence(runner, order, expected_rate):
    results = [runner(order, refinement) for refinement in (1, 2, 3)]
    errors = [result["err"] for result in results]

    assert errors[0] > errors[1] > errors[2]
    assert convergence_rate(results) >= 0.75 * expected_rate
