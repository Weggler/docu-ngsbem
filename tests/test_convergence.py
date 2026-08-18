from math import log

import pytest

from convergence_timing.Laplace_DtN_Convergence import run_laplace_dtn
from convergence_timing.Laplace_NtD_Convergence import run_laplace_ntd
from convergence_timing.Maxwell_Mie_Convergence import run_maxwell_mie


def convergence_rate(results):
    """For fixed order on a two dimensional surface mesh, ndof ~ h**(-2).
    Thus error ~ h**p is equivalent to error ~ ndof**(-p/2).
    """

    coarse, fine = results[0], results[-1]
    return log(coarse["err"] / fine["err"]) / log(fine["ndof"] / coarse["ndof"])


@pytest.mark.parametrize(
    ("runner", "order", "expected_rate"),
    [
        pytest.param(run_laplace_dtn, 1, 0.5, id="laplace-dtn"),
        pytest.param(run_laplace_ntd, 1, 0.5, id="laplace-ntd"),
        pytest.param(run_maxwell_mie, 1, 1.0, id="maxwell-mie"),
        pytest.param(run_laplace_dtn, 3, 1.5, id="laplace-dtn-high-order"),
        pytest.param(run_laplace_ntd, 3, 1.5, id="laplace-ntd-high-order"),
        pytest.param(run_maxwell_mie, 3, 2.0, id="maxwell-mie-high-order"),
    ],
)
def test_convergence(runner, order, expected_rate):
    results = [runner(order, refinement) for refinement in (1, 2, 3)]
    errors = [result["err"] for result in results]

    assert errors[0] > errors[1] > errors[2]
    assert convergence_rate(results) >= 0.75 * expected_rate
