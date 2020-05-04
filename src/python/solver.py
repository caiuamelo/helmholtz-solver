from callbacks import dispatch_callbacks
from complex_linear_system import convert_complex_linear_system_to_real,\
    extract_complex_solution


def solve_2d_hellmholtz(
    mesh,
    omega,
    callbacks=None,
):
    if callbacks is None:
        callbacks = {}

    Ke, f = _prepare_linear_system(mesh, omega)
    Ke, f = convert_complex_linear_system_to_real(Ke, f)
    P = extract_complex_solution(_solve_linear_system(Ke, f))

    dispatch_callbacks(callbacks, 'on_after_solve_2d_hellmholtz', P=P, mesh=mesh)


def _prepare_linear_system(mesh, omega):
    # TODO: Implement
    return None, None


def _solve_linear_system(Ke, f):
    # TODO: Implement
    return None
