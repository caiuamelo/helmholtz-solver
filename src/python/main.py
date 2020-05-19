from solver import solve_2d_hellmholtz
from mesh import build_mesh
import plotting

if __name__ == "__main__":
    mesh = build_mesh(
        size_x=200.0,
        size_y=200.0,
        nx=80,
        ny=80
    )
    solve_2d_hellmholtz(
        mesh,
        omega=150.0,
        callbacks=plotting.get_callbacks()
    )
