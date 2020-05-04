from solver import solve_2d_hellmholtz
from mesh import build_mesh
import plotting

if __name__ == "__main__":
    mesh = build_mesh()
    solve_2d_hellmholtz(
        mesh,
        omega=300.0,
        callbacks=plotting.get_callbacks()
    )
