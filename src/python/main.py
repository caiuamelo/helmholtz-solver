from solver import solve_2d_hellmholtz
from mesh import build_mesh
import plotting
import time


if __name__ == "__main__":
    import globals as g
    g.set_cpp(True)

    s = time.time()
    mesh = build_mesh(size_x=200.0, size_y=200.0, nx=200, ny=200)
    solve_2d_hellmholtz(mesh, omega=150.0)
    print(time.time() - s)
