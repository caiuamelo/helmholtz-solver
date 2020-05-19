from matplotlib import pyplot as plt


def plot_results(P, mesh, *args, **kwargs):
    P = P.real
    P = P.reshape((mesh.nx + 1, mesh.ny + 1))
    plt.imshow(P)
    plt.show()


def get_callbacks():
    return {
        "on_after_solve_2d_hellmholtz": [plot_results,],
    }
