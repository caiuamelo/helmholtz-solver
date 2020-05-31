#ifndef _FORWARD_SOLVER_H
#define _FORWARD_SOLVER_H

#include <Eigen/Eigen>

namespace fwi_ls {

Eigen::Array<std::complex<double>, 4, 4> build_local_Ke(
    Eigen::Ref<Eigen::Array<double, 4, 2> const> const& element_points,
    double omega,
    double mu,
    double eta
);

Eigen::Array<double, 4, 1> build_local_f(
    Eigen::Array<double, 4, 2> element_points,
    double S_e
);

}

#endif