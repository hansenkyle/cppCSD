#include "solver.h"

#include <array>
#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <vector>


Solver::Kernel::Kernel() {
    M << 2.0, 1.0, 1.0, 2.0;
    M*= (1.0/6);
    L << 0.5, 0.5, -0.5, -0.5;
    Lb << -1, 0, 0, -1;

    A = Eigen::Matrix4d::Zero();
    b = Eigen::Vector4d::Zero();
}