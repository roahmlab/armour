#ifndef HEADER_H
#define HEADER_H

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include "device_launch_parameters.h"
#include <cstdio>
#include <omp.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <fstream>
#include <cstring>
#include <boost/numeric/interval.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <cassert>
#include <vector>
#include <cstdint>

#define WARNING_PRINT printf

// intervals
namespace bn = boost::numeric;
namespace bi = bn::interval_lib;

using Interval = bn::interval<
        double, 
        bi::policies<
            bi::save_state<bi::rounded_transc_std<double> >,
            bi::checking_base<double>
        > 
    >;

// interval matrices
namespace Eigen {
    namespace internal {
        template<typename X, typename S, typename P>
        struct is_convertible<X, bn::interval<S,P> > {
            enum { value = is_convertible<X,S>::value };
        };

        template<typename S, typename P1, typename P2>
        struct is_convertible<bn::interval<S,P1>, bn::interval<S,P2> > {
            enum { value = true };
        };
    }
}

typedef Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic> MatrixXInt;

using std::vector;
using std::cout;
using std::endl;
using std::sort;
using std::swap;
using std::min;
using std::max;

#endif