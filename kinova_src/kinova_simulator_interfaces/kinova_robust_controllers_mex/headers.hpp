#ifndef HEADERS
#define HEADERS

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <boost/numeric/interval.hpp>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <chrono>

#define _USE_MATH_DEFINES

namespace bn = boost::numeric;
namespace bi = bn::interval_lib;

using Interval = bn::interval<
        double, 
        bi::policies<
            bi::save_state<bi::rounded_transc_std<double> >,
            bi::checking_base<double>
        > 
    >;

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

typedef Eigen::Matrix<Interval, Eigen::Dynamic, 1> VectorXint;

using std::vector;
using std::string;
using std::map;
using std::stringstream;
using std::ostream;
using std::ifstream;

using std::max;
using std::min;

using std::cout;
using std::endl;

#endif