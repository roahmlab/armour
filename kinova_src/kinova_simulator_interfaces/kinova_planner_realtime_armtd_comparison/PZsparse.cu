#ifndef PZ_SPARSE_CPP
#define PZ_SPARSE_CPP

#include "PZsparse.h"

/*
Helper functions
*/

double getCenter(const Interval& a) {
    return (a.lower() + a.upper()) * 0.5;
}

double getRadius(const Interval& a) {
    return (a.upper() - a.lower()) * 0.5;
}

Eigen::MatrixXd getCenter(const MatrixXInt& a) {
    Eigen::MatrixXd res(a.rows(), a.cols());

    for (uint i = 0; i < a.rows(); i++) {
        for (uint j = 0; j < a.cols(); j++) {
            res(i, j) = getCenter(a(i, j));
        }
    }

    return res;
}

Eigen::MatrixXd getRadius(const MatrixXInt& a) {
    Eigen::MatrixXd res(a.rows(), a.cols());

    for (uint i = 0; i < a.rows(); i++) {
        for (uint j = 0; j < a.cols(); j++) {
            res(i, j) = getRadius(a(i, j));
        }
    }

    return res;
}

bool Monomial_sorter_degree(Monomial const& lhs, Monomial const& rhs) {
    return lhs.degree < rhs.degree;
}

/*
Initialization
*/

PZsparse::PZsparse(uint NRows_inp, uint NCols_inp) {
    NRows = NRows_inp;
    NCols = NCols_inp;
    center = Eigen::MatrixXd::Zero(NRows, NCols);
    independent = Eigen::MatrixXd::Zero(NRows, NCols);
}

PZsparse::PZsparse(const PZsparse& pz_inp) {
    NRows = pz_inp.NRows;
    NCols = pz_inp.NCols;
    center = pz_inp.center;
    polynomial = pz_inp.polynomial;
    independent = pz_inp.independent;
}

// 1x1 PZ
PZsparse::PZsparse(double center_inp) {
    NRows = 1;
    NCols = 1;
    center.resize(NRows, NCols);
    center(0) = center_inp;
    independent = Eigen::MatrixXd::Zero(NRows, NCols);
}

// NxM PZ
PZsparse::PZsparse(const Eigen::MatrixXd& center_inp) {
    NRows = center_inp.rows();
    NCols = center_inp.cols();
    center = center_inp;
    independent = Eigen::MatrixXd::Zero(NRows, NCols);
}

// 1x1 PZ
// PZsparse::PZsparse(double center_inp, double uncertainty_percent) {
//     NRows = 1;
//     NCols = 1;
//     center.resize(NRows, NCols);
//     center(0) = center_inp;
//     independent.resize(NRows, NCols);
//     independent(0) = uncertainty_percent * fabs(center_inp);
// }

// NxM PZ
PZsparse::PZsparse(const Eigen::MatrixXd& center_inp, double uncertainty_percent) {
    NRows = center_inp.rows();
    NCols = center_inp.cols();
    center = center_inp;
    independent = uncertainty_percent * center_inp.cwiseAbs();
}

// // 1x1 PZ
// PZsparse::PZsparse(Interval interval_inp) {
//     NRows = 1;
//     NCols = 1;
//     center.resize(NRows, NCols);
//     center(0) = getCenter(interval_inp);
//     independent.resize(NRows, NCols);
//     independent(0) = getRadius(interval_inp);
// }

PZsparse::PZsparse(double center_inp, Interval independent_inp) {
    NRows = 1;
    NCols = 1;
    center.resize(NRows, NCols);
    center(0) = center_inp + getCenter(independent_inp);
    independent.resize(NRows, NCols);
    independent(0) = getRadius(independent_inp);
}

// 1x1 PZ
PZsparse::PZsparse(double center_inp, double* coeff_inp, uint64_t degree_inp[][NUM_FACTORS * 6], uint num_monomials) {
    NRows = 1;
    NCols = 1;
    
    center.resize(NRows, NCols);
    center(0) = center_inp;

    polynomial.reserve(num_monomials);

    for (uint i = 0; i < num_monomials; i++) {
        polynomial.emplace_back(coeff_inp[i], convertDegreeToHash(degree_inp[i]));
    }

    independent = Eigen::MatrixXd::Zero(NRows, NCols);

    simplify();
}

// 1x1 PZ
PZsparse::PZsparse(double center_inp, double* coeff_inp, uint64_t degree_inp[][NUM_FACTORS * 6], uint num_monomials, Interval independent_inp) {
    NRows = 1;
    NCols = 1;
    
    center.resize(NRows, NCols);
    center(0) = center_inp;

    polynomial.reserve(num_monomials);

    for (uint i = 0; i < num_monomials; i++) {
        polynomial.emplace_back(coeff_inp[i], convertDegreeToHash(degree_inp[i]));
    }

    // assume independent_inp is centered at 0
    independent.resize(NRows, NCols);
    independent(0) = getRadius(independent_inp);

    simplify();
}

// 3x3 PZ
PZsparse::PZsparse(const double roll, const double pitch, const double yaw) {
    NRows = 3;
    NCols = 3;
    center.resize(NRows, NCols);
    
    center(0,0) = cos(pitch)*cos(yaw);
    center(0,1) = -cos(pitch)*sin(yaw);
    center(0,2) = sin(pitch);
    center(1,0) = cos(roll)*sin(yaw) + cos(yaw)*sin(pitch)*sin(roll);
    center(1,1) = cos(roll)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw);
    center(1,2) = -cos(pitch)*sin(roll);
    center(2,0) = sin(roll)*sin(yaw) - cos(roll)*cos(yaw)*sin(pitch);
    center(2,1) = cos(yaw)*sin(roll) + cos(roll)*sin(pitch)*sin(yaw);
    center(2,2) = cos(pitch)*cos(roll);

    independent = Eigen::MatrixXd::Zero(NRows, NCols);
}

// 3x3 PZ
PZsparse::PZsparse(double cos_center_inp, double* cos_coeff_inp, uint64_t cos_degree_inp[][NUM_FACTORS * 6], uint cos_num_monomials,
                   double sin_center_inp, double* sin_coeff_inp, uint64_t sin_degree_inp[][NUM_FACTORS * 6], uint sin_num_monomials,
                   const uint axis) {
    NRows = 3;
    NCols = 3;
    
    makeRotationMatrix(center, cos_center_inp, sin_center_inp, axis);

    polynomial.reserve(cos_num_monomials + sin_num_monomials);

    Eigen::MatrixXd coeff_temp;
    for (uint i = 0; i < cos_num_monomials; i++) {
        makeRotationMatrix(coeff_temp, cos_coeff_inp[i], 0, axis, true);
        polynomial.emplace_back(coeff_temp, convertDegreeToHash(cos_degree_inp[i]));
    }

    for (uint i = 0; i < sin_num_monomials; i++) {
        makeRotationMatrix(coeff_temp, 0, sin_coeff_inp[i], axis, true);
        polynomial.emplace_back(coeff_temp, convertDegreeToHash(sin_degree_inp[i]));
    }

    // assume independent_inp is centered at 0
    // makeRotationMatrix(independent, getRadius(cos_independent_inp), getRadius(sin_independent_inp), axis);
    independent = Eigen::MatrixXd::Zero(3, 3);

    simplify();
}

/*
Internal functions
*/

void PZsparse::makeRotationMatrix(Eigen::MatrixXd& R, const double cosElt, const double sinElt, const uint axis, bool startFromZero) {
    if (startFromZero) {
        R = Eigen::MatrixXd::Zero(3,3);
    }
    else {
        R = Eigen::MatrixXd::Identity(3,3);
    }

    const double negSinElt = -1.0 * sinElt;

    switch (axis) {
        case 0: // fixed joints
            // NOTE:
            // This is just a simplified implementation!!!
            // The rotation matrix is just identity matrix for all fixed joints of Fetch
            // don't do anything
            return;
        case 1: // rx
            R(1,1) = cosElt;
            R(1,2) = negSinElt;
            R(2,1) = sinElt;
            R(2,2) = cosElt;
            break;
        case 2: // ry
            R(0,0) = cosElt;
            R(0,2) = sinElt;
            R(2,0) = negSinElt;
            R(2,2) = cosElt;
            break;
        case 3: // rz
            R(0,0) = cosElt;
            R(0,1) = negSinElt;
            R(1,0) = sinElt;
            R(1,1) = cosElt;
            break;
        default:
            WARNING_PRINT("Undefined axis");
            throw -1;
    }
}

bool PZsparse::internalCheck() const {
    if (center.rows() != NRows) {
        WARNING_PRINT("PZsparse error: center matrix number of rows not consistent!");
        return false;
    }
    if (center.cols() != NCols) {
        WARNING_PRINT("PZsparse error: center matrix number of columns not consistent!");
        return false;
    }
    if (independent.rows() != NRows) {
        WARNING_PRINT("PZsparse error: independent generator matrix number of rows not consistent!");
        return false;
    }
    if (independent.cols() != NCols) {
        WARNING_PRINT("PZsparse error: independent generator matrix number of columns not consistent!");
        return false;
    }
    for (uint i = 0; i < independent.rows(); i++) {
        for (uint j = 0; j < independent.cols(); j++) {
            if (independent(i, j) < 0) {
                WARNING_PRINT("PZsparse error: independent generator matrix has negative entry!");
                return false;
            }
        }
    }
    return true;
}

void PZsparse::simplify() {
    assert(internalCheck());

    sort(polynomial.begin(), polynomial.end(), Monomial_sorter_degree);

    Eigen::MatrixXd reduce_amount(NRows, NCols); 
    reduce_amount.setZero();

    vector<Monomial> polynomial_new;
    polynomial_new.reserve(polynomial.size());

    uint i = 0;
    while (i < polynomial.size()) {
        uint j;
        uint64_t degree = polynomial[i].degree;

        for (j = i + 1; j < polynomial.size(); j++) {
            if (polynomial[j].degree != degree) {
                break;
            }

            polynomial[i].coeff += polynomial[j].coeff;
        }

        Eigen::MatrixXd temp = polynomial[i].coeff;
        if (temp.norm() <= SIMPLIFY_THRESHOLD) {
            reduce_amount += temp.cwiseAbs();
        }
        else {
            polynomial_new.emplace_back(polynomial[i]);
        }

        i = j;
    }

    // for (uint i = 0; i < polynomial.size(); i++) {
    //     if (polynomial[i].coeff.norm() == 0) {
    //         continue;
    //     }

    //     uint64_t degree = polynomial[i].degree;

    //     for (uint j = 0; j < polynomial.size(); j++) {
    //         if (j == i || polynomial[j].coeff.norm() == 0) {
    //             continue;
    //         }
    //         if (polynomial[j].degree == degree) {
    //             polynomial[i].coeff += polynomial[j].coeff;
    //             polynomial[j].coeff.setZero();
    //         }
    //     }

    //     Eigen::MatrixXd temp = polynomial[i].coeff;
    //     if (temp.norm() <= SIMPLIFY_THRESHOLD) {
    //         reduce_amount += temp.cwiseAbs();
    //     }
    //     else {
    //         polynomial_new.emplace_back(polynomial[i]);
    //     }
    // }

    polynomial = polynomial_new;

    if (reduce_amount.norm() != 0) {
        independent = independent + reduce_amount;
    }
}

void PZsparse::reduce() {
    assert(internalCheck());

    vector<Monomial> polynomial_new;
    polynomial_new.reserve(polynomial.size());

    for (auto it : polynomial) {
        if (it.degree < max_hash_dependent_k_only) { // only dependent on k
            polynomial_new.emplace_back(it.coeff, it.degree);
        }
        else {
            independent += it.coeff.cwiseAbs();
        }
    }

    polynomial = polynomial_new;
}

Eigen::MatrixXd PZsparse::reduce_link_PZ() {
    assert(internalCheck());
    assert(NRows == 3 && NCols == 1);

    Eigen::MatrixXd link_independent_generators(3, 3 + 3);
    link_independent_generators.setZero();

    vector<Monomial> polynomial_new;
    polynomial_new.reserve(polynomial.size());

    int j = 0;

    for (auto it : polynomial) {
        if (it.degree < max_hash_dependent_k_only) { // only dependent on k
            polynomial_new.emplace_back(it.coeff, it.degree);
        }
        else if (it.degree < max_hash_dependent_k_links_only && (it.degree & dependent_k_mask) == 0) { // only dependent on link x, y, z generators
            assert(j < 3);
            link_independent_generators.col(j++) = it.coeff;
        }
        else {
            independent += it.coeff.cwiseAbs();
        }
    }

    polynomial = polynomial_new;

    link_independent_generators(0, 3) = independent(0);
    link_independent_generators(1, 4) = independent(1);
    link_independent_generators(2, 5) = independent(2);

    return link_independent_generators;
}

MatrixXInt PZsparse::slice(const double* factor) {
    assert(internalCheck());
    MatrixXInt res(NRows, NCols);
    Eigen::MatrixXd res_center = center;
    Eigen::MatrixXd res_radius = independent;

    for (auto it : polynomial) {
        Eigen::MatrixXd resTemp = it.coeff;

        if (it.degree < (1 << (2 * NUM_FACTORS))) { // only dependent on k
            convertHashToDegree(it.degree);

            for (uint j = 0; j < NUM_FACTORS; j++) {
                resTemp *= pow(factor[j], degreeArray[j]);
            }

            res_center += resTemp;
        }
        else { // this line should never be triggered if you run reduce first
            res_radius += resTemp.cwiseAbs();
        }
    }

    for (uint i = 0; i < NRows; i++) {
        for (uint j = 0; j < NCols; j++) {
            res(i,j) = Interval(res_center(i,j) - res_radius(i,j),
                                res_center(i,j) + res_radius(i,j));
        }
    }

    return res;
}

void PZsparse::slice(Eigen::MatrixXd* gradient, const double* factor) {
    assert(internalCheck());

    for (uint k = 0; k < NUM_FACTORS; k++) {
        gradient[k] = Eigen::MatrixXd::Zero(NRows, NCols);
    }

    Eigen::Array<Eigen::MatrixXd, NUM_FACTORS, 1> resTemp;

    for (auto it : polynomial) {
        if (it.degree <= (1 << (2 * NUM_FACTORS))) { // only dependent on k
            for (uint k = 0; k < NUM_FACTORS; k++) {
                resTemp[k] = it.coeff;
            }

            convertHashToDegree(it.degree);

            for (uint j = 0; j < NUM_FACTORS; j++) {
                for (uint k = 0; k < NUM_FACTORS; k++) {
                    if (j == k) { // differentiate this!
                        if (degreeArray[j] == 0) { // monomial unrelated to k
                            resTemp[k] = Eigen::MatrixXd::Zero(NRows, NCols);
                        }
                        else {
                            resTemp[k] *= degreeArray[j] * pow(factor[j], degreeArray[j] - 1);
                        }
                    }
                    else {
                        resTemp[k] *= pow(factor[j], degreeArray[j]);
                    }
                }
            }

            for (uint k = 0; k < NUM_FACTORS; k++) {
                gradient[k] += resTemp[k];
            }
        }
    }
}

void PZsparse::slice(Eigen::Vector3d* gradient, const double* factor) {
    assert(NRows == 3 && NCols == 1);
    assert(internalCheck());

    for (uint k = 0; k < NUM_FACTORS; k++) {
        gradient[k].setZero();
    }

    Eigen::Array<Eigen::Vector3d, NUM_FACTORS, 1> resTemp;

    for (auto it : polynomial) {
        if (it.degree <= (1 << (2 * NUM_FACTORS))) { // only dependent on k
            for (uint k = 0; k < NUM_FACTORS; k++) {
                resTemp[k] = it.coeff;
            }

            convertHashToDegree(it.degree);

            for (uint j = 0; j < NUM_FACTORS; j++) {
                for (uint k = 0; k < NUM_FACTORS; k++) {
                    if (j == k) { // differentiate this!
                        if (degreeArray[j] == 0) { // monomial unrelated to k
                            resTemp[k].setZero();
                        }
                        else {
                            resTemp[k] *= degreeArray[j] * pow(factor[j], degreeArray[j] - 1);
                        }
                    }
                    else {
                        resTemp[k] *= pow(factor[j], degreeArray[j]);
                    }
                }
            }

            for (uint k = 0; k < NUM_FACTORS; k++) {
                gradient[k] += resTemp[k];
            }
        }
    }
}

void PZsparse::slice(double* gradient, const double* factor) {
    assert(internalCheck());
    assert(NRows == 1 && NCols == 1);

    memset(gradient, 0, NUM_FACTORS * sizeof(double));

    double resTemp[NUM_FACTORS] = {0};

    for (auto it : polynomial) {
        if (it.degree <= (1 << (2 * NUM_FACTORS))) { // only dependent on k
            for (uint k = 0; k < NUM_FACTORS; k++) {
                resTemp[k] = it.coeff(0);
            }

            convertHashToDegree(it.degree);

            for (uint j = 0; j < NUM_FACTORS; j++) {
                for (uint k = 0; k < NUM_FACTORS; k++) {
                    if (j == k) { // differentiate this!
                        if (degreeArray[j] == 0) { // monomial unrelated to k
                            resTemp[k] = 0;
                        }
                        else {
                            resTemp[k] *= degreeArray[j] * pow(factor[j], degreeArray[j] - 1);
                        }
                    }
                    else {
                        resTemp[k] *= pow(factor[j], degreeArray[j]);
                    }
                }
            }

            for (uint k = 0; k < NUM_FACTORS; k++) {
                gradient[k] += resTemp[k];
            }
        }
    }
}

MatrixXInt PZsparse::toInterval() {
    assert(internalCheck());

    MatrixXInt res(NRows, NCols);
    Eigen::MatrixXd res_center = center;
    Eigen::MatrixXd res_radius = independent;

    for (auto it : polynomial) {
        res_radius += it.coeff.cwiseAbs();
    }

    for (uint i = 0; i < NRows; i++) {
        for (uint j = 0; j < NCols; j++) {
            res(i,j) = Interval(res_center(i,j) - res_radius(i,j),
                                res_center(i,j) + res_radius(i,j));
        }
    }

    return res;
}

void PZsparse::convertHashToDegree(uint64_t degree) {
    for (uint64_t i = 0; i < NUM_FACTORS * 6; i++) {
        degreeArray[i] = degree & DEGREE_MASK[i];    
        degree >>= MOVE_BIT_INC[i];
    }

    return;
}

uint64_t convertDegreeToHash(const uint64_t* degreeArray) {
    uint64_t degree = 0;
    uint64_t move_bit = 0;

    for (uint64_t i = 0; i < NUM_FACTORS * 6; i++) {
        if (degreeArray[i] > 1) {
            WARNING_PRINT("degree can not be larger than 1!");
            throw;
        }

        degree += (degreeArray[i] << move_bit);

        move_bit += MOVE_BIT_INC[i];
    }

    return degree;
}

std::ostream& operator<<(std::ostream& os, PZsparse& a) {
    // if (independent_only) {
    //     Interval temp = center + independent;
    //     cout << "[ " << temp.lower() << ", " << temp.upper() << " ]\n\n";
    //     return;
    // }

    os << a.center << " +...\n";

    for (auto it : a.polynomial) {
        os << '(' << it.coeff << ')';
        
        a.convertHashToDegree(it.degree);
        
        os << " * k^(";
        for (uint j = 0; j < NUM_FACTORS; j++) {
            os << a.degreeArray[j];
        }
        os << ") ";

        os << " * qde^(";
        for (uint j = 0; j < NUM_FACTORS; j++) {
            os << a.degreeArray[j + NUM_FACTORS * 1];
        }
        os << ") ";

        os << " * qdae^(";
        for (uint j = 0; j < NUM_FACTORS; j++) {
            os << a.degreeArray[j + NUM_FACTORS * 2];
        }
        os << ") ";

        os << " * qddae^(";
        for (uint j = 0; j < NUM_FACTORS; j++) {
            os << a.degreeArray[j + NUM_FACTORS * 3];
        }
        os << ") ";

        os << " * cosqe^(";
        for (uint j = 0; j < NUM_FACTORS; j++) {
            os << a.degreeArray[j + NUM_FACTORS * 4];
        }
        os << ") ";

        os << " * sinqe^(";
        for (uint j = 0; j < NUM_FACTORS; j++) {
            os << a.degreeArray[j + NUM_FACTORS * 5];
        }
        os << ") ";

        os << " +...\n";
    }

    os << "[ " << -a.independent << ", " << a.independent << " ]\n\n";

    return os;
}

std::ostream& operator<<(std::ostream& os, const MatrixXInt& a) {
    for (uint i = 0; i < a.rows(); i++) {
        for (uint j = 0; j < a.cols(); j++) {
            os << "[ " << a(i,j).lower() << ", " << a(i,j).upper() << "] ";
        }
        os << '\n';
    }

    return os;
}

/*
Arithmetic
*/

PZsparse PZsparse::operator() (int row_id, int col_id) const {
    assert(internalCheck());
    assert(row_id < NRows);
    assert(col_id < NCols);

    PZsparse res(1, 1);
    
    res.center = center.block(row_id, col_id, 1, 1);

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(it.coeff.block(row_id, col_id, 1, 1), it.degree);
    }

    res.independent = independent.block(row_id, col_id, 1, 1);

    return res;

}

PZsparse PZsparse::operator=(const double a) {
    NRows = 1;
    NCols = 1;
    center.resize(NRows, NCols);
    center(0) = a;
    polynomial.clear();
    independent = Eigen::MatrixXd::Zero(NRows, NCols);
    return *this;
}

// PZsparse PZsparse::operator=(const Interval& a) {
//     center = getCenter(a);
//     polynomial.clear();
//     independent = a - center;
//     return *this;
// }

PZsparse PZsparse::operator=(const PZsparse& a) {
    NRows = a.NRows;
    NCols = a.NCols;
    center = a.center;
    polynomial = a.polynomial;
    independent = a.independent;
    return *this;
}

PZsparse PZsparse::operator-() {
    assert(internalCheck());

    PZsparse res(NRows, NCols);
    
    res.center = -center;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(-it.coeff, it.degree);
    }

    res.independent = -independent;

    return res;
}

PZsparse PZsparse::operator+(const PZsparse& a) {
    assert(internalCheck());
    assert(a.NRows == NRows || a.NCols == NCols); // check if they are add-able

    PZsparse res(NRows, NCols);

    res.center = center + a.center;

    res.polynomial.reserve(polynomial.size() + a.polynomial.size());

    res.polynomial.insert(res.polynomial.end(), polynomial.begin(), polynomial.end());

    for (auto it : a.polynomial) {
        res.polynomial.push_back(it);
    }

    res.independent = independent + a.independent;

    res.simplify();
    
    return res;
}

PZsparse PZsparse::operator+(const double a) {
    assert(internalCheck());

    PZsparse res = *this;

    for (uint i = 0; i < NRows; i++) {
        for (uint j = 0; j < NCols; j++) {
            res.center(i,j) += a;
        }
    }

    return res;
}

PZsparse operator+(const double a, const PZsparse& b) {
    assert(b.internalCheck());

    PZsparse res = b;

    for (uint i = 0; i < b.NRows; i++) {
        for (uint j = 0; j < b.NCols; j++) {
            res.center(i,j) += a;
        }
    }

    return res;
}

PZsparse PZsparse::operator+=(const PZsparse& a) {
    assert(internalCheck());
    assert(a.NRows == NRows || a.NCols == NCols); // check if they are add-able

    center += a.center;

    polynomial.reserve(polynomial.size() + a.polynomial.size());

    for (auto it : a.polynomial) {
        polynomial.push_back(it);
    }

    independent += a.independent;

    simplify();
    
    return *this;
}

PZsparse PZsparse::operator-(const PZsparse& a) {
    assert(internalCheck());    
    assert(a.NRows == NRows || a.NCols == NCols); // check if they are add-able

    PZsparse res(NRows, NCols);

    res.center = center - a.center;

    res.polynomial.reserve(polynomial.size() + a.polynomial.size());

    res.polynomial.insert(res.polynomial.end(), polynomial.begin(), polynomial.end());

    for (auto it : a.polynomial) {
        res.polynomial.emplace_back(-it.coeff, it.degree);
    }

    res.independent = independent + a.independent;

    res.simplify();
    
    return res;
}

PZsparse PZsparse::operator-(const double a) {
    assert(internalCheck());

    PZsparse res = *this;

    for (uint i = 0; i < NRows; i++) {
        for (uint j = 0; j < NCols; j++) {
            res.center(i,j) -= a;
        }
    }

    return res;
}

PZsparse operator-(const double a, const PZsparse& b) {
    assert(b.internalCheck());

    PZsparse res = b;

    for (uint i = 0; i < b.NRows; i++) {
        for (uint j = 0; j < b.NCols; j++) {
            res.center(i,j) -= a;
        }
    }

    return res;
}

PZsparse PZsparse::operator*(const PZsparse& a) {
    assert(internalCheck());
    assert(NCols == a.NRows || (NRows == 1 && NCols == 1) || (a.NRows == 1 && a.NCols == 1));

    PZsparse res;

    if (NRows == 1 && NCols == 1) {
        res.NRows = a.NRows;
        res.NCols = a.NCols;
    }
    else if (a.NRows == 1 && a.NCols == 1) {
        res.NRows = NRows;
        res.NCols = NCols;
    }
    else {
        res.NRows = NRows;
        res.NCols = a.NCols;
    }

    // center * center
    if (NRows == 1 && NCols == 1) {
        res.center = center(0) * a.center;
    }
    else if (a.NRows == 1 && a.NCols == 1) {
        res.center = center * a.center(0);
    }
    else {
        res.center = center * a.center;
    }

    res.polynomial.reserve(polynomial.size() + a.polynomial.size() + polynomial.size() * a.polynomial.size());
    // a.center * polynomial
    for (auto it : polynomial) {
        if (NRows == 1 && NCols == 1) {
            res.polynomial.emplace_back(it.coeff(0) * a.center, it.degree);
        }
        else if (a.NRows == 1 && a.NCols == 1) {
            res.polynomial.emplace_back(it.coeff * a.center(0), it.degree);
        }
        else {
            res.polynomial.emplace_back(it.coeff * a.center, it.degree);
        }
    }

    // center * a.polynomial
    for (auto it : a.polynomial) {
        if (NRows == 1 && NCols == 1) {
            res.polynomial.emplace_back(center(0) * it.coeff, it.degree);
        }
        else if (a.NRows == 1 && a.NCols == 1) {
            res.polynomial.emplace_back(center * it.coeff(0), it.degree);
        }
        else {
            res.polynomial.emplace_back(center * it.coeff, it.degree);
        }
    }

    // polynomial * a.polynomial (degree for each factor shouldn't be larger than 1)
    // Eigen::MatrixXd reduce_amount_1 = Eigen::MatrixXd::Zero(NRows, a.NCols);

    for (auto it1 : polynomial) {
        for (auto it2 : a.polynomial) {
            Eigen::MatrixXd multiply_coeff = it1.coeff * it2.coeff;

            if (NRows == 1 && NCols == 1) {
                multiply_coeff = it1.coeff(0) * it2.coeff;
            }
            else if (a.NRows == 1 && a.NCols == 1) {
                multiply_coeff = it1.coeff * it2.coeff(0);
            }
            else {
                multiply_coeff = it1.coeff * it2.coeff;
            }

            // Do not have to check carry
            // if we already know the maximum degree in the polynomial
            res.polynomial.emplace_back(multiply_coeff, it1.degree + it2.degree);
        }
    }

    // a.independent * (center + polynomial)
    Eigen::MatrixXd reduce_amount_2 = center.cwiseAbs();

    for (auto it : polynomial) {
        reduce_amount_2 += it.coeff.cwiseAbs();
    }

    if (NRows == 1 && NCols == 1) {
        reduce_amount_2 = reduce_amount_2(0) * a.independent;
    }
    else if (a.NRows == 1 && a.NCols == 1) {
        reduce_amount_2 *= a.independent(0);
    }
    else {
        reduce_amount_2 *= a.independent;
    }
    
    // independent * (a.center + a.polynomial)
    Eigen::MatrixXd reduce_amount_3 = a.center.cwiseAbs();

    for (auto it : a.polynomial) {
        reduce_amount_3 += it.coeff.cwiseAbs();
    }
    
    if (NRows == 1 && NCols == 1) {
        reduce_amount_3 = independent(0) * reduce_amount_3;
    }
    else if (a.NRows == 1 && a.NCols == 1) {
        reduce_amount_3 = independent * reduce_amount_3(0);
    }
    else {
        reduce_amount_3 = independent * reduce_amount_3;
    }

    // independent * a.independent + add reduced intervals
    Eigen::MatrixXd reduce_amount = reduce_amount_2 + reduce_amount_3;

    if (NRows == 1 && NCols == 1) {
        res.independent = independent(0) * a.independent + reduce_amount;
    }
    else if (a.NRows == 1 && a.NCols == 1) {
        res.independent = independent * a.independent(0) + reduce_amount;
    }
    else {
        res.independent = independent * a.independent + reduce_amount;
    }

    res.simplify();
    
    return res;
}

PZsparse PZsparse::operator*(const double a) {
    assert(internalCheck());

    PZsparse res(NRows, NCols);

    res.center = center * a;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(a * it.coeff, it.degree);
    }

    res.independent = independent * fabs(a);

    return res;
}

PZsparse operator*(const double a, const PZsparse& b) {
    assert(b.internalCheck());

    PZsparse res(b.NRows, b.NCols);

    res.center = b.center * a;

    res.polynomial.reserve(b.polynomial.size());

    for (auto it : b.polynomial) {
        res.polynomial.emplace_back(a * it.coeff, it.degree);
    }

    res.independent = b.independent * fabs(a);

    return res;
}

PZsparse PZsparse::operator/(const double a) {
    assert(internalCheck());

    PZsparse res(NRows, NCols);

    res.center = center / a;

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(it.coeff / a, it.degree);
    }

    res.independent = independent / fabs(a);

    return res;
}

PZsparse PZsparse::transpose() {
    assert(internalCheck());

    PZsparse res(NCols, NRows);

    res.center = center.transpose();

    res.polynomial.reserve(polynomial.size());

    for (auto it : polynomial) {
        res.polynomial.emplace_back(it.coeff.transpose(), it.degree);
    }

    res.independent = independent.transpose();

    return res;
}

void PZsparse::addOneDimPZ(const PZsparse& a, uint row_id, uint col_id) {
    assert(internalCheck());
    assert(a.NRows == 1 && a.NCols == 1);
    assert(row_id < NRows);
    assert(col_id < NCols);

    center(row_id, col_id) += a.center(0);

    for (auto it : a.polynomial) {
        Eigen::MatrixXd temp_coeff = Eigen::MatrixXd::Zero(NRows, NCols);
        temp_coeff(row_id, col_id) = it.coeff(0);
        polynomial.emplace_back(temp_coeff, it.degree);
    }

    independent(row_id, col_id) += a.independent(0);

    simplify();
}

PZsparse stack(const PZsparseArray& a) {
    assert(a.cols() == 1);

    for (uint i = 0; i < a.rows(); i++) {
        assert(a(i, 0).NRows == 1 && a(i, 0).NCols == 1);
    }

    PZsparse res(a.rows(), 1);

    for (uint i = 0; i < a.rows(); i++) {
        res.center(i, 0) = a(i, 0).center(0);
    }

    res.polynomial.reserve(a.rows() * a(0, 0).polynomial.size());
    for (uint i = 0; i < a.rows(); i++) {
        for (auto it : a(i, 0).polynomial) {
            Eigen::MatrixXd temp_coeff = Eigen::MatrixXd::Zero(a.rows(), 1);
            temp_coeff(i) = it.coeff(0);
            res.polynomial.emplace_back(temp_coeff, it.degree);
        }
    }

    for (uint i = 0; i < a.rows(); i++) {
        res.independent(i, 0) = a(i, 0).independent(0);
    }

    res.simplify();

    return res;
}

PZsparse cross(const Eigen::MatrixXd& a, const PZsparse& b) {
    assert(a.rows() == 3 && a.cols() == 1 && b.NRows == 3 && b.NCols == 1);

    PZsparseArray res(3, 1);

    PZsparse b0 = b(0, 0);
    PZsparse b1 = b(1, 0);
    PZsparse b2 = b(2, 0);

    res(0, 0) = a(1, 0) * b2 - a(2, 0) * b1;
    res(1, 0) = a(2, 0) * b0 - a(0, 0) * b2;
    res(2, 0) = a(0, 0) * b1 - a(1, 0) * b0;

    return stack(res);
}

PZsparse cross(const PZsparse& a, const PZsparse& b) {
    assert(a.NRows == 3 && a.NCols == 1 && b.NRows == 3 && b.NCols == 1);

    PZsparseArray res(3, 1);

    PZsparse a0 = a(0, 0);
    PZsparse a1 = a(1, 0);
    PZsparse a2 = a(2, 0);
    PZsparse b0 = b(0, 0);
    PZsparse b1 = b(1, 0);
    PZsparse b2 = b(2, 0);

    res(0, 0) = a1 * b2 - a2 * b1;
    res(1, 0) = a2 * b0 - a0 * b2;
    res(2, 0) = a0 * b1 - a1 * b0;

    return stack(res);
}

PZsparse cross(const PZsparse& a, const Eigen::MatrixXd& b) {
    assert(a.NRows == 3 && a.NCols == 1 && b.rows() == 3 && b.cols() == 1);

    PZsparseArray res(3, 1);

    PZsparse a0 = a(0, 0);
    PZsparse a1 = a(1, 0);
    PZsparse a2 = a(2, 0);

    res(0, 0) = a1 * b(2, 0) - a2 * b(1, 0);
    res(1, 0) = a2 * b(0, 0) - a0 * b(2, 0);
    res(2, 0) = a0 * b(1, 0) - a1 * b(0, 0);

    return stack(res);
}

#endif