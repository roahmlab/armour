#ifndef COLLISION_CHECKING_H
#define COLLISION_CHECKING_H

#include "Trajectory.h"

#define BUFFER_OBSTACLE_GENERATOR_NUM (MAX_OBSTACLE_GENERATOR_NUM + 6)
#define COMB_NUM BUFFER_OBSTACLE_GENERATOR_NUM * (BUFFER_OBSTACLE_GENERATOR_NUM - 1) / 2

__constant__ double dev_obstacles[MAX_OBSTACLE_NUM * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3];

__constant__ uint dev_combA[COMB_NUM];
__constant__ uint dev_combB[COMB_NUM];

class Obstacles {
public:
    int num_obstacles = 0;
    const double* obstacles = nullptr;

    Eigen::Vector3d* dev_A = nullptr;
    double* d = nullptr;
    double* dev_d = nullptr;
    double* delta = nullptr;
    double* dev_delta = nullptr;

	Eigen::Vector3d* dev_buffered_c = nullptr;
	Eigen::Matrix<double, 3, BUFFER_OBSTACLE_GENERATOR_NUM>* dev_buffered_G = nullptr;

    Eigen::Matrix<double, 3, 3 + 3>* dev_link_independent_generators = nullptr;

    Eigen::Vector3d* dev_link_sliced_center = nullptr;
    Eigen::Vector3d* dev_dk_link_sliced_center = nullptr;

    double* dev_link_c = nullptr;
    double* dev_grad_link_c = nullptr;

	Obstacles(const double* obstacles_inp, const int num_obstacles_inp);

	~Obstacles();

	void initializeHyperPlane(const Eigen::Matrix<double, 3, 3 + 3>* link_independent_generators);

	void linkFRSConstraints(Eigen::Vector3d* link_sliced_center,
                            Eigen::Vector3d* dk_link_sliced_center,
							double* link_c,
							double* grad_link_c = nullptr);
};

__global__ void bufferObstaclesKernel(const Eigen::Matrix<double, 3, 3 + 3>* link_independent_generators, 
                                      Eigen::Vector3d* buffered_c, 
                                      Eigen::Matrix<double, 3, BUFFER_OBSTACLE_GENERATOR_NUM>* buffered_G, 
                                      int link_id);

__global__ void polytope_PH(Eigen::Vector3d* buffered_c, 
                            Eigen::Matrix<double, 3, BUFFER_OBSTACLE_GENERATOR_NUM>* buffered_G, 
                            Eigen::Vector3d* A, 
                            double* d, 
                            double* delta,
                            int link_id);

__global__ void checkCollisionKernel(Eigen::Vector3d* A, 
                                     double* d, 
                                     double* delta, 
									 Eigen::Vector3d* link_sliced_center,
									 Eigen::Vector3d* dk_link_sliced_center,
									 int link_id, 
                                     double* link_c, 
                                     double* grad_link_c);

#endif