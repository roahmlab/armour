#ifndef COLLISION_CHECKING_CU
#define COLLISION_CHECKING_CU

#include "CollisionChecking.h"

Obstacles::Obstacles(const double* obstacles_inp, const int num_obstacles_inp) {
    obstacles = obstacles_inp;
    num_obstacles = num_obstacles_inp;

    if (num_obstacles > MAX_OBSTACLE_NUM) {
        WARNING_PRINT("Number of obstacles larger than MAX_OBSTACLE_NUM !\n");
        throw;
    }

    // build pre-defined hyper-plane equations for collision checking between links and obstacles
    if (obstacles != nullptr && num_obstacles > 0) {
        cudaMalloc((void**)&dev_A, NUM_JOINTS * NUM_TIME_STEPS * MAX_OBSTACLE_NUM * COMB_NUM * sizeof(Eigen::Vector3d));
        cudaMalloc((void**)&dev_d, NUM_JOINTS * NUM_TIME_STEPS * MAX_OBSTACLE_NUM * COMB_NUM * sizeof(double));
		cudaMalloc((void**)&dev_delta, NUM_JOINTS * NUM_TIME_STEPS * MAX_OBSTACLE_NUM * COMB_NUM * sizeof(double));

        // initialize constant memory
        // obstacle data
        cudaMemcpyToSymbol(dev_obstacles, obstacles, num_obstacles * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3 * sizeof(double));

        // combination index
        uint combA[COMB_NUM], combB[COMB_NUM];
        uint a_id = 0, b_id = 1;
        for (uint i = 0; i < COMB_NUM; i++) {
            combA[i] = a_id;
            combB[i] = b_id;

            if (b_id < BUFFER_OBSTACLE_GENERATOR_NUM - 1) {
                b_id++;
            }
            else {
                a_id++;
                b_id = a_id + 1;
            }
        }

        cudaMemcpyToSymbol(dev_combA, combA, COMB_NUM * sizeof(uint));
        cudaMemcpyToSymbol(dev_combB, combB, COMB_NUM * sizeof(uint));

        cudaMalloc((void**)&dev_buffered_c, NUM_TIME_STEPS * MAX_OBSTACLE_NUM * sizeof(Eigen::Vector3d));
        cudaMalloc((void**)&dev_buffered_G, NUM_TIME_STEPS * MAX_OBSTACLE_NUM * sizeof(Eigen::Matrix<double, 3, BUFFER_OBSTACLE_GENERATOR_NUM>));

		cudaMalloc((void**)&dev_link_independent_generators, NUM_TIME_STEPS * NUM_JOINTS * sizeof(Eigen::Matrix<double, 3, 3 + 3>));

        cudaMalloc((void**)&dev_link_sliced_center, NUM_TIME_STEPS * NUM_JOINTS * sizeof(Eigen::Vector3d));
        cudaMalloc((void**)&dev_dk_link_sliced_center, NUM_TIME_STEPS * NUM_JOINTS * NUM_FACTORS * sizeof(Eigen::Vector3d));

        cudaMalloc((void**)&dev_link_c, NUM_TIME_STEPS * MAX_OBSTACLE_NUM * sizeof(double));
        cudaMalloc((void**)&dev_grad_link_c, NUM_TIME_STEPS * MAX_OBSTACLE_NUM * NUM_FACTORS * sizeof(double));
    }
}

Obstacles::~Obstacles() {
    cudaFree(dev_A);
    cudaFree(dev_d);
    cudaFree(dev_delta);

    cudaFree(dev_buffered_c);
    cudaFree(dev_buffered_G);

	cudaFree(dev_link_independent_generators);

    cudaFree(dev_link_sliced_center);
    cudaFree(dev_dk_link_sliced_center);

    cudaFree(dev_link_c);
    cudaFree(dev_grad_link_c);
}

void Obstacles::initializeHyperPlane(const Eigen::Matrix<double, 3, 3 + 3>* link_independent_generators) {
    if (num_obstacles == 0) return;

    cudaMemcpy(dev_link_independent_generators, link_independent_generators, NUM_TIME_STEPS * NUM_JOINTS * sizeof(Eigen::Matrix<double, 3, 3 + 3>), cudaMemcpyHostToDevice);

    for (uint i = 0; i < NUM_JOINTS; i++) {
        dim3 block1(3, num_obstacles);
        bufferObstaclesKernel << < NUM_TIME_STEPS, block1 >> > (dev_link_independent_generators, dev_buffered_c, dev_buffered_G, i);

        dim3 grid1(NUM_TIME_STEPS, num_obstacles);
        polytope_PH << < grid1, COMB_NUM >> > (dev_buffered_c, dev_buffered_G, dev_A, dev_d, dev_delta, i);
    }

    cudaDeviceSynchronize();
}

void Obstacles::linkFRSConstraints(Eigen::Vector3d* link_sliced_center,
                                   Eigen::Vector3d* dk_link_sliced_center,
                                   double* link_c,
                                   double* grad_link_c) {
    if (num_obstacles == 0) return;
                                    
    bool ifComputeGradient = true;

    cudaMemcpy(dev_link_sliced_center, link_sliced_center, NUM_TIME_STEPS * NUM_JOINTS * sizeof(Eigen::Vector3d), cudaMemcpyHostToDevice);

    if (dk_link_sliced_center == nullptr || grad_link_c == nullptr) {
        ifComputeGradient = false;
    }
    else {
        cudaMemcpy(dev_dk_link_sliced_center, dk_link_sliced_center, NUM_TIME_STEPS * NUM_JOINTS * NUM_FACTORS * sizeof(Eigen::Vector3d), cudaMemcpyHostToDevice);
    }

    for (uint i = 0; i < NUM_JOINTS; i++) {
        dim3 grid1(NUM_TIME_STEPS, num_obstacles);
        
        if (ifComputeGradient) {
            checkCollisionKernel << < grid1, COMB_NUM >> > (dev_A, dev_d, dev_delta,
                                                            dev_link_sliced_center,
                                                            dev_dk_link_sliced_center,
                                                            i, 
                                                            dev_link_c, 
															dev_grad_link_c);
        }
        else {
            checkCollisionKernel << < grid1, COMB_NUM >> > (dev_A, dev_d, dev_delta,
                                                            dev_link_sliced_center,
                                                            nullptr,
                                                            i, 
                                                            dev_link_c,
															nullptr);
        }
            
        if (link_c != nullptr) {
            cudaMemcpy(link_c + i * NUM_TIME_STEPS * num_obstacles, dev_link_c, NUM_TIME_STEPS * num_obstacles * sizeof(double), cudaMemcpyDeviceToHost);
        }
        if (grad_link_c != nullptr) {
            cudaMemcpy(grad_link_c + i * NUM_TIME_STEPS * num_obstacles * NUM_FACTORS, dev_grad_link_c, NUM_TIME_STEPS * num_obstacles * NUM_FACTORS * sizeof(double), cudaMemcpyDeviceToHost);
        }
    }
}

__global__ void bufferObstaclesKernel(const Eigen::Matrix<double, 3, 3 + 3>* link_independent_generators, 
                                      Eigen::Vector3d* buffered_c, 
                                      Eigen::Matrix<double, 3, BUFFER_OBSTACLE_GENERATOR_NUM>* buffered_G, 
                                      int link_id) {
	uint time_id = blockIdx.x;
	uint obs_id = threadIdx.y;
	uint p_id = threadIdx.x;
	uint num_obstacles = blockDim.y;

	__shared__ double independent_generators[3][3 + 3]; // equivalent to an Eigen matrix

	if (obs_id == 0) {
		for (uint i = 0; i < 3 + 3; i++) {
			independent_generators[p_id][i] = link_independent_generators[time_id * NUM_JOINTS + link_id](p_id, i);
		}	
	}

	__syncthreads();

	// copy obstacle center to buffered obstacle center
	buffered_c[time_id * num_obstacles + obs_id](p_id) = dev_obstacles[obs_id * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3 + p_id];

	// copy obstacle generator to buffered obstacle generator
	for (uint i = 0; i < MAX_OBSTACLE_GENERATOR_NUM; i++) {
		buffered_G[time_id * num_obstacles + obs_id](p_id, i) = dev_obstacles[(obs_id * (MAX_OBSTACLE_GENERATOR_NUM + 1) + i + 1) * 3 + p_id];
	}

	// copy link independent generators to buffered obstacle generator
	for (uint i = 0; i < 3 + 3; i++) {
		buffered_G[time_id * num_obstacles + obs_id](p_id, i + MAX_OBSTACLE_GENERATOR_NUM) = independent_generators[p_id][i];
	}
}

__global__ void polytope_PH(Eigen::Vector3d* buffered_c, 
                            Eigen::Matrix<double, 3, BUFFER_OBSTACLE_GENERATOR_NUM>* buffered_G, 
                            Eigen::Vector3d* A, 
                            double* d, 
                            double* delta,
							int link_id) {
	uint time_id = blockIdx.x;
	uint obs_id = blockIdx.y;
	uint num_obstacles = gridDim.y;
	uint p_id = threadIdx.x;

	uint a_id = dev_combA[p_id];
	uint b_id = dev_combB[p_id];

	__shared__ double G[BUFFER_OBSTACLE_GENERATOR_NUM][3];
	__shared__ double c[3];

	if (p_id < BUFFER_OBSTACLE_GENERATOR_NUM) {
		G[p_id][0] = buffered_G[time_id * num_obstacles + obs_id](0, p_id);
		G[p_id][1] = buffered_G[time_id * num_obstacles + obs_id](1, p_id);
		G[p_id][2] = buffered_G[time_id * num_obstacles + obs_id](2, p_id);
	}

	if (p_id < 3) {
		c[p_id] = buffered_c[time_id * num_obstacles + obs_id](p_id);
	}

	__syncthreads();

	Eigen::Vector3d generator_cross;
	generator_cross(0) = G[a_id][1] * G[b_id][2] - G[a_id][2] * G[b_id][1];
    generator_cross(1) = G[a_id][2] * G[b_id][0] - G[a_id][0] * G[b_id][2];
    generator_cross(2) = G[a_id][0] * G[b_id][1] - G[a_id][1] * G[b_id][0];

	double norm_cross = generator_cross.norm();

	Eigen::Vector3d C;

	if (norm_cross > 0) {
		C = generator_cross / norm_cross;
	}
	else {
		C.setZero();
	}

	// A = C
	A[(((time_id * NUM_JOINTS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id)] = C;

	// d = C * c
	d[((time_id * NUM_JOINTS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id] = C[0] * c[0] + C[1] * c[1] + C[2] * c[2];

	// delta = sum(abs(C * G))
	double delta_res = 0.0;

	for (uint j = 0; j < BUFFER_OBSTACLE_GENERATOR_NUM; j++) {
		delta_res += fabs(C[0] * G[j][0] + C[1] * G[j][1] + C[2] * G[j][2]);
	}

	delta[((time_id * NUM_JOINTS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id] = delta_res;
}

__global__ void checkCollisionKernel(Eigen::Vector3d* A, 
                            		 double* d, 
                            		 double* delta, 
									 Eigen::Vector3d* link_sliced_center,
									 Eigen::Vector3d* dk_link_sliced_center,
									 int link_id, 
                                     double* link_c, 
									 double* grad_link_c) {
	uint time_id = blockIdx.x;
	uint obs_id = blockIdx.y;
	uint num_obstacles = gridDim.y;
	uint p_id = threadIdx.x;

	__shared__ double pos_res[COMB_NUM];
	__shared__ double neg_res[COMB_NUM];

	Eigen::Vector3d A_elt = A[(((time_id * NUM_JOINTS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id)];
	double d_elt = d[(((time_id * NUM_JOINTS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id)];
	double delta_elt = delta[(((time_id * NUM_JOINTS + link_id) * num_obstacles + obs_id) * COMB_NUM + p_id)];

	Eigen::Vector3d center_elt = link_sliced_center[time_id * NUM_JOINTS + link_id];

	if (A_elt.norm() > 0) {
		pos_res[p_id] = A_elt.dot(center_elt) - (d_elt + delta_elt);
		neg_res[p_id] = -A_elt.dot(center_elt) - (-d_elt + delta_elt);
	}
	else {
		pos_res[p_id] = -100000000;
		neg_res[p_id] = -100000000;
	}

	__syncthreads();

	if (p_id == 0) {
		double max_elt = -100000000;
		uint max_id = 0;
		bool pos_or_neg = false; // false -> pos, true -> neg

		for (uint i = 0; i < COMB_NUM; i++) {
			if (pos_res[i] > max_elt) {
				max_elt = pos_res[i];
				max_id = i;
				pos_or_neg = false;
			}

			if (neg_res[i] > max_elt) {
				max_elt = neg_res[i];
				max_id = i;
				pos_or_neg = true;
			}
		}

		if (link_c != nullptr) {
			link_c[time_id * num_obstacles + obs_id] = -max_elt;
		}

		if (dk_link_sliced_center != nullptr && grad_link_c != nullptr) {
			Eigen::Vector3d max_A_elt = A[(((time_id * NUM_JOINTS + link_id) * num_obstacles + obs_id) * COMB_NUM + max_id)];

			for (uint k = 0; k < NUM_FACTORS; k++) {
				if (!pos_or_neg) {
					grad_link_c[(time_id * num_obstacles + obs_id) * NUM_FACTORS + k] = -max_A_elt.dot(dk_link_sliced_center[(time_id * NUM_JOINTS + link_id) * NUM_FACTORS + k]);
				}
				else {
					grad_link_c[(time_id * num_obstacles + obs_id) * NUM_FACTORS + k] = max_A_elt.dot(dk_link_sliced_center[(time_id * NUM_JOINTS + link_id) * NUM_FACTORS + k]);
				}
			}
		}
	}
}

#endif