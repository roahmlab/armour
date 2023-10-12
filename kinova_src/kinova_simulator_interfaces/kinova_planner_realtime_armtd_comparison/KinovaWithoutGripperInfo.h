#ifndef KINOVA_INFO_H
#define KINOVA_INFO_H

#include "Headers.h"

// Hard code robot physical properties here
// Kinova 7 dof

// number of joints
#define NUM_JOINTS 7

// number of factors (variables) in polynomial / number of actuated motors
// The last one joint (end effector) of Kinova is fixed
#define NUM_FACTORS 7

// 1,2,3 -> x,y,z, negative sign means rotate in reverse direction, 0 means fixed joint
const int axes[NUM_JOINTS] = { 3,3,3,3,3,3,3 }; 

// joint position translation element w.r.t previous joint frame, same as xyz in urdf
const double trans[(NUM_JOINTS + 1) * 3] = {  0,          0,     0.15643,
											0,   0.005375,    -0.12838,
											0,   -0.21038,   -0.006375,
											0,   0.006375,    -0.21038,
											0,   -0.20843,   -0.006375,
											0, 0.00017505,    -0.10593,
											0,   -0.10593, -0.00017505,
											0,          0,           0 };


// joint position rotation element w.r.t previous joint frame, same as rpy in urdf
const double rots[NUM_JOINTS * 3] = { M_PI,       0, 0,
                                    M_PI * 0.5, 0, 0,
								   -M_PI * 0.5, 0, 0,
								    M_PI * 0.5, 0, 0,
								   -M_PI * 0.5, 0, 0,
								    M_PI * 0.5, 0, 0,
								   -M_PI * 0.5, 0, 0 };

// link mass
const double mass[NUM_JOINTS] = { 1.3773, 1.1636, 1.1636, 0.9302, 0.6781, 0.6781, 0.5 }; 
const double mass_uncertainty = 0.03;

// link center of mass
const double com[NUM_JOINTS * 3] = {-0.000023, -0.010364,  -0.07336,
									-0.000044,  -0.09958, -0.013278,
									-0.000044, -0.006641, -0.117892,
									-0.000018, -0.075478, -0.015006,
									 0.000001, -0.009432, -0.063883,
									 0.000001, -0.045483,  -0.00965,
									 0.000281,  0.011402, -0.029798 };
const double com_uncertainty = 0.0;

// link inertia
const double inertia[NUM_JOINTS * 9] = {  0.00457, 0.000001, 0.000002, 0.000001, 0.004831, 0.000448, 0.000002, 0.000448, 0.001409,
										0.011088, 0.000005, 0, 0.000005, 0.001072, -0.000691, 0, -0.000691, 0.011255,
										0.010932, 0, -0.000007, 0, 0.011127, 0.000606, -0.000007, 0.000606, 0.001043,
										0.008147, -0.000001, 0, -0.000001, 0.000631, -0.0005, 0, -0.0005, 0.008316,
										0.001596, 0, 0, 0, 0.001607, 0.000256, 0, 0.000256, 0.000399,
										0.001641, 0, 0, 0, 0.00041, -0.000278, 0, -0.000278, 0.001641,
										0.000587, 0.000003, 0.000003, 0.000003, 0.000369, -0.000118, 0.000003, -0.000118, 0.000609 };
const double inertia_uncertainty = 0.03;

// joint friction
// const double friction[NUM_JOINTS] = {0.5217383101288284, 0.5769579059927288, 0.4213397946418778, 0.4945515376566732, 0.1611070502661354, 0.1333491185514130, 0.1434440181717370};
const double friction[NUM_JOINTS] = {0.0}; // disable friction in Matlab simulation for now

// joint damping
// const double damping[NUM_JOINTS] = {10.5, 7.4064845817230722, 9.9727633408172860, 8.2667950822503915, 8.8572249026528151, 8.7110831569332845, 8.8881903638306934};
const double damping[NUM_JOINTS] = {0.0}; // disable damping in Matlab simulation for now

// joint armature / motor transmission inertia
const double armature[NUM_JOINTS] = {8.03, 11.9962024615303644, 9.0025427861751517, 11.5806439316706360, 8.4665040917914123, 8.8537069373742430, 8.8587303664685315};

// other robot info
// 1000.0 means they are continuous joints and have no limits
const double state_limits_lb[NUM_FACTORS] = { -1000.0,  -2.41,  -1000.0,  -2.66,  -1000.0,  -2.23,  -1000.0 }; // rad
const double state_limits_ub[NUM_FACTORS] = {  1000.0,   2.41,   1000.0,   2.66,   1000.0,   2.23,   1000.0 }; // rad

const double speed_limits[NUM_FACTORS] = { 1.3963, 1.3963, 1.3963, 1.3963, 1.2218, 1.2218, 1.2218 }; // rad/s

const double torque_limits[NUM_FACTORS] = { 56.7, 56.7, 56.7, 56.7, 29.4, 29.4, 29.4 }; // N*m These are WARNING torque limits

const double gravity = 9.81;

// link zonotope information
const double link_zonotope_center[NUM_JOINTS][3] = {{ 0.000000, -0.001297, -0.088375 },
													{ 0.000000, -0.089400, -0.007877 },
													{ 0.000000, -0.001502, -0.129375 },
													{ 0.000000, -0.087450, -0.013648 },
													{ 0.000001, -0.009023, -0.071752 },
													{ 0.000000, -0.041661, -0.009251 },
													{ 0.000000, -0.018585, -0.033462 }};

const double link_zonotope_generators[NUM_JOINTS][3] = {{ 0.046358, 0.047354, 0.086000 },
														{ 0.046000, 0.135400, 0.047501 },
														{ 0.046000, 0.047501, 0.127000 },
														{ 0.046000, 0.133450, 0.042293 },
														{ 0.034999, 0.044023, 0.069252 },
														{ 0.035000, 0.076739, 0.044076 },
														{ 0.045500, 0.056085, 0.030963 }};

// ultimate bound
const double alpha = 10.0;
const double V_m = 1e-2;
const double M_max = 15.79635774;
const double M_min = 5.095620491878957;
const double eps = sqrt(2 * V_m / M_min);
const double K = 5.0;
const double qe = eps / K;
const double qde = 2 * eps;
const double qdae = eps;
const double qddae = 2 * K * eps;

#endif
