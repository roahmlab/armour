# ARMOUR Online Planning Implementation In C++ & CUDA

## Prerequisites
* [CUDA](https://developer.nvidia.com/cuda-downloads)
* [Ipopt](https://coin-or.github.io/Ipopt/INSTALL.html)
* [libboost-dev](https://www.boost.org/)

## Getting Started
1. Define your robot physical properties in a header file "Robot"Info.h. Check KinovaInfo.h and FetchInfo.h as examples (will include a matlab or python script that automatically reads the urdf and prints the data to file).

2. In Parameters.h, include your "Robot"Info.h.

3. Note that you have to configure the Ipopt path if it is installed in a different path on your computer! (introduce some other fancy ways to compile code)

4. In 'armour_main.cpp', change the variable 'pathname' to 'WHERE-YOU-INSTALL-ARMTD/armtd-dev/cuda-dev/PZsparse-Bernstein/results/'.

5. Run './compile.sh' in the terminal to compile the C++ program.

## Files
### armour_main.cpp
The main program that computes all polynomial zonotopes (forward kinematics & input constraints) and solve the optimization problem using Ipopt.

### PZsparse.h
This header file defines a polynomial zonotope class. It also defines a 3x1 vector polynomial zonotope class and a 3x3 matrix polynomial zonotope class, together with all necessary arithmetic.

### Trajectory.h
This header file should define a trajectory class that can compute the desired trajectory given parameter _k_ and time _t_. It should also define the function that computes a polynomial zonotope that over-approximates the desired trajectory over a given time interval.

The header file in this repo right now gives a desired trajectory represented using a 5-degree Bezier curve.

### Dynamics.h
This header file defines a class that helps compute the polynomial zonotope for forward kinematics and RNEA (recursive newton euler algorithm).

### CollisionChecking.h
This header file defines a 3D obstacle class that can perform collision checking between obstacles and any point in gpu.

### NLPclass.h
This header file defines an NLP class required by Ipopt. 

## Parameters
All hyperparameters are defined in Parameters.h and all robot physical properties should be defined in "Robot"Info.h

* NUM_JOINTS: number of joints (including fixed joints).

* NUM_FACTORS: number of actuated&parameterized joints, Note that we assume that all actuated joints are sequentially located at the beginning of the kinematics chain.

* SIMPLIFY_THRESHOLD: If the coefficient of a monomial in the polynomial zonotope is smaller than this threshold, then it will be reduced and over-approximated into the independent interval of the polynomial zonotope.

* NUM TIME STEPS: number of time intervals split in the planning horizon (This must be an even number for the current 5-degree Bezier curve class!)

* k range[NUM_FACTORS]: The range of the destination away from the initial position at the end of the desired trajectory (check ARMOUR paper)

* MAX OBSTACLE NUM: maximum number of obstacles for memory preallocation (number of obstacles should NEVER be larger than this!)

* MAX OBSTACLE GENERATOR NUM: number of generators for all obstacle zonotope

* NUM THREADS: number of threads to parallel in cpu for polynomial zonotope computation (note this should NEVER be larger than the number of processors on your computer)

* COLLISION AVOIDANCE CONSTRAINT VIOLATION THRESHOLD: threshold for collision avoidance constraint considered to be violated (unit: meter)

* TORQUE INPUT CONSTRAINT VIOLATION THRESHOLD: threshold for input constraint considered to be violated (unit: Newton * meter)

* IPOPT_OPTIMIZATION_TOLERANCE: Ipopt option. Click [here](https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_tol) for more info

* IPOPT_MAX_CPU_TIME: Ipopt option. Click [here](https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_max_cpu_time) for more info

* IPOPT_PRINT_LEVEL: Ipopt option. Click [here](https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_print_level) for more info

* IPOPT_MU_STRATEGY: Ipopt option. Click [here](https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_mu_strategy) for more info

* IPOPT_LINEAR_SOLVER: Ipopt option. Click [here](https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_linear_solver) for more info

(You can add more Ipopt options in armour_main.cpp)

## Usage
Interface with Matlab is implemented by calling a compiled C++ program in Matlab through function 'system'.
It allows users to run operating system commands and get the output in Matlab.
In that way, Matlab will not be affected even if there is error in the C++ program.
The communication between Matlab and the C++ program is initiated through reading and writing different files since the communication in Matlab simulation does not require synchronization.
Since the Matlab scripts may call the C++ program in different paths, we use absolute path for all the input and output files here.
Make sure you change the "pathname" defined in armour_main.cpp when you run the Matlab tests in a different computer.

### Input
Input is stored as a text file in 'armtd-dev/cuda-dev/PZsparse-Bernstein/results/armour_main.in'

* Initial joint position: NUM_FACTORS floating numbers

* Initial joint velocity: NUM_FACTORS floating numbers

* Initial joint acceleration: NUM_FACTORS floating numbers

* Desired joint position: NUM_FACTORS floating numbers

* Number of obstacles: 1 integer number

* Obstacle zonotopes: (Number of obstacles) * (MAX OBSTACLE GENERATOR NUM + 1) * 3 floating numbers

* t_plan: 1 floating number range from 0 to 1 (unit: second), minimize the distance between the desired joint position and the joint trajectory at t_plan. For example, in the middle of the planning, this should be 0.5 in order to achieve receding horizon planning. At the end of the planning (the last HLP waypoint has been provided), this should be 1 because we want the robot to stop at the final goal position.

### Output
Output is stored as multiple text files in 'armtd-dev/cuda-dev/PZsparse-Bernstein/results/armour_main_xxx.out'

* armour_main.out stores k_opt, which contains NUM_FACTORS floating numbers

* armour_main_joint_position_center.out stores the center of all joint FRS in all time intervals, when tracking desired trajectories that corresponds to k_opt

* armour_main_joint_position_radius.out stores the radius of all joint FRS in all time intervals, when tracking desired trajectories that corresponds to k_opt

* armour_main_control_input_radius.out stores the radius of the control input PZ

* armour_main_constraints.out stores all constraint values that corresponds to k_opt (the first NUM TIME STEPS * NUM FACTORS entries are just the center of the control input PZ, so together with armour_main_control_input_radius.out you can get the control input PZ corresponds to the optimal parameter that ipopt finds)

