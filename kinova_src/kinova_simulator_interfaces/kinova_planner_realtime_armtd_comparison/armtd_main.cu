#include "NLPclass.h"
#include "BufferPath.h"

const std::string inputfilename = pathname + "armtd.in";
const std::string outputfilename1 = pathname + "armtd.out";
const std::string outputfilename2 = pathname + "armtd_joint_position_center.out";
const std::string outputfilename3 = pathname + "armtd_joint_position_radius.out";
const std::string outputfilename4 = pathname + "armtd_constraints.out";

int main() {
/*
Section I:
    Parse input
    There is no check and warning, so be careful!
*/
    // Here is an example of required input
    // double q0[NUM_FACTORS] = {0.6543, -0.0876, -0.4837, -1.2278, -1.5735, -1.0720, 0};
    // double qd0[NUM_FACTORS] = {0, 0, 0, 0, 0, 0, 0};
    // double qdd0[NUM_FACTORS] = {0, 0, 0, 0, 0, 0, 0};
    // double q_des[NUM_FACTORS] = {0.6831, 0.009488, -0.2471, -0.9777, -1.414, -0.9958, 0};

    // const int num_obstacles = 10;
    // const double obstacles[num_obstacles * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3] = {-0.28239,  -0.33281, 0.88069, 0.069825, 0, 0, 0,  0.09508, 0, 0, 0, 0.016624,
    //                                                                             -0.19033,  0.035391,  1.3032,  0.11024, 0, 0, 0, 0.025188, 0, 0, 0, 0.014342,
    //                                                                             0.67593, -0.085841, 0.43572,  0.17408, 0, 0, 0,  0.07951, 0, 0, 0,  0.18012,
    //                                                                             0.75382,   0.51895,  0.4731, 0.030969, 0, 0, 0,  0.22312, 0, 0, 0,  0.22981,
    //                                                                             0.75382,   0.51895,  0.4731, 0.030969, 0, 0, 0,  0.22312, 0, 0, 0,  0.22981,
    //                                                                             -0.28239,  -0.33281, 0.88069, 0.069825, 0, 0, 0,  0.09508, 0, 0, 0, 0.016624,
    //                                                                             -0.19033,  0.035391,  1.3032,  0.11024, 0, 0, 0, 0.025188, 0, 0, 0, 0.014342,
    //                                                                             0.67593, -0.085841, 0.43572,  0.17408, 0, 0, 0,  0.07951, 0, 0, 0,  0.18012,
    //                                                                             0.75382,   0.51895,  0.4731, 0.030969, 0, 0, 0,  0.22312, 0, 0, 0,  0.22981,
    //                                                                             0.75382,   0.51895,  0.4731, 0.030969, 0, 0, 0,  0.22312, 0, 0, 0,  0.22981};

    // declare this first and make sure we always have a new output
    std::ofstream outputstream1(outputfilename1);

    double q0[NUM_FACTORS] = {0.0};
    double qd0[NUM_FACTORS] = {0.0};
    double q_des[NUM_FACTORS] = {0.0};

    double c_cos_q_des[NUM_FACTORS * NUM_TIME_STEPS];
    double g_cos_q_des[NUM_FACTORS * NUM_TIME_STEPS];
    double r_cos_q_des[NUM_FACTORS * NUM_TIME_STEPS];

    double c_sin_q_des[NUM_FACTORS * NUM_TIME_STEPS];
    double g_sin_q_des[NUM_FACTORS * NUM_TIME_STEPS];
    double r_sin_q_des[NUM_FACTORS * NUM_TIME_STEPS];

    double k_range[NUM_FACTORS];

    int num_obstacles = 0;
    double obstacles[MAX_OBSTACLE_NUM * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3] = {0.0};

    std::ifstream inputstream(inputfilename);
    if (!inputstream.is_open()) {
        WARNING_PRINT("        CUDA & C++: Error reading input files !\n");
        outputstream1 << -1;
        outputstream1.close();
        throw;
    }
    for (int i = 0; i < NUM_FACTORS; i++) {
        inputstream >> q0[i];
    }
    for (int i = 0; i < NUM_FACTORS; i++) {
        inputstream >> qd0[i];
    }
    for (int i = 0; i < NUM_FACTORS; i++) {
        inputstream >> q_des[i];
    }
    for (int i = 0; i < NUM_FACTORS; i++) {
        for (int j = 0; j < NUM_TIME_STEPS; j++) {
            inputstream >> c_cos_q_des[i * NUM_TIME_STEPS + j];
        }
        for (int j = 0; j < NUM_TIME_STEPS; j++) {
            inputstream >> g_cos_q_des[i * NUM_TIME_STEPS + j];
        }
        for (int j = 0; j < NUM_TIME_STEPS; j++) {
            inputstream >> r_cos_q_des[i * NUM_TIME_STEPS + j];
        }
        for (int j = 0; j < NUM_TIME_STEPS; j++) {
            inputstream >> c_sin_q_des[i * NUM_TIME_STEPS + j];
        }
        for (int j = 0; j < NUM_TIME_STEPS; j++) {
            inputstream >> g_sin_q_des[i * NUM_TIME_STEPS + j];
        }
        for (int j = 0; j < NUM_TIME_STEPS; j++) {
            inputstream >> r_sin_q_des[i * NUM_TIME_STEPS + j];
        }
        inputstream >> k_range[i];
    }
    inputstream >> num_obstacles;
    if (num_obstacles > MAX_OBSTACLE_NUM || num_obstacles < 0) {
        WARNING_PRINT("Number of obstacles larger than MAX_OBSTACLE_NUM !\n");
        outputstream1 << -1;
        outputstream1.close();
        throw;
    }
    if (num_obstacles > 0) {
        for (int i = 0; i < num_obstacles * (MAX_OBSTACLE_GENERATOR_NUM + 1) * 3; i++) {
            inputstream >> obstacles[i];
        }
    }

    inputstream.close();
     
/*
Section II:
    Initialize all polynomial zonotopes
*/
    PZsparse p[NUM_TIME_STEPS * NUM_FACTORS * 3];
    double jointPositionRadius[NUM_TIME_STEPS * NUM_FACTORS * 3] = {0.0};

    ConstantAccelerationCurve traj(q0, qd0, 
                                   c_cos_q_des, g_cos_q_des, r_cos_q_des,
                                   c_sin_q_des, g_sin_q_des, r_sin_q_des,
                                   k_range);

    Obstacles O(obstacles, num_obstacles);  

    KinematicsDynamics kd(&traj);
    Eigen::Matrix<double, 3, 3 + 3> link_independent_generators[NUM_TIME_STEPS * NUM_JOINTS];

    // do not count time for memory allocation
    auto start1 = std::chrono::high_resolution_clock::now();                                                         

    omp_set_num_threads(NUM_THREADS);

    int openmp_t_ind; // loop index

    try {
        #pragma omp parallel for shared(traj) private(openmp_t_ind) schedule(dynamic, 1)
        for(openmp_t_ind = 0; openmp_t_ind < NUM_TIME_STEPS; openmp_t_ind++) {
            traj.makePolyZono(openmp_t_ind);
        }
    }
    catch (int errorCode) {
        WARNING_PRINT("        CUDA & C++: Error creating JRS! Check previous error message!");
        return -1;
    }

    try {
        #pragma omp parallel for shared(kd, link_independent_generators) private(openmp_t_ind) schedule(dynamic)
        for(openmp_t_ind = 0; openmp_t_ind < NUM_TIME_STEPS; openmp_t_ind++) {
            // compute link PZs through forward kinematics
            kd.fk(openmp_t_ind);

            // reduce non-only-k-dependent generators so that slice takes less time
            for (int i = 0; i < NUM_JOINTS; i++) {
                link_independent_generators[openmp_t_ind * NUM_JOINTS + i] = kd.links(i, openmp_t_ind).reduce_link_PZ();
            }
        }
    }
    catch (int errorCode) {
        WARNING_PRINT("        CUDA & C++: Error computing link PZs and nominal torque PZs! Check previous error message!");
        return -1;
    }

    try {
        O.initializeHyperPlane(link_independent_generators);
    }
    catch (int errorCode) {
        WARNING_PRINT("        CUDA & C++: Error initializing collision checking hyperplanes! Check previous error message!");
        return -1;
    }

    auto stop1 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - start1);
    cout << "        CUDA & C++: Time taken by generating trajectory & forward kinematics: " << duration1.count() << " milliseconds" << endl;

/*
Section III:
    Solve the optimization problem using IPOPT
*/
    auto start2 = std::chrono::high_resolution_clock::now();

    SmartPtr<armtd_NLP> mynlp = new armtd_NLP();
	mynlp->set_parameters(q_des, &traj, &kd, &O);

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", IPOPT_OPTIMIZATION_TOLERANCE);
	app->Options()->SetNumericValue("max_wall_time", IPOPT_MAX_WALL_TIME);
	app->Options()->SetIntegerValue("print_level", IPOPT_PRINT_LEVEL);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", IPOPT_LINEAR_SOLVER);
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-6);

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded ) {
		WARNING_PRINT("Error during initialization!");
        outputstream1 << -1;
        outputstream1.close();
        throw;
    }

    // Ask Ipopt to solve the problem
	status = app->OptimizeTNLP(mynlp);

    if (status == Maximum_CpuTime_Exceeded) {
        cout << "        CUDA & C++: Ipopt maximum CPU time exceeded!\n";
    }

    auto stop2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - start2);
    cout << "        CUDA & C++: Time taken by Ipopt: " << duration2.count() << " milliseconds" << endl;

/*
Section IV:
    Prepare output
*/
    // output k_opt
    // set precision to 10 decimal digits
    outputstream1 << std::setprecision(10);

    if (mynlp->feasible) {
        for (int i = 0; i < NUM_FACTORS; i++) {
            outputstream1 << mynlp->solution[i] << '\n';
        }
    }
    else {
        outputstream1 << -1 << '\n';
    }

    // output time cost (in milliseconds) in C++
    outputstream1 << duration1.count() + duration2.count();

    outputstream1.close();

    // output FRS and other information, you can comment them if they are unnecessary
    std::ofstream outputstream2(outputfilename2);
    outputstream2 << std::setprecision(10);
    for (int i = 0; i < NUM_TIME_STEPS; i++) {
        for (int j = 0; j < NUM_JOINTS; j++) {
            for (int l = 0; l < 3; l++) {
                outputstream2 << mynlp->link_sliced_center[i * NUM_JOINTS + j](l) << ' ';
            }
            outputstream2 << '\n';
        }
    }
    outputstream2.close();

    std::ofstream outputstream3(outputfilename3);
    outputstream3 << std::setprecision(10);
    for (int i = 0; i < NUM_TIME_STEPS; i++) {
        for (int j = 0; j < NUM_JOINTS; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3 + 3; l++) {
                    outputstream3 << link_independent_generators[i * NUM_JOINTS + j](k, l) << ' ';
                }
                outputstream3 << '\n';
            }
        }
    }
    outputstream3.close();

    std::ofstream outputstream4(outputfilename4);
    outputstream4 << std::setprecision(6);
    for (int i = 0; i < mynlp->constraint_number; i++) {
        outputstream4 << mynlp->g_copy[i] << '\n';
    }
    outputstream4.close();

    return 0;
}
