# g++ -o test PZ_tests.cu Dynamics.cu Trajectory.cu PZsparse.cu -fopenmp -std=c++14 -lstdc++ -ldl -lm -lc -lgomp -O2
nvcc -o test PZ_tests.cu CollisionChecking.cu Dynamics.cu Trajectory.cu PZsparse.cu -Xcompiler -fopenmp -std=c++14 -O2 -I/usr/local/include -I/usr/local/include/coin-or -L/usr/local/lib -L/usr/lib -lipopt -lquadmath -lstdc++ -ldl -lm -lc -lgomp
