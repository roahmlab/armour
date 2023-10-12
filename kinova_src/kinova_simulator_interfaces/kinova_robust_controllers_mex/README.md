# Mexed Version of Robust Controller

## Prerequisite
- Ubuntu 20.04 or 22.04
- libeigen3-dev (3.3.7)
- libboost-dev (1.71)

## Usage

We use kinova (without any gripper) as an example to show how to create the robust controller mex instance for any arbitrary robot given the URDF file.

### Compile
 - Run 'compile.m' script to compile.

### Input and Output
 - In Matlab, the kinova_controller can be called using the following example:
```
[u, tau, v] = kinova_controller(Kr, alpha, V_max, r_norm_threshold, q, qd, q_des, qd_des, qdd_des);
```
 - Add more explanations in the future
