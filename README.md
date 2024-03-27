![](https://github.com/roahmlab/armour/blob/main/assets/armour_logo.png?raw=true)

# Autonomous Robust Manipulation via Optimization with Uncertainty-aware Reachability
**Authors:** Jonathan Michaux (jmichaux@umich.edu), Patrick Holmes (pdholmes@umich.edu), Bohao Zhang (jimzhang@umich.edu), Che Chen (cctom@umich.edu), Baiyue Wang (baiyuew@umich.edu), Shrey Sahgal (shreyps@umich.edu), Tiancheng Zhang (zhangtc@umich.edu), Sidhartha Dey (sid.dey@agilityrobotics.com), Shreyas Kousik (skousik@gatech.edu), and Ram Vasudevan (ramv@umich.edu). 

- This work is supported by the Ford Motor Company via the Ford-UM Alliance under award N022977, National Science Foundation Career Award 1751093 and by the Office of Naval Research under Award Number N00014-18-1-2575.
- `ARMOUR` was developed in [Robotics and Optimization for Analysis of Human Motion (ROAHM) Lab](http://www.roahmlab.com/) at University of Michigan - Ann Arbor.
- Please check out our [project page](https://roahmlab.github.io/armour/)!

## Introduction
Robotic manipulators have the potential to assist humans in a wide variety of collaborative settings, such as manufacturing, package delivery, and in-home care.
However, such settings are typically constrained and uncertain; nevertheless, the robot must operate in a safety-critical fashion. 
This makes it challenging to directly apply high-torque manipulators that can ignore uncertainty due to their own mass and the mass of manipulated objects.
Instead, it is necessary to develop motion planning and control strategies that can operate safely by accounting for these types of uncertainty in real-time.
In this context, safety means avoiding collisions while obeying joint position, velocity, and torque limits.
To address the safety challenge, this paper proposes **Autonomous Robust Manipulation via Optimization with Uncertainty-aware Reachability**, a method for guaranteed-safe, real-time manipulator motion planning and control.
An overview of this method is given in figure below.

<img height="360" src="/assets/armour_summary_figure.png"/>

<!---<img height="270" src="/assets/armour_method_figure.pdf"/>-->

## Dependency
The repo has been verified on MATLAB R>=2021b and Ubuntu >= 20.04

This repo depends on the following repos:
 - [CORA 2021](https://tumcps.github.io/CORA/pages/archive/v2021/index.html)
 
You need to download this repo and add to your MATLAB path.

This repo assumes that you have installed the following libraries:

 - libboost-dev
 - libeigen3-dev (3.3.7)
 - libipopt
 - libcoinhsl
 
#### Install Boost C++ library
Simply run the following command:

     sudo apt install libboost-dev 

#### Install Eigen3
In this work, we use Boost C++ library for interval arithmetic computation. 
We further put intervals into Eigen library to create interval matrix and interval vector. 
However, currently we only know how to do this for Eigen 3.3.7.
We have not found any solutions to dump Boost intervals to the latest version of Eigen yet.
If you are working in Ubuntu 20.04, the default Eigen library version should be Eigen 3.3.7, so you can install the correct version of Eigen simply by running the following command:

     sudo apt install libeigen3-dev 

If you are working in later version of Ubuntu, you would have to manually install Eigen 3.3.7.
We provide a possible way to do this in the instructions below.
Download `eigen-3.3.7` by following [this link](https://gitlab.com/libeigen/eigen/-/releases/3.3.7).

     cd ~/Downloads
     tar -xvzf eigen-3.3.7.tar.gz
     mv eigen-3.3.7 /your/favorite/path/
     cd /your/favorite/path/eigen-3.3.7
     mkdir build && cd $_
     cmake ..
     sudo make
     sudo make install
 
#### Install Ipopt and HSL (TODO: Maybe add more introduction, but will be fairly long)
libipopt and libcoinhsl could be very annoying to install and to work with MATLAB. 
Suppose libipopt and libcoinhsl are both installed in /usr/local/lib.
You need to add that path to both user's environmental variable 'LD_LIBRARY_PATH' and MATLAB's environment variable 'LD_LIBRARY_PATH'
Check [here](https://www.mathworks.com/help/matlab/matlab_external/set-run-time-library-path-on-linux-systems.html) and [here](https://stackoverflow.com/questions/13428910/how-to-set-the-environmental-variable-ld-library-path-in-linux) for more information.

## Building
Run 
 - initialize.m
 - kinova_src/initialize.m
 
in MATLAB before you run any other scripts!

## Usage
Before running any scripts, make sure you run the initalization scripts successfully and put the Ipopt libraries in the proper path.

All of our results in the paper is developed based on [Kinova Gen3](https://www.kinovarobotics.com/product/gen3-robots). 
All of the related test scripts are included in [kinova_src](https://github.com/roahmlab/armour/tree/main/kinova_src).
Check the [README](https://github.com/roahmlab/armour/blob/main/kinova_src/README.md) in that folder for more information.

## License

`ARMOUR` is released under a [GNU license](https://github.com/roahmlab/armour/blob/main/LICENSE). For a list of all code/library dependencies, please check dependency section. For a closed-source version of `ARMOUR` for commercial purpose, please contact the authors. 

An overview of the theoretical and implementation details has been published in arxiv. 
If you use `ARMOUR` in an academic work, please cite using the following BibTex entry:
```
@article{article,
author = {Michaux, Jonathan and Holmes, Patrick and Zhang, Bohao and Chen, Che and Wang, Baiyue and Sahgal, Shrey and Zhang, Tiancheng and Dey, Sidhartha and Kousik, Shreyas and Vasudevan, Ram},
year = {2023},
month = {01},
pages = {},
title = {Can't Touch This: Real-Time, Safe Motion Planning and Control for Manipulators Under Uncertainty}
doi={10.48550/arXiv.2301.13308}}
```

