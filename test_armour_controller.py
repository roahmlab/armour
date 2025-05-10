import pinocchio as pin
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import scipy.io
import pybullet as p
import time
import matplotlib.pyplot as plt

from kinova_dynamics import set_position, integrate, armour_Bezier_trajectory

# import sys
# sys.path.append('/workspaces/armour/RobustMovers/build/lib') # be careful about the path
# import kinova_controller_armour_nanobind as controller_armour
import RobustMovers.build.lib.kinova_controller_armour_nanobind as controller_armour

# Kinova URDF model
model_path = "RAPTOR/Robots/kinova-gen3/kinova.urdf"
config_path = "kinova_model_parameters_with_armature.yaml"

# Ratio of uncertainties in the inertial parameters
# make sure these are consistent with kinova_model_parameters.yaml
mass_model_uncertainty = 0.05 # 5%
com_model_uncertainty = 0.05 # 5%
inertia_model_uncertainty = 0.05 # 5%

if __name__ == "__main__":
    # Load the URDF model
    model = pin.buildModelFromUrdf(model_path)
    
    # Add uncertainties to the model
    for i in range(model.nv):
        mass_uncertainty_range = (2 * mass_model_uncertainty * np.random.rand() - mass_model_uncertainty)
        model.inertias[i].mass *= (1 + mass_uncertainty_range)
        com_uncertainty_range = (2 * com_model_uncertainty * np.random.rand() - com_model_uncertainty)
        model.inertias[i].lever *= (1 + com_uncertainty_range)
        inertia_uncertainty_range = (2 * inertia_model_uncertainty * np.random.rand() - inertia_model_uncertainty)
        model.inertias[i].inertia *= (1 + inertia_uncertainty_range)
        
    # Set the armature
    model.armature = np.array([8.03, 11.9962024615303644, 9.0025427861751517, 11.5806439316706360, 8.4665040917914123, 8.8537069373742430, 8.8587303664685315])
    
    # Controller parameters
    Kr = 5 * np.ones(model.nv)
    V_max = 0.02
    alpha = 20
    r_norm_threshold = 1e-12
    controller_instance = controller_armour.kinova_controller_armour_pybindwrapper(
        model_path, config_path, Kr, V_max, alpha, r_norm_threshold
    )
    
    def controller(q, q_d, qd, qd_d, qd_dd):
        return controller_instance.update(q, q_d, qd, qd_d, qd_dd)
    
    # Set the initial condition (to be aligned with the beginning of the desired trajectory)
    q0 = np.random.rand(7) * 2 * np.pi
    q_d0 = np.random.rand(7) * 0.5 * np.pi
    q_dd0 = np.random.rand(7) * 0.2 * np.pi
    q1 = q0 + 0.2
    T = 4
    
    def desired_trajectory(t):
        return armour_Bezier_trajectory(q0, q_d0, q_dd0, q1, t, T)
    
    # Set the simulation time (4 seconds, 1000 steps for visualization)
    ts_sim = np.linspace(0, T, 1000)

    # Integrate the dynamics using pinocchio
    pos, vel, tau = integrate(model, ts_sim, np.concatenate([q0, q_d0]), desired_trajectory, controller)
    
    # Recover the desired trajectory
    desired_pos = np.zeros((len(ts_sim), 7))
    desired_vel = np.zeros((len(ts_sim), 7))
    for i, t in enumerate(ts_sim):
        desired_pos[i], desired_vel[i], _ = desired_trajectory(t)
        
    tracking_error = pos - desired_pos
    tracking_error_vel = vel - desired_vel
    
    # Initialize the PyBullet simulation
    p.connect(p.GUI)
    p.setGravity(0, 0, -9.81)
    robot = p.loadURDF(model_path, useFixedBase=True, basePosition=[0,0,0], baseOrientation=[0,0,0,1])
    num_joints = p.getNumJoints(robot)
    
    # Playback the simulation
    for q in pos:
        set_position(robot, q)
        time.sleep(0.01)
        
    # input("Press Enter to continue...")
    p.disconnect()
    
    # Plot the tracking error
    plt.figure()
    for i in range(7):
        plt.subplot(2, 4, i+1)
        plt.plot(ts_sim, tracking_error[:,i] * 180 / np.pi, label=f"Joint {i}")
        plt.xlabel("Time [s]")
        plt.ylabel("Tracking error [deg]")
        plt.legend()
    plt.savefig("tracking_error.png")
    
    # Plot the control input
    plt.figure()
    for i in range(7):
        plt.subplot(2, 4, i+1)
        plt.plot(ts_sim, tau[:,i], label=f"Joint {i}")
        plt.xlabel("Time [s]")
        plt.ylabel("Control input [Nm]")
        plt.legend()
    plt.savefig("control_input.png")

