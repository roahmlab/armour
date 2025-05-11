import numpy as np
import pickle
import matplotlib.pyplot as plt
import pybullet as p
import time
import pinocchio as pin
import scipy.io as sio
import argparse
import os

from kinova_dynamics import goal_distance, set_position, integrate, armour_Bezier_trajectory

import RAPTOR.build.lib.armour_nanobind as armour
import RobustMovers.build.lib.kinova_controller_armour_nanobind as controller_armour

### initializations
parser = argparse.ArgumentParser(description="Receive an integer.")
parser.add_argument('--world_idx', type=int, required=True, help='An integer value')
args = parser.parse_args()

# Load the URDF file
urdf_filename = "RAPTOR/Robots/kinova-gen3/kinova.urdf"

p.connect(p.DIRECT)
p.setGravity(0, 0, -9.81)
robot = p.loadURDF(urdf_filename, useFixedBase=True)
num_joints = p.getNumJoints(robot)

model = pin.buildModelFromUrdf(urdf_filename)
model.armature = np.array([8.03, 11.9962024615303644, 9.0025427861751517, 11.5806439316706360, 8.4665040917914123, 8.8537069373742430, 8.8587303664685315])
model.gravity.linear = np.array([0, 0, -9.81])

# Initialize the Armour planner
display_info = False
planner_config_filename = "kinova_planner_parameters.yaml"
planner = armour.ArmourPybindWrapper(urdf_filename, planner_config_filename, display_info)

# Initialize the Armour controller (make sure these are consistent with kinova_planner_parameters.yaml)
controller_config_path = "kinova_model_parameters_with_armature.yaml"
Kr = 5 * np.ones(model.nv)
V_max = 0.01
alpha = 20
r_norm_threshold = 1e-12
controller_instance = controller_armour.kinova_controller_armour_pybindwrapper(
    urdf_filename, controller_config_path, Kr, V_max, alpha, r_norm_threshold, "ignore"
)

def controller(q, q_d, qd, qd_d, qd_dd):
    return controller_instance.update(q, q_d, qd, qd_d, qd_dd)

# Load the world file
dir_path = "/workspaces/armour/saved_worlds/random/"
world_files = []
for filename in os.listdir(dir_path):
    file_path = os.path.join(dir_path, filename)
    if os.path.isfile(file_path) and filename.endswith('.csv'):
        world_files.append(file_path)
    
world_idx = args.world_idx
print(world_files[world_idx])
world_info = np.loadtxt(world_files[world_idx], delimiter=',')

start = world_info[0, :7]
goal = world_info[1, :7]
set_position(robot, start)

# Initialize the obstacles
obstacles = np.zeros((world_info.shape[0]-3, 9))
for i in range(3, world_info.shape[0]):
    obstacles[i-3, :3] = world_info[i, :3]
    obstacles[i-3, 3:6] = 0.0
    obstacles[i-3, 6:9] = world_info[i, 3:6]

# create and visualize the obstacles in pybullet
for i, obs in enumerate(obstacles):
    position = obs[0:3]
    orientation = p.getQuaternionFromEuler(obs[3:6])
    size = obs[6:9]
    
    boxCollisionShapeId = p.createCollisionShape(shapeType=p.GEOM_BOX,
                                                 halfExtents=size*0.5)
    boxVisualShapeId = p.createVisualShape(shapeType=p.GEOM_BOX,
                                           halfExtents=size*0.5,
                                           rgbaColor=[1,0,0,0.2])
    box_id = p.createMultiBody(baseMass=0,
                               baseVisualShapeIndex=boxVisualShapeId,
                               baseCollisionShapeIndex=boxCollisionShapeId,
                               basePosition=position,
                               baseOrientation=orientation)

# set up the local planner
duration = 2.0 # duration of the trajectory
t_plan = 0.5 # minimize the distance between the middle of the trajectory and the local goal

lookahead_distance = 0.2 # distance to the local goal

planner.set_obstacles(obstacles)

# ipopt settings
planner.set_ipopt_parameters(
    1e-5,          # ipopt_tol
    1e-6,          # ipopt_constr_viol_tol
    10.0,          # ipopt_obj_scaling_factor
    0.7,           # ipopt_max_wall_time
    0,             # ipopt_print_level
    "monotone",    # mu_strategy
    "ma86",        # ipopt_linear_solver
    False          # ipopt_gradient_check
)

waypoints = []
solutions = []

positions = np.array([])
velocities = np.array([])
torques = np.array([])

print("start", start)

q0 = start
q_d0 = np.zeros(7)
q_dd0 = np.zeros(7)

status = "not_reached"

for iter in range(400):
    # target information
    if np.linalg.norm(goal - q0) < lookahead_distance:
        q_des = goal
        t_plan = 1.0 # minimize the distance between the end of the trajectory and the final goal
    else:
        q_des = q0 + lookahead_distance * (goal - q0) / np.linalg.norm(goal - q0)
        t_plan = 0.5
    waypoints.append(q_des)
    
    k_center = np.zeros(7)
    k_range = np.pi/24 * np.ones(7)
    
    # set up the local planner
    planner.set_trajectory_parameters(
        q0, \
        q_d0, \
        q_dd0, \
        k_center, \
        k_range, \
        duration, \
        q_des, \
        t_plan
    )

    # optimize the trajectory
    k, ifFeasible = planner.optimize()
    
    print(k)
    
    if ifFeasible:
        solutions.append(k)
        
        q_next = q0 + k_center + k * k_range
        
        def desired_trajectory(t):
            return armour_Bezier_trajectory(q0, q_d0, q_dd0, q_next, t, duration)
        
        if np.linalg.norm(goal - q0) < lookahead_distance:
            local_duration = duration
        else:
            local_duration = 0.5 * duration
        
        ts_sim = np.linspace(0, local_duration, 1000)
        try:
            pos, vel, tau = integrate(model, ts_sim, np.concatenate([q0, q_d0]), desired_trajectory, controller)
        except Exception as e:
            status = 'torque limit exceeded'
            break
        
        # visualize the trajectory and check collision in pybullet
        for q in pos:
            set_position(robot, q)
            contacts = p.getContactPoints()
            if len(contacts) > 0:
                status = "crashed"
                break
            
        if status == "crashed":
            break
        
        q0_prev = q0
        q_d0_prev = q_d0
        q_dd0_prev = q_dd0
        q0, q_d0, q_dd0 = desired_trajectory(local_duration)
        
        if positions.shape[0] == 0:
            positions = pos
            velocities = vel
            torques = tau
        else:
            positions = np.vstack((positions, pos))
            velocities = np.vstack((velocities, vel))
            torques = np.vstack((torques, tau))
    else:
        q_next = q0_prev + k_center + solutions[-1] * k_range
        
        def desired_trajectory(t):
            return armour_Bezier_trajectory(q0_prev, q_d0_prev, q_dd0_prev, q_next, t, duration)
        
        ts_sim = np.linspace(0.5 * duration, duration, 1000)
        try:
            pos, vel, tau = integrate(model, ts_sim, np.concatenate([q0, q_d0]), desired_trajectory, controller)
        except Exception as e:
            status = 'torque limit exceeded'
            break
        
        # visualize the trajectory and check collision in pybullet
        for q in pos:
            set_position(robot, q)
            contacts = p.getContactPoints()
            if len(contacts) > 0:
                status = "crashed"
                break
            
        if status == "crashed":
            break
        
        q0, q_d0, q_dd0 = desired_trajectory(duration)
        
        if positions.shape[0] == 0:
            positions = pos
            velocities = vel
            torques = tau
        else:
            positions = np.vstack((positions, pos))
            velocities = np.vstack((velocities, vel))
            torques = np.vstack((torques, tau))

    if goal_distance(q0 - goal) < 0.05:
        print("Goal reached")
        status = "success"
        break
    else:
        print(np.linalg.norm(q0 - goal))
        print(np.linalg.norm(q0 - q_des))

sio.savemat(f'results/armour_trajectories_{world_idx}.mat', 
            {'status': status,
             'positions': positions, 
             'velocities': velocities,
             'torques': torques})