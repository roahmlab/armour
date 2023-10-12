# Making an Arm Agent

In this tutorial, we'll construct arm agents. First, we'll make a 2-D, 3-link, 3-DOF arm. Then we'll make a 3-D, 3-link, 6-DOF arm.

## 1. Making a 2-D Arm Agent

We will create our arm as a MATLAB class. First, let's create a skeleton for this class:

```matlab
classdef my_arm < robot_arm_agent
	methods
		function A = my_arm(varargin)
			% we will fill this part in
		end
end
```

Save this in a file named `my_arm.m`.

### 1.1. Required Properties

To create an arm, we have to specify, at a minimum, the following properties:

- `n_links_and_joints`: the number of links and joints of the arm; for these agents, we are using the convention that every joint is succeeded by a link (and every link, except the baselink, is preceded by a joint), hence the number of links (besides the baselink) and joints is the same
- `dimension`: a scalar value set to either 2 or 3, which is the dimension of the arm's workspace
- `link_sizes`: a 2-by-N or 3-by-N array, where the number of rows is the dimension of the arm's workspace, and N is the number of links/joints of the arm
- `joint_locations`: a 4-by-N or 6-by-N array, where the number of of rows is twice the dimension of the arm's workspace and each column is the location of the joint in its predecessor and successor links' local coordinate frames; the first 2 (or 3) rows are the location of the joint in its predecessor link's frame, and the last 2 (or 3) rows are the location of the joint in its successor link's frame
- `joint_state_limits`: a 2-by-N array, where each column corresponds to the state of one joint (each joint's state is a real number in radians); the first row is the minimum value and the second row is the maximum value
- `joint_speed_limits`: a 2-by-N array just like `joint_limits`, but for the joint speeds as opposed to joint positions
- `joint_input_limits`: a 2-by-N array just like `joint_limits`, but for the joint inputs (which are commanded accelerations) as opposed to the joint positions

All other properties can be auto-filled by the `robot_arm_agent` superclass. For a rundown of these other properties, take a look at Section 3 below.

### 1.2 Suggested Properties

The following properties should _probably_ be set by the user, because the default values may not produce a nice arm:

- `link_shapes`: this is a 1-by-N cell array where each cell contains the strings 'box' or 'oval' (for a 2-D arm) or 'cuboid' or 'ellipsoid' for a 3-D arm; this obviously specifies the shape, and defaults to boxes or cuboids for the respective dimensions
- `joint_axes`: this is a 2-by-N or 3-by-N array where each column is the axis of the corresponding joint in that joint's predecessor link's local coordinate frame; by default, these are all the z-axis, which is fine for a 2-D arm, but makes for a very boring 3-D arm

Now let's make an arm!



### 1.2. Making the 2-D Arm

Now we can fill in the arm class by specifying the values of the properties listed above in the arm's constructor method:

```matlab
classdef my_arm < robot_arm_agent
	methods
		function A = my_arm(varargin)
			% we will fill this part in
		end
end
```





## 2. Making a 3-D Arm Agent

Coming soon!



## 3. Arm Agent Properties

Coming soon!



## 4. Arm Agent Methods

Coming soon!