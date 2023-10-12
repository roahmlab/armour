function A = reset_agent_with_chomp_traj(A,traj_filename)
traj = textread(traj_filename) ;
traj = traj(:,1:6)' ;
N = size(traj,2) ;
dt = 0.03 ;
t = linspace(0,dt*N,N) ;
A.time = t ;
A.state = zeros(A.n_states,N) ;
A.state(A.joint_state_indices,:) = traj ;
end