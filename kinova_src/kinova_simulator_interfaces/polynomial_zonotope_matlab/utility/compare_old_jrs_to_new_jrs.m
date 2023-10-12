% patrick 20220330
% just a quick script to plot a comparison of the old offline JRSs (found
% through CORA reachability) and new JRS (just from poly zono arithmetic)

figure(1); clf; hold on;

load('/Users/pdholmes/Documents/MATLAB/dynamics-dev/dynamics/rnea/polynomial_zonotope/jrs_saved/JRS_1.210.mat')
taylor_degree = 1; % used to take sine/cosine taylor expansion
joint_axes = [0; 0; 1];
traj_type = 'bernstein';
q = pi/2;
dq = pi/6;
ddq = pi/12;
[Q_des, Qd_des, Qdd_des, Q, Qd, Qd_a, Qdd_a, R_des, R_t_des, R, R_t, id, id_names] = create_jrs_online(q, dq, ddq, joint_axes, taylor_degree, traj_type);

% for i = 1:100
% p = plot(JRS{i}, [1, 2], 'b', 'Filled', true);
% p.FaceAlpha = 0.05;
% p.EdgeAlpha = 0.2;
% end

for i = 1:5:100
% pz = exactCartProd(cos(Q_des{i}{1}, deg), sin(Q_des{i}{1}, deg));
pz = R_des{i}{1}*[1; 0; 0];
p2 = plot(pz, [1, 2], 'g', 'Filled', true);
p2.FaceAlpha = 0.05;
p2.EdgeAlpha = 0.2;
end
% 
% axis equal;

% just compare final time step
figure(2); clf; hold on;
% plot(JRS{100}, [1, 2], 'b');
% pz = exactCartProd(cos(Q_des{100}{1}, deg), sin(Q_des{100}{1}, deg));
pz = R_des{100}{1}*[1; 0; 0];
plot(pz, [1, 2], 'g');


% plot JRSs to see if they make sense...
t = 0.005:0.01:0.995;
figure(3); clf; hold on;
figure(4); clf; hold on;
figure(5); clf; hold on;

for i = 1:length(t)
	figure(3);
	plot(t(i), Q_des{i}{1}.c + abs(Q_des{i}{1}.G) + abs(Q_des{i}{1}.Grest), 'k.')
	plot(t(i), Q_des{i}{1}.c - abs(Q_des{i}{1}.G) - abs(Q_des{i}{1}.Grest), 'k.')
	plot(t(i), Q{i}{1}.c + abs(Q{i}{1}.G) + abs(Q{i}{1}.Grest), 'r.')
	plot(t(i), Q{i}{1}.c - abs(Q{i}{1}.G) - abs(Q{i}{1}.Grest), 'r.')
    
    Q_tmp = getSubset(Q_des{i}{1}, 1, 0.5);
	plot(t(i), Q_tmp.c + abs(Q_tmp.Grest), 'g.')
    plot(t(i), Q_tmp.c - abs(Q_tmp.Grest), 'g.')

	figure(4);
	plot(t(i), Qd_des{i}{1}.c + abs(Qd_des{i}{1}.G) + abs(Qd_des{i}{1}.Grest), 'k.')
	plot(t(i), Qd_des{i}{1}.c - abs(Qd_des{i}{1}.G) - abs(Qd_des{i}{1}.Grest), 'k.')
	plot(t(i), Qd{i}{1}.c + abs(Qd{i}{1}.G) + abs(Qd{i}{1}.Grest), 'r.')
	plot(t(i), Qd{i}{1}.c - abs(Qd{i}{1}.G) - abs(Qd{i}{1}.Grest), 'r.')	
	plot(t(i), Qd_a{i}{1}.c + abs(Qd_a{i}{1}.G) + abs(Qd_a{i}{1}.Grest), 'b.')
	plot(t(i), Qd_a{i}{1}.c - abs(Qd_a{i}{1}.G) - abs(Qd_a{i}{1}.Grest), 'b.')	

	figure(5);
	plot(t(i), Qdd_des{i}{1}.c + abs(Qdd_des{i}{1}.G) + abs(Qdd_des{i}{1}.Grest), 'k.')
	plot(t(i), Qdd_des{i}{1}.c - abs(Qdd_des{i}{1}.G) - abs(Qdd_des{i}{1}.Grest), 'k.')
	plot(t(i), Qdd_a{i}{1}.c + abs(Qdd_a{i}{1}.G) + abs(Qdd_a{i}{1}.Grest), 'b.')
	plot(t(i), Qdd_a{i}{1}.c - abs(Qdd_a{i}{1}.G) - abs(Qdd_a{i}{1}.Grest), 'b.')	
end

