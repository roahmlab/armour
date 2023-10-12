% script for demonstrating taking cos/sin of a zonotope
% using polynomial zonotopes, and plotting the output.

clear; clc;
figure(1); clf; hold on;

%%% define PZ version of a regular zonotope
c = pi/4;
g1 = 3*pi/8;
g2 = pi/16;
expMat = eye(2);
id = [1; 2];
theta = polyZonotope_ROAHM(c, [g1, g2], [], expMat, id);

%%% degree of taylor expansion for cos/sin
degree = 6;

%%% take taylor expansions
cos_theta = cos(theta, degree);
sin_theta = sin(theta, degree);

%%% reduce output if desired:
n_reduce = 40;
cos_theta = reduce(cos_theta, 'girard', n_reduce);
sin_theta = reduce(sin_theta, 'girard', n_reduce);

%%% combine and plot
cos_sin_theta = exactCartProd(cos_theta, sin_theta);
p = plot(cos_sin_theta, [1, 2], 'FaceColor', [0.4 0.4 0.4], 'EdgeColor',...
    [0 0 0], 'Filled', true, 'LineWidth', 2);
p.FaceAlpha = 0.1;
p.EdgeAlpha = 0.5;

%%% take a slice and plot
slice_pt = -0.75;
cos_sin_theta_slice = getSubset(cos_sin_theta, 1, slice_pt);
p2 = plot(cos_sin_theta_slice, [1, 2], 'FaceColor', [0 1 0], 'EdgeColor',...
    [0 0.8 0], 'Filled', true, 'LineWidth', 2);
p.FaceAlpha = 0.2;
p.EdgeAlpha = 0.5;

%%% format
axis equal; axis square;
grid on;
xlabel('cos(theta)');
ylabel('sin(theta)');
plot(0, 0, 'k.', 'MarkerSize', 20);

