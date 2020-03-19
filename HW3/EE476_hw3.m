%Read in the initial data for stars as a data table
data = readtable("stars.txt");
% Make the column headers look nice
data.Properties.VariableNames = {'ID', 'R', 'i', 'Omega', 'phi', 'theta'};
data.i = deg2rad(data.i);
data.Omega = deg2rad(data.Omega);
data.phi = deg2rad(data.phi);
data.theta = deg2rad(data.theta);
% Using the initial data, set up the initial state for the state variables
% equations utilized from the second part of instructions, t=0
% x1 = x position, x2 = x velocity, x3 = y position, x4 = y velocity
% x5 = z position, x6 = z velocity
initial = table(...
    ... % x1(0)
    data.R .* (cos(data.phi) .* cos(data.Omega) - ...
    sin(data.phi) .* cos(data.i) .* sin(data.Omega)), ...
    ... % x2(0)
    orbit_speed(data.R) .* (-sin(data.phi) .* cos(data.Omega) - ...
    cos(data.phi) .* cos(data.i) .* sin(data.Omega)), ...
    ... % x3(0)
    data.R .* (cos(data.phi) .* sin(data.Omega) + ...
    sin(data.phi) .* cos(data.i) .* cos(data.Omega)), ...
    ... % x4(0)
    orbit_speed(data.R) .* (-sin(data.phi) .* sin(data.Omega) + ...
    cos(data.phi) .* cos(data.i) .* cos(data.Omega)), ...
    ... % x5(0)
    data.R .* (sin(data.phi) .* sin(data.i)), ...
    ... % x6(0)
    orbit_speed(data.R) .* (cos(data.phi) .* sin(data.i)), ...
    'VariableNames', {'x1i', 'x2i' , 'x3i', 'x4i', 'x5i', 'x6i'});

% x, dx/dt, y, dy/dt, z, dz/dt
ss_initial = [-1.12,-211.48,-3.10,129.12,-3.40,-48.23];
% time (Myrs)
t = linspace(0, 90, 5000);
% For debugging: planet ID
id = 5;

% **********
% * Solve space ship motion
% **********
[ss_t, ss_sol] = ode23(@(t, x) ssmodel(x), t, ss_initial);
ss_x = ss_sol(:, 1);
ss_y = ss_sol(:, 3);
ss_z = ss_sol(:, 5);


% **********
% * Solve stars' motion
% **********
% n, defined in instructions
n = orbit_speed(data.R)./data.R;

x = (data.R) .* (cos(n.*t + data.phi) .* cos(data.Omega) - ...
    sin(n.*t + data.phi) .* cos(data.i) .* sin(data.Omega));
y = data.R .* (cos(n.*t + data.phi) .* sin(data.Omega) + ...
    sin(n.*t + data.phi) .* cos(data.i) .* cos(data.Omega));
z = data.R .* (sin(n.*t + data.phi) .* sin(data.i));

% loc = table(...
%     ... % x
%     data.R(id) .* (cos(n(id)*t + data.phi(id)) .* cos(data.Omega(id)) - ...
%     sin(n(id)*t + data.phi(id)) .* cos(data.i(id)) .* sin(data.Omega(id))), ...
%     ... % y
%     data.R(id) .* (cos(n(id)*t + data.phi(id)) .* sin(data.Omega(id)) + ...
%     sin(n(id)*t + data.phi(id)) .* cos(data.i(id)) .* cos(data.Omega(id))), ...
%     'VariableNames', {'x1', 'y1'});

%% Plot
figure(1)
clf
xlim([-40 40])
ylim([-40 40])
hold on

% F(100) = struct('cdata',[],'colormap',[]);
num_stars = length(n); % Set the number of stars to plot
for i = 1:100
    if i == 1
        ax_star = plot3(x(1:num_stars, i), y(1:num_stars, i), ...
                        z(1:num_stars, i), '.', 'MarkerSize', 1, ...
                        'Color', 'white'); % [0.3010, 0.7450, 0.9330]
        ax_ss = plot3(ss_x(i), ss_y(i), ss_z(i), 'r.', 'MarkerSize',30);
        h = gca;
        set(h, 'Color', [0.30, 0.30, 0.30], 'XColor', [1 1 1], ...
            'YColor', [1 1 1]);
        set(gcf, 'Color', 'k');
        axis equal
        set(ax_star, 'XLimInclude', 'off', 'ZLimInclude', 'off');
        set(ax_ss, 'XLimInclude', 'off', 'ZLimInclude', 'off');
        xlabel('X', 'fontsize', 20);
        ylabel('Y', 'FontSize', 20);
        zlabel('Z', 'FontSize', 20);
        title('Propegation of stars and spaceship', 'FontSize', 20, ...
              'Color', 'white');
        legend('\color{white}Star', '\color{white}Space Ship', ...
               'FontSize', 18, 'Color', 'k');
        rotate3d on
%         view(60,20);
    end
    pause(0.1)
    ax_star.XData = x(1:num_stars, i);
    ax_star.YData = y(1:num_stars, i);
    ax_star.ZData = z(1:num_stars, i);
    ax_ss.XData = ss_x(i);
    ax_ss.YData = ss_y(i);
    ax_ss.ZData = ss_z(i);
%     F(i) = getframe(gcf);
end


 
 % defined in instructions
 function vc = orbit_speed(r)
     k = [...
         0.00287729,   0.0023821,  -0.0010625,   0.000198502, ...
         -1.88428e-05, 9.70521e-07, -2.70559e-08, 3.7516e-10, -1.94316e-12
         ];
     vc = 1 ...
          ./ (k(9)*r.^8 + k(8)*r.^7 + k(7)*r.^6 + k(6)*r.^5 ...
             + k(5)*r.^4 + k(4)*r.^3 + k(3)*r.^2 + k(2)*r.^1 + k(1));
 end
 
 % defined in instructions
 function f = acceleration(r)
     f = orbit_speed(r)^2 / r;
 end
 
 % state space model according to the velocity and acceleration equations
 % x(1) = x
 % x(2) = dx/dt
 % x(3) = y
 % x(4) = dy/dt
 % x(5) = z
 % x(6) = dz/dt
 function dx = ssmodel(x)
     dx = zeros(size(x));
     r = sqrt(x(1)^2 + x(3)^2 + x(5)^2);
     f = acceleration(r); 
     dx(1) =  x(2);
     dx(2) = -x(1) / r * f;
     dx(3) =  x(4);
     dx(4) = -x(3) / r * f;
     dx(5) =  x(6);
     dx(6) = -x(5) / r * f;
 end