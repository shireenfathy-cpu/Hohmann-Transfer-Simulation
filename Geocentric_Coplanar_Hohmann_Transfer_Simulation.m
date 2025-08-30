% By Shireen Fathy
% Project: Geocentric Coplanar Hohmann Transfer Simulation

%% ===== Earth Constants =====
R_E = 6378.137;       % Earth radius [km]
mu  = 398600;         % Gravitational parameter [km^3/s^2]

fprintf('Geocentric Hohmann Transfer Simulation\n');

%% ===== User Inputs =====
h_initial = input('Enter initial orbit altitude [km]: ');
if h_initial < 160
    error('Unstable orbit! Try again.');
end
r_initial = R_E + h_initial;        % Initial orbit radius

h_final = input('Enter final orbit altitude [km]: ');
if h_final + R_E <= r_initial
    error('Final orbit smaller than initial! Try again.');
end
r_final = R_E + h_final;            % Final orbit radius

%% ===== Circular Orbit Velocities =====
v_initial = sqrt(mu / r_initial);   % Initial orbit velocity
v_final   = sqrt(mu / r_final);     % Final orbit velocity

%% ===== Transfer Orbit Calculations =====
a_transfer = (r_initial + r_final)/2;                % Semi-major axis
e_transfer = (r_final - r_initial)/(r_final + r_initial);  % Eccentricity
r_perigee  = a_transfer * (1 - e_transfer);         % Perigee distance
r_apogee   = a_transfer * (1 + e_transfer);         % Apogee distance
v_perigee  = sqrt((2*mu*r_final)/(r_initial*(r_initial+r_final)));  % Velocity at perigee
v_apogee   = sqrt((2*mu*r_initial)/(r_final*(r_initial+r_final)));  % Velocity at apogee
time_of_flight = (pi * sqrt(a_transfer^3 / mu)) / 3600; % Time of flight [hours]

%% ===== Delta V Calculations =====
dV1 = v_perigee - v_initial;      % First burn
dV2 = v_final - v_apogee;         % Second burn
dV_total = dV1 + dV2;             % Total delta V

%% ===== Display Results =====
fprintf('\n===== Results =====\n');

fprintf('Initial Orbit:\n  Altitude = %.2f km\n  Radius   = %.2f km\n  Velocity = %.4f km/s\n\n', ...
    h_initial, r_initial, v_initial);

fprintf('Final Orbit:\n  Altitude = %.2f km\n  Radius   = %.2f km\n  Velocity = %.4f km/s\n\n', ...
    h_final, r_final, v_final);

fprintf('Transfer Orbit:\n  Semi-major axis = %.2f km\n  Eccentricity = %.4f\n', ...
    a_transfer, e_transfer);
fprintf('  Perigee = %.2f km\n  Apogee  = %.2f km\n', r_perigee, r_apogee);
fprintf('  Velocity at Perigee = %.4f km/s\n  Velocity at Apogee = %.4f km/s\n', v_perigee, v_apogee);
fprintf('  Time of Flight = %.4f hours\n\n', time_of_flight);

fprintf('Delta V:\n  Delta V1 = %.4f km/s\n  Delta V2 = %.4f km/s\n  Total Delta V = %.4f km/s\n', ...
    dV1, dV2, dV_total);

%% ===== 3D Plot =====
figure; 
hold on; 
grid on; 
axis equal; 
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Geocentric Hohmann Transfer');

% Plot Earth
ax = gca;
plotEarthSphere(ax,'km');  % keep function name fixed

%% ===== Circular Orbits =====
n_points = 1000;
theta = linspace(0,2*pi,n_points);

% Initial orbit (green)
r_orbit_init = ones(1,n_points)*r_initial;
[X_init,Y_init] = pol2cart(theta,r_orbit_init);
Z_init = zeros(1,n_points);
plot3(X_init,Y_init,Z_init,'-c','LineWidth',1.5); 

% Final orbit (orange)
r_orbit_final = ones(1,n_points)*r_final;
[X_final,Y_final] = pol2cart(theta,r_orbit_final);
Z_final = zeros(1,n_points);
plot3(X_final,Y_final,Z_final,'-y','LineWidth',1.5); 

%% ===== Transfer Orbit Animation =====
coeT = [a_transfer e_transfer 0 0 0 0];
[rT, vT] = coe2eci(mu, coeT);

T_period = 2*pi*a_transfer*sqrt(a_transfer/mu);
sim_time = -(0.5*T_period/360);

% Initialize spacecraft marker (red)
sc_marker = plot3(0,0,0,'ro','MarkerSize',6,'MarkerFaceColor','b');

rx = zeros(1,361); ry = zeros(1,361); rz = zeros(1,361);

for i = 1:361
    sim_time = sim_time + (0.5*T_period/360);
    [r,~] = propagateTwoBody(mu, sim_time, rT, vT);

    rx(i) = r(1); ry(i) = r(2); rz(i) = r(3);

    % Update spacecraft marker position
    set(sc_marker,'XData',rx(i),'YData',ry(i),'ZData',rz(i));

    % Draw transfer trajectory (red line)
    if i>1
        plot3(rx(i-1:i), ry(i-1:i), rz(i-1:i), '-m','LineWidth',1.5);
    end

    pause(0.005);
end

% Final spacecraft position marker (black)
plot3(rx(end), ry(end), rz(end), 'ok','MarkerSize',7,'MarkerFaceColor','b'); 


%% ===== STK Integration =====
app = actxserver('STK11.application');
app.Visible = 1;
root = app.Personality2;

try
    root.CloseScenario();
catch
end

root.NewScenario('Hohmann_Scenario_MATLAB');
scenario = root.CurrentScenario;

inclination = 28.5;

% Initial satellite
sat_init = scenario.Children.New('eSatellite','Sat_Init');
try; 
sat_init.SetPropagatorType('ePropagatorTwoBody');
 end
sat_init.Propagator.InitialState.Representation.AssignClassical('eCoordinateSystemICRF', r_initial,0,inclination,0,0,0);
sat_init.Propagator.Propagate();

% Transfer satellite
sat_trans = scenario.Children.New('eSatellite','Sat_Trans');
try; 
sat_trans.SetPropagatorType('ePropagatorTwoBody'); 
end
sat_trans.Propagator.InitialState.Representation.AssignClassical('eCoordinateSystemICRF', a_transfer,e_transfer,inclination,0,0,0);
sat_trans.Propagator.Propagate();

% Final satellite
sat_final = scenario.Children.New('eSatellite','Sat_Final');
try; 
sat_final.SetPropagatorType('ePropagatorTwoBody');
 end
sat_final.Propagator.InitialState.Representation.AssignClassical('eCoordinateSystemICRF', r_final,0,inclination,0,0,0);
sat_final.Propagator.Propagate();

disp('Finished: Initial, Transfer, and Final orbits created and propagated in STK.');
