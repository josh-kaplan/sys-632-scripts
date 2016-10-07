%% Homework 3
%
% Josh Kaplan
% contact@joshkaplan.org
%
% Problem Statement:
% Assume the YachtSat spacecraft will be deployed into a circular, 
% 300km altitude, parking orbit with an inclination of 57 deg. and must 
% maneuver to a circular, 800km altitude, mission orbit with an 
% inclination of 70 deg.
%
close all; clear all; clc;

%% Setup

% Earth characteristics
mu = 3.986004418e5;     % µ = G*m       [km^3*s^-2]
r_e = 6378;             % radius        [km]

% Parking Orbit
alt_1 = 300;            % altitude      [km]
r_1 = alt_1 + r_e;      % radius        [km]
i_1 = 57;               % inclination   [degrees]

% Mission Orbit
alt_2 = 800;            % altitude      [km]
r_2 = alt_2 + r_e;      % radius        [km]
i_2 = 70;               % inclination   [degrees]

%% Problem 3.1 
% Determine the velocity, specific mechanical energy and period of the 
% parking orbit.
v_1 = sqrt(mu/r_1);
e_1 = -mu/(2*r_1);
P_1 = 2*pi*sqrt(r_1^3/mu);

fprintf('\n---------- 3.1 ----------\n')
fprintf('vel.   = %10.3f km/s\n', v_1)
fprintf('sp. e. = %10.3f km^2*s-2\n', e_1)
fprintf('period = %10.3f s\n', P_1)

%% Problem 3.2
% Determine the velocity, specific mechanical energy and period of the 
% mission orbit.
v_2 = sqrt(mu/r_2);
e_2 = -mu/(2*r_2);
P_2 = 2*pi*sqrt(r_2^3/mu);

fprintf('\n---------- 3.2 ----------\n')
fprintf('vel.   = %10.3f km/s\n', v_2)
fprintf('sp. e. = %10.3f km^2*s-2\n', e_2)
fprintf('period = %10.3f s\n', P_2)


%% Problem 3.3
% Calculate the minimum ?V necessary to move the spacecraft from its 
% parking orbit to the mission orbit. Describe the process in words.

% Transfer orbit
a_t = (r_1 + r_2)/2;            % semimajor axis        [km]
r_p = r_1;                      % radius at perigee     [km]
r_a = r_2;                      % radius at apogee      [km]

% SMAD, p. 147
v_p = sqrt(mu*(2/r_1 - 1/a_t)); % velocity at perigee   [km/s]
v_a = sqrt(mu*(2/r_2 - 1/a_t)); % velocity at apogee    [km/s]

% delta_V_1 - Parking Orbit to transfer orbit
dV_1 = v_p - v_1;  

% delta_V_2 - Transfer orbit to mission-sized orbit
dV_2 = v_2 - v_a;

% Plane change -> dV = 2*V_i*sin(ø/2)
% ?V_p - the plane change delta-V, occurs at apogee of tranfer orbit
%
% I'm not sure where this burn can actually occur. It depends on RAAN
% I think. But we'll ignore it for now and just say it occurs at tranfer
% apogeee.
dV_p = 2*v_a*sin((i_2 - i_1)/2);

dV_total = dV_1 + dV_2 + dV_p;

fprintf('\n---------- 3.3 ----------\n')
fprintf('Total delta-V = %6.3f km/s\n', dV_total);


%% Problem 3.4
% Mission planners have decided to budget for three 90 deg re-phasing 
% maneuvers during the mission life time of each spacecraft. Each maneuver 
% must be performed in one day. How much total ?V should be budgeted for 
% all of these maneuvers combined?
fprintf('\n---------- 3.4 ----------\n')

rephasing_angle = 90;

fprintf('Mission Orbit:\n')
fprintf('  Period          %8.0f s\n', P_2);
fprintf('  Radius          %8.0f km\n', r_2);
fprintf('  Velocity        %8.0f km/s\n', v_2);

% Determine smallest possible rephasing orbit
r_phasing_apo = r_2;
r_phasing_per = r_e + 100;
a_phasing = (r_phasing_apo + r_phasing_per) / 2;

P_phasing = 2*pi*sqrt(a_phasing^3 / mu);

v_phasing_apo = sqrt(mu*(2/r_phasing_apo - 1/a_phasing));
v_phasing_per = sqrt(mu*(2/r_phasing_per - 1/a_phasing));

delta_V_phasing = abs(v_2 - v_phasing_apo);

drift_rate = 1080 * delta_V_phasing / v_2;

% fprintf('Smallest Possible Phasing orbit:\n')
% fprintf('  Period          %8.0f s\n', P_phasing);
% fprintf('  Semimajor Axis  %8.0f km\n', a_phasing);
% fprintf('  Rad @ apogee    %8.0f km\n', r_phasing_apo);
% fprintf('  Rad @ perigee   %8.0f km\n', r_phasing_per);
% fprintf('  Vel @ apogee    %8.3f km/s\n', v_phasing_apo);
% fprintf('  Vel @ perigee   %8.3f km/s\n', v_phasing_per);
% fprintf('  Drift rate      %8.2f deg/orbit\n', drift_rate);

% Iterate through smaller possible orbits to solve
lo = P_phasing;
hi = P_2;
mid = (hi+lo) / 2;
p = mid;

% This the tolerance for how close to a whole number or 
% orbits we need to be.
TOL = .01;

% Maximum number of iterations
MAX_ITER = round(p/2);     

% Our starting best delta-V.
best_dV = delta_V_phasing;
best_p = p;

possible_solutions = 0;

fprintf('Solving phasing maneuver ... ')
for i=1:MAX_ITER
    clear r_phasing_per;
    clear a_phasing; 
    clear P_phasing;
    clear v_phasing_apo;
    clear v_phasing_per;
    clear delta_V_phasing;
    clear drift_rate;

    % Determine orbit geometry
    a_phasing = nthroot(mu*(p / (2*pi))^2, 3);
    
    % Determine delta-V to get into phasing orbit 
    v_phasing_apo = sqrt(mu*(2/r_phasing_apo - 1/a_phasing));
    delta_V_phasing = abs(v_2 - v_phasing_apo);

    % Calculate drift rate in degrees per orbit
    drift_rate = 1080 * delta_V_phasing / v_2;
    n_orbits = rephasing_angle / drift_rate;
    
    % If maneuver is longer than a day, decrease our period estimate
    if n_orbits * p > 60*60*24
        p = round(p - 1);
    % Otherwise, increase it (larger period means less delta-V)
    else
        % If it meets our tolerance and is smaller current optimal 
        % delta-V, keep track of it.
        if abs(n_orbits - round(n_orbits)) < TOL && delta_V_phasing < best_dV
            best_dV = abs(v_phasing_apo - v_2);
            best_p = p;
            possible_solutions = possible_solutions + 1;
            %fprintf('Possible Solution: i=%i, p=%.0f\n', i, p);
        end
        p = round(p + 1);
    end
end
fprintf('done.\n')
fprintf('Found %i possible solutions in %d iterations.\n', possible_solutions, i)

p = best_p;
% Determine orbit geometry
a_phasing = nthroot(mu*(p / (2*pi))^2, 3);
r_phasing_per = 2*a_phasing - r_phasing_apo;

% Determine orbital velocities
v_phasing_apo = sqrt(mu*(2/r_phasing_apo - 1/a_phasing));
v_phasing_per = sqrt(mu*(2/r_phasing_per - 1/a_phasing));

% Determine delta-V to get into phasing orbit 
delta_V_phasing = abs(v_2 - v_phasing_apo);

% Calculate drift rate in degrees per orbit
drift_rate = 1080 * delta_V_phasing / v_2;
n_orbits = rephasing_angle / drift_rate;

fprintf('Phasing orbit:\n')
fprintf('  Period          %8.0f s\n', p);
fprintf('  Semimajor Axis  %8.0f km\n', a_phasing);
fprintf('  Rad @ apogee    %8.0f km\n', r_phasing_apo);
fprintf('  Rad @ perigee   %8.0f km\n', r_phasing_per);
fprintf('  Vel @ apogee    %8.3f km/s\n', v_phasing_apo);
fprintf('  Vel @ perigee   %8.3f km/s\n', v_phasing_per);
fprintf('  delta-V_phase   %8.3f km/s\n', delta_V_phasing)
fprintf('  delta-V_total   %8.3f km/s\n', 2*delta_V_phasing)
fprintf('  Drift rate      %8.2f deg/orbit\n', drift_rate);
fprintf('  No. of Orbits   %8.2f orbits\n', n_orbits);
fprintf('  Drift           %8.2f deg\n', drift_rate*n_orbits)
fprintf('  Drift Time      %8.2f hours\n', p*n_orbits / (60*60))

delta_V_3phasingmaneuvers = 3 * 2*delta_V_phasing;
fprintf('\n  dV, 3 maneuvers   %6.3f km/s\n', delta_V_3phasingmaneuvers)

%% Problem 3.5 
% If the nominal ballistic coefficient (BC) for each spacecraft is 75kg/m2 
% and the worst case atmospheric density (?) estimate is 4.39*10-14 kg/m3
% at 800 km altitude how much spacecraft ?V should be budgeted for drag 
% compensation during a 6 year mission lifetime?

BC = 75;
rho = 4.39e-14;

% BC = m / (C_d * A)
% Drag ?V = ?*(1/BC)*?*a*v in m/s per orbit
delta_V_drag = pi * (1/BC) * rho * r_2 * v_2; 

% Total ?V due to drag over lifetime
life = 6 * 365.25 * 25 * 60 * 60;
lifetime_orbits = life / P_2;
delta_V_drag = delta_V_drag * lifetime_orbits;

% Convert to km/s for consistency
delta_V_drag = delta_V_drag * 1000;

fprintf('\n---------- 3.5 ----------\n')
fprintf('Total delta-V due to drag %8.3f km/s\n', delta_V_drag);

%% Problem 3.6
% International orbital debris mitigation guidelines recommend each 
% spacecraft be de-orbited from its mission orbit to the surface of the 
% earth at end-of-life (EOL). How much ?V will this de-orbit maneuver 
% require?

r_deorbit_per = r_e + 100; % perigee of deorbit trajectory [km]
r_deorbit_apo = r_2;
a_deorbit = (r_deorbit_apo + r_deorbit_per) / 2;

v_deorbit_apo = sqrt(mu*(2/r_deorbit_apo - 1/a_deorbit));
delta_V_deorbit = abs(v_2 - v_deorbit_apo);

fprintf('\n---------- 3.6 ----------\n')
fprintf('Total delta-V to deorbit %8.3f km/s\n', delta_V_deorbit);


%% Problem 3.7
% 3.7 Based on the above calculations, what would be the total mission ?V 
% budget for each spacecraft without margin? With margin?

total_delta_V = dV_total + delta_V_3phasingmaneuvers + delta_V_drag + delta_V_deorbit;

fprintf('\n---------- 3.7 ----------\n')
fprintf('delta-V budget:\n')
fprintf('  parking to mission orbit %8.3f km/s\n', dV_total);
fprintf('  3 re-phasing maneuvers   %8.3f km/s\n', delta_V_3phasingmaneuvers);
fprintf('  drag compensation        %8.3f km/s\n', delta_V_drag);
fprintf('  deorbit maneuver         %8.3f km/s\n', delta_V_deorbit);
fprintf('------------------------------------------\n')
fprintf('  TOTAL                    %8.3f km/s\n', total_delta_V);

