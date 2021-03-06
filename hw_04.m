%% Homework 4
%
% Josh Kaplan
% contact@joshkaplan.org
%
% As a general rule of thumb for imaging missions, a minimum 60 deg 
% elevation angle (?) will be assumed for this payload.
close all; clear all; clc;

%% Setup 
mu  = 3.986004418e5;    % km^3*s-2
r_e = 6378;             % radius of earth   [km]
alt = 800;              % orbit alt         [km]
r   = r_e + alt;        % orbit radius      [km]
e   = 60;               % elevation angle   [deg]

%% Problem 4.1 
% For the assumed mission orbit in Question 3, and the assumed minimum 
% elevation angle, what is the spacecraft sensor field-of-regard (FOR)?

lambda_o = acosd(r_e / r);                  % SMAD, eq.  5-24
rho = 90 - lambda_o;                        % SMAD, fig. 5-13
nadir_angle = asind(cosd(e) * sind(rho));   % SMAD, eq. 5-25
FOR = 2 * nadir_angle;

fprintf('\n---------- 4.1 ----------\n')
fprintf('Field of Regard (FOR)  %6.2f deg\n', FOR)

%% Problem 4.2 
% What swath width (SW) (in km) is generated by the above FOR? 
lambda = 90 - nadir_angle - e;    % swath width   [deg]
swath_width = sind(2*lambda)*r_e;  % swath width   [km] 

fprintf('\n---------- 4.2 ----------\n')
fprintf('Swath Width            %6.2f deg\n', lambda)
fprintf('Swath Width            %6.0f km\n', swath_width)

%% Problem 4.3
% What is the slant range to the edge of the FOR?

D = r_e * sind(lambda) / sind(nadir_angle);

fprintf('\n---------- 4.3 ----------\n')
fprintf('Slant Range            %6.2f km\n', D)

%% Problem 4.4 
% Can a one spacecraft constellation with the above described FOR and 
% mission orbit achieve full earth coverage twice each day in daylight? 
% Why or why not?

P = 2*pi*sqrt(r^3/mu);
ACR = (4*pi/P) * sind(lambda); 

coverage_per_orbit = ACR*r_e^2 * P;
orbits_per_day = 24*60*60 / P;
coverage_per_day = coverage_per_orbit * orbits_per_day;

fprintf('\n---------- 4.4 ----------\n')
fprintf('ACR                    %.3e km^2/s\n', ACR*r_e^2)
fprintf('Area per Orbit         %.3e km^2/s\n', coverage_per_orbit)
fprintf('Coverage per Day       %.3e km^2/s\n', coverage_per_day)
