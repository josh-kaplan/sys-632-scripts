%% Homework 3.3 Plots
close all; clear all; clc;

% Orbit 1 - 300km altitude, circular orbit
a    = 6678;            % semimajor axis
x1_0 = 0;               % center coordinates
y1_0 = 0;
t    = 0:pi/100:2*pi;
x1   = x1_0 + a*cos(t);
y1   = y1_0 + a*sin(t);

% Orbit 2 - 800km altitude, circular orbit
a    = 7178;            % semimajor axis
x2_0 = 0;               % center coordinates
y2_0 = 0;
t    = 0:pi/100:2*pi;
x2   = x2_0 + a*cos(t);
y2   = y2_0 + a*sin(t);

% Transfer orbit
a    = (7178 + 6678) / 2;       % semimajor axis
e    = 1 - (6678 / a);          % eccentricity
b    = a*(1-e^2)/(1+cosd(90));  % semiminor axis
xt_0 = 0;                       % x, center of ellipse
yt_0 = (7178+6678)/2 - 6678;              % y, center of ellipse
t    = -pi/2:pi/100:pi/2;         % time step
xt   = xt_0+a*cos(t);
yt   = yt_0+b*sin(t);

red = [244/255, 67/255, 54/255];
blue = [33/255, 150/255, 243/255];
%green = [118/255, 255/255, 3/255];
green = [76/255, 175/255, 80/255];

hold on
plot(0,0,'w.')                  % center point
plot(x1, y1, 'Color', blue)     % parking orbit
plot(x2, y2, 'Color', red)      % mission orbit
plot(xt, yt, 'g.')              % transfer orbit
hold off
set(gca,'Color',[0.1 0.1 0.13]);
%set(gca,'Color',[0.95 0.95 0.95]);
%axis([-7200,7200,-7200,7200])
axis off
