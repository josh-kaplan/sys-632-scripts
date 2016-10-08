%% Homework 3.3 Plots, Polar
close all; clear all; clc;

% Orbit 1 - 300km altitude, circular orbit
a    = 6678;            % semimajor axis
e    = 0;
t1   = 0:pi/100:2*pi;
N    = length(t1);
r1   = zeros(0, N);
for i=1:length(t1)
    r1(i) = a*(1-e^2)/(1+e*cos(t1(i)));
end

% Orbit 2 - 800km altitude, circular orbit
a    = 7178;            % semimajor axis
e    = 0;
t2   = 0:pi/100:2*pi;
N    = length(t2);
r2   = zeros(N, 0);
for i=1:length(t2)
    r2(i) = a*(1-e^2)/(1+e*cos(t2(i)));
end

% Transfer orbit
a    = (7178 + 6678) / 2;       % semimajor axis
e    = 1 - (6678 / a);          % eccentricity
t3   = 0:pi/100:pi;
N    = length(t3);
r3   = zeros(N, 0);
for i=1:length(t3)
    r3(i) = a*(1-e^2)/(1+e*cos(t3(i)));
end

red = [244/255, 67/255, 54/255];
blue = [33/255, 150/255, 243/255];
%green = [118/255, 255/255, 3/255];
green = [76/255, 175/255, 80/255];

polarplot(0,0, 'ko')
hold on
polarplot(t1, r1,'Color', blue)
polarplot(t2, r2, 'Color', red)
polarplot(t3, r3, 'Color', green, 'LineStyle', 'none', 'Marker', '.')
hold off
axis off
