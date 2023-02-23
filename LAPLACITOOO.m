% Define initial conditions
t = 0:0.1:1000;
a = [421700, 671100, 1070400]; % semi-major axes in meters
e = [0.004, 0.009, 0.001]; % eccentricities
m = [8.93e22, 4.8e22, 1.48e23]; % masses in kg
G = 6.674e-11; % gravitational constant

% Calculate orbital periods and angular frequencies
P = 2*pi*sqrt(a.^3./(G*(m(1)+m(2)+m(3))));
omega = 2*pi./P;

% Define initial angles
theta = [0, 0, 0];

% Define function to calculate derivatives of theta (im using ode45)
dtheta_dt = @(t,theta) omega - [1, -2, 1].*omega.*cos([theta(1)-theta(2), theta(2)-theta(3), theta(1)-theta(3)]);

% Integrate the equations of motion
[t,theta] = ode45(dtheta_dt, t, theta);

% Convert angles to degrees
theta = rad2deg(theta);

% Plot the results
figure;
plot(t,theta(:,1),'g-',t,theta(:,2),'b-',t,theta(:,3),'r-');
xlabel('Time (s)');
ylabel('Angle (degrees)');
legend('Io','Europa','Ganymede');
title('Laplace Resonance Model');
