


Y0 = [100; 4; 1; 0.1; 0.05; 0.08; 5; 0; 10; 0; 0; -200];

U = [10; 2; 0; 5; 2; 20];


tspan = [0 10];


options = odeset('RelTol', 1e-3, 'AbsTol', 1e-5);
[t, y] = ode15s(@(t, y) F16_deq(t, y, U), tspan, Y0, options);



figure;
plot3(y(:,10), y(:,11), y(:,12));
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
title('F-16 3D Trajectory');
grid on;

% Linear Velocities Plot
figure;
plot(t, y(:,1), t, y(:,2), t, y(:,3));
xlabel('Time (s)');
ylabel('Linear Velocity (m/s)');
title('Linear Velocities');
legend('u', 'v', 'w');
grid on;

% Angular Velocities Plot
figure;
plot(t, y(:,4), t, y(:,5), t, y(:,6));
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Angular Velocities');
legend('p', 'q', 'r');
grid on;

% Euler Angles Plot
figure;
plot(t, y(:,7), t, y(:,8), t, y(:,9));
xlabel('Time (s)');
ylabel('Euler Angles (deg)');
title('Euler Angles');
legend('phi', 'theta', 'psi');
grid on;