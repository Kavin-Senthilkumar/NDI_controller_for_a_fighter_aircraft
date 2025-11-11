


Y0 = [100; 4; 1; 0.1; 0.05; 0.08; 5; 0; 10; 0; 0; -200];

U = [10; 2; 0; 5; 2; 20];


%tspan = [0 10];

sim('F16.slx');

%options = odeset('RelTol', 1e-3, 'AbsTol', 1e-5);
%[t, y] = ode15s(@(t, y) F16_deq(t, y, U), tspan, Y0, options);

t = out.Y.time;
y1 = out.Y.Data(2,:);



plot(t,y1);