%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code by Lakshay%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all
clc


gravity = [-3.7114 0 0]';

mass_dry = 1505; 

mass_wet = 1905; 

alph = 4.53e-4; 

rho_1 = 4972; 

rho_2 = 13260; 

%Convex Linear Operators
E = [eye(3) zeros(3,4)];
F = [0 0 0 0 0 0 1];
E_u = [eye(3) zeros(3,1)];
E_v = [zeros(3,3) eye(3) zeros(3,1)];

e_1 = [1 0 0]';
e_2 = [0 1 0]';
e_3 = [0 0 1]';

gamma = pi*(4/180); 

S = [e_2';e_3'];

c = e_1/tan(gamma);

N = 55; %discrete time steps


position_0 = [1500 500 2000]'; 
velocity_0 = [-75 0 100]';

%optimal final time
tf_opt = 78.4; 
%time increment
time_step = tf_opt/N;

%Define the State transition matrices 
A_c = [zeros(3,3) eye(3) zeros(3,1);zeros(4,7)];
B_c = [[zeros(3,3);eye(3);0 0 0] [0 0 0 0 0 0 -alph]'];

A = expm(A_c*time_step); %continuous time A matrix

B = A*(eye(7)*time_step - A_c*time_step^2/2)*B_c; %continuous time B matrix

Lambda_k = zeros(N*7,4);
Lambda_k(1:7,1:4) = B;
for k = 2:N
    
    Lambda_k((k-1)*7+1:k*7,:) = A*Lambda_k((k-2)*7+1:(k-1)*7,:) + B; 
end

Psi_k = zeros(N*7,4*(N+1));
for k = 2:(N)
    %Next time step of gravities effect is dependent on the state
    %transition matrix and the previous time step
    Psi_k((k-1)*7+1:k*7,:) = A*Psi_k((k-2)*7+1:(k-1)*7,:);
    Psi_k((k-1)*7+1:k*7,((k*4-7):(4*k-4))) = B;
end

% Mass after the change of variables
z0 = log(mass_wet-alph*rho_2*time_step*(0:N)');

% Initial state vector
y0 = [position_0; velocity_0; log(mass_wet)];

s(1:N,7) = 0;
for i = 1:N
    s(i,:) = (7*i-6):(7*i);
end

cvx_begin
    variable eta((N+1)*4)
    variable y(N*7)

    % Objective function
    minimize(norm(y(end-6:end-4),2))

    subject to

    % Convexified thrust constraint
    for k = 0:N
        norm(E_u*eta(4*k+1:4*k+4), 2) <= eta(4*k+4);
    end

    % Thrust constraint 1
    eta(4) <= rho_2*exp(-z0(1)).*(1-(F*y0-z0(1)));
    rho_1*exp(-z0(1))*(1-(F*y0-z0(1))+0.5*(F*y0-z0(1)).^2) <= eta(4);

    for k = 1:N
        % Cone constraints
        norm(S*E*(y(s(k, :))-y(s(N, :))), 2)-c'*(E*(y(s(k, :)))) <= 0;

        % Thrust constraints
        eta(4*(k)+4) <= rho_2*exp(-z0(k+1)).*(1-(F*y(s(k, :))-z0(k+1)));

        rho_1*exp(-z0(k+1))*(1-(F*y(s(k, :))-z0(k+1))+...
            0.5*(F*y(s(k, :))-z0(k+1)).^2) <= eta(4*k+4);

        % System dynamics constraint
        y(s(k, :)) == A^k*y0+Lambda_k(s(k, :), :)*[gravity; 0]+...
                        Psi_k(s(k, :), :)*eta;
     end

    % Fuel mass constraint
    y(end) >= log(mass_dry);

    % Final height is 0 constraint
    y(end-6) == 0;

    % Final velocity constraint
    for i = 1:3
        y(end-i) == 0;
    end
cvx_end

% Converting output into manageable format
dist(1:3, N+1) = 0;
dist(1:3, 1) = position_0;
mass(1) = mass_wet;
for i = 1:N
    dist(1:3, i+1) = y((7*i-6):(7*i-4));
    mass(i+1) = y(7*i);
end



plot(dist(2,:),dist(3,:))
grid on 
xlabel('y (meters)')
ylabel('z (meters)')
title('Optimal Trajectory')
figure

plot(dist(3,:),dist(1,:))
grid on 
xlabel('z (meters)')
ylabel('x (meters)')
title('Optimal Trajectory')
figure

plot3(dist(3, :), dist(2, :), dist(1, :), '-')
grid on
xlabel('z (meters)')
ylabel('y (meters)')
zlabel('x (meters)')
title('3D Trajectory')
figure


plot(0:time_step:tf_opt,dist)
grid on
legend('x','y','z');
xlabel('Time (seconds)')
ylabel('Position (meters)')
title('Trajectory vs. Time')
figure


plot(0:time_step:tf_opt,exp(mass))
grid on
xlabel('Time (seconds)')
ylabel('Mass (kg)')
title('Mass vs. Time')
grid on