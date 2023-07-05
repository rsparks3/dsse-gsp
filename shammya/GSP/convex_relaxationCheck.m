clc;
clear all;
close all;
rng(1)   
define_constants;
INCLUDE_LOAD = 0;
casestr = 'case85.m';
mpc_case = loadcase(casestr);

mpc_case_results = runpf(mpc_case);
[Y_bus,~,~] = makeYbus(mpc_case);
Y_L = mpc_case.bus(:,[PD,QD])./mpc_case.baseMVA;
Y_L = Y_L(:,1)-1i*Y_L(:,2);
Y_DD = full(Y_bus)+INCLUDE_LOAD*diag(Y_L); % diagonal matrix 


%%
nb = length(mpc_case.bus);
number_of_sensors = nb;
number_of_eigenvectors = 40;
  % Get the number of buses
M =zeros(number_of_sensors,nb); % Initialize the selection matrix with zeros
sensor_locations = 1:nb;
for r = 1:number_of_sensors
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end

s_noise = 0*random_distribution(nb); % A customized function to generate random values 
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
% Add complex noise
v_noise = 0*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;


% Select measurements according to sensor locations
sm= s(sensor_locations);
vm = vcomplex(sensor_locations);

%%
cvx_solver mosek
cvx_begin 
    variable W(nb,nb) complex semidefinite
    minimize (1*norm(sm - M*conj(diag(Y_DD*W)), 2) + 1*norm(vm.*vm-diag(M*W*M'),2) ) ; 
    subject to
        W == W';

cvx_end
%%


f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);

G= graph(f,t);

voltage_angle= zeros(nb,1);
for i = 2:nb
    path = shortestpath(G,1,i);
    angle_sum = 0;
    for j = 1:length(path)-1
        angle_sum = angle_sum - angle(W(path(j),path(j+1)));
    end
    voltage_angle(i,1) = angle_sum;
    
end

figure()
subplot(211)
plot(abs(vm),'linewidth',1.5)
hold on 
grid on 
plot(sqrt(abs(diag(W))),'linewidth',1.5)
legend('Original','Convex Relaxation')
xlabel('Bus Number')
xlim([1 nb])
title('Voltage Magnitude')
ylabel('Voltage Magnitude (pu)')
hold off


subplot(212)
plot(mpc_case_results.bus(:,VA),'linewidth',1.5)
hold on 
plot(voltage_angle*180/pi,'linewidth',1.5)
grid on
hold off
xlabel('Bus Number')
ylabel('Voltage Angle (degree)')
xlim([1 nb])
xlabel('Bus Number')
legend('Original','Convex Relaxation')
title('Voltage Angle')









