clc;
clear all;
close all;

define_constants
casestr = 'case85';
opt = mpoption('verbose',0,'out.all',0);
mpc_case = loadcase(casestr);
ADD_NOISE = 1;
nb = length(mpc_case.bus);
[Y_bus,Y_f,Y_t] = makeYbus(mpc_case);
Y_DD = full(Y_bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= graph(f,t,edge_weights);

sensor_locations = 1:nb;


sensor_locations = sensor_locations';
sensor_locations = unique(sensor_locations);
number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end
mpc_case_results = runpf(mpc_case);
s_noise = ADD_NOISE*random_distribution(nb); % A customized function to generate random values 
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
% Add complex noise
v_noise = ADD_NOISE*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;

c = Y_DD*vcomplex;
c_noise = ADD_NOISE*random_distribution(nb);
c = c + c_noise + 1i*c_noise;

cm = abs(c(sensor_locations));
sm= s(sensor_locations);
vm = abs(vcomplex(sensor_locations));

cvx_clear
% cvx_precision high
% cvx_solver mosek
cvx_begin 
    variable W(nb,nb) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb)
    variable a_est(nb)
    minimize (1*norm(sm - M*s_est, 2) + 5*norm(vm.*vm-M*(v_est),2)  +  5*norm(cm.*cm-M*(a_est),2))

    subject to
    %         sm == M*s_est;
    W == W';
    s_est == conj(diag(Y_DD*W));
    v_est == diag(W);
    a_est == diag(Y_DD*W*Y_DD');
cvx_end
voltage_angle= zeros(nb,1);
for i = 2:nb
    path = shortestpath(G,1,i);
    angle_sum = 0;
    for j = 1:length(path)-1
        angle_sum = angle_sum - angle(W(path(j),path(j+1)));
    end
    voltage_angle(i,1) = angle_sum;
    
end

voltage_mag = sqrt(abs(diag(W)));
figure(2) 
subplot(2,1,1)
h = plot(G);
plot(abs(vcomplex),'linewidth',1.25)
hold on 
grid on 
plot(voltage_mag,'linewidth',1.25)
legend('true','Relaxation')
ylabel('Voltage magnitude (pu)')
xlabel('Bus Numbers')
xlim([1 nb])
title('Reconstruction of Vm');
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
hold off

subplot(2,1,2)
plot(angle(vcomplex)*180/pi,'linewidth',1.25)
hold on 
grid on 
plot(voltage_angle*180/pi,'linewidth',1.25)
legend('true','Relaxation')
ylabel('Voltage Angle (Degree)')
xlabel('Bus Numbers')
xlim([1 nb])
title('Reconstruction of Va');
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
hold off

