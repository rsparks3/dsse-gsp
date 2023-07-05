clc;
clear all;
close all;
rng(1)   
define_constants;
INCLUDE_LOAD = 0;
casestr = 'case85.m';
mpc_case = loadcase(casestr); % load case to get back to the original network
nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
G= graph(f,t);

%% Calcualtion for Y matrix
[Y_bus,~,~] = makeYbus(mpc_case);
Y_L = mpc_case.bus(:,[PD,QD])./mpc_case.baseMVA;
Y_L = Y_L(:,1)-1i*Y_L(:,2);
Y_DD = full(Y_bus)+INCLUDE_LOAD*diag(Y_L); % diagonal matrix 

[U_Y, Lambda_Y] = eig(Y_DD);%./mpcint.baseMVA);
[val_Y, ind_Y]=sort(diag(Lambda_Y),'ascend');
U_Y_arr = U_Y(:,ind_Y');

% change into complex orthogonal matrices:
r_vect = zeros(nb,1);
for i=1:nb
        r_vect(i) = 1./sqrt(transpose(U_Y_arr(:,i))*U_Y_arr(:,i));
end
    
U_Y_new = U_Y_arr*diag(r_vect);
fprintf('Max Difference in estimate for Y:%f\n',max(abs(U_Y_new * diag(val_Y) * transpose(U_Y_new) - Y_DD),[],'all'));

%% 


mpc_case_results = runpf(mpc_case); % run the power flow
number_of_eigenvectors = 8;
U_k = U_Y_new (:,1:number_of_eigenvectors);

% number_of_eigenvectors = ceil(1*nb);
  % Get the number of buses
 % Initialize the selection matrix with zeros
[sensor_locations,~] = sensor_placement_matpower(mpc_case);

sensor_locations = [sensor_locations 13 9 73 48 67];

sensor_locations = sensor_locations(sensor_locations~=60);
sensor_locations = sensor_locations(sensor_locations~=74);
% sensor_locations = sensor_locations(sensor_locations~=75);
sensor_locations = sensor_locations(sensor_locations~=72);
sensor_locations = sensor_locations(sensor_locations~=79);
sensor_locations = sensor_locations(sensor_locations~=31);
sensor_locations = sensor_locations(sensor_locations~=43);
sensor_locations = sensor_locations(sensor_locations~=55);
sensor_locations = sensor_locations(sensor_locations~=56);
sensor_locations = sensor_locations(sensor_locations~=85);


% [~,sensor_locations] = sensor_placement_matpower(mpc_case);
% [~,sensor_locations]= optimal_placement_greedy_parfor(U_k,36,[1]);
% sensor_locations = 1:nb;
sensor_locations = sensor_locations';
number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end


% Choosing the first selected eigenvectors

s_noise = 0*random_distribution(nb); % A customized function to generate random values 
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
% Add complex noise
v_noise = 0*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;

%%
% Select measurements according to sensor locations
sm= s(sensor_locations);
vm = abs(vcomplex(sensor_locations));
w_3 = 0;
% 
% cvx_solver mosek
% cvx_begin 
%     variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
%     variable s_est(nb) complex
%     variable v(nb) 
%     minimize ( 1*norm(s_est - conj(diag(Y_DD*U_k*tildeW*U_k')), 2) + 7.5*norm(v-diag(U_k*tildeW*transpose(U_k)),2) )  
%     
%     subject to
% %         sm == M*s_est;
%         sm == M*s_est;
%         tildeW == tildeW'; 
%         vm == M*v;
%         v(1) == 1;
%         v(2:end) >= 0;
%         sum(real(s_est)) <=0.05; 
%         sum(imag(s_est)) <=0.05; 
% %         real(s_est(1)) + sum(real(s_est(2:end))) <= 0.05;
% %         imag(s_est(1)) + sum(imag(s_est(2:end))) <= 0.05;
% %         real(s_est(1)) >= 0;
% %         imag(s_est(1)) >= 0;
% %         real(s_est(2:end)) <= 0;
% %         imag(s_est(2:end)) <= 0;
% 
%         %         abs(s_est(1)) - sum(abs(s_est(2:end))) <=0.05; 
%         
% 
% cvx_end


cvx_solver mosek
cvx_begin 
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    minimize ( 100*norm(sm - M*s_est, 2) + 750*norm(vm.*vm-diag(M*U_k*tildeW*transpose(U_k)*M'),2) + 100*norm(sum(real(s_est)),1) + 100*norm(sum(imag(s_est)),1))  
    
    subject to
%         sm == M*s_est;
        tildeW == tildeW'; 
        s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
%         sum(real(s_est)) <=0.05; 
%         sum(imag(s_est)) <=0.05; 
%         real(s_est(1)) + sum(real(s_est(2:end))) <= 0.05;
%         imag(s_est(1)) + sum(imag(s_est(2:end))) <= 0.05;
%         real(s_est(1)) >= 0;
%         imag(s_est(1)) >= 0;
%         real(s_est(2:end)) <= 0;
%         imag(s_est(2:end)) <= 0;

        %         abs(s_est(1)) - sum(abs(s_est(2:end))) <=0.05; 
        

cvx_end

%%
W = U_k * tildeW * U_k';
voltage_angle= zeros(nb,1);
for i = 2:nb
    path = shortestpath(G,1,i);
    angle_sum = 0;
    for j = 1:length(path)-1
        angle_sum = angle_sum - angle(W(path(j),path(j+1)));
    end
    voltage_angle(i,1) = angle_sum;
    
end

ypred = load('ypred.mat').loc;
ytest = load('ytest.mat').loc;
vm_predict = sqrt(abs(diag(W)));

figure 
subplot(2,2,[1 3])
h = plot(G);
title(strcat('Number of Sensors:', string(number_of_sensors)))
highlight(h,sensor_locations,'Nodecolor','k');
subplot(222)
plot(abs(vcomplex),'linewidth',1.25)
hold on 
grid on 
plot(vm_predict,'linewidth',1.25)
hold on 
% plot(ypred(1:85),'linewidth',1.25);
% plot(abs(v_recon_1))
% plot(abs(v_recon_2))
legend('true','Relaxation','ANN')
ylabel('Voltage magnitude (pu)')
xlabel('Bus Numbers')
xlim([1 nb])
title(strcat('Reconstruction of Vm with w_3 = ', string(w_3)))
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
hold off





subplot(224)
plot(angle(vcomplex)*180/pi,'linewidth',1.5)
hold on 
plot(voltage_angle*180/pi,'linewidth',1.5)
hold on 
% plot(ypred(86:end),'linewidth',1.25);
grid on
hold off
xlabel('Bus Number')
ylabel('Voltage Angle (degree)')
xlim([1 nb])
xlabel('Bus Number')
xlim([1 nb])
legend('Original','Convex Relaxation','ANN')
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
title(strcat('Reconstruction of Va with w_3 = ', string(w_3)))

max_angle_diff = norm(angle(vcomplex) - voltage_angle,2);
max_voltage_diff = norm(abs(vcomplex) - vm_predict,2);
fprintf('Maximum Angle Difference: %f \n',max_angle_diff);
fprintf('Maximum Magnitude Difference: %f \n',max_voltage_diff);

%%

v_predict = vm_predict.*exp(1j*voltage_angle); % Predicting the complex voltages
s_predict = diag(v_predict) * conj(Y_bus*v_predict);
s_test = diag(vcomplex) * conj(Y_bus*vcomplex);

figure

subplot(211)
plot(real(s_test(sensor_locations)));
hold on 
grid on 
plot(real(s_predict(sensor_locations)))
% xlim([1 nb]) 
xlabel('Bus Numbers')
ylabel('Real Power Injection')
legend('Original Value at S','Predicted Value')
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;

subplot(212)
plot(imag(s_test(sensor_locations)));
hold on 
grid on 
plot(imag(s_predict(sensor_locations)))
% xlim([1 nb])
xlabel('Bus Numbers')
ylabel('Reactive Power Injection')
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
legend('Original Value','Predicted Value')


%%
non_sensor_locations = setdiff(1:85,sensor_locations);

figure

subplot(211)
plot(real(s_test(non_sensor_locations)));
hold on 
grid on 
plot(real(s_predict(non_sensor_locations)))
% xlim([1 nb]) 
xlabel('Bus Numbers')
ylabel('Real Power Injection')
legend('Original Value at non-S','Predicted Value')
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;

subplot(212)
plot(imag(s_test(non_sensor_locations)));
hold on 
grid on 
plot(imag(s_predict(non_sensor_locations)))
% xlim([1 nb])
xlabel('Bus Numbers')
ylabel('Reactive Power Injection')
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
legend('Original Value','Predicted Value')

%%
figure
subplot(211)
plot(real(sm),'linewidth',1.5);
hold on
plot(real(M*s_est),'linewidth',2);
legend('sm','Ms-est')

subplot(212)
plot(imag(sm),'linewidth',1.5);
hold on
plot(imag(M*s_est),'linewidth',2);
legend('sm','Ms-est')


%%
[U_w,S_w,~] = svd(W);
v_tilde = sqrt(S_w(1,1)) * U_w(:,1); 

figure 
plot(abs(vcomplex));
hold on 
plot(abs(v_tilde))

figure 
plot(angle(vcomplex));
hold on 
plot(angle(v_tilde))

[U_w,S_w] = eigs(W,1);
% v_tilde = sqrt(S_w(1,1)) * U_w(:,1); 
v_tilde = sqrt(S_w) * U_w(:,1); 

figure 
plot(abs(vcomplex));
hold on 
plot(abs(v_tilde))

figure 
plot(angle(vcomplex));
hold on 
plot(angle(v_tilde))

















