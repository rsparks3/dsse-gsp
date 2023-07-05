% clc;
clear all;
% close all;
% rng(1)   
define_constants;
INCLUDE_LOAD = 0;
ADD_NOISE = 1;
casestr = 'case85';
mpc_case = loadcase(casestr); % load case to get back to the original network
nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= graph(f,t,edge_weights);
F= full(sparse(1:nb-1,f,1,nb-1,nb));
T= full(sparse(1:nb-1,t,1,nb-1,nb));
%% Calcualtion for Y matrix
[Y_bus,Y_f,Y_t] = makeYbus(mpc_case);
Y_DD = full(Y_bus);

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
mpc_case_results = runpf(mpc_case); % run the power flow
%% 

number_of_eigenvectors = 24;
U_V_new = load(strcat(casestr,'_U_V','.mat')).U_V_new;
U_k = U_Y_new (:,1:number_of_eigenvectors);

[~,sensor_locations] = sensor_placement_matpower(mpc_case);
ammeter = load('ammeter.mat');

sensor_to_exclude=[76 56 43 22 85 77 22 47 11 13 80 69 40 45 38 67 37 23 66 84 82 78 15 11 54 55 51 18 20 71 75 16 39 34 36 42 74 79 72 59 11 62];
for s_e = 1:length(sensor_to_exclude)
    sensor_locations = sensor_locations(sensor_locations~=sensor_to_exclude(s_e));
end
sensor_locations = [sensor_locations [32 12 34 5 68 27 54 9 48 64 16 8 25 58 22 19 22 3]];


% [~,sensor_locations] = sensor_placement_matpower(mpc_case);

% sensor_locations = 1:nb;
% sensor_locations = sensor_locations(1:ceil(length(sensor_locations)/1));
% deg_ranks = centrality(G,'degree');
% sensor_locations = find(deg_ranks>=3);
% deg_ranks_closeness = centrality(G,'eigenvector','Importance',G.Edges.Weight);
% sensor_locations = find(deg_ranks_closeness>0.00001);
% [sorted_rank, sort_order] = sort(deg_ranks_closeness,'descend');
% sensor_locations = sort_order(1:19);
% [~,sensor_locations]= optimal_placement_greedy_parfor(U_k,30,[2]);
sensor_locations = sensor_locations';
sensor_locations = unique(sensor_locations);
number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end

% Choosing the first selected eigenvectors

s_noise = ADD_NOISE*random_distribution(nb); % A customized function to generate random values 
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
% Add complex noise
v_noise = ADD_NOISE*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;

c = Y_f*vcomplex;
c_noise = ADD_NOISE*random_distribution(nb-1);
c = c + c_noise + 1i*c_noise;
ammeter_locations  = load('ammeter.mat').ammeter_location;
ammeter_locations = ammeter_locations(ammeter_locations~=45);
ammeter_locations = ammeter_locations(ammeter_locations~=46);
ammeter_locations = ammeter_locations(ammeter_locations~=13);
ammeter_locations = ammeter_locations(ammeter_locations~=61);
ammeter_locations = ammeter_locations(ammeter_locations~=78);
ammeter_locations = ammeter_locations(ammeter_locations~=81);
ammeter_locations = ammeter_locations(ammeter_locations~=76);
ammeter_locations = ammeter_locations(ammeter_locations~=74);
ammeter_locations = ammeter_locations(ammeter_locations~=11);
ammeter_locations = ammeter_locations(ammeter_locations~=31);
ammeter_locations = [ammeter_locations 48 21 43 53 79 22 39];
cm = abs(c(ammeter_locations));
M_a =zeros(length(ammeter_locations),nb-1);
for r = 1:length(ammeter_locations)
    M_a (r,ammeter_locations(r)) = 1; % Put 1 in sensor locations
end

% Select measurements according to sensor locations
sm= s(sensor_locations);
vm = abs(vcomplex(sensor_locations));
w_3 = 0;

figure(1)
h = plot(G);
title(strcat('Number of Sensors:', string(number_of_sensors)))
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)
cvx_clear
cvx_solver mosek
cvx_precision high
cvx_begin quiet
    cvx_precision high
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb)
    variable a_est(nb-1) 
    minimize (1*norm(sm - M*s_est, 2) + 30*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M_a*a_est,2) )  
    
    subject to
%         sm == M*s_est;
        tildeW == tildeW'; 
        s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
        v_est == diag(U_k*tildeW*(U_k)');
        sum(real(s_est)) <=0.1; 
        sum(imag(s_est)) <=0.1; 
        a_est == diag(Y_f*U_k*tildeW*U_k'*Y_f');
        U_k(1,:)*(tildeW*U_k(1,:)') == 1;
cvx_end


[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
max_angle_diff = norm(angle(vcomplex) - voltage_angle,2);
max_voltage_diff = norm(abs(vcomplex) - voltage_mag,2);
fprintf('Maximum Angle Diff : %f & Magnitude Difference: %f and optval:%f\n',max_angle_diff,max_voltage_diff,cvx_optval);

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

%%

% cvx_clear
% cvx_begin quiet
%     cvx_precision high
%     variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
%     variable s_est(nb) complex
%     variable v_est(nb)
%     variable a_est(nb-1) 
%     minimize (norm(sm - M*s_est, 2) + 30*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M_a*a_est,2) + 0.0*trace(U_k*tildeW*U_k'))  
%     
%     subject to
% %         sm == M*s_est;
%         tildeW == tildeW'; 
%         s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
%         v_est == diag(U_k*tildeW*(U_k)');
%         sum(real(s_est)) <=0.1; 
%         sum(imag(s_est)) <=0.1; 
%         a_est == diag(Y_f*U_k*tildeW*U_k'*Y_f');
%         U_k(1,:)*(tildeW*U_k(1,:)') == 1;
% cvx_end
% 
% 
% [V,D] = eigs(U_k*tildeW*U_k',1,'largestabs');
% delta = 1-D/trace(U_k*tildeW*U_k');
% w = 0;
% cvx_clear
% cvx_solver mosek
% cvx_precision high
% cvx_begin quiet
%     cvx_precision high
%     variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
%     variable s_est(nb) complex
%     variable v_est(nb)
%     variable a_est(nb-1) 
%     minimize (norm(sm - M*s_est, 2) + 30*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M_a*a_est,2) + 0.0*trace(U_k*tildeW*U_k'))  
%     
%     subject to
% %         sm == M*s_est;
%         tildeW == tildeW'; 
%         s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
%         v_est == diag(U_k*tildeW*(U_k)');
%         sum(real(s_est)) <=0.1; 
%         sum(imag(s_est)) <=0.1; 
%         a_est == diag(Y_f*U_k*tildeW*U_k'*Y_f');
%         U_k(1,:)*(tildeW*U_k(1,:)') == 1;
%         V'*U_k*tildeW*U_k'*V >= 1*trace(U_k*tildeW*U_k');
% cvx_end
% M = 24;
% cvx_clear
% cvx_begin
%     variable location(nb,1)
%     minimize (-lambda_min(transpose((U_k)) * diag(location)* (U_k) ) )
%     subject to
%         0 <= location <= 1;
%         ones(1,nb) * location == M ; 
% cvx_end

c = T'*Y_t*vcomplex - F'*Y_f*vcomplex;
c_noise = ADD_NOISE*random_distribution(nb);
c = c + c_noise + 1i*c_noise;
cm = abs(c(sensor_locations));

figure(1)
h = plot(G);
title(strcat('Number of Sensors:', string(number_of_sensors)))
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
% highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)
cvx_clear
cvx_solver mosek
cvx_precision high
cvx_begin 
    cvx_precision high
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb)
    variable a_est(nb) 
    minimize (1*norm(sm - M*s_est, 2) + 30*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M*a_est,2) )  
    
    subject to
%         sm == M*s_est;
        tildeW == tildeW'; 
        s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
        v_est == diag(U_k*tildeW*(U_k)');
        sum(real(s_est)) <=0.1; 
        sum(imag(s_est)) <=0.1; 
        a_est == diag((T'*Y_t-F'*Y_f)*U_k*tildeW*U_k'*(T'*Y_t-F'*Y_f)');
        U_k(1,:)*(tildeW*U_k(1,:)') == 1;
cvx_end


[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
max_angle_diff = norm(angle(vcomplex) - voltage_angle,2);
max_voltage_diff = norm(abs(vcomplex) - voltage_mag,2);
fprintf('Maximum Angle Diff : %f & Magnitude Difference: %f and optval:%f\n',max_angle_diff,max_voltage_diff,cvx_optval);

figure(3) 
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