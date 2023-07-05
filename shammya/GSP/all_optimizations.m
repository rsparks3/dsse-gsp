% clc;
clear all;
% close all;
rng(1)   
define_constants;
INCLUDE_LOAD = 0;
ADD_NOISE = 1;
casestr = 'case85';
mpc_case = loadcase(casestr); % load case to get back to the original network
nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
G= graph(f,t);
%% Calcualtion for Y matrix
[Y_bus,Y_f,~] = makeYbus(mpc_case);
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

%% 
mpc_case_results = runpf(mpc_case); % run the power flow
number_of_eigenvectors = 18;
U_V_new = load(strcat(casestr,'_U_V','.mat')).U_V_new;
U_k = U_Y_new (:,1:number_of_eigenvectors);

[~,sensor_locations] = sensor_placement_matpower(mpc_case);
ammeter = load('ammeter.mat');
sensor_locations = [sensor_locations];
% sensor_locations = sensor_locations(sensor_locations~=60);
% sensor_locations = sensor_locations(sensor_locations~=74);
% % sensor_locations = sensor_locations(sensor_locations~=75);
% sensor_locations = sensor_locations(sensor_locations~=72);
% sensor_locations = sensor_locations(sensor_locations~=79);
% sensor_locations = sensor_locations(sensor_locations~=31);
% sensor_locations = sensor_locations(sensor_locations~=43);
% sensor_locations = sensor_locations(sensor_locations~=55);
% sensor_locations = sensor_locations(sensor_locations~=56);
% sensor_locations = sensor_locations(sensor_locations~=85);


% [~,sensor_locations] = sensor_placement_matpower(mpc_case);
% [~,sensor_locations]= optimal_placement_greedy_parfor(U_k,36,[1]);
% sensor_locations = 1:nb;
% sensor_locations = sensor_locations(1:ceil(length(sensor_locations)/2));
sensor_locations = sensor_locations';
number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end
figure
h = plot(G);
title(strcat('Number of Sensors:', string(number_of_sensors)))
highlight(h,sensor_locations,'Nodecolor','k');
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
cm = abs(c(ammeter_locations));
M_a =zeros(length(ammeter_locations),nb-1);
for r = 1:length(ammeter_locations)
    M_a (r,ammeter_locations(r)) = 1; % Put 1 in sensor locations
end

% Select measurements according to sensor locations
sm= s(sensor_locations);
vm = abs(vcomplex(sensor_locations));
w_3 = 0;

%% Formulation 1 
cvx_clear
cvx_solver mosek
cvx_begin quiet
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb)
    variable a_est(nb-1) complex
    minimize ( 100*norm(sm - M*s_est, 2) + 550*norm(vm.*vm-M*(v_est),2) + norm(cm.*cm-M_a*a_est,2) + ...
                        100*norm(sum(real(s_est)),1) + 100*norm(sum(imag(s_est)),1))  
    
    subject to
%         sm == M*s_est;
        tildeW == tildeW'; 
        s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
        v_est == diag(U_k*tildeW*(U_k)');
        a_est == diag(Y_f*U_k*tildeW*U_k'*Y_f');
        U_k(1,:)*(tildeW*U_k(1,:)') == 1;
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


[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
max_angle_diff = norm(angle(vcomplex) - voltage_angle,2);
max_voltage_diff = norm(abs(vcomplex) - voltage_mag,2);
fprintf('Maximum Angle Doff : %f & Magnitude Difference: %f in Formulation 1\n',max_angle_diff,max_voltage_diff);


%% Formulation 2

cvx_clear
cvx_solver mosek
cvx_begin quiet
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    variable a_est(nb-1) complex
    minimize ( 100*norm(sm - M*s_est, 2) + 1900*norm(vm.*vm-diag(M*U_k*tildeW*(U_k)'*M'),2) + norm(cm.*cm-M_a*a_est,2) + ...
                        00*norm(sum(real(s_est)),1) + 00*norm(sum(imag(s_est)),1))  
    
    subject to
%         sm == M*s_est;
        tildeW == tildeW'; 
        s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
        sum(real(s_est)) <=0.05; 
        sum(imag(s_est)) <=0.05; 
        a_est == diag(Y_f*U_k*tildeW*U_k'*Y_f');
        U_k(1,:)*(tildeW*U_k(1,:)') == 1;
%         real(s_est(1)) + sum(real(s_est(2:end))) <= 0.05;
%         imag(s_est(1)) + sum(imag(s_est(2:end))) <= 0.05;
%         real(s_est(1)) >= 0;
%         imag(s_est(1)) >= 0;
%         real(s_est(2:end)) <= 0;
%         imag(s_est(2:end)) <= 0;

        %         abs(s_est(1)) - sum(abs(s_est(2:end))) <=0.05; 
        

cvx_end


[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
max_angle_diff = norm(angle(vcomplex) - voltage_angle,2);
max_voltage_diff = norm(abs(vcomplex) - voltage_mag,2);
fprintf('Maximum Angle Doff : %f & Magnitude Difference: %f in Formulation 2\n',max_angle_diff,max_voltage_diff);




%% Formulation 3
cvx_clear
cvx_solver mosek
cvx_begin quiet
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb)
    variable a_est(nb-1) 
    minimize ( norm(sm - M*s_est, 2) + 25*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M_a*a_est,2))  
    
    subject to
%         sm == M*s_est;
        tildeW == tildeW'; 
        s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
        v_est == diag(U_k*tildeW*(U_k)');
        sum(real(s_est)) <=0.1; 
        sum(imag(s_est)) <=0.1; 
        a_est == diag(Y_f*U_k*tildeW*U_k'*Y_f');
        U_k(1,:)*(tildeW*U_k(1,:)') == 1;
%         real(s_est(1)) + sum(real(s_est(2:end))) <= 0.05;
%         imag(s_est(1)) + sum(imag(s_est(2:end))) <= 0.05;
%         real(s_est(1)) >= 0;
%         imag(s_est(1)) >= 0;
%         real(s_est(2:end)) <= 0;
%         imag(s_est(2:end)) <= 0;

        %         abs(s_est(1)) - sum(abs(s_est(2:end))) <=0.05; 
        

cvx_end


[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
max_angle_diff = norm(angle(vcomplex) - voltage_angle,2);
max_voltage_diff = norm(abs(vcomplex) - voltage_mag,2);
fprintf('Maximum Angle Diff : %f & Magnitude Difference: %f in Formulation 3\n',max_angle_diff,max_voltage_diff);

%%  Formulation 4
cvx_clear
cvx_solver mosek
cvx_begin quiet
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    variable a_est(nb-1) complex
    minimize ( 100*norm(sm - M*s_est, 2) + 1900*norm(vm.*vm-diag(M*U_k*tildeW*(U_k)'*M'),2) + norm(cm.*cm-M_a*a_est,2) + ...
                        100*norm(sum(real(s_est)),1) + 100*norm(sum(imag(s_est)),1))  
    
    subject to
%         sm == M*s_est;
        tildeW == tildeW'; 
        s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
        a_est == diag(Y_f*U_k*tildeW*U_k'*Y_f');
        U_k(1,:)*(tildeW*U_k(1,:)') == 1;
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


[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
max_angle_diff = norm(angle(vcomplex) - voltage_angle,2);
max_voltage_diff = norm(abs(vcomplex) - voltage_mag,2);
fprintf('Maximum Angle Doff : %f & Magnitude Difference: %f in Formulation 4\n',max_angle_diff,max_voltage_diff);



%% Formulation 5: adding current measurments 
cvx_clear
cvx_solver mosek
cvx_begin quiet
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb) complex
    variable a_est(nb-1) complex
    minimize ( 100*norm(sm - M*s_est, 2) + 2500*norm(vm.*vm-M*(v_est),2) + 300*norm(cm.*cm-M_a*a_est,2) + ...
                        100*norm(sum(real(s_est)),1) + 100*norm(sum(imag(s_est)),1))  
    
    subject to
%         sm == M*s_est;
        tildeW == tildeW'; 
        s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
        a_est == diag(Y_f*U_k*tildeW*U_k'*Y_f');
        U_k(1,:)*(tildeW*U_k(1,:)') == 1;
        v_est == diag(U_k*tildeW*(U_k)');
       

cvx_end
[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
max_angle_diff = norm(angle(vcomplex) - voltage_angle,2);
max_voltage_diff = norm(abs(vcomplex) - voltage_mag,2);
fprintf('Maximum Angle Diff : %f & Magnitude Difference: %f in Formulation 5\n',max_angle_diff,max_voltage_diff);


%%

figure 
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


