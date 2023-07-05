clc;
clear all;
close all;
rng(1)   
define_constants;
INCLUDE_LOAD = 0;
casestr = 'case85.m';
mpc_case = loadcase(casestr);
nb = length(mpc_case.bus);

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


% number_of_eigenvectors = ceil(1*nb);
  % Get the number of buses
 % Initialize the selection matrix with zeros
[sensor_locations,~] = sensor_placement_matpower(mpc_case);
% [~,sensor_locations] = sensor_placement_matpower(mpc_case);
sensor_locations = [sensor_locations 2];
% Exclude the slack bus measurement
sensor_locations = sensor_locations(sensor_locations~=1);
% sensor_locations = sensor_locations(sensor_locations~=74);
% % sensor_locations = sensor_locations(sensor_locations~=75);
% sensor_locations = sensor_locations(sensor_locations~=72);
% sensor_locations = sensor_locations(sensor_locations~=79);
% sensor_locations = sensor_locations(sensor_locations~=31);
% sensor_locations = sensor_locations(sensor_locations~=43);
% sensor_locations = sensor_locations(sensor_locations~=55);
% sensor_locations = sensor_locations(sensor_locations~=56);
% sensor_locations = sensor_locations(sensor_locations~=85);
% sensor_locations = 1:80;
sensor_locations = unique(sensor_locations);
number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end

mpc_case = loadcase(casestr); % load case to get back to the original network
mpc_case_results = runpf(mpc_case); % run the power flow
% Choosing the first selected eigenvectors

number_of_eigenvectors = 8; %chosen arbitarily
U_k = U_Y_new (:,1:number_of_eigenvectors);
s_noise = 0;%*random_distribution(nb); % A customized function to generate random values 
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
% Add complex noise
v_noise = 0;%*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;
vm = v(sensor_locations);
v_mag = abs(vm);


f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
% figure
G= graph(f,t);


%%
% voltage_error = [];
% for gamma = 0.001:0.001:1
%     cvx_clear
%     cvx_begin quiet
%         variable v_new(nb) complex
%         minimize (norm(vm-M*v_new,2) + gamma*norm((eye(nb,nb) - (U_k*U_k'))*v_new,2))
%         subject to 
%             v_new(1) == 1;
%     cvx_end
% voltage_error = [voltage_error norm(v_new-vcomplex,2)];
% end
% [~,min_gamma_loc] = min(voltage_error);
% gamma = 0.001:0.001:1;
% min_gamma = gamma(min_gamma_loc);

min_gamma = 0.87;
cvx_clear
cvx_solver SDPT3
cvx_begin
    variable v_new_cvx(nb) complex
    minimize (norm(vm-M*v_new_cvx,2) + min_gamma*norm((eye(nb,nb) - (U_k*U_k'))*v_new_cvx,2))
    subject to 
        v_new_cvx(1) == 1;
cvx_end

%%
% gamma_range = 10000:1:30000;
% norm_error = zeros(1,length(gamma_range));
% gamma = 0.0001;
% 
% % 
% for i = 1:length(gamma_range)
%     v_gamma = (transpose(M)*M + gamma_range(i)*(U_k*U_k') - gamma_range(i) * eye(size(nb,nb)))\rhs;
%     norm_error(1,i) = norm(v_gamma-vcomplex,2);
% end
%%
% gamma = 0.0175;
% v_new = (U_k)*pinv(M*U_k)*v(sensor_locations);
% v_gamma = inv(transpose(M)*M + gamma*(U_k*U_k') - gamma * eye(size(nb,nb)) )*rhs;




%% Using the gamma formulation
pi_1 = eye(size(nb,nb)) - (U_k*U_k') ;

pi_plus = pi_1(:,1);
pi_minus = pi_1(:,2:end);
M_bar = M(:,2:end);

cvx_clear
%cvx_solver SDPT3
cvx_begin
    variable v_new_cvx_1(nb-1) complex
    minimize ((norm(vm-M_bar*v_new_cvx_1,2)) + min_gamma*(norm(pi_plus + pi_minus*v_new_cvx_1,2)))
cvx_end
v_new_cvx_1 = [1;v_new_cvx_1];
gamma=min_gamma;
%gamma = 10;
condition_check = transpose(M_bar)*M_bar + gamma* transpose(pi_minus) * pi_minus ;
rcond(condition_check)
rhs = transpose(M_bar) * vm - gamma* transpose(pi_minus) * pi_plus;
v_bar = pinv(condition_check)* rhs;
v_bar = [1;v_bar];



% condition_check_l = pinv(transpose(M)*M + gamma * transpose(pi_1)*pi_1);
% rhs_l = transpose(M)*vm;
% % Lagrange multiplier for equality constraint-affects everyone..
% RR = condition_check_l*rhs_l; r_eq = (RR(1)-1)./(condition_check_l(1,1));
% 
% v_bar_new = condition_check_l*(rhs_l - (r_eq)*[1;zeros(nb-1,1)]);


%% Using the pseudo inverse formula after modifying for excluding the slack bus
v_new = (U_k(2:end,:))*pinv(M_bar*U_k(2:end,:))*v(sensor_locations);
v_new = [1;v_new];
%%
subplot(2,1,1)
plot(abs(vcomplex),'linewidth',1.5)
hold on 
grid on 
plot((abs(v_new_cvx_1)),'linewidth',1.5)
hold on 
plot((abs(v_bar)),'linewidth',1.5)
hold on 
plot((abs(v_new)),'linewidth',1.5)
xlim([1 nb])
legend('true','cvx-gamma','gradient-gamma','pseudo-inverse')
xlabel('Bus Numbers')
title('Voltage Magnitude Plot Using Projection')


subplot(2,1,2)
plot(angle(vcomplex),'linewidth',1.5)
hold on 
grid on 
plot((angle(v_new_cvx_1)),'linewidth',1.5)
hold on 
plot((angle(v_bar)),'linewidth',1.5)
hold on 
plot((angle(v_new)),'linewidth',1.5)
hold off
legend('true','cvx-gamma','gradient-gamma','pseudo-inverse')
xlim([1 nb])
xlabel('Bus Numbers')
title('Voltage Angle Plot Using Projection')


%%
% cvx_clear
% cvx_solver mosek
% cvx_begin
%     variable W(nb,nb) complex semidefinite
%     minimize (norm(conj(s(sensor_locations))- diag(M*Y_DD*W*M'),2) + norm(v_mag-diag(M*W*M'),2) + norm(pi_1*W*pi_1,2))
%     subject to 
%         W == W';
%         W(1,1) == 1;
%     
% cvx_end
% 
% figure 
% plot(abs(vcomplex),'linewidth',1.5)
% hold on 
% plot(abs(sqrt(diag(W))),'linewidth',1.5)
% hold off 









