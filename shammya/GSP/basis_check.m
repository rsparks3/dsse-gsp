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
% 
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
number_of_eigenvectors = number_of_sensors;
number_of_eigenvectors = 8;
U_k = U_Y_new (:,1:number_of_eigenvectors);
s_noise = 0*random_distribution(nb); % A customized function to generate random values 
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
% Add complex noise
v_noise = 0*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;
vm = v(sensor_locations);
rhs = transpose(M) * vm;

f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
figure
G= graph(f,t);
subplot(2,2,[1 3])
h = plot(G);
highlight(h,sensor_locations,'Nodecolor','k');
title(strcat('Number of Sensors:', string(number_of_sensors), ' for case: ',casestr))

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
gamma = 0.0175;
v_new = (U_k)*pinv(M*U_k)*v(sensor_locations);
v_gamma = inv(transpose(M)*M + gamma*(U_k*U_k') - gamma * eye(size(nb,nb)) )*rhs;


subplot(2,2,2)
plot(abs(vcomplex))
hold on 
grid on 
plot((abs(v_new)))
hold on 
plot((abs(v_gamma)))
hold off
xlim([1 nb])
% plot(abs(v_recon_1))
% plot(abs(v_recon_2))
legend('true','Relaxation','Gamma')
xlabel('Bus Numbers')
title('Voltage Magnitude Plot Using Projection')


subplot(2,2,4)
plot(angle(vcomplex))
hold on 
grid on 
plot((angle(v_new)))
hold on 
plot((angle(v_gamma)))
hold off
% plot(abs(v_recon_1))
% plot(abs(v_recon_2))
legend('true','Relaxation','Gamma')
xlim([1 nb])
xlabel('Bus Numbers')
title('Voltage Angle Plot Using Projection')

%%

pi_1 = eye(size(nb,nb)) - (U_k*U_k') ;
pi_plus = pi_1(:,1);
pi_minus = pi_1(:,2:end);
M_bar = M(:,2:end);

gamma=0.0175;
condition_check = transpose(M_bar)*M_bar - gamma* transpose(pi_minus) * pi_minus ;
rcond(condition_check)
rhs = transpose(M_bar) * vm + gamma* transpose(pi_minus) * pi_plus;
v_bar = inv(condition_check) * rhs;

%%
voltage_error = [];
for gamma = 0.001:0.001:1

    cvx_clear
    % cvx_solver mosek
    cvx_begin quiet
        variable v_new(nb) complex
        minimize (norm(vm-M*v_new,2) + gamma*norm((eye(nb,nb) - (U_k*U_k'))*v_new,2))
        subject to 
            v_new(1) == 1;
        cvx_end
voltage_error = [voltage_error norm(v_new-vcomplex,2)];
end
disp(min(voltage_error))


gamma = 0.0175
    cvx_clear
    % cvx_solver mosek
    cvx_begin quiet
        variable v_new(nb) complex
        minimize (norm(vm-M*v_new,2) + gamma*norm((eye(nb,nb) - (U_k*U_k'))*v_new,2))
        subject to 
            v_new(1) == 1;
        cvx_end

subplot(2,2,2)
plot(abs(vcomplex))
hold on 
grid on 
plot((abs(v_new)))
hold on 
hold off
xlim([1 nb])
% plot(abs(v_recon_1))
% plot(abs(v_recon_2))
legend('true','Relaxation','Gamma')
xlabel('Bus Numbers')
title('Voltage Magnitude Plot Using Projection')


subplot(2,2,4)
plot(angle(vcomplex))
hold on 
grid on 
plot((angle(v_new)))
hold on 
hold off
% plot(abs(v_recon_1))
% plot(abs(v_recon_2))
legend('true','Relaxation','Gamma')
xlim([1 nb])
xlabel('Bus Numbers')
title('Voltage Angle Plot Using Projection')







