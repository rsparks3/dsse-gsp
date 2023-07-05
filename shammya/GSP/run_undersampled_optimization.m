function [] = run_undersampled_optimization(G,casestr,sensor_locations,ammeter_locations)
define_constants;
% rng(1);
ADD_NOISE = 1;
mpc_case = loadcase(casestr); % load case to get back to the original network
nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);

%% Calcualtion for Y matrix
[Y_bus,Y_f,~] = makeYbus(mpc_case);
Y_DD = full(Y_bus);

[U_Y, Lambda_Y] = eig(Y_DD);%./mpcint.baseMVA);
[~, ind_Y]=sort(diag(Lambda_Y),'ascend');
U_Y_arr = U_Y(:,ind_Y');

% change into complex orthogonal matrices:
r_vect = zeros(nb,1);
for i=1:nb
        r_vect(i) = 1./sqrt(transpose(U_Y_arr(:,i))*U_Y_arr(:,i));
end
    
U_Y_new = U_Y_arr*diag(r_vect);
mpc_case_results = runpf(mpc_case); % run the power flow
%% 
sensor_locations = [44 46 83 19 41 61 53 50 73 65 70 14 29 33 26 7];
number_of_eigenvectors = 23;
U_k = U_Y_new (:,1:number_of_eigenvectors);

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
cm = abs(c(ammeter_locations));
M_a =zeros(length(ammeter_locations),nb-1);
for r = 1:length(ammeter_locations)
    M_a (r,ammeter_locations(r)) = 1; % Put 1 in sensor locations
end

% Select measurements according to sensor locations
sm= s(sensor_locations);
vm = abs(v(sensor_locations));


figure(101)
h = plot(G);
title(strcat('Number of Sensors:', string(number_of_sensors)))
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)
cvx_clear
cvx_solver mosek
cvx_begin quiet
    cvx_precision high
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb)
    variable a_est(nb-1) 
    minimize (norm(sm - M*s_est, 2) + 5*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M_a*a_est,2) + 0.0*trace(U_k*tildeW*U_k'))  
    
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

figure(102) 
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




end

