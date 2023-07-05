function [] = run_undersampled_optimization_buscurrent(G,casestr,sensor_locations)
define_constants;
% rng(1);
ADD_NOISE = 1;

opt = mpoption('verbose',0,'out.all',0);
mpc_case = loadcase(casestr); % load case to get back to the original network
nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
F= full(sparse(1:G.numedges,f,1,G.numedges,G.numnodes));
T= full(sparse(1:G.numedges,t,1,G.numedges,G.numnodes));
degree_vector = degree(G,1:G.numnodes);
degree_vector(degree_vector>1) = 2 ;
%% Calcualtion for Y matrix
[Y_bus,Y_f,Y_t] = makeYbus(mpc_case);
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
mpc_case_results = runpf(mpc_case,opt); % run the power flow
%% 

number_of_eigenvectors = 24;
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

% c = (F'*Y_f-T'*Y_t)./degree_vector*vcomplex;
CM = Y_bus;
c = CM*vcomplex;
c_noise = ADD_NOISE*random_distribution(nb);
c = c + c_noise + 1i*c_noise;
cm = abs(c(sensor_locations));

% Select measurements according to sensor locations
sm= s(sensor_locations);
vm = abs(v(sensor_locations));


figure(101)
h = plot(G);
title(strcat('Number of Sensors:', string(number_of_sensors)))
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
% highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)
cvx_clear
cvx_solver mosek
cvx_begin quiet
    cvx_precision high
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb)
    variable a_est(nb) 
    minimize (norm(sm - M*s_est, 2) + 15*norm(vm.*vm-M*(v_est),2) + 10*norm(cm.*cm-M*a_est,2))  
    
    subject to
%         sm == M*s_est;
        tildeW == tildeW'; 
        s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
        v_est == diag(U_k*tildeW*(U_k)');
        sum(real(s_est)) <=0.1; 
        sum(imag(s_est)) <=0.1; 
        a_est == diag(CM*U_k*tildeW*U_k'*CM');
        U_k(1,:)*(tildeW*U_k(1,:)') == 1;
cvx_end


[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
rms_angle_diff = rms(angle(vcomplex) - voltage_angle);
rms_voltage_diff = rms(abs(vcomplex) - voltage_mag);
fprintf('RMS Angle Diff : %f & Magnitude Difference: %f and optval:%f\n',rms_angle_diff,rms_voltage_diff,cvx_optval);

[max_angle_diff,bus_angle] = max(abs(angle(vcomplex) - voltage_angle));
max_voltage_diff = max(abs(abs(vcomplex) - voltage_mag));
fprintf('Maximum Angle Diff : %f & Magnitude Difference: %f and optval:%f\n',180/pi*max_angle_diff,max_voltage_diff,cvx_optval);



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

