% Ths code intends to do the community detection using the voltage
% covariace matrix as outlined in the paper "ON MODELING VOLTAGE PHASOR
% MEASUREMENTS AS GRAPH SIGNALS". 
% The work in the paper is done for tranmissions systems, this code
% investigates whether the same approach will work for distribution systems
% too

clc;
clear all;
close all;
rng('default')    % Force the random generator 
define_constants;
ADD_NOISE = 0 ; 
casestr = 'case85';
sensor_per_cluster = 4;
mpc_case = loadcase(casestr); % load case to get back to the original network
nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= digraph(f,t,edge_weights);
ncluster = 8;
C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6],'m',[.8 .2 .9],[.1 .4 .7]}; % Cell array of colors.
F= full(sparse(1:G.numedges,f,1,G.numedges,G.numnodes));
T= full(sparse(1:G.numedges,t,1,G.numedges,G.numnodes));

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
Yuk = Y_bus * U_Y_new;
% Yik=(F'*Y_f-T'*Y_t)./degree_vector'*U_Y_new;
Yik=(T'*Y_t+F'*Y_f)*U_Y_new;

%%
% load the covariance matrices that have been already generated
PMU_cov = load(strcat(casestr,'_covariance','.mat')).covariance;
[PMUV, PMUD]  = eigs(PMU_cov, ncluster, 'largestabs');
% cluster as two different dimensions...
[PMUidx, PMUC,sumP] = kmeans([real(PMUV), imag(PMUV)], ncluster);
ammeter_locations = [];
sensor_locations = [];
Adjacency_clusters_cov=zeros(length(PMUidx));
for i=1:length(PMUidx)
  edges = find(PMUidx==PMUidx(i)); %from the same group
  Adjacency_clusters_cov(i,edges')=1;
end
n_eigs = 24;
U_k = U_Y_new(:,1:n_eigs);
figure(4)
subplot(121)
h = plot(G);
title('Clustering using the covariance matrix as the Laplacian matrix')
for i = 1:ncluster
    nodes = find(PMUidx == i);
    [~,sl,al]= greedy_placement(U_k,sensor_per_cluster,[],nodes',G,f,t,3,Yuk(:,1:n_eigs),Yik(:,1:n_eigs));
    sensor_locations = [sensor_locations sl'];
    ammeter_locations = [ammeter_locations al];
    highlight(h,nodes,'NodeColor',C{i},'MarkerSize',6);
end
subplot(122)
h = plot(G);
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)
%%
G= graph(f,t,edge_weights);
%%
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


% highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)
%%
number_of_eigenvectors = n_eigs;
cvx_clear
cvx_solver mosek
cvx_begin quiet
cvx_precision high
variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
variable s_est(nb) complex
variable v_est(nb)
variable a_est(nb)
minimize (1*norm(sm - M*s_est, 2) + 10*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M*a_est,2))

subject to
    tildeW == tildeW';
    s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
    v_est == diag(U_k*tildeW*(U_k)');
    sum(real(s_est)) <=0.05;
    sum(imag(s_est)) <=0.05;
    a_est == diag(Y_DD*U_k*tildeW*U_k'*Y_DD');
    U_k(1,:)*(tildeW*U_k(1,:)') == 1;
cvx_end


[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
rmse_angle  = rms(angle(vcomplex)-voltage_angle) ;
rmse_voltage = rms(abs(vcomplex) - voltage_mag);
fprintf('optval:%f RMS Angle:%f, RMS Voltage: %f\n',cvx_optval,rmse_angle,rmse_voltage);

figure
subplot(211)
plot(voltage_mag)
hold on
plot(abs(vcomplex))
legend('Reconstructed','Original')
subplot(212)
plot(voltage_angle)
hold on
plot(angle(vcomplex))
legend('Reconstructed','Original')