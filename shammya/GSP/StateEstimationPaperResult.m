clc;
clear all;
close all;
define_constants;
casestr = 'case85';
ADD_NOISE = 01;
opt = mpoption('verbose',0,'out.all',0);
mpc_case = loadcase(casestr); % load case to get back to the original network
base_case = mpc_case;
bmult=1;
PD_case = base_case.bus(:,PD);
QD_case = base_case.bus(:,QD);

nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= graph(f,t,edge_weights);

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
mpc_case_results = runpf(mpc_case,opt); % run the power flow

number_of_eigenvectors = nb ;
U_k_gft = U_Y_new(:,1:number_of_eigenvectors);
[U_k,~,U_k_perp ] = gram_schemidt(U_k_gft);

Yuk = Y_bus * U_k;
  

SENSOR_PER_CLUSTER = ceil(nb*bmult);
[~,sensor_locations,al]= greedy_placement(U_k,SENSOR_PER_CLUSTER,[],1:nb,G,f,t,4,Yuk,[]);

% sensor_locations = randperm(85,ceil(nb*0.75));

number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end

%%
s_noise = ADD_NOISE*random_distribution(nb); % A customized function to generate random values
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
v_noise = ADD_NOISE*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;

c = Y_DD*vcomplex;
c_noise = ADD_NOISE*random_distribution(nb);
c = c + c_noise + 1i*c_noise;
cm = abs(c(sensor_locations));
sm= conj(s(sensor_locations));
vm = abs(v(sensor_locations));

id = eye(nb,nb);
HV = [];

for i = 1:nb
    HV=[HV;transpose(kron(id(:,i),id(:,i)))];
end
HV= HV(sensor_locations',:);

HS = [];
for i = 1:nb
    HS=[HS;kron(transpose(id(:,i)),Y_DD(i,:))];
end
HS= HS(sensor_locations',:);

HC = [];
for i = 1:nb
    HC=[HC;kron(conj(Y_DD(i,:)),Y_DD(i,:))];
end
HC= HC(sensor_locations',:);

xm=[vm.*vm;sm;cm.*cm];
H_M = [HV;HS;HC];

%%
weights_s = zeros(length(sensor_locations),1);
weights_c = zeros(length(sensor_locations),1);
for i = 1:length(sensor_locations)
    weights_s(i) = norm(Y_DD(i,i));
    weights_c(i) = (norm(Y_DD(i,i))*norm(Y_DD(i,i)));
end
weights = [ones(number_of_sensors,1);min(weights_s)*ones(number_of_sensors,1);min(weights_c)*ones(number_of_sensors,1)];
% weights = ones(3*number_of_sensors,1);

cvx_clear
cvx_solver Mosek
cvx_precision low
cvx_begin quiet
    variable W(nb,nb) complex semidefinite
    variable ep(size(xm))
    variable ep1(size(vm))
    variable ep2(size(sm))
    variable ep3(size(cm))
%     minimize (norm(ep1)+1*norm(ep2)+1*norm(ep3)) 
%     minimize norm(ep)
    minimize norm(weights.*(xm - H_M*vec(W))) 
    subject to
        H_M(1,:)*vec(W) == 1;
%         ep1 == vm.*vm - HV*vec(W);
%         ep2 == (sm) - HS*vec(W);
%         ep3 == cm.*cm - HC*vec(W);
%         ep == xm - H_M*vec(W);
     
cvx_end

if strcmpi('Solved',cvx_status) || strcmpi('Inaccurate/Solved',cvx_status)
    [voltage_mag2,voltage_angle2] = calculate_voltage_phasor(W,[],G);
    figure
    subplot(211)
    plot(voltage_mag2)
    hold on
    plot(abs(vcomplex))
    legend('Reconstructed','Original')
    subplot(212)
    plot(voltage_angle2)
    hold on
    plot(angle(vcomplex))
    legend('Reconstructed','Original')
end