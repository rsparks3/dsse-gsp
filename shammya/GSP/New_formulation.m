clc;
clear all;
close all;
define_constants;
casestr = 'case85';
ADD_NOISE = 0;
opt = mpoption('verbose',0,'out.all',0);
mpc_case = loadcase(casestr); % load case to get back to the original network
base_case = mpc_case;

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

sensor_locations = 1:nb;
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
sm= s(sensor_locations);
vm = abs(v(sensor_locations));

%% Sanity Check of the equations 
number_of_eigenvectors = 15;
U_k_gft = U_Y_new (:,1:number_of_eigenvectors);

% if you want the SVD idea uncomment these 3 lines
% [U_k_p,S_u_k,~] = svd(U_k_gft);
% U_k = U_k_p(:,1:number_of_eigenvectors);
% U_k_perp = U_k_p(:,number_of_eigenvectors+1:nb);



[U_k,~,U_k_perp ] = gram_schemidt(U_k_gft);


v_tilde = (U_k)'*v;
var_epsilon = (U_k_perp)'*v;

epsilon_k = kron(conj(v),v)- kron(conj(U_k),U_k) * kron(conj(v_tilde),v_tilde);
fprintf('Validating Proposition 1 \n')
fprintf('Norm of Epsilon by original definition: %f \n', norm(epsilon_k)^2);
prop_1 = 2*norm(v_tilde)^2*norm(var_epsilon)^2 + norm(var_epsilon)^4;
fprintf('RHS from Proposition 1: %f \n', prop_1 );



%% New Formulation
id = eye(nb,nb);
H1 = [];
for i = 1:nb
    H1=[H1;kron(transpose(id(:,i)),Y_DD(i,:))];
end

H2 = [];
for i = 1:nb
    H2=[H2;transpose(kron(id(:,i),id(:,i)))];
end


H3 = [];
for i = 1:nb
    H3=[H3;kron(conj(Y_DD(i,:)),Y_DD(i,:))];
end
H_M = [H1;H2;H3];
fprintf('Condition Number for H_M : %f \n',cond(H_M));
U_k_kron = kron(conj(U_k),U_k);
H_tilde = H_M*U_k_kron;
% Considering All measurments
measurement = [conj(sm);vm.*vm;cm.*cm];
% Calculate the right hand side of equation 30
H_M_inv = pinv(H_M);
U_H = U_k_kron'*H_M_inv;

%% Solving Equation 30
cvx_clear
cvx_solver SeDuMi
cvx_precision low
cvx_begin
   variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite 
   variable x_m_p(number_of_eigenvectors*number_of_eigenvectors,1)
   minimize norm(x_m_p-vec(tildeW))
   subject to
    x_m_p == U_H * measurement;
cvx_end

if strcmpi('Solved',cvx_status)
    [voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
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
    rmse_angle  = rms(angle(vcomplex)-voltage_angle) ;
    rmse_voltage = rms(abs(vcomplex) - voltage_mag);
    fprintf('optval:%f RMS Angle:%f, RMS Voltage: %f\n',cvx_optval,rmse_angle,rmse_voltage);
    
end
%% Solving Equation 27
cvx_clear
cvx_solver Mosek_3
cvx_precision low
cvx_begin
    variable W1(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable ep1(number_of_sensors*3,1)

    minimize (norm(measurement- H_tilde*vec(W1))) 
    subject to
     
cvx_end
[voltage_mag2,voltage_angle2] = calculate_voltage_phasor(W1,U_k,G);
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
rmse_angle  = rms(angle(vcomplex)-voltage_angle2) ;
rmse_voltage = rms(abs(vcomplex) - voltage_mag2);
fprintf('optval:%f RMS Angle:%f, RMS Voltage: %f\n',cvx_optval,rmse_angle,rmse_voltage);


%%
U_k_perp = U_k_perp(:,1:number_of_eigenvectors);
A = H_M*[kron(conj(U_k),U_k_perp) kron(conj(U_k_perp),U_k) kron(U_k_perp,conj(U_k_perp))] ;  
[V,S_v,~]= svd(A);
V_k = V(:,end-number_of_eigenvectors*number_of_eigenvectors+1:end);

cvx_clear
cvx_solver SeDuMi
cvx_precision low
cvx_begin
   variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite 
   variable x_m_p(number_of_eigenvectors*number_of_eigenvectors,1)
   minimize norm(x_m_p-vec(tildeW))
   subject to
    x_m_p == V_k' * measurement;
cvx_end

if strcmpi('Solved',cvx_status)
    [voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
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
    rmse_angle  = rms(angle(vcomplex)-voltage_angle) ;
    rmse_voltage = rms(abs(vcomplex) - voltage_mag);
    fprintf('optval:%f RMS Angle:%f, RMS Voltage: %f\n',cvx_optval,rmse_angle,rmse_voltage);
    
end
