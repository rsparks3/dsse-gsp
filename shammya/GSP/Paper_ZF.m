% This code generates the plots for the matpower 85 bus test case
% considering sensors at all the buses, this code tries to justify that low
% rank approximation can help use to enable a state estimation formulation
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

Yddh= Y_DD';

%% Sanity Check of the equations 
id = eye(nb,nb);
H1 = [];

for i = 1:nb
    H1=[H1;transpose(kron(id(:,i),id(:,i)))];
end

w = kron(conj(v),v);
vm_kron = H1*w;

H2 = [];
for i = 1:nb
    H2=[H2;kron(transpose(id(:,i)),Y_DD(i,:))];
end

sm_kron = conj(H2*w);

H3 = [];
for i = 1:nb
    H3=[H3;kron(conj(Y_DD(i,:)),Y_DD(i,:))];
end

cm_kron = H3*w;

W = v*v';
diag_v = diag(v*v');



cvx_clear
cvx_solver SeDuMi
cvx_precision low
cvx_begin
    variable w_cvx(nb*nb,1)
    variable ep1(nb,1)
    variable ep2(nb,1)
    variable ep3(nb,1)
    minimize (norm(ep1)+norm(ep2)+norm(ep3)) 
    subject to
        ep1 == vm.*vm - H1*w_cvx;
        ep2 == conj(sm) - H2*w_cvx;
        ep3 == cm.*cm - H3*w_cvx;
        
cvx_end
tildeW = reshape(w_cvx,[nb,nb]);
[voltage_mag,voltage_angle] = calculate_voltage_phasor(transpose(tildeW),U_k,G);
%%
id = eye(nb,nb);
term1 = transpose(kron(id(:,1),id(:,1)));
H = [];

for i = 1:nb
    H=[H;kron(transpose(id(:,i)),Y_DD(i,:))];
end

for i = 1:nb
    H=[H;transpose(kron(id(:,i),id(:,i)))];
end

for i = 1:nb
%     H=[H;kron(conj(Y_DD(i,:)),Y_DD(i,:))];
    H=[H;kron(transpose(Y_DD(:,i)),Y_DD(i,:))];
end

E = [];Ebar = [];
for i = 1:nb
    for j = i:nb
%         if j ~=i
            E = [E;transpose(kron(id(:,j),id(:,i)))];
            Ebar = [Ebar;-1*transpose(kron(id(:,i),id(:,j)))];
%         end
    end
end

% H_sdp = [term1 zeros(size(term1)); H zeros(size(H)); E Ebar];
H_sdp = [H zeros(size(H)); E Ebar];
pinv_hsdp = pinv(H_sdp);
all_meas = [sm;vm;cm;zeros(nb/2*(nb-1),1)];
ww = pinv_hsdp*all_meas;

w = ww(1:nb*nb);
wconj = ww(nb*nb+1:end);
tildeW = reshape(w,[nb,nb]);
[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
%%
number_of_eigenvectors = nb;
U_k = U_Y_new (:,1:number_of_eigenvectors);
Uk_UkH =  zeros(number_of_eigenvectors,number_of_eigenvectors*number_of_eigenvectors);
YUk_UkH =  zeros(number_of_eigenvectors,number_of_eigenvectors*number_of_eigenvectors);
UkH_YkH =  zeros(number_of_eigenvectors,number_of_eigenvectors*number_of_eigenvectors);
A=U_k;
C=U_k';
for k = 1:nb
   Uk_UkH(k,:) =  transpose(vec(C(:,k)*A(k,:)));
end
A=Y_DD*U_k;
for k = 1:nb
   YUk_UkH(k,:) =  transpose(vec(C(:,k)*A(k,:)));
end

C = U_k'*Y_DD';
for k = 1:nb
   UkH_YkH(k,:) =  transpose(vec(C(:,k)*A(k,:)));
end

tic
cvx_clear
cvx_solver Mosek_2
cvx_precision low
cvx_begin 
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb)
    variable a_est(nb)
    minimize 1

    subject to
        tildeW == tildeW';
%         U_k(1,:)*(tildeW*U_k(1,:)') == 1;
% The main three equations
%         vm.*vm == diag(tildeW);
%         sm == conj(diag(Y_DD*tildeW));
%         cm.*cm == diag(Y_DD*tildeW*Y_DD');
% THe equations with Tilde W
%         sm == conj(diag(Y_DD*U_k*tildeW*U_k'));
%         vm.*vm == diag(U_k*tildeW*(U_k)');
%         cm.*cm  == diag(Y_DD*U_k*tildeW*U_k'*Y_DD');  
        
        vm.*vm == Uk_UkH*vec(transpose(tildeW)); 
        sm == conj( YUk_UkH*vec(transpose(tildeW)) ); 
        cm.*cm == UkH_YkH*vec(transpose(tildeW)); 
        

cvx_end
toc
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

% csvwrite('Paper_Figure1.csv',[transpose(1:nb) voltage_mag voltage_angle abs(vcomplex) angle(vcomplex)])

% v_new = (U_k)*pinv(M*U_k)*vcomplex(sensor_locations);
% figure(2)
% plot(angle(v_new))
% hold on 
% plot(angle(vcomplex))






