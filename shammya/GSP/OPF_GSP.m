% clc;
clear all;
close all;
% rng(1)
define_constants;
%% OPF Case
opt = mpoption('verbose',0,'out.all',0);
opf_opt = mpoption();
opf_opt.opf.ac.solver='MIPS';
opf_mpc_case = loadcase('case85_opf.m');

opf_mpc_case.gencost = [
	2	0	0	3	0.11	0	0;
	2	0	0	3	0.085	0	0;
];

opf_case_results = runopf(opf_mpc_case,opf_opt);
%%



INCLUDE_LOAD = 0;
ADD_NOISE = 1;
casestr = 'case85';

mpc_case = loadcase(casestr); % load case to get back to the original network
base_case = mpc_case;
nb = length(mpc_case.bus);
nbr = length(mpc_case.branch);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
P_demand = mpc_case.bus(:,PD) / mpc_case.baseMVA;
Q_demand = mpc_case.bus(:,QD) / mpc_case.baseMVA;
mpc_case.bus(:,BS) = 0;
Q_shunt = mpc_case.bus(:,BS) / mpc_case.baseMVA;

G= graph(f,t);

[Y_bus,~,~] = makeYbus(mpc_case);
[mpc_case_results,success] = runpf(mpc_case,opt);
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s_complex = diag(vcomplex) * conj(Y_bus*vcomplex); % Calculating power S = V* conjugate(Y*V) & add noise
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

%%
number_of_eigenvectors = nb-75;
%U_V_new = load(strcat(casestr,'_U_V','.mat')).U_V_new;
U_k = U_Y_new (:,1:number_of_eigenvectors);
slack_bus = find(mpc_case.bus(:,BUS_TYPE)==3);
generation_bus = [1 5];
slack_bus_per_unit = mpc_case.bus(slack_bus,VM);
% slack_bus_per_unit = 1.093;
mult = zeros(nb,1);
mult(generation_bus) = 1;
C = [0.11;0.085];
% P_gen=zeros(nb,1);
% Q_gen=zeros(nb,1);

cvx_clear
cvx_solver SDPT3
cvx_precision low
cvx_begin quiet
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable pg(nb) 
    variable qg(nb) 
    variable e(nb)
%     minimize sum(C.*pg(generation_bus).*pg(generation_bus)+C.*qg(generation_bus).*qg(generation_bus))+...
%                     1*norm(pg.*mult-P_demand-real(conj(diag(Y_DD*U_k*tildeW*U_k'))))+ 1*norm(qg.*mult-Q_demand-imag(conj(diag(Y_DD*U_k*tildeW*U_k'))))
    minimize sum(C.*pg(generation_bus).*pg(generation_bus)+0*C.*qg(generation_bus).*qg(generation_bus))
    subject to
       0.9^2 <= diag(U_k * tildeW *(U_k)') <= 1.1^2;
       U_k(1,:)*(tildeW*U_k(1,:)') == slack_bus_per_unit^2;
       tildeW == tildeW';
%        pg(setdiff(1:nb,generation_bus))==0;
%        qg(setdiff(1:nb,generation_bus))==0;
%        -P_demand (2:end)  <= real(s(2:end)) <= -P_demand (2:end) ;
%        -Q_demand (2:end)  <= imag(s(2:end)) <= -Q_demand (2:end) ;
        sum(pg(generation_bus)) >= sum(P_demand);
%         sum(qg(generation_bus)) >= sum(Q_demand);
        norm(pg.*mult-P_demand-real(conj(diag(Y_DD*U_k*tildeW*U_k')))) <= 0.3
        norm(qg.*mult-Q_demand-imag(conj(diag(Y_DD*U_k*tildeW*U_k')))) <= 0.3
%        (pg.*mult-P_demand) + 1i*(qg.*mult + Q_shunt-Q_demand) == (conj(diag(Y_DD*U_k*tildeW*U_k')))
%        s_complex(2:end) == (conj(diag(Y_DD(2:end,2:end)*U_k(2:end,:)*tildeW*U_k(2:end,:)')))

         
cvx_end
[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
% 
% if strcmpi(cvx_status,'solved')
%     figure(1)
%     plot(abs(opf_case_results.bus(:,VM)),'linewidth',1.5)
%     hold on 
%     plot(abs(voltage_mag),'linewidth',1.5)
%     legend('original','reconstructed')
% else
%     fprintf('CVX Statues:%s\n',cvx_status)
% end

%% Without GSP
cvx_clear
cvx_solver SDPT3
cvx_precision low
cvx_begin 
    variable W(nb,nb) complex semidefinite
    variable pg_convex(nb) 
    variable qg_convex(nb) 
    minimize sum(C.*pg_convex(generation_bus).*pg_convex(generation_bus))
    subject to
        0.9^2 <= diag(W) <= 1.1^2;
        W(1,1) == slack_bus_per_unit^2;
        W == W';
        (pg_convex.*mult-P_demand) + 1i*(qg_convex.*mult -Q_demand) == (conj(diag(Y_DD*W)))
cvx_end
[voltage_mag_convex,voltage_angle_convex] = calculate_voltage_phasor(W,U_k,G);

if strcmpi(cvx_status,'solved')
    figure(4)
    plot(abs(voltage_mag_convex),'linewidth',1.5)
    hold on 
    plot(abs(voltage_mag),'linewidth',1.5)
    legend('no GSP','GSP')
else
    fprintf('CVX Statues:%s\n',cvx_status)
end


%% comparison of voltages
mpc_case_modified = opf_mpc_case;
mpc_case_modified.gen(:,PG) = pg(generation_bus);
mpc_case_modified_results = runpf(mpc_case_modified,opt);
vcomplex_modified = mpc_case_modified_results.bus(:,VM).*exp(1j*mpc_case_modified_results.bus(:,VA)*pi/180); % Preparing the complex voltages













































































%% Branch Flow Model

F=sparse(1:nbr,f,1,nbr,nb);
T=sparse(1:nbr,t,1,nbr,nb);
M=F-T;
r=mpc_case.branch(:,BR_R);
x=mpc_case.branch(:,BR_X);
Rline=r.*speye(nbr,nbr);
Xline=x.*speye(nbr,nbr);
TFT=T*F';
I=speye(size(TFT));
M=M(:,2:end);
R=2*(M\Rline)*((I-TFT)\T);
X=2*(M\Xline)*((I-TFT)\T);
V=[slack_bus_per_unit^2 ;slack_bus_per_unit^2+R*P_demand+X*Q_demand];
L = M\(2*diag(r)*inv(I-TFT)*diag(r) + 2*diag(x)*inv(I-TFT)*diag(x) - diag(r.*r+x.*r));
figure(2)
plot(abs(vcomplex),'linewidth',1.5)
hold on 
plot(abs(V),'linewidth',1.5)
legend('original','reconstructed')


cvx_clear
cvx_solver SDPT3
cvx_precision low
cvx_begin 
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_br(nbr)
    variable c_br(nbr)
    variable pg(nb) 
    variable qg(nb) 
    minimize 10*pg(generation_bus).*pg(generation_bus)+10*qg(generation_bus).*qg(generation_bus)
    subject to 
        0.8^2 <= diag(U_k * tildeW *(U_k)') <= 1.1^2;
        U_k(1,:)*(tildeW*U_k(1,:)') == slack_bus_per_unit^2;
        diag(U_k(2:end,:)*(tildeW*U_k(2:end,:)')) == slack_bus_per_unit^2 + R*(P_demand-pg.*mult)+X*(Q_demand-qg.*mult);
cvx_end

[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);

if strcmpi(cvx_status,'solved')
    figure(3)
    plot(abs(vcomplex),'linewidth',1.5)
    hold on 
    plot(abs(voltage_mag),'linewidth',1.5)
    legend('original','reconstructed')
else
    fprintf('CVX Statues:%s\n',cvx_status)
end