% clc;
clear all;
close all;
% rng(1)
define_constants;
bus_number = 300;
[n,e] = single_feeder_gen(bus_number);

mpc_case = matpower_fmt(n,e);
mpc_case = parallel_branch_join(mpc_case);

%% OPF Case
opt = mpoption('verbose',0,'out.all',0);
opf_opt = mpoption();
opf_opt.opf.ac.solver='MIPS';
opf_mpc_case = mpc_case;
gen_row = opf_mpc_case.gen;
gen_row = [gen_row;gen_row];
gen_row(2,1) = 15;
opf_mpc_case.gen = gen_row;

opf_mpc_case.gencost = [
    2	0	0	3	0.11	0	0;
    2	0	0	3	0.085	0	0;
    ];

[opf_case_results,success_opf] = runopf(opf_mpc_case,opf_opt);
if success_opf==0
    error('MATPOWER OPF FAILED.')
end
%%
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
number_of_eigenvectors = 35;
%U_V_new = load(strcat(casestr,'_U_V','.mat')).U_V_new;
U_k = U_Y_new (:,1:number_of_eigenvectors);
slack_bus = find(mpc_case.bus(:,BUS_TYPE)==3);
generation_bus = [1 5];
slack_bus_per_unit = mpc_case.bus(slack_bus,VM);
% slack_bus_per_unit = 1.093;
mult = zeros(nb,1);
mult(generation_bus) = 1;
Cs = [0.11;0.085];
% P_gen=zeros(nb,1);
% Q_gen=zeros(nb,1);

%% Covex Relaxation Code With GSP
Uk_UkH =  zeros(number_of_eigenvectors,number_of_eigenvectors^2);
YUk_UkH =  zeros(number_of_eigenvectors,number_of_eigenvectors^2);
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
%%
cvx_clear
cvx_solver SDPT3
cvx_precision low
cvx_begin
variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
variable pg(nb)
variable qg(nb)
expression delta_p 
expression delta_q
%     variable e(nb)
minimize sum(Cs.*pg(generation_bus).*pg(generation_bus))
%     minimize sum(C.*pg(generation_bus).*pg(generation_bus)+0*C.*qg(generation_bus).*qg(generation_bus))
subject to
    for k = 1:nb
        0.9^2 <= real(Uk_UkH(k,:)*vec(transpose(tildeW))) <= 1.1^2;
    end
    U_k(1,:)*(tildeW*U_k(1,:)') == slack_bus_per_unit^2;
    tildeW == tildeW';
    sum(pg(generation_bus)) >= sum(P_demand);
    sum(qg(generation_bus)) >= sum(Q_demand);
    delta_p =0;
    delta_q =0;
    for k = 1:nb
        delta_p = delta_p + (pg(k)*mult(k)-P_demand(k) - real(conj(YUk_UkH(k,:)*vec(transpose(tildeW)))))^2;
        delta_q = delta_q + (qg(k)*mult(k)-Q_demand(k) - imag(conj(YUk_UkH(k,:)*vec(transpose(tildeW)))))^2;
    end
    delta_p <= 0.5;
    delta_q <= 0.5;
cvx_end
[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);



%% Without GSP
% cvx_clear
% cvx_solver SDPT3
% cvx_precision low
% cvx_begin 
%    variable W(nb,nb) complex semidefinite
%    variable pg_convex(nb) 
%    variable qg_convex(nb) 
%    minimize sum(Cs.*pg_convex(generation_bus).*pg_convex(generation_bus))
%    subject to
%        0.9^2 <= diag(W) <= 1.1^2;
%        W(1,1) == slack_bus_per_unit^2;
%        W == W';
%        (pg_convex.*mult-P_demand) + 1i*(qg_convex.*mult -Q_demand) == (conj(diag(Y_DD*W)));
% cvx_end
% [voltage_mag_convex,voltage_angle_convex] = calculate_voltage_phasor(W,U_k,G);
