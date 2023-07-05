% clc;
clear all;
close all;
rng(1)
define_constants;
bus_number = 600;
opt = mpoption('verbose',0,'out.all',0);
[n,e] = single_feeder_gen(bus_number);

mpc_case = matpower_fmt(n,e);
mpc_case = parallel_branch_join(mpc_case);

generation_bus= randperm(bus_number,int16(bus_number/10));

if sum(ismember(generation_bus,1))
else
   generation_bus(1) = 1 ;
end
generation_bus = sort(generation_bus,'ascend');

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
mpc_case.branch(1,TAP) = 1;

[Y_bus,~,~] = makeYbus(mpc_case);
[mpc_case_results,success] = runpf(mpc_case);
pgmax = sum(mpc_case.bus(:,PD))/mpc_case.baseMVA/length(generation_bus);
%% OPF Case

opf_opt = mpoption();
opf_opt.opf.ac.solver='MIPS';
opf_mpc_case = mpc_case;
gen_row_base = opf_mpc_case.gen;
gencost_base = [2	0	0	3	0.11	0	0];
gencost = gencost_base;
gen_row = gen_row_base;
for r = 2:length(generation_bus)
    gen_row = [gen_row;gen_row_base];
    gencost = [gencost;gencost_base];
end


gen_row(:,1) = generation_bus';
gen_row (:,PMAX) = 10*pgmax*mpc_case.baseMVA;
gen_row (:,PMIN) = 0;

opf_mpc_case.gen = gen_row;
b = 0.011; a = 0.85;
cost = (b-a).*rand(length(generation_bus),1) + a;
gencost(:,5) = cost;


opf_mpc_case.gencost = gencost;

[opf_case_results,success_opf] = runopf(opf_mpc_case,opf_opt);
if success_opf==0
    error('MATPOWER OPF FAILED.')
end
%%

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
number_of_eigenvectors = 250;
%U_V_new = load(strcat(casestr,'_U_V','.mat')).U_V_new;
U_k = U_Y_new (:,1:number_of_eigenvectors);
slack_bus = find(mpc_case.bus(:,BUS_TYPE)==3);
% generation_bus = [1 5];
slack_bus_per_unit = mpc_case.bus(slack_bus,VM);
% slack_bus_per_unit = 1.093;
mult = zeros(nb,1);
mult(generation_bus) = 1;
Cs = cost;
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
disp('Matrix Multiplication Completed.')
non_generation_bus = transpose(setdiff(1:nb,generation_bus));
eps_p_limit = real(sum(get_losses(mpc_case_results)))/mpc_case.baseMVA;
eps_q_limit = imag(sum(get_losses(mpc_case_results)))/mpc_case.baseMVA;
%%
cvx_clear
cvx_solver SDPT3
cvx_precision low
cvx_begin
variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
variable pg(length(generation_bus)) nonnegative
variable qg(length(generation_bus)) 
% expression delta_p(nb) 
% expression delta_q(nb)
variable delta_p(nb) 
variable delta_q(nb)
% minimize (sum(Cs.*pg.*pg)+1000*delta_p+1000*delta_q)
minimize (sum(Cs.*pg.*pg)+100000*sum_square(delta_p)+100000*sum_square(delta_q))
% minimize (sum(Cs.*pg.*pg))

subject to
    for k = 1:nb
        0.95^2 <= real(Uk_UkH(k,:)*vec(transpose(tildeW))) <= 1.05^2;
    end
	pg(2:end) <= pgmax;
	-0.75*pgmax <= qg(2:end) <= 0.75*pgmax;
%     U_k(1,:)*(tildeW*U_k(1,:)') == slack_bus_per_unit^2;
    tildeW == tildeW';
    sum(pg) >= sum(P_demand);
    sum(qg) >= sum(Q_demand);
    for k = 1:nb
         if ismember(k,generation_bus)
            delta_p(k) ==  pg(generation_bus==k)-P_demand(k) - real(YUk_UkH(k,:)*vec(transpose(tildeW)));
            delta_q(k) ==  qg(generation_bus==k)-Q_demand(k) + imag((YUk_UkH(k,:)*vec(transpose(tildeW))));
        else
            delta_p (k) == -P_demand(k) - real((YUk_UkH(k,:)*vec(transpose(tildeW))));
            delta_q (k) ==  -Q_demand(k) + imag((YUk_UkH(k,:)*vec(transpose(tildeW))));
         end
    end
%      for k = 1:nb
%          if ismember(k,generation_bus)
%             0 ==  pg(generation_bus==k)-P_demand(k) - real(YUk_UkH(k,:)*vec(transpose(tildeW)));
%             0 ==  qg(generation_bus==k)-Q_demand(k) + imag((YUk_UkH(k,:)*vec(transpose(tildeW))));
%         else
%             0 == -P_demand(k) - real((YUk_UkH(k,:)*vec(transpose(tildeW))));
%             0 ==  -Q_demand(k) + imag((YUk_UkH(k,:)*vec(transpose(tildeW))));
%          end
%     end
%     sum(delta_p) <= eps_p_limit;
%     sum(delta_q) <= eps_q_limit;
    

cvx_end
if strcmpi(cvx_status,'solved')
    [voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
    save(strcat('OPF_',string(bus_number),'.mat'),'bus_number','pg','qg','generation_bus','Cs','number_of_eigenvectors','mpc_case','voltage_mag','voltage_angle')
end





%% Without GSP
% cvx_clear
% cvx_solver SDPT3
% cvx_precision low
% cvx_begin 
%    variable W(nb,nb) complex semidefinite
%    variable pg_convex(nb) nonnegative
%    variable qg_convex(nb) 
%    minimize sum(Cs.*pg_convex(generation_bus).*pg_convex(generation_bus))
%    subject to
%        0.95^2 <= diag(W) <= 1.05^2;
% %        W(1,1) == slack_bus_per_unit^2;
%        W == W';
%        sum(pg_convex(generation_bus)) >= sum(P_demand);
%        sum(qg_convex(generation_bus)) >= sum(Q_demand);
%        (pg_convex.*mult-P_demand) + 1i*(qg_convex.*mult -Q_demand) == (conj(diag(Y_DD*W)));
%        for k = 2:length(generation_bus)
%             pg_convex(generation_bus(k)) <= pgmax;
%            -0.75*pgmax <= qg_convex(generation_bus(k)) <= 0.75*pgmax;
% 
%        end
% cvx_end
% [voltage_mag_convex,voltage_angle_convex] = calculate_voltage_phasor(W,U_k,G);
% pg_error = rms(pg_convex(generation_bus)-pg);
% fprintf('Pg Error:%f\n',rms(pg_convex(generation_bus)-pg))
% save(strcat('OPF_',string(bus_number),'.mat'),'bus_number','pg','qg','generation_bus','Cs','number_of_eigenvectors','mpc_case','voltage_mag','voltage_angle','pg_error');
% 
% 
% 











%%

%     for k = 1:nb
%          if ismember(k,generation_bus)
%             delta_p = delta_p + (pg(generation_bus==k)-P_demand(k) - real(conj(YUk_UkH(k,:)*vec(transpose(tildeW)))))^2;
%             delta_q = delta_q + (qg(generation_bus==k)-Q_demand(k) - imag(conj(YUk_UkH(k,:)*vec(transpose(tildeW)))))^2;
%         else
%             delta_p = delta_p + (-P_demand(k) - real(conj(YUk_UkH(k,:)*vec(transpose(tildeW)))))^2;
%             delta_q = delta_q + (-Q_demand(k) - imag(conj(YUk_UkH(k,:)*vec(transpose(tildeW)))))^2;
%          end
%     end