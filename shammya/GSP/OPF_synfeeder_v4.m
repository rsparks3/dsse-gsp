clc;
clear all;
close all;
rng(1)
define_constants;
bus_number = 400;
number_of_eigenvectors = 100;
gsp_approach = 1;

opt = mpoption('verbose',0,'out.all',0);
[n,e] = single_feeder_gen(bus_number,bus_number/20,0);

mpc_case = matpower_fmt(n,e);
mpc_case = parallel_branch_join(mpc_case);

generation_bus= randperm(bus_number,int16(bus_number/10));
% generation_bus (bus_number/10) = 100;

if sum(ismember(generation_bus,1))
else
   generation_bus(1) = 1 ;
end
generation_bus = sort(generation_bus,'ascend');

base_case = mpc_case;
base_case.bus(:,BS) = 0;


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
% mpc_case.branch(:,BR_B) = 0;
[Y_bus,Y_f,~] = makeYbus(mpc_case);
[mpc_case_results,success] = runpf(mpc_case);
pgmax = sum(mpc_case.bus(:,PD))/mpc_case.baseMVA/length(generation_bus);
power_flow_limit = 1.5*mpc_case_results.branch(:,PF)/mpc_case.baseMVA;
% pgmax(1) = 999;
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
gen_row (1,PMAX) = 99;

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
mpc_case_nocap = mpc_case;
mpc_case_nocap.branch(:,BR_B) = 0;
[Y_bus_nocap,~,~] = makeYbus(mpc_case_nocap);
Y_DD = full(Y_bus);
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages

idx = sub2ind(size(Y_DD),mpc_case.branch(:,1),mpc_case.branch(:,2));
line_y = -Y_DD(idx);
line_y = conj(line_y);
vv = vcomplex(mpc_case.branch(:,1))*(vcomplex(mpc_case.branch(:,1))-vcomplex(mpc_case.branch(:,2)))';
power_flow = mpc_case_results.baseMVA* diag(vv*diag((line_y)));

s_complex = diag(vcomplex) * conj(Y_bus*vcomplex); % Calculating power S = V* conjugate(Y*V) & add noise


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

%U_V_new = load(strcat(casestr,'_U_V','.mat')).U_V_new;
U_k = U_Y_new (:,1:number_of_eigenvectors);
slack_bus = find(mpc_case.bus(:,BUS_TYPE)==3);
% generation_bus = [1 5];
slack_bus_per_unit = mpc_case.bus(slack_bus,VM);
% slack_bus_per_unit = 1.093;
mult = zeros(nb,1);
mult(generation_bus) = 1;
Cs = cost;
[U_k_ortho,~,U_k_ortho_perp ] = gram_schemidt(U_k);
U_k = U_k_ortho;
% P_gen=zeros(nb,1);
% Q_gen=zeros(nb,1);

%% Covex Relaxation Code With GSP
A=U_k;
C=U_k';


H_f = [];
H_t = [];
for i = 1:length(f)
    H_f = [H_f;kron(transpose(C(:,f(i))),A(f(i),:))];
    H_t = [H_t;kron(transpose(C(:,t(i))),A(f(i),:))];
end
H_ft = H_f-H_t;


non_generation_bus = transpose(setdiff(1:nb,generation_bus));
eps_p_limit = real(sum(get_losses(mpc_case_results)))/mpc_case.baseMVA;
eps_q_limit = imag(sum(get_losses(mpc_case_results)))/mpc_case.baseMVA;
pgmax_vec = zeros(nb,1);
pgmax_vec(generation_bus) = 10*pgmax/100*randi([100 150],length(generation_bus),1);
id = speye(nb,nb);
[r_y,c_y] = size(Y_DD);
[r_i,c_i] = size (id);
H2 = zeros(nb,c_y*c_i);
for i = 1:nb
    H2(i,:)=kron(transpose(id(:,i)),Y_DD(i,:));
end
H2 = sparse(H2);

H1 = zeros(nb,c_i*c_i);

for i = 1:nb
    H1(i,:)=transpose(kron(id(:,i),id(:,i)));
end
H1 = sparse(H1);

H2_nocap = zeros(nb,c_y*c_i);
for i = 1:nb
    H2_nocap(i,:)=kron(transpose(id(:,i)),Y_bus_nocap(i,:));
end
H2_nocap = sparse(H2_nocap);

U_k_kron = kron(conj(U_k),U_k);
H1_tilde = H1*U_k_kron;
H2_tilde = H2*U_k_kron;
H2_nocap_tilde = H2_nocap*U_k_kron;

disp('Matrix Multiplication Completed.')
%%
gsp_approach = 1;
if gsp_approach == 1
    tstart_gsp = tic;
    cvx_clear
    cvx_solver Mosek_3
    cvx_precision low
    cvx_begin
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable pg(nb) nonnegative
    variable qg(nb) 
    variable delta_p(nb) 
    variable delta_q(nb)
    minimize (sum(Cs.*pg(generation_bus).*pg(generation_bus))+1*sum_square(delta_p)+1*sum_square(delta_q))
    subject to
        0.95^2 <= H1_tilde*vec(tildeW) <= 1.05^2;
        pg(2:end) <= pgmax_vec(2:end);
        -0.75*pgmax_vec(2:end) <= qg(2:end) <= 0.75*pgmax_vec(2:end);

%         tildeW == tildeW';
%         sum(pg(generation_bus)) >= sum(P_demand);
%         sum(qg(generation_bus)) >= sum(Q_demand);
        sum(pg(generation_bus)) - sum(P_demand) - 1i*(sum(qg(generation_bus)) -sum(Q_demand)) == sum(sum(H2_nocap_tilde*vec(tildeW)));
        (pg - P_demand) - 1i*(qg - Q_demand) - H2_tilde*vec(tildeW) == delta_p + 1i*delta_q;
%         pg - P_demand - real(H2_tilde*vec(tildeW)) == 0;
%         qg - Q_demand + imag(H2_tilde*vec(tildeW)) == 0;
%         for br = 1:length(mpc_case.branch)
%             -power_flow_limit(br) <= real((U_k(f(br),:)*(tildeW*U_k(f(br),:)')-U_k(f(br),:)*(tildeW*U_k(t(br),:)')  )*conj(line_y(br))) <= power_flow_limit(br);
%         end
        -power_flow_limit <= real( (H_ft*vec(tildeW)) .*line_y) <= power_flow_limit;

    cvx_end
    if strcmpi(cvx_status,'solved')
        telapsed_gsp = toc(tstart_gsp)
        [voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
        save(strcat('OPF_GSP_',string(bus_number),'.mat'),'bus_number','pg','qg','generation_bus','Cs','number_of_eigenvectors','mpc_case','voltage_mag','voltage_angle','opf_case_results','telapsed_gsp','delta_p','delta_q')
    end
end

%% Without GSP

tstart_no_gsp = tic;
gsp_approach = 0;
if gsp_approach == 0
    cvx_clear
    cvx_solver Mosek_3
    cvx_precision low
    cvx_begin quiet
       variable W(nb,nb) complex semidefinite
       variable pg_convex(nb) nonnegative
       variable qg_convex(nb) 
       minimize sum(Cs.*pg_convex(generation_bus).*pg_convex(generation_bus))
       subject to
           0.95^2 <= diag(W) <= 1.05^2;
           (pg_convex-P_demand)  == real(conj(diag(Y_DD*W)));
           (qg_convex-Q_demand) == imag(conj(diag(Y_DD*W)));
           pg_convex(2:end) <= pgmax_vec(2:end);
           -0.75*pgmax_vec(2:end) <= qg_convex(2:end) <= 0.75*pgmax_vec(2:end);
           for br = 1:length(mpc_case.branch)
               ff = mpc_case.branch(br,1);
               tt = mpc_case.branch(br,2);
               -power_flow_limit(br) <= real((W(ff,ff)-W(ff,tt))*line_y(br)) <= power_flow_limit(br);
           end
    cvx_end
    if strcmpi(cvx_status,'solved')
        [voltage_mag_convex,voltage_angle_convex] = calculate_voltage_phasor(W,U_k,G);
        telapsed = toc(tstart_no_gsp)
%         pg_error = rms(pg_convex-pg);
%         fprintf('Pg Error:%f\n',rms(pg_convex(generation_bus)-pg(generation_bus)))
%         save(strcat('OPF_NOGSP_',string(bus_number),'.mat'),'bus_number','telapsed','pg_convex','qg_convex','generation_bus','Cs','number_of_eigenvectors','mpc_case','voltage_mag_convex','voltage_angle_convex','opf_case_results');
    end
end



%% rersolve the case
mpc_case_new=opf_mpc_case;
% mpc_case_new.gen(2:end,PG)= mpc_case_new.baseMVA*pg(generation_bus(2:end));
mpc_case_new.gen = mpc_case_new.gen(1,:);
mpc_case_new.bus(generation_bus(2:end),PD) = mpc_case_new.bus(generation_bus(2:end),PD) - mpc_case_new.baseMVA*pg(generation_bus(2:end));
mpc_case_new.bus(generation_bus(2:end),QD) = mpc_case_new.bus(generation_bus(2:end),QD) - mpc_case_new.baseMVA*qg(generation_bus(2:end));
mpc_case_new_results = runpf(mpc_case_new);
new_voltage = mpc_case_new_results.bus(:,VM);
[loss,fchg] = get_losses(mpc_case_new_results);
vcomplex = mpc_case_new_results.bus(:,VM).*exp(1j*mpc_case_new_results.bus(:,VA)*pi/180); % Preparing the complex voltages

s_complex = diag(vcomplex) * conj(Y_bus*vcomplex); % Calculating power S = V* conjugate(Y*V) & add noise
sum(loss)-1i*sum(fchg)
% mean(abs(new_voltage-voltage_mag))












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