clc;
clear all;
close all;
rng(1)
define_constants;
bus_number = 400;
number_of_eigenvectors = 100;
gsp_approach = 1;
slack_multi = 0;

opt = mpoption('verbose',0,'out.all',0);
[n,e] = single_feeder_gen(bus_number,bus_number/20,0);

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
power_flow_limit = 1.25*mpc_case_results.branch(:,PF)/mpc_case.baseMVA;
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

vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s_complex = diag(vcomplex) * conj(Y_bus*vcomplex); % Calculating power S = V* conjugate(Y*V) & add noise
Y_DD = full(Y_bus);
idx = sub2ind(size(Y_DD),mpc_case.branch(:,1),mpc_case.branch(:,2));
line_y = -Y_DD(idx);
line_y = conj(line_y);
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

H_f = [];
H_t = [];
H_tt = [];
for i = 1:length(f)
    H_f = [H_f;kron(transpose(C(:,f(i))),A(f(i),:))];
    H_t = [H_t;kron(transpose(C(:,t(i))),A(f(i),:))];
    H_tt = [H_tt;kron(transpose(C(:,t(i))),A(t(i),:))];
end
H_ft = H_f-H_t;


disp('Matrix Multiplication Completed.')
non_generation_bus = transpose(setdiff(1:nb,generation_bus));
eps_p_limit = real(sum(get_losses(mpc_case_results)))/mpc_case.baseMVA;
eps_q_limit = imag(sum(get_losses(mpc_case_results)))/mpc_case.baseMVA;
pgmax_vec = zeros(nb,1);
pgmax_vec(generation_bus) = 10*pgmax/100*randi([100 150],length(generation_bus),1);
Y_DD = sparse(Y_DD);
%%
br_b = 0.5*mpc_case.branch(:,BR_B);
if gsp_approach == 1
    tstart_gsp = tic;
    cvx_clear
    cvx_solver Mosek_3
    cvx_precision low
    cvx_begin quiet
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable pg(nb) nonnegative
    variable qg(nb) 
    variable delta_p(nb) 
    variable delta_q(nb)
    minimize (sum(Cs.*pg(generation_bus).*pg(generation_bus))+0*sum_square(delta_p)+0*sum_square(delta_q))

    subject to
        0.95^2 <= real(Uk_UkH*vec(transpose(tildeW))) <= 1.05^2;
        pg(2:end) <= slack_multi*pgmax_vec(2:end);
        -0.75*pgmax_vec(2:end)*slack_multi <= qg(2:end) <= 0.75*pgmax_vec(2:end)*slack_multi;
        sum(pg(generation_bus)) - sum(P_demand) - 1i*(sum(qg(generation_bus)) -sum(Q_demand)) == sum(sum(YUk_UkH*vec(transpose(tildeW))));
        pg - P_demand - real(YUk_UkH*vec(transpose(tildeW))) == delta_p;
        qg - Q_demand + imag(YUk_UkH*vec(transpose(tildeW))) == delta_q;
        -power_flow_limit <= real( (H_ft*vec(tildeW)) .*line_y) <= power_flow_limit;
    cvx_end
    if strcmpi(cvx_status,'solved')
        telapsed_gsp = toc(tstart_gsp)
        [voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
        if slack_multi == 1
%         save(strcat('OPF_GSP_',string(bus_number),'.mat'),'bus_number','pg','qg','generation_bus','Cs','number_of_eigenvectors','mpc_case','voltage_mag','voltage_angle','opf_case_results','telapsed_gsp')
        end
    end
end

%% re-solve the case
mpc_case_new_gsp=opf_mpc_case;
% mpc_case_new.gen(2:end,PG)= mpc_case_new.baseMVA*pg(generation_bus(2:end));
mpc_case_new_gsp.gen = mpc_case_new_gsp.gen(1,:);
mpc_case_new_gsp.gen(slack_bus,VG) = voltage_mag(1);
mpc_case_new_gsp.bus(generation_bus(2:end),PD) = mpc_case_new_gsp.bus(generation_bus(2:end),PD) - mpc_case_new_gsp.baseMVA*pg(generation_bus(2:end));
mpc_case_new_gsp.bus(generation_bus(2:end),QD) = mpc_case_new_gsp.bus(generation_bus(2:end),QD) - mpc_case_new_gsp.baseMVA*qg(generation_bus(2:end));
mpc_case_new_results_gsp = runpf(mpc_case_new_gsp,opt);
new_voltage = mpc_case_new_results_gsp.bus(:,VM);
[loss,fchg,tchg] = get_losses(mpc_case_new_results_gsp);
vcomplex = mpc_case_new_results_gsp.bus(:,VM).*exp(1j*mpc_case_new_results_gsp.bus(:,VA)*pi/180); % Preparing the complex voltages

s_complex = diag(vcomplex) * conj(Y_bus*vcomplex); % Calculating power S = V* conjugate(Y*V) & add noise

flow_gsp = real((H_ft*vec(tildeW)) .*line_y)*mpc_case_new_gsp.baseMVA;
flow_gsp_matpower = mpc_case_new_results_gsp.branch(:,PF);
fprintf('Voltage Magnitude Diff:%f\n',max(abs(voltage_mag-new_voltage)));
fprintf('Voltage Angle Diff:%f\n',max(abs(voltage_angle-angle(vcomplex))));
fprintf('Difference in Real Flow in MW:%f\n',max(abs(flow_gsp-flow_gsp_matpower)))

figure(1)
subplot(311)
plot(new_voltage)
hold on 
plot(voltage_mag)
legend('matpower','voltage gsp')

subplot(312)
plot(angle(vcomplex))
hold on 
plot(voltage_angle)
legend('matpower','voltage gsp')

subplot(313)
plot(flow_gsp_matpower)
hold on 
plot(flow_gsp)
legend('matpower','voltage gsp')

