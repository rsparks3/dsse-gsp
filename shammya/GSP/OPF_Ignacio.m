clc;
clear all;
close all;
% rng(1)
%% Define the global constants
SYN_FEEDER = true; % Whether the code is expecting a synthetic feeder or a general Matpower case, for general matpower  case, set to false, otherwise true 
opt = mpoption('verbose',0,'out.all',0);
define_constants; % Define the matpower constants
MATPOWER_CASE_FILE = 'case85'; % MATPOWER case file name, ignored if SYN_FEEDER = 1
BUS_NUMBER = 100; % Size of the synthetic feeder, ignored if SYN_FEEDER= 0
number_of_eigenvectors = 70; 
ORTHORNOMAL_DECOMPOSITION = 1; % 1 = Gram Schemidt, 2 = SVD
GSP = true; % if false, uses the basic convex relaxation equation and avoids GSP integration
cvx_solver SDPT3; % Choose appropriate solver for CVX
%% Load the appropriate case file
if SYN_FEEDER == true
    [n,e] = single_feeder_gen(BUS_NUMBER);
    mpc_case = matpower_fmt(n,e);
    mpc_case = parallel_branch_join(mpc_case);
    mpc_case.branch(:,BR_B) = 0; % Forcing the case BRANCH SUSCEPTENCE to 0
    mpc_case.bus(:,BS) = 0; % Shunt Reactance 0
    mpc_case.branch(1,TAP) = 1; % Force the TAP to be 1
else
    mpc_case = loadcase(MATPOWER_CASE_FILE);
    mpc_case.branch(:,BR_B) = 0; % Forcing the case BRANCH SUSCEPTENCE to 0
    mpc_case.bus(:,BS) = 0; % Shunt Reactance 0
%     mpc_case.branch(1,TAP) = 1; % Force the TAP to be 1
end

%% Only considering slack bus as the generation bus
generation_bus= [1];
% generation_bus (bus_number/10) = 100;

if sum(ismember(generation_bus,1))
else
   generation_bus(1) = 1 ;
end
generation_bus = sort(generation_bus,'ascend');


nb = length(mpc_case.bus);
nbr = length(mpc_case.branch);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
P_demand = mpc_case.bus(:,PD) / mpc_case.baseMVA;
Q_demand = mpc_case.bus(:,QD) / mpc_case.baseMVA;
slack_bus = mpc_case.bus(:,BUS_TYPE)==3;

G= graph(f,t);
[Y_bus,~,~] = makeYbus(mpc_case);
[mpc_case_results,success] = runpf(mpc_case,opt);
pgmax = sum(mpc_case.bus(:,PD))/mpc_case.baseMVA/length(generation_bus);
power_flow_limit = 1.25*mpc_case_results.branch(:,PF)/mpc_case.baseMVA;

%% OPF By Matpower 
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

%% Check All the Values
Y_DD = full(Y_bus);
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
idx = sub2ind(size(Y_DD),mpc_case.branch(:,1),mpc_case.branch(:,2));
line_y = conj(-Y_DD(idx));
vv = vcomplex(mpc_case.branch(:,1))*(vcomplex(mpc_case.branch(:,1))-vcomplex(mpc_case.branch(:,2)))';
power_flow = mpc_case_results.baseMVA* diag(vv*diag((line_y)));

s_complex = diag(vcomplex) * conj(Y_bus*vcomplex); % Calculating power S = V* conjugate(Y*V) & add noise



%%
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
U_k_gft = U_Y_new (:,1:number_of_eigenvectors);
if ORTHORNOMAL_DECOMPOSITION  == 1
    [U_k,~,~] = gram_schemidt(U_k_gft);
elseif ORTHORNOMAL_DECOMPOSITION  == 2
    [U_k,~,~] = svd(U_k_gft);
    U_k = U_k (:,1:number_of_eigenvectors);
else
    error('Orthonormal Decomposition Mode not Understood, Use 1 for Gram Schemidt, 2 for SVD')
end

mult = zeros(nb,1);
mult(generation_bus) = 1;

% Swapping so that utility has the highest cost
Cs = cost;
[temp,loc] = max(Cs);
Cs(loc)= Cs(1);
Cs(1)= temp;



%% Covex Relaxation Code With GSP
A=U_k;
C=U_k';

H_f = [];
H_t = [];
H_tt = [];
for i = 1:length(f)
    H_f = [H_f;kron(transpose(C(:,f(i))),A(f(i),:))];
    H_t = [H_t;kron(transpose(C(:,t(i))),A(f(i),:))];
    H_tt = [H_tt;kron(transpose(C(:,t(i))),A(t(i),:))];
end
H_ft = H_f-H_t;

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

U_k_kron = kron(conj(U_k),U_k);
H1_tilde = H1*U_k_kron;
H2_tilde = H2*U_k_kron;

disp('Matrix Multiplication Completed.')
%%
close all 
if GSP == true
    tstart_gsp = tic;
    cvx_clear
    cvx_solver SDPT3
    cvx_precision low
    cvx_begin 
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable pg(nb) nonnegative
    variable qg(nb) 
    variable delta_p(nb)  
    variable delta_q(nb)
%     minimize (0.01*norm(delta_p)+0.01*norm(delta_q) )
%     minimize norm(H2_tilde*vec(tildeW))
    minimize 1
%     minimize (sum(Cs.*pg(generation_bus).*pg(generation_bus)))
    subject to
        0.75^2 <= H1_tilde*vec(tildeW) <= 1.05^2;
        pg(2:end) == 0;
        qg(2:end) == 0;
        pg(1) >= sum(P_demand);
        qg(1) >= sum(Q_demand);
%         (pg - P_demand) - 1i*(qg - Q_demand) == H2_tilde*vec(tildeW);
        (pg - P_demand) - 1i*(qg - Q_demand) - H2_tilde*vec(tildeW) == delta_p + 1i*delta_q;
%         (pg - P_demand) == real( H2_tilde*vec(tildeW));
%         (qg - Q_demand) == imag(conj(H2_tilde*vec(tildeW)));
        U_k(1,:)*(tildeW*U_k(1,:)') == 1;


    cvx_end
    if strcmpi(cvx_status,'solved')
        telapsed_gsp = toc(tstart_gsp)
        [voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
        mpc_case_new_gsp=opf_mpc_case;
        % mpc_case_new.gen(2:end,PG)= mpc_case_new.baseMVA*pg(generation_bus(2:end));
        mpc_case_new_gsp.gen(slack_bus,VG) = voltage_mag(1);
        mpc_case_new_results_gsp = runpf(mpc_case_new_gsp,opt);
        new_voltage = mpc_case_new_results_gsp.bus(:,VM);
        [loss,fchg,tchg] = get_losses(mpc_case_new_results_gsp);
        vcomplex = mpc_case_new_results_gsp.bus(:,VM).*exp(1j*mpc_case_new_results_gsp.bus(:,VA)*pi/180); % Preparing the complex voltages


        flow_gsp = real((H_ft*vec(tildeW)) .*line_y)*mpc_case_new_gsp.baseMVA;
        flow_gsp_matpower = mpc_case_new_results_gsp.branch(:,PF);
        fprintf('Voltage Magnitude Diff:%f\n',max(abs(voltage_mag-new_voltage)));
        fprintf('Voltage Angle Diff:%f\n',max(abs(voltage_angle-angle(vcomplex))));
        fprintf('Difference in Real Flow in MW:%f\n',max(abs( abs(flow_gsp)-abs(flow_gsp_matpower))))

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
        plot(abs(flow_gsp_matpower))
        hold on 
        plot(abs(flow_gsp))
        legend('matpower','voltage gsp')
     end
end
