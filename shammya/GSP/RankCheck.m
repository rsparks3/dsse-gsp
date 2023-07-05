clc;
clear all;
close all;

define_constants;

mpc_case = loadcase('case85.m');
base_case = mpc_case;
mpc_case_results = runpf(mpc_case);
number_of_scenarios = 500;
V = zeros(length(mpc_case.bus),number_of_scenarios);
Vangle = zeros(length(mpc_case.bus),number_of_scenarios);
PD_case = zeros(length(mpc_case.bus),number_of_scenarios);
QD_case = zeros(length(mpc_case.bus),number_of_scenarios);

V(:,1) = mpc_case_results.bus(:,VM);
Vangle(:,1) =  mpc_case_results.bus(:,VA);
PD_case(:,1) = base_case.bus(:,PD);
QD_case(:,1) = base_case.bus(:,QD);
avg_pd = mean(base_case.bus(:,PD));
avg_qd = mean(base_case.bus(:,QD));
i = 1;
while i<number_of_scenarios
    mult = ones(1,length(mpc_case.bus));
%     Randomly choose the number of loads to be changed
    number_location = randi([1, length(mpc_case.bus)],1,1);
%     Randomly choose the locations of those loads
    locations = randperm(length(mpc_case.bus),number_location);
%     Get Multipliers between -0.5 and 1.5 and place them in the multiplier
%     matrix
    mult(1,locations) = (1.5+1.5)*rand(length(locations),1)-1.5;
%     Change real and reactive power 
    mpc_case.bus(:,PD) = base_case.bus(:,PD).*mult' ; 
    mpc_case.bus(:,QD) = base_case.bus(:,QD).*mult' ;
    for l = 1:length(locations)
        if mpc_case.bus(locations(l),PD) == 0
            mpc_case.bus(locations(l),PD) = avg_pd*mult(locations(l));
        end
        if mpc_case.bus(locations(l),QD) == 0
            mpc_case.bus(locations(l),QD) = avg_qd*mult(locations(l));
        end
        
    end
    
    [mpc_case_results,success] = runpf(mpc_case);
    if success == 1
        i = i+1;
        V(:,i) = mpc_case_results.bus(:,VM);
        Vangle(:,i) = mpc_case_results.bus(:,VA);
        PD_case(:,i) = mpc_case.bus(:,PD);
        QD_case(:,i) = mpc_case.bus(:,QD);
    end
end
j=1i;
Vcomplex = zeros(length(mpc_case.bus),number_of_scenarios);
for r = 1:length(mpc_case.bus)
    for c = 1:number_of_scenarios
        Vcomplex (r,c) = V(r,c) * exp(j*Vangle(r,c)/180);
    end
    
end


covaraince = Vcomplex * Vcomplex' / number_of_scenarios;
[E_l,D,E_R] = eig(covaraince);
eigen = log10(abs(sort(diag(D),'descend')));
subplot(321)
plot(eigen,'-x','linewidth',1.5)
ylabel('Sorted Eigen Value')
xlabel('Number of Eigen Values')
title('EigenValue plot for the Covariance Matrix')

subplot(322)
imagesc(V)
title('Plot for Voltage')

[U,S,V_1] = svd(Vcomplex);
subplot(323)
singular_value = log10(abs(sort(diag(S),'descend')));
plot(singular_value,'-o','linewidth',1.5)
ylabel('sorted Singular Values')
title('Singular Values')

subplot(324)
imagesc(Vcomplex*Vcomplex')
title('Plot for V*V^t')


subplot(325)
imagesc(PD_case + abs(min(PD_case)))
title('Plot for PD')

subplot(326)
imagesc(QD_case + abs(min(QD_case)))
title('Plot for QD')



f = mpc_case.branch(:,1);
t = mpc_case.branch(:,2);

[nmap, rmap, fnew, tnew] = NodeRelabling(f,t,1);
V_ordered = zeros(length(mpc_case.bus),number_of_scenarios);
Vangleordered = zeros(length(mpc_case.bus),number_of_scenarios);
for i = 1:length(mpc_case.bus(:,VM))
    V_ordered(nmap(i),:) = Vcomplex(i,:); 
    Vangleordered(nmap(i),:) = Vangle(i,:); 
end

covaraince = V_ordered * V_ordered' / number_of_scenarios;
figure
[E_l,D,E_R] = eig(covaraince);
eigen = sort(diag(D),'descend');
subplot(321)
plot(eigen,'-x','linewidth',1.5)
ylabel('Sorted Eigen Value')
xlabel('Number of Eigen Values')
title('EigenValue plot for the Covariance Matrix')

subplot(322)
imagesc(V_ordered)
title('Plot for Ordered Voltage')

[U,S,V_1] = svd(V_ordered);
subplot(323)
plot(sort(diag(S),'descend'),'-o','linewidth',1.5)
ylabel('sorted Singular Values')
title('Singular Values')

subplot(325)
imagesc(PD_case + abs(min(PD_case)))
title('Plot for PD')

subplot(326)
imagesc(QD_case + abs(min(QD_case)))
title('Plot for QD')

% figure 
% subplot(121)
% plot(graph(mpc_case.branch(:,1),mpc_case.branch(:,2)))
% title('original graph')
% 
% subplot(122)
% plot(graph(fnew,tnew))
% title('Relabeled graph')



[Y_bus,~,~] = makeYbus(mpc_case);
%%
% genbus = mpcint.gen(:,GEN_BUS); % this is the bus number for each generator
% % you might want to check that the generator is on:
% gen_status = mpcint.gen(:,GEN_STATUS) > 0;
% 
% % now you have only active generators:
% genbus = genbus(gen_status);
% 
% % now make sure there are no repeat indices:
% genbus = unique(genbus);
% % update ng
% ngold = ng;
% ng = length(genbus);
%% let's make a boolean mask for buses



%% change the U_D to complex orthogonal case. 
% 
% load('Gen_admittances_vector.mat');
% Y_gen = -1i*diag(Gen_admittances_vector(bool_mask)); % diagonal matrix
% 
% load('Delta_diag.mat');
% 
% load('M_gen_vec.mat');
% 
% M_all = diag(M_gen); % diagonal matrix

%% load 'impedances' (?)
nb = size(mpc_case.bus,1); %this is the number of buses
nl = size(mpc_case.branch,1); % this is the number of branches
ng = size(mpc_case.gen,1); % this is the number of generators (some might be switched off)
Y_L = mpc_case.bus(:,[PD,QD])./mpc_case.baseMVA;
Y_L = Y_L(:,1)-1i*Y_L(:,2);

Y_DD = full(Y_bus)+diag(Y_L); % diagonal matrix 

%% diagonal matrix \Delta
% Delta_matrix = diag([diag(Y_gen)+1j*Y_sh(bool_mask,2); Y_L(~bool_mask)+1j*Y_sh(~bool_mask,2)]);
% Delta_matrix = diag([diag(Y_gen)+1j*Y_sh(bool_mask,2);zeros(nb-ng,1)]);
%% Y+\Delta
[U_Y, Lambda_Y] = eig(Y_DD);%./mpcint.baseMVA);
[val_Y, ind_Y]=sort(diag(Lambda_Y),'ascend');
U_Y_arr = U_Y(:,ind_Y');

%% change into complex orthogonal matrices:
r_vect = zeros(nb,1);
for i=1:nb
        r_vect(i) = 1./sqrt(transpose(U_Y_arr(:,i))*U_Y_arr(:,i));
end
    
U_new = U_Y_arr*diag(r_vect);
[E_l_Y,D_Y,E_R_Y] = eig(full(Y_bus));
eigen_Y = log10(abs(sort(diag(D_Y),'ascend')));

plot(eigen_Y,'-r')
hold on 
plot(eigen,'-b')
legend('eigenY','eigenV')
[coeff,score,latent] = pca(abs(V_ordered)');




