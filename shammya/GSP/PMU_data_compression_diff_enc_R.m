%% compression of phasor data: 
clear all
%close all
clc;
%%
addpath(genpath('matpower'));
define_constants; % this creates a bunch of variables to index into mpc matrices.
fname = 'case_ACTIVSg2000';
mpc = loadcase(fname);
% for some of the commands to work you need to have all bus numbers
% consecutive starting at 1, etc. ext2int does that for you.
mpcint = ext2int(mpc);
%% some useful numbers to have around
nb = size(mpc.bus,1); %this is the number of buses
nl = size(mpc.branch,1); % this is the number of branches
ng = size(mpc.gen,1); % this is the number of generators (some might be switched off)
%% This is the laplacian
B = makeBdc(mpcint);
%% This is not Laplacian but is only approximately so
Y_bus = makeYbus(mpcint);
%%
genbus = mpcint.gen(:,GEN_BUS); % this is the bus number for each generator
% you might want to check that the generator is on:
gen_status = mpcint.gen(:,GEN_STATUS) > 0;

% now you have only active generators:
genbus = genbus(gen_status);

% now make sure there are no repeat indices:
genbus = unique(genbus);
% update ng
ngold = ng;
ng = length(genbus);
%% let's make a boolean mask for buses
bool_mask = false(nb,1);
bool_mask(genbus) = true;
%% shunt elements
Y_sh  = mpcint.bus(:,[GS,BS])./mpcint.baseMVA ;
%% remove shunt elements from diagonal of Ybus and add it to the load values. 
Y_bus = Y_bus - diag(1j*Y_sh(:,2));

%% separate into different blocks
L_11 = Y_bus(bool_mask, bool_mask); 
L_12 = Y_bus(bool_mask, ~bool_mask);
L_21 = Y_bus(~bool_mask, bool_mask);
L_22 = Y_bus(~bool_mask, ~bool_mask);
Bord = [L_11,L_12;
    L_21,L_22];

Bord = full(Bord);


%% change the U_D to complex orthogonal case. 

load('Gen_admittances_vector.mat');
Y_gen = -1i*diag(Gen_admittances_vector(bool_mask)); % diagonal matrix

load('Delta_diag.mat');

load('M_gen_vec.mat');

M_all = diag(M_gen); % diagonal matrix

%% load 'impedances' (?)
Y_L = mpcint.bus(:,[PD,QD])./mpcint.baseMVA;
Y_L = Y_L(:,1)-1i*Y_L(:,2);

Y_DD = diag(Y_L(~bool_mask)); % diagonal matrix 

%% diagonal matrix \Delta
Delta_matrix = diag([diag(Y_gen)+1j*Y_sh(bool_mask,2); Y_L(~bool_mask)+1j*Y_sh(~bool_mask,2)]);
Delta_matrix = diag([diag(Y_gen)+1j*Y_sh(bool_mask,2);zeros(nb-ng,1)]);
%% Y+\Delta
[U_Y, Lambda_Y] = eig(full(Bord+Delta_matrix));%./mpcint.baseMVA);
[val_Y, ind_Y]=sort(diag(Lambda_Y),'ascend');
U_Y_arr = U_Y(:,ind_Y');

%% change into complex orthogonal matrices:
r_vect = zeros(nb,1);
for i=1:nb
        r_vect(i) = 1./sqrt(transpose(U_Y_arr(:,i))*U_Y_arr(:,i));
end
    
    U_new = U_Y_arr*diag(r_vect);

%% data:
%load('ACTIVSg_downsampled.mat');
load('ACTIVSg_5k_samples.mat');
va = unwrap(angle_chopped*pi/180,[],2);
voltage_phasors  = v_chopped.*exp(1i*va);
%GFT_U_D = (U_Y_arr)\voltage_phasors;
voltage_phasors = voltage_phasors(:,1:1000);
%% input e_t: 
E = inv(Y_gen)*diag(1j*Y_sh(bool_mask,2))*voltage_phasors(1:ng,:) + voltage_phasors(1:ng,:)+ inv(Y_gen)*(L_11*voltage_phasors(1:ng,:)+ L_12*voltage_phasors(ng+1:end,:));
ip = [Y_gen*E; zeros(nb-ng,1000)];

%nu_v = abs(E);
%% correct Y_gen
% to_be_corrected = mean(nu_v,2);
% to_be_corrected(to_be_corrected<1)=1;
% Y_gen_new = Y_gen*diag(to_be_corrected);
% 
% E_new = voltage_phasors(1:ng,:)+ inv(Y_gen_new)*(L_11*voltage_phasors(1:ng,:)+ L_12*voltage_phasors(ng+1:end,:));
% ip_new = [Y_gen_new*E_new; zeros(nb-ng,300)];
% 
% Delta_matrix_new = diag([diag(Y_gen_new)+1j*Y_sh(bool_mask,2); Y_L(~bool_mask)+1j*Y_sh(~bool_mask,2)]);

% [U_Y_new, Lambda_Y_new] = eig(full(Bord+Delta_matrix_new));%./mpcint.baseMVA);
% [val_Y_new, ind_Y_new]=sort(diag(Lambda_Y_new),'ascend');
% U_Y_arr_new = U_Y_new(:,ind_Y_new');



%% Y_red matrix:
Q_red = Y_gen-Y_gen*inv(L_11+Y_gen - L_12*inv(L_22+(Y_DD))*transpose(L_12))*Y_gen;
Gamma = -imag(Q_red);
Q_tilde = ( Gamma - diag(Gamma*ones(length(Gamma),1)));

%% new S_red used:
S_red = diag(1./sqrt(M_gen))*Q_tilde*diag(1./sqrt(M_gen));
[U_red,Lam_red]=eig(S_red);
[sort_val_red,sort_idx_red]=sort(diag(Lam_red),'ascend');
U_red = U_red(:,sort_idx_red);


[U_P,Lam_P]=eig(inv(diag(M_gen))*Q_tilde);
[sort_val_P,sort_idx_P]=sort(diag(Lam_P),'ascend');
U_P = U_P(:,sort_idx_P); 



%% try the reconstruction here: 
 K = 2000;
%ip = full(Bord+Delta_matrix)*voltage_phasors;
ip = U_new(:,1:K)*inv(diag(val_Y(1:K)))*transpose(U_new(1:ng,:))\voltage_phasors;
%U_D_inv  = inv(U_Y_arr_new);
Bases_A = U_new(:,1:K)*inv(diag(val_Y(1:K)))*transpose(U_new(1:ng,:));%+(0-1i*(0.1))*eye(K)

volt_new = Bases_A*ip;
% volt_new  = full(Bord+Delta_matrix)\ip;
figure;
plot(real(voltage_phasors(:,1)),'-*');
hold on
% plot(real(volt_recons(:,1)),'-o'); 
 plot(real(volt_new(:,1)),'-o'); 

%% input angles:
% amplitude = abs([ip(1:ng,:)]); 
% amplitude_nu = log(amplitude);
% phase_new = angle(ip(1:ng,:));

E_g = inv(Y_gen)*ip;

x_t = diag(sqrt(M_gen))*log(E_g);
increasing_order_GFT = U_red\x_t;
%increasing_order_GFT = (U_P)\log(abs(ip(1:ng,:)));
%increasing_order_GFT = (U_P)\(angle(ip(1:ng,:)));

figure;
subplot(2,1,1);
plot(sort_val_red,real(increasing_order_GFT(:,1)),'-*');
subplot(2,1,2);
plot(sort_val_red,imag(increasing_order_GFT(:,1)),'-*');

%% parameters..:
% h1_new  = 2 - (1/10)^2*(sort_val_P)-(1/10)*0.1;
% h2_new = -1+(1/10)*0.1;

h1_all_amp = []; h2_all_amp=[];

collect_sigma_amp = [];
collect_sigma_old=[];
collect_func = [];
%omega = -pi:pi/100:pi;
omega_1 = 0:pi/100:2*pi;
z = exp(-1i.*omega_1);
T_end = 700;
%% approx:
figure;
for i=1:length(sort_val_red)
    % y = iddata(transpose(increasing_order_GFT(i,1:200)),[],(1/10));
    % sys = armax(y,[2 0]);
    y = transpose(increasing_order_GFT(i,3:T_end));
    A_matrix = [transpose(increasing_order_GFT(i,1:T_end-2)), transpose(increasing_order_GFT(i,2:T_end-1))];
    hh = A_matrix\y;
    
    yy = A_matrix * hh;

    h1_all_amp =[h1_all_amp;hh(2)];
    h2_all_amp = [h2_all_amp;hh(1)];
    
%     y = transpose(increasing_order_GFT(i,3:200));
%     A_matrix = [transpose(increasing_order_GFT(i,1:198)), transpose(increasing_order_GFT(i,2:199))];
%    yy = A_matrix * [h2_new;h1_new(i)];
    sig_new  = norm(yy-y,'fro').^2;
    collect_sigma_amp=[collect_sigma_amp;sig_new];
%    

%    yy_old = A_matrix * [h2_new;h1_new(i)];
%    sig_old = norm(yy_old-y,'fro').^2;
%    collect_sigma_old = [collect_sigma_old;sig_old];
   
   func = (sig_new./length(yy))./abs(1-hh(2)*z-hh(1)*z.^2).^2;
   collect_func = [collect_func;10*log10(func)];
   %plot(omega./(2*pi),10*log10(func./max(func)));
   %hold on
%    figure;
%    plot(imag(y),'-*');
%    hold on
%    plot(imag(yy),'-o');
%    plot(yy_old,'-o');
%    close all;

% figure;
% plot(y-yy,'-*');
% close all;

end
xlabel('f');
ylabel('$\phi(\omega)$');
% figure;
% subplot(2,1,1);
% plot(collect_sigma,'-*');
% subplot(2,1,2);
% plot((M_gen),'-o')

%% plot3d
figure;
mesh(omega_1'./(2*pi),sort_val_red,collect_func);
xlabel('\omega/(2\pi)');
ylabel('\lambda_{red}');
zlabel('10log(\phi(\omega))');
title('2-D frequency response in time and graph frequency domains');
%%
figure;
subplot(2,1,1);
plot(sort_val_red,real(h1_all_amp),'*');
subplot(2,1,2);
plot(sort_val_red,imag(h1_all_amp),'*');

figure;
subplot(2,1,1);
plot(sort_val_red,real(h2_all_amp),'*');
subplot(2,1,2);
plot(sort_val_red,imag(h2_all_amp),'*');

%% do the polynomial fitting
[a_1_real,S1_real] = polyfit(sort_val_red, real(h1_all_amp),1);
a_1_real_fit  = polyval(a_1_real,sort_val_red,[]);

figure;
plot(sort_val_red,real(h1_all_amp),'o');
hold on
plot(sort_val_red,a_1_real_fit,'-*');

a_1_imag = polyfit(sort_val_red, imag(h1_all_amp),1);
a_1_imag_fit  = polyval(a_1_imag,sort_val_red,[]);


a_2_real = polyfit(sort_val_red, real(h2_all_amp),1);
a_2_real_fit  = polyval(a_2_real,sort_val_red,[]);

a_2_imag = polyfit(sort_val_red, imag(h2_all_amp),1);
a_2_imag_fit  = polyval(a_2_imag,sort_val_red,[]);

%% look at the fit:
hh_new = [a_2_real_fit+1i*a_2_imag_fit, a_1_real_fit+1i*a_1_imag_fit];

%% tested norm:
norm_actual = []; norm_fit=[]; sigma_new = [];
TT = 200;
T_end = 700;
collect_func_fit=[];
for i=1:length(a_1_real_fit)
y = transpose(increasing_order_GFT(i,T_end+3:T_end+TT));
A_matrix = [transpose(increasing_order_GFT(i,T_end+1:T_end-2+TT)), transpose(increasing_order_GFT(i,2+T_end:T_end-1+TT))];
hh = A_matrix\y;
    
 yy = A_matrix * hh;
 norm_actual = [norm_actual; (1/200)*norm(y-yy,2).^2];
 yy_new = A_matrix*transpose(hh_new(i,:));
 norm_fit = [norm_fit; (1/200)*norm(y-yy_new,2).^2];
 
 func = (sig_new./length(yy))./abs(1-hh_new(i,2)*z-hh_new(i,1)*z.^2).^2;
   collect_func_fit = [collect_func_fit;10*log10(func)];
%  figure;
%  plot(real(yy),'-*');
%  hold on
%  %plot(imag(yy),'-o');
%  plot(real(yy_new),'-o');
    
end
figure;
plot(abs(norm_actual-norm_fit),'-*');
%% poles:
sys = tf([1,0,0],[1, -h1_all_amp(10), -h2_all_amp(10) ],0.1);
figure;
pzmap(sys);
sys2 = tf([1,0,0],[1, -hh_new(10,2), -hh_new(10,1) ],0.1);
figure;
pzmap(sys2);

all_poles_old=zeros(length(hh_new),2);
all_poles = zeros(length(hh_new),2);
for i=1:length(hh_new)
    sys2 = tf([1,0,0],[1, -h1_all_amp(i), -h2_all_amp(i) ],0.1);
    all_poles_old(i,:)  = transpose(pole(sys2));
    sys2 = tf([1,0,0],[1, -hh_new(i,2), -hh_new(i,1) ],0.1);
    all_poles(i,:)  = transpose(pole(sys2));
end

figure;
plot3(sort_val_red,real(all_poles_old(:,2)),imag(all_poles_old(:,2)),'o');
hold on
plot3(sort_val_red,real(all_poles(:,2)),imag(all_poles(:,2)),'*-');
grid on

figure;
plot3(sort_val_red,real(all_poles_old(:,1)),imag(all_poles_old(:,1)),'o');
hold on
plot3(sort_val_red,real(all_poles(:,1)),imag(all_poles(:,1)),'*-');
grid on


