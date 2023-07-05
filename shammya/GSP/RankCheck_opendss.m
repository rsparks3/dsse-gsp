clc;
clear all;
close all;
[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText = DSSCircObj.Text;
DSSText.command = 'clear';
DSSText.command = 'Compile C:\models\123Bus\IEEE123Master.DSS';
DSSText.Command = 'solve';
number_of_scenarios = 500;
Buses = getBusInfo(DSSCircObj);
base_results = Buses;
Base_Loads = getLoadInfo(DSSCircObj);
PD_case = zeros(length(Base_Loads),500);
QD_case = zeros(length(Base_Loads),500);
all_loads = struct();
BusVoltage = zeros(3*length(Buses),500);
BusVoltage_A = zeros(3*length(Buses),500);
BusVoltage_B = zeros(3*length(Buses),500);
BusVoltage_C = zeros(3*length(Buses),500);
for i = 1:500
    number_location = randi([1, length(Base_Loads)],1,1);
    locations = randperm(length(Base_Loads),number_location);
    mult = (1.5+0.5)*rand(length(locations),1)-0.5;
    for l = 1:length(locations)
        DSSCircuit.Loads.Name = Base_Loads(l).name;
        DSSCircuit.Loads.kW = mult(l)*Base_Loads(l).kW;
        DSSCircuit.Loads.kvar = mult(l)*Base_Loads(l).kvar;
    end
    load_sum = 0;
    loads = getLoadInfo(DSSCircObj);
    for ll = 1: length(loads)
        PD_case(ll,i) = loads(ll).kW;
        QD_case(ll,i) = loads(ll).kvar;
        load_sum = load_sum + loads(ll).kW;
    end
    fprintf('Total Load %f for iteration %d\n', load_sum,i);
    
    DSSCircuit.Solution.Solve();
    
    if (DSSCircuit.Solution.Converged == 1)
        Buses = getBusInfo(DSSCircObj);
        counter = 1;
        counter_A = 1;
        counter_B = 1;
        counter_C = 1;
        for b= 1:length(Buses)
            phase_voltage = Buses(b).phaseVoltagesPU .* exp(1i*Buses(b).voltagebaseAngle);
            for p = 1:3
                BusVoltage(counter,i) = phase_voltage(p);
                counter = counter+1;
            end
            phase_nodes = Buses(b).nodes;
            for p = 1:length(phase_nodes)
                if (phase_nodes(p)) == 1
                    BusVoltage_A(counter_A,i) = phase_voltage(phase_nodes(p));
                    counter_A = counter_A + 1 ;
                elseif (phase_nodes(p)) == 2
                    BusVoltage_B(counter_B,i) = phase_voltage(phase_nodes(p));
                    counter_B = counter_B + 1 ;
                elseif (phase_nodes(p)) == 3
                    BusVoltage_C(counter_C,i) = phase_voltage(phase_nodes(p));
                    counter_C = counter_C + 1 ;
                end
            end
        end
    end
    
    
    
end
BusVoltage_A = BusVoltage_A(1:counter_A-1,:);
BusVoltage_B = BusVoltage_B(1:counter_B-1,:);
BusVoltage_C = BusVoltage_C(1:counter_C-1,:);
%% Clear the zero rows
missing_rows = [];
for j= 1:counter-1
    if (abs(sum(BusVoltage(j,:))))==0
        missing_rows = [missing_rows j];
    end
end
BusVoltage(missing_rows,:) = [];
%%

V= BusVoltage;
covaraince = V * V' / number_of_scenarios;
[E_l,D,E_R] = eig(covaraince);
eigen = log10(abs(sort(diag(D),'descend')));
subplot(321)
plot(eigen,'-x','linewidth',1.5)
ylabel('Sorted Eigen Value')
xlabel('Number of Eigen Values')
title('EigenValue plot for the Covariance Matrix')

subplot(322)
imagesc(abs(V))
title('Plot for Voltage')

[U,S,V_1] = svd(V);
subplot(323)
plot(log10(sort(diag(S),'descend')),'-o','linewidth',1.5)
ylabel('sorted Singular Values')
title('Singular Values')

subplot(324)
imagesc(abs(V*V'))
title('Plot for V*V^t')


subplot(325)
imagesc(PD_case + abs(min(PD_case)))
title('Plot for PD')

subplot(326)
imagesc(QD_case + abs(min(QD_case)))
title('Plot for QD')


%%
base_y = DSSCircuit.SystemY;
nb = sqrt(length(base_y)/2);
Y = zeros(nb,nb);
k= 1;

for i = 1:nb
    for j = 1:nb
        Y(i,j)= base_y(1,k)+1i*base_y(1,k+1);
        k = k+2;
    end
end

bus_a = [];
bus_b = [] ;
bus_c = [];
y_node_order = DSSCircuit.YNodeOrder;
ordered_bus_names = cell(numel(y_node_order),1);
for i = 1:numel(y_node_order)
    bus_name = y_node_order{i};
    bus = split(bus_name,'.');
    if strcmp(bus{2},'1')
        bus_a = [bus_a i];
    elseif strcmp(bus{2},'2')
        bus_b = [bus_b i];
    elseif strcmp(bus{2},'3')
        bus_c = [bus_c i];
    end
    ordered_bus_names{i} = bus{1};
end
ordered_bus_names = unique(ordered_bus_names,'stable');
Y_a = Y(bus_a',:);
Y_a = Y_a(:,bus_a);

Y_b = Y(bus_b',:);
Y_b = Y_b(:,bus_b);

Y_c = Y(bus_c',:);
Y_c = Y_c(:,bus_c);

%%
Y_DD = Y; % diagonal matrix

[U_Y, Lambda_Y] = eig(Y_DD);%./mpcint.baseMVA);
[val_Y, ind_Y]=sort(diag(Lambda_Y),'ascend');
U_Y_arr = U_Y(:,ind_Y');

% change into complex orthogonal matrices:
r_vect = zeros(length(Y_DD),1);
for i=1:length(Y_DD)
    r_vect(i) = 1./sqrt(transpose(U_Y_arr(:,i))*U_Y_arr(:,i));
end

U_Y_new = U_Y_arr*diag(r_vect);
fprintf('Max Difference in estimate for Y:%f\n',max(abs(U_Y_new * diag(val_Y) * transpose(U_Y_new) - Y_DD),[],'all'));
figure(1)
plot(log10(abs(val_Y)))


[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText = DSSCircObj.Text;
DSSText.command = 'clear';
DSSText.command = 'Compile C:\models\123Bus\IEEE123Master.DSS';
DSSText.Command = 'solve';
Buses = getBusInfo(DSSCircObj,ordered_bus_names);

voltage_magnitude = [];
full_voltage_magnitude = [];
for i = 1:length(Buses)
    v = Buses(i).phaseVoltages.* exp(1i*Buses(i).voltagebaseAngle);
    voltage_magnitude = [voltage_magnitude; (v(Buses(i).nodes))'];
    full_voltage_magnitude = [full_voltage_magnitude; v' ];
end

number_of_eigenvectors = 250;
U_k = U_Y_new (:,1:number_of_eigenvectors);
v_new = (U_k)*pinv(U_k)*voltage_magnitude;

figure(2)
plot(abs(v_new))
hold on
plot(abs(voltage_magnitude))

%% Filling the Y matrix
% a = [1 2 3 4 5 6];
% pos = [ 2 4];
% b = zeros(length(a)+length(pos),1);
% c2 = 1;
% c1 = 1;
% while (c2 <= length(b))
%     if ~ismember(c2,pos)
%         b(c2) = a(c1);
%         c1 = c1+1;
%     end
%     c2 = c2 + 1 ;
% end
total_nodes = [1 2 3];
insert_location = [];
total_length = numel(ordered_bus_names)*3;
for i = 1:numel(ordered_bus_names)
    loc = i*3-2:i*3;
    diff = loc(setdiff(total_nodes,Buses(i).nodes));
    insert_location = [insert_location diff];
end

full_Y = zeros(total_length,length(Y));
for i = 1:length(Y)
    a = Y(:,i);
    b = zeros(length(a)+length(insert_location),1);
    c2 = 1;
    c1 = 1;
    while (c2 <= length(b))
        if ~ismember(c2,insert_location)
            b(c2) = a(c1);
            c1 = c1+1;
        end
        c2 = c2 + 1 ;
    end
    full_Y(:,i)=b;
%     full_Y(i,:)=b';
end
c1 = 1 ;
c2 = 1 ;
Complete_Y = zeros(total_length,total_length);
while (c2 <= total_length)
    if ~ ismember(c2,insert_location)
        Complete_Y(:,c2) = full_Y(:,c1);
        c1 = c1+1;
    end
    c2 = c2+1;
end

Y_DD = Complete_Y; % diagonal matrix

[U_Y, Lambda_Y] = eig(Y_DD);%./mpcint.baseMVA);
[val_Y, ind_Y]=sort(diag(Lambda_Y),'ascend');
U_Y_arr = U_Y(:,ind_Y');

% change into complex orthogonal matrices:
r_vect = zeros(length(Y_DD),1);
for i=1:length(Y_DD)
    r_vect(i) = 1./sqrt(transpose(U_Y_arr(:,i))*U_Y_arr(:,i));
end

U_Y_new = U_Y_arr*diag(r_vect);

fprintf('Max Difference in estimate for Y:%f\n',max(abs(U_Y_new * diag(val_Y) * transpose(U_Y_new) - Y_DD),[],'all'));
figure(4)
plot(log10(abs(val_Y)),'linewidth',1.5)
title('Eigen Values Sorted')


U_k = U_Y_new (:,1:360) ;
full_v_new = (U_k)*pinv(U_k)*full_voltage_magnitude;

figure(3)
plot(abs(full_v_new),'linewidth',2)
hold on
plot(abs(full_voltage_magnitude),'linewidth',2)
title('interpolation with all measurements ')
legend('interpolated','Original')


%% just Place Sensors
three_phase_sensor_locations = randperm(numel(ordered_bus_names));
three_phase_sensor_locations = three_phase_sensor_locations(1:70);
sensor_locations = [];
for i = 1:length(three_phase_sensor_locations)
    sensor_locations = [sensor_locations three_phase_sensor_locations(i)*3-2:three_phase_sensor_locations(i)*3];
end
U_k = U_Y_new (:,1:380) ;
number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,total_length);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end
full_v_new_with_sensor = (U_k)*pinv(M*U_k)*full_voltage_magnitude(sensor_locations);

figure(5)
plot(abs(full_v_new_with_sensor),'linewidth',2)
hold on
plot(abs(full_voltage_magnitude),'linewidth',2)
legend('interpolated','Original')
title('Interpolation with Selected Sensors ')