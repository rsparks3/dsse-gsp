clc;
clear all;
close all;
eigenvectors_reduction = 00;
BusSize = 123;
ADD_LOAD = 00;
[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText.command = 'clear';
%%
[file_directory,slack_bus_name]= ChooseOpenDSSModel(BusSize);
DSSText.command = file_directory;
DSSText.Command = 'vsource.source.enabled=no';
DSSText.command = 'batchedit load..* enabled=no';
DSSText.Command = 'solve';

base_y = DSSCircuit.SystemY;
Ymatrix=reshape(base_y,[length(DSSCircuit.AllNodeNames)*2,length(DSSCircuit.AllNodeNames)]); % "2" is because the real/imaginary parts are separate
Ymatrix=Ymatrix';
Y_dss = Ymatrix(:,1:2:end) + sqrt(-1)*Ymatrix(:,2:2:end);

ynode_order = DSSCircuit.YNodeOrder;
[sorted_y_node,sort_index] =  sort(ynode_order);

sorted_y_dss= Y_dss(sort_index,:);
sorted_y_dss = sorted_y_dss(:,sort_index);

bus_names = cell(size(sorted_y_node));
for j = 1:length(sorted_y_node)
    name = split(sorted_y_node(j),'.');
    bus_names{j} = name{1};
end

bus_name_loc=containers.Map();
unique_bus_names = unique(bus_names,'stable');

for i = 1:length(unique_bus_names)
    bus_name_loc(unique_bus_names{i})=find(strcmpi(bus_names,unique_bus_names{i}));
end


%%
% Recompile the model to get a valid power flow
DSSText.command = file_directory;
DSSText.Command = 'solve';
% Retrieve all the voltages from the OpenDSS result
AllBusVoltage = DSSCircuit.AllbusVolts;
vcomplex = AllBusVoltage(1:2:end) + 1i * AllBusVoltage(2:2:end) ; % unit: Volt
vcomplex = transpose(vcomplex);
% sort the voltages the same way we sorted the y bus

sorted_vcomplex=vcomplex(sort_index);
% Calculate the original bus injection
s_original = diag(vcomplex) * conj(Y_dss*vcomplex)/1000;
% calculate the sorted bus injection values
s_sorted = diag(sorted_vcomplex) * conj(sorted_y_dss*sorted_vcomplex)/1000;
% Run a sanity check to ensure that the power calculations are correct
disp(max(abs(s_sorted-s_original(sort_index))));






% Lets run the interpolation equation to see whether the sorted y_dss works or not
Y_DD = sorted_y_dss; % diagonal matrix

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
BusInfo = getBusInfo(DSSCircObj,unique_bus_names);
phase_distribution = zeros(length(unique_bus_names),3);
base_voltage = 1000*[BusInfo(:).kVBase]';

for i = 1:length(unique_bus_names)
    ph = BusInfo(i).nodes;
    mapping = bus_name_loc(upper(unique_bus_names{i}));
    phase_distribution(i,ph) = mapping;
end
%%
number_of_eigenvectors = length(Y_DD)-eigenvectors_reduction;
% vcomplex = full_voltage_magnitude;
U_k = U_Y_new (:,1:number_of_eigenvectors) ;
v_new = U_k*transpose(U_k)*sorted_vcomplex;
figure(3)
% subplot(311)
% loc = find(abs(full_voltage_magnitude(loc))>0);
subplot(211)
plot(abs(v_new),'linewidth',1.5)
hold on
plot(abs(sorted_vcomplex),'linewidth',1.5)
legend('Interpolated','Original')
subplot(212)
plot(angle(v_new),'linewidth',1.5)
hold on
plot(angle(sorted_vcomplex),'linewidth',1.5)
legend('Interpolated','Original')




%% Sensor placement

Yuk = sorted_y_dss * U_Y_new;
Yuk = Yuk(:,1:number_of_eigenvectors);
M_s = [];
M= 92;
M_tilde = M-length(M_s); % total number to be found
[N,~]=size(U_k);
i=1;
j_tilde = 1:length(sorted_y_dss);
j_tilde(M_s)=[];
while i <= M_tilde
    sigma_min_uk=zeros(length(j_tilde),1);
    sigma_min_yuk=zeros(length(j_tilde),1);
    for j=1: length(j_tilde)
        A = U_k([M_s;j_tilde(j)],:);
        [~,sigma_min_uk(j),~,flag] = svds(A,1,'smallest');
        if(flag)
            continue;
        end
        B = Yuk([M_s;j_tilde(j)],:);
        [~,sigma_min_yuk(j),~,flag] = svds(B,1,'smallest');
        if(flag)
            continue;
        end
    end
%         [~,min_ind]=max(abs(sigma_min_uk)+abs(sigma_min_yuk));
    [~,min_ind]= max(min([sigma_min_uk sigma_min_yuk],[],2));
    M_s = [M_s;j_tilde(min_ind)];
    j_tilde(min_ind)=[];
    i=i+1;
end
P_S = eye(N);
P_S=P_S(:,M_s');
bus_name_sensor = unique(bus_names(M_s));

sensor_locations = zeros(length(bus_names),1);

for i = 1:length(bus_name_sensor)
    x = find(strcmp(bus_names,bus_name_sensor(i)));
    sensor_locations(x,1)=1;
end
%% Construct the graph for node relabeling
from_bus_nodes = [];
to_bus_nodes = [];
lines = getLineInfo(DSSCircObj);
transformers = getTransformerInfo(DSSCircObj);
bus_names = cell(length(lines)*2 + length(transformers)*2,1);
counter = 1;
for i = 1:length(lines)
    name= split(lines(i).bus1,'.');
    
    bus_names{counter} = name{1};
    counter = counter +1 ;
    name = split(lines(i).bus2,'.');
    bus_names{counter} =name{1};
    counter = counter +1 ;
end

for i = 1:length(transformers)
    name = split(transformers(i).bus1,'.');
    bus_names{counter} = name{1};
    counter = counter +1 ;
    name = split(transformers(i).bus2,'.');
    bus_names{counter} = name{1};
    counter = counter +1 ;
end


bus_info=getBusInfo(DSSCircObj,unique_bus_names);
inner_map = containers.Map(unique_bus_names,1:length(unique_bus_names));
% looping through the lines for finding all the from and nodes
for i = 1:length(lines)
    name= split(lines(i).bus1,'.');
    from_bus_loc = inner_map(upper(name{1}));
    from_bus_nodes = [from_bus_nodes from_bus_loc];
    name= split(lines(i).bus2,'.');
    to_bus_loc = inner_map(upper(name{1}));
    to_bus_nodes = [to_bus_nodes to_bus_loc];
end
% Append the from and to nodes of transformers to the corresponding list
for i = 1:length(transformers)
    name= split(transformers(i).bus1,'.');
    from_bus_loc = inner_map(upper(name{1}));
    from_bus_nodes = [from_bus_nodes from_bus_loc];
    name= split(transformers(i).bus2,'.');
    to_bus_loc = inner_map(upper(name{1}));
    to_bus_nodes = [to_bus_nodes to_bus_loc];
end
%% Process the graph 
remove_list=[];
for i = 1:length(from_bus_nodes)
    if ~ismember(i,remove_list)
        f=from_bus_nodes(i);
        t= to_bus_nodes(i);
        for j = 1:length(from_bus_nodes)
            if j ~=i
                if (from_bus_nodes(j)==f && to_bus_nodes(j)==t)
                    remove_list = [remove_list j];
                end
            end
        end
    end
end

from_bus_nodes(remove_list) = [];
to_bus_nodes(remove_list) = [];

[nmap, rmap, fnew, tnew] = NodeRelabling(from_bus_nodes',to_bus_nodes',inner_map(upper(slack_bus_name)));

%%
node_ordering_opendss = DSSCircuit.YNodeOrder;
node_ordering_opendss_names = cell(size(node_ordering_opendss));
G = digraph(fnew,tnew);

for i = 1:length(node_ordering_opendss)
    name= split(node_ordering_opendss(i),'.');
    node_ordering_opendss_names{i} = name{1};
end
opendss_bus_names = DSSCircuit.AllBusNames;
testkeys = keys(inner_map);

nb = length(sorted_y_dss);
number_of_sensors = sum(sensor_locations);
sensor_locations= find(sensor_locations==1);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end







%% Covex Relaxation Code
% clear tildeW

sm = s_sorted*1000;
sm = sm(sensor_locations);
vm = abs(sorted_vcomplex);
vm=vm(sensor_locations);
c = sorted_y_dss*sorted_vcomplex;
cm = abs(c);
cm = cm(sensor_locations);

save('C:\Users\asus\Dropbox (ASU)\Shammya\Projects\Neptune 2.0\GSP\37bus.mat')



