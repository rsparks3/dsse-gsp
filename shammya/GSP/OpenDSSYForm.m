clc;
clear ;
close all;
number_of_eigenvectors = 60;
ADD_LOAD = 00;
[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText.command = 'clear';
%%
DSSText.command = 'Compile C:\models\123Busm\IEEE123Master.DSS';
% DSSText.command = 'Compile C:\models\37Bus\ieee37.DSS';
% DSSText.command = 'Compile C:\models\37Bus\ieee37.DSS';
% DSSText.command = 'Compile C:\models\34Bus\ieee34Mod2.dss';
% DSSText.command = 'Compile C:\models\13Bus\IEEE13Nodeckt.dss';
% DSSText.command = 'Compile C:\models\test3bus\3bus.dss';
% DSSText.command = 'Compile C:\models\8500-Node\Master-unbal.dss';
DSSText.Command = 'vsource.source.enabled=no';
DSSText.command = 'batchedit load..* enabled=no';
DSSText.Command = 'solve';

base_y = DSSCircuit.SystemY;
Ymatrix=reshape(base_y,[length(DSSCircuit.AllNodeNames)*2,length(DSSCircuit.AllNodeNames)]); % "2" is because the real/imaginary parts are separate
Ymatrix=Ymatrix';
Y_dss = Ymatrix(:,1:2:end) + sqrt(-1)*Ymatrix(:,2:2:end);

%%
DSSText.command = 'Compile C:\models\123Busm\IEEE123Master.DSS';
% DSSText.command = 'Compile C:\models\37Bus\ieee37.DSS';
% DSSText.command = 'Compile C:\models\37Bus\ieee37.DSS';
% DSSText.command = 'Compile C:\models\34Bus\ieee34Mod2.dss';
% DSSText.command = 'Compile C:\models\13Bus\IEEE13Nodeckt.dss';
% DSSText.command = 'Compile C:\models\test3bus\3bus.dss';
% DSSText.command = 'Compile C:\models\8500-Node\Master-unbal.dss';

DSSText.Command = 'vsource.source.enabled=no';
DSSText.command = 'batchedit load..* enabled=yes';
DSSText.Command = 'solve';

base_y = DSSCircuit.SystemY;
Ymatrix=reshape(base_y,[length(DSSCircuit.AllNodeNames)*2,length(DSSCircuit.AllNodeNames)]); % "2" is because the real/imaginary parts are separate
Ymatrix=Ymatrix';
Y_dss_load = Ymatrix(:,1:2:end) + sqrt(-1)*Ymatrix(:,2:2:end);





%%
DSSText.command = 'Compile C:\models\123Busm\IEEE123Master.DSS';
% DSSText.command = 'Compile C:\models\37Bus\ieee37.DSS';
% DSSText.command = 'Compile C:\models\37Bus\ieee37.DSS';
% DSSText.command = 'Compile C:\models\34Bus\ieee34Mod2.dss';
% DSSText.command = 'Compile C:\models\13Bus\IEEE13Nodeckt.dss';
% DSSText.command = 'Compile C:\models\test3bus\3bus.dss';
% DSSText.command = 'Compile C:\models\8500-Node\Master-unbal.dss';
DSSText.Command = 'solve';
lines = getLineInfo(DSSCircObj);
transformers = getTransformerInfo(DSSCircObj);
capacitors = getCapacitorInfo(DSSCircObj);

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

bus_names = unique(bus_names,'stable');
inner_map = containers.Map(bus_names,1:length(bus_names));
Y = zeros(length(bus_names)*3,length(bus_names)*3);

% Phase Check
for i = 1:length(lines)
    name= split(lines(i).bus1,'.');
    
end

%%
from_bus_nodes = [];
to_bus_nodes = [];
for i = 1:length(lines)
    name= split(lines(i).bus1,'.');
    from_bus_loc = inner_map(name{1});
    from_bus_nodes = [from_bus_nodes from_bus_loc];
    name= split(lines(i).bus2,'.');
    to_bus_loc = inner_map(name{1});
    to_bus_nodes = [to_bus_nodes to_bus_loc];
    line_y = lines(i).Yprim;
    numphases = lines(i).numPhases;
    phases = lines(i).phases;
    line_y = transpose(reshape(line_y,[numphases*2*2,numphases*2]));
    line_y = line_y(1:numphases,1:numphases*2);
    y_matrix = [];
    for j=1:2:numphases*2
        y_matrix = [y_matrix line_y(:,j)+1i*line_y(:,j+1)];
    end
    temp = zeros(3,3);
    if sum(phases) == 3
        temp = y_matrix;
    elseif sum(phases) == 1
        ph = find(lines(i).phases==1);
        if length(y_matrix) ~= 1
            fprintf('Issues for line:%d\n',i);
        else
            temp(ph,ph)= y_matrix;
        end
    elseif sum(phases) == 2
        ph = find(lines(i).phases==1);
        temp(ph(1),ph(1)) = y_matrix(1,1);
        temp(ph(1),ph(2)) = y_matrix(1,2);
        temp(ph(2),ph(1)) = y_matrix(2,1);
        temp(ph(2),ph(2)) = y_matrix(2,2);
    end
    
    r_c = 1;
    for k1= from_bus_loc*3-2:from_bus_loc*3
        Y(k1,to_bus_loc*3-2:to_bus_loc*3) = Y(k1,to_bus_loc*3-2:to_bus_loc*3)-temp(r_c,:);
        r_c = r_c+1;
    end
    r_c = 1;
    for k1= to_bus_loc*3-2:to_bus_loc*3
        Y(k1,from_bus_loc*3-2:from_bus_loc*3) = Y(k1,from_bus_loc*3-2:from_bus_loc*3)-transpose(temp(:,r_c));
        r_c = r_c+1;
    end
    r_c = 1;
    for k1= from_bus_loc*3-2:from_bus_loc*3
        Y(k1,from_bus_loc*3-2:from_bus_loc*3) = Y(k1,from_bus_loc*3-2:from_bus_loc*3) +  temp(r_c,:);
        r_c = r_c+1;
    end
    r_c = 1;
    for k1= to_bus_loc*3-2:to_bus_loc*3
        Y(k1,to_bus_loc*3-2:to_bus_loc*3) = Y(k1,to_bus_loc*3-2:to_bus_loc*3) +  temp(r_c,:);
        r_c = r_c+1;
    end
end

for i = 1:length(transformers)
    name= split(transformers(i).bus1,'.');
    from_bus_loc = inner_map(name{1});
    from_bus_nodes = [from_bus_nodes from_bus_loc];
    name= split(transformers(i).bus2,'.');
    to_bus_loc = inner_map(name{1});
    to_bus_nodes = [to_bus_nodes to_bus_loc];
    line_y = transformers(i).Yprim;
    numphases = transformers(i).numPhases+1;
    phases = transformers(i).phases;
    line_y = transpose(reshape(line_y,[numphases*2*2,numphases*2]));
    numphases = numphases-1;
    [~,c] = size(line_y);
    transformer_y=[];
    for j = 1:2:c
        transformer_y=[transformer_y line_y(:,j)+1i*line_y(:,j+1)];
    end
    
    %     line_y = line_y(1:numphases,1:numphases*2);
    
    if sum(phases) == 3
        transformer_y([numphases+1 (numphases+1)*2],:)=[];
        transformer_y(:,[numphases+1 (numphases+1)*2])=[];
        [r,c] = size(transformer_y);
        
        from_from_y_matrix = transformer_y(1:r/2,1:c/2);
        to_to_y_matrix = transformer_y(numphases+1:r,numphases+1:c);
        from_to_y_matrix = transformer_y(numphases+1:r,1:c/2);
    end
    
    if sum(phases)==1
        ph = find(phases==1);
        from_from_y_matrix = zeros(3,3);
        from_from_y_matrix(ph,ph) = transformer_y(1,1);
        to_to_y_matrix = zeros(3,3);
        to_to_y_matrix(ph,ph) = transformer_y(3,3);
        from_to_y_matrix = zeros(3,3);
        from_to_y_matrix(ph,ph) = transformer_y(1,3);
    end
    
    
    
    
    %     y_matrix = [];
    %     for j=1:2:numphases*2
    %         y_matrix = [y_matrix line_y(:,j)+1i*line_y(:,j+1)];
    %     end
    %     temp = zeros(3,3);
    %     if numphases == 3
    %         temp = y_matrix;
    %     elseif sum(phases) == 1
    %         ph = find(transformers(i).phases==1);
    %         if length(y_matrix) ~= 1
    %             fprintf('Issues for line:%d\n',i);
    %         else
    %             temp(ph,ph)= y_matrix;
    %         end
    %     elseif sum(phases) == 2
    %         ph = find(lines(i).phases==1);
    %         temp(ph(1),ph(1)) = y_matrix(1,1);
    %         temp(ph(1),ph(2)) = y_matrix(1,2);
    %         temp(ph(2),ph(1)) = y_matrix(2,1);
    %         temp(ph(2),ph(2)) = y_matrix(2,2);
    %     end
    r_c = 1;
    for k1= from_bus_loc*3-2:from_bus_loc*3
        Y(k1,to_bus_loc*3-2:to_bus_loc*3) = Y(k1,to_bus_loc*3-2:to_bus_loc*3) + from_to_y_matrix(r_c,:);
        r_c = r_c+1;
    end
    r_c = 1;
    for k1= to_bus_loc*3-2:to_bus_loc*3
        Y(k1,from_bus_loc*3-2:from_bus_loc*3) = Y(k1,from_bus_loc*3-2:from_bus_loc*3)+ transpose(from_to_y_matrix(:,r_c));
        r_c = r_c+1;
    end
    r_c = 1;
    for k1= from_bus_loc*3-2:from_bus_loc*3
        Y(k1,from_bus_loc*3-2:from_bus_loc*3) = Y(k1,from_bus_loc*3-2:from_bus_loc*3) +  from_from_y_matrix(r_c,:);
        r_c = r_c+1;
    end
    r_c = 1;
    for k1= to_bus_loc*3-2:to_bus_loc*3
        Y(k1,to_bus_loc*3-2:to_bus_loc*3) = Y(k1,to_bus_loc*3-2:to_bus_loc*3) +  to_to_y_matrix(r_c,:);
        r_c = r_c+1;
    end
end


for i = 1:length(capacitors)
    name = split(capacitors(i).busName,'.');
    from_bus_loc = inner_map(name{1});
    line_y = capacitors(i).Yprim;
    nphases = capacitors(i).numPhases;
    phases = capacitors(i).phases;
    temp = zeros(3,3);
    if nphases == 1
        ph = find(phases==1);
        temp(ph,ph) = line_y(1,1)+1i*line_y(1,2);
    elseif nphases==3
        line_y = transpose(reshape(line_y,[nphases*2*2,nphases*2]));
        line_y = line_y(1:nphases,1:nphases*2);
        y_matrix = [];
        for j=1:2:nphases*2
            y_matrix = [y_matrix line_y(:,j)+1i*line_y(:,j+1)];
        end
        temp=y_matrix;
    end
    
    r_c = 1;
    for k1= from_bus_loc*3-2:from_bus_loc*3
        Y(k1,from_bus_loc*3-2:from_bus_loc*3) = Y(k1,from_bus_loc*3-2:from_bus_loc*3) +  temp(r_c,:);
        r_c = r_c+1;
    end
    
    
end
%%
if ADD_LOAD == 1
    loads = getLoadInfo(DSSCircObj);
    for i = 1:length(loads)
        name = split(loads(i).busName,'.');
        from_bus_loc = inner_map(name{1});
        if loads(i).isDelta
            if loads(i).numPhases == 1
                line_y = loads(i).Yprim;
                phase = find(loads(i).phases==1);
                temp = zeros(3,3);
                temp(phase(1),phase(1)) = line_y(1,1)+1i*line_y(1,2);
                temp(phase(2),phase(2)) = line_y(1,1)+1i*line_y(1,2);
                temp(phase(1),phase(2)) = line_y(1,3)+1i*line_y(1,4);
                temp(phase(2),phase(1)) = line_y(1,3)+1i*line_y(1,4);
            elseif loads(i).numPhases == 3
                line_y = loads(i).Yprim;
                temp = transpose(reshape(line_y,[6,3]));
                temp = [temp(:,1)+1i*temp(:,2) temp(:,3)+1i*temp(:,4) temp(:,5)+1i*temp(:,6)];
            end
            
        else
            line_y = loads(i).Yprim;
            if loads(i).numPhases == 3
                temp = zeros(3,3);
                y = line_y(1,1)+1i*line_y(1,2);
                temp(1,1) = y;
                temp(2,2) = y;
                temp(3,3) = y;
            elseif loads(i).numPhases == 1
                temp = zeros(3,3);
                phase = find(loads(i).phases==1);
                temp(phase,phase)= line_y(1,1)+1i*line_y(1,2);
                
            end
            
        end
        r_c = 1;
        for k1= from_bus_loc*3-2:from_bus_loc*3
                Y(k1,from_bus_loc*3-2:from_bus_loc*3) = Y(k1,from_bus_loc*3-2:from_bus_loc*3) +  temp(r_c,:);
                r_c = r_c+1;
        end
    end
end


%%
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

[nmap, rmap, fnew, tnew] = NodeRelabling(from_bus_nodes',to_bus_nodes',inner_map('sourcebus'))

node_ordering_opendss = DSSCircuit.YNodeOrder;
node_ordering_opendss_names = cell(size(node_ordering_opendss));
G = digraph(fnew,tnew);

for i = 1:length(node_ordering_opendss)
    name= split(node_ordering_opendss(i),'.');
    node_ordering_opendss_names{i} = name{1};
end
opendss_bus_names = DSSCircuit.AllBusNames;
testkeys = keys(inner_map);

%%
Buses = getBusInfo(DSSCircObj,bus_names);


full_voltage_magnitude = [];
for i = 1:length(Buses)
    v = (Buses(i).phaseVoltages.* exp(1i*Buses(i).voltagebaseAngle));%/(1000*Buses(i).kVBase);
    full_voltage_magnitude = [full_voltage_magnitude; v' ];
end

AllBusVoltage = DSSCircuit.AllbusVolts;
vcomplex = AllBusVoltage(1:2:end) + 1i * AllBusVoltage(2:2:end) ; % unit: Volt
vcomplex = transpose(vcomplex);


% s_v1 = diag(full_voltage_magnitude) * conj(Y*full_voltage_magnitude)/1000;
s_v2 = diag(vcomplex) * conj(Y_dss*vcomplex)/1000;
% s_v3 = diag(vcomplex) * conj(Y*vcomplex)/1000;
%%
% bus_a = 1:3:length(Y);
% bus_b = 2:3:length(Y);
% bus_c = 3:3:length(Y);
%
% % Y_a =zeros(size(Y));
% %
% % for i = 1:length(bus_a)
% %     Y_a(:,bus_a(i)) = Y(:,bus_a(i));
% %     Y_a(bus_a(i),:) = Y (bus_a(i),:);
% % end
%
% Y_a = Y(bus_a',:);
% Y_a = Y_a(:,bus_a);
%
% Y_b = Y(bus_b',:);
% Y_b = Y_b(:,bus_b);
%
% Y_c = Y(bus_c',:);
% Y_c = Y_c(:,bus_c);
% phase_a_voltage = full_voltage_magnitude(bus_a);
% phase_b_voltage = full_voltage_magnitude(bus_b);
% phase_c_voltage = full_voltage_magnitude(bus_c);

%%
% loc = find(abs(full_voltage_magnitude)>0);
% Y(:,setdiff(1:length(Y),loc))=[];
% Y(setdiff(1:length(Y),loc),:)=[];
% Y=Y/(4.16*4.16/100);
%%
Y_DD = Y_dss; % diagonal matrix

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




%%
number_of_eigenvectors = length(Y_DD)-45;
% vcomplex = full_voltage_magnitude;
U_k = U_Y_new (:,1:number_of_eigenvectors) ;
v_new = U_k*transpose(U_k)*vcomplex;
figure(3)
% subplot(311)
% loc = find(abs(full_voltage_magnitude(loc))>0);
subplot(211)
plot(abs(v_new),'linewidth',1.5)
hold on
plot(abs(vcomplex),'linewidth',1.5)
legend('Interpolated','Original')
subplot(212)
plot(angle(v_new),'linewidth',1.5)
hold on
plot(angle(vcomplex),'linewidth',1.5)
legend('Interpolated','Original')

ten_bus_matrix_map = [1 1 1 3 3 4 2 7 7 5];
bus_numbering = [1 2 6 7 8 9 3 4 5 10];
%% Covex Relaxation Code
clear tildeW
s = diag(vcomplex) * conj(Y_dss*vcomplex)/1000;
s_1 = diag(full_voltage_magnitude) * conj(Y*full_voltage_magnitude)/1000;
sm= s*1000;
vm = abs(vcomplex);
c = Y_dss*vcomplex;
cm = abs(c);
cvx_clear
cvx_solver SDPT3
cvx_precision low
cvx_begin 

    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(size(sm)) complex
    variable v_est(size(vm))
    variable a_est(size(cm))
    minimize (15*norm(sm - s_est, 2) + 5*norm(vm.*vm-(v_est),2) + 0*norm(cm.*cm-a_est,2) )

    subject to
    %         sm == M*s_est;
    tildeW == tildeW';
    s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
    v_est == diag(U_k*tildeW*(U_k)');
%     0 <= sum(real(s_est)) <=150;
%     0 <= sum(imag(s_est)) <=310;
%     a_est == diag(Y_DD*U_k*tildeW*U_k'*Y_DD');
%     U_k(1,:)*(tildeW*U_k(1,:)') == 1;
    U_k(1:3,:)*(tildeW*U_k(1:3,:)') == vcomplex(1:3)*vcomplex(1:3)';
cvx_end
W = U_k * tildeW *(U_k)';


%% Voltage Reconstruction for the 37 bus test case
all_voltage = zeros(size(vcomplex));
all_voltage(1:3,1) = vcomplex(1:3);
for node = 2:length(rmap)
    to = node;
    from = predecessors(G,to);
    to = (find(nmap==to));
    testind = cellfun(@(x)isequal(x,to),values(inner_map));
    msg_key = testkeys(testind);
    to= find(strcmp(opendss_bus_names,msg_key));
    
    from = (find(nmap==from)); 
    testind = cellfun(@(x)isequal(x,from),values(inner_map));
    msg_key = testkeys(testind);
    from= find(strcmp(opendss_bus_names,msg_key));
    source_voltage = all_voltage(from*3-2:from*3,1);
    all_voltage(to*3-2:to*3,1) = transpose(1/trace(source_voltage*source_voltage')*W(to*3-2:to*3,from*3-2:from*3)*source_voltage);
    
end

figure(6)
subplot(211)
plot(abs(all_voltage(4:end)),'linewidth',1.5);
hold on 
plot(abs(vcomplex(4:end)),'linewidth',1.5);
legend('Estimated','Actual');

subplot(212)
plot(angle(all_voltage(4:end)),'linewidth',1.5);
hold on 
plot(angle(vcomplex(4:end)),'linewidth',1.5);
legend('Estimated','Actual');


%%
% V1 = vcomplex(1:3);
% W21 = W(4:6,1:3);
% V2 = 1/trace(V1*V1')*W21*V1;
% V3 = 1/trace(V1*V1')*W(16:18,1:3)*V1;
% V4 = 1/trace(V3*V3')*W(19:21,16:18)*V3;
% 
% voltr = transpose(reshape(vcomplex,[3,10]));
% voltr = voltr(bus_numbering',:);
% 
% all_voltage = zeros(10,3);
% all_voltage(1,:) = transpose(vcomplex(1:3));
% for i = 2:10
%     row_pos = ten_bus_matrix_map(i);
%     bus_number = bus_numbering(i);
%     source_voltage = transpose(all_voltage(row_pos,:));
%     all_voltage(i,:) = transpose(1/trace(source_voltage*source_voltage')*W(bus_number*3-2:bus_number*3,bus_numbering(row_pos)*3-2:bus_numbering(row_pos)*3)*source_voltage);
% end

max(max(abs(all_voltage-voltr)))

figure(5)
subplot(321)
plot(abs(all_voltage(:,1)),'linewidth',1.5)
hold on 
plot(abs(voltr(:,1)),'linewidth',1.5)

subplot(322)
plot(angle(all_voltage(:,1)),'linewidth',1.5)
hold on 
plot(angle(voltr(:,1)),'linewidth',1.5)

subplot(323)
plot(abs(all_voltage(:,2)),'linewidth',1.5)
hold on 
plot(abs(voltr(:,2)),'linewidth',1.5)

subplot(324)
plot(angle(all_voltage(:,2)),'linewidth',1.5)
hold on 
plot(angle(voltr(:,2)),'linewidth',1.5)

subplot(325)
plot(abs(all_voltage(:,3)),'linewidth',1.5)
hold on 
plot(abs(voltr(:,3)),'linewidth',1.5)

subplot(326)
plot(angle(all_voltage(:,3)),'linewidth',1.5)
hold on 
plot(angle(voltr(:,3)),'linewidth',1.5)


 
