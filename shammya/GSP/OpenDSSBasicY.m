clc;
clear ;
close all;
[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText.command = 'clear';
DSSText.command = 'Compile C:\models\123Busm\IEEE123Master.DSS';
% DSSText.command = 'Compile C:\models\37Bus\ieee37.DSS';
% DSSText.command = 'Compile C:\models\37Bus\ieee37.DSS';
% DSSText.command = 'Compile C:\models\34Bus\ieee34Mod2.dss';
% DSSText.command = 'Compile C:\models\13Busm\IEEE13Nodeckt.dss';
% DSSText.command = 'Compile C:\models\test3bus\3bus.dss';
% DSSText.Command = 'solve';
ADD_LOAD = 0;

%%

DSSText.Command = 'vsource.source.enabled=no';
DSSText.command = 'batchedit load..* enabled=no';
DSSText.Command = 'solve';

base_y = DSSCircuit.SystemY;
Ymatrix=reshape(base_y,[length(DSSCircuit.AllNodeNames)*2,length(DSSCircuit.AllNodeNames)]); % "2" is because the real/imaginary parts are separate
Ymatrix=Ymatrix';
Y_dss = Ymatrix(:,1:2:end) + sqrt(-1)*Ymatrix(:,2:2:end);


DSSText.Command = 'vsource.source.enabled=no';
DSSText.command = 'batchedit load..* enabled=yes';
DSSText.Command = 'solve';

base_y = DSSCircuit.SystemY;
Ymatrix=reshape(base_y,[length(DSSCircuit.AllNodeNames)*2,length(DSSCircuit.AllNodeNames)]); % "2" is because the real/imaginary parts are separate
Ymatrix=Ymatrix';
Y_dss_load = Ymatrix(:,1:2:end) + sqrt(-1)*Ymatrix(:,2:2:end);

DSSText.Command = 'vsource.source.enabled=yes';
DSSText.command = 'batchedit load..* enabled=yes';
DSSText.Command = 'solve';

AllBusVoltage = DSSCircuit.AllbusVolts;
vcomplex = AllBusVoltage(1:2:end) + 1i * AllBusVoltage(2:2:end) ; % unit: Volt
vcomplex = transpose(vcomplex);
s = diag(vcomplex) * conj(Y_dss*vcomplex)/1000;


Y_DD = Y_dss_load; % diagonal matrix

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


number_of_eigenvectors = 150;
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


% y_node_order = DSSCircuit.YNodeOrder;
% ordered_bus_names = cell(numel(y_node_order),1);
% 
% for i = 1:numel(y_node_order)
%     bus_name = y_node_order{i};
%     bus = split(bus_name,'.');
% %     if strcmp(bus{2},'1')
% %         bus_a = [bus_a i];
% %     elseif strcmp(bus{2},'2')
% %         bus_b = [bus_b i];
% %     elseif strcmp(bus{2},'3')
% %         bus_c = [bus_c i];
% %     end
%     ordered_bus_names{i} = bus{1};
% end
% ordered_bus_names = unique(ordered_bus_names,'stable');
% buses = getBusInfo(DSSCircObj,ordered_bus_names);
% vcomplex = [buses.phaseVoltages].*exp(1j*[buses.voltagebaseAngle]);
