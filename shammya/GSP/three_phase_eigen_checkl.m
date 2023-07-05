clc;
clear ;
close all;
[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText.command = 'clear';
% DSSText.command = 'Compile C:\models\123Bus\IEEE123Master.DSS';
% DSSText.command = 'Compile C:\models\13Busm\IEEE13Nodeckt.dss';
DSSText.command = 'Compile C:\models\test3bus\3bus.dss';
DSSText.Command = 'solve';


lines = getLineInfo(DSSCircObj);
transformers = getTransformerInfo(DSSCircObj);
bus_names = cell(length(lines)*2 + length(transformers)*2,1);
counter = 1;
for i = 1:length(lines)
    phases = lines(i).phases;
    numphases = lines(i).numPhases;
    line_y = lines(i).Yprim;
    if sum(phases) == 3
        line_y = transpose(reshape(line_y,[numphases*2*2,numphases*2]));
        line_y = line_y(1:numphases,1:numphases*2);
        y_matrix = [];
        for j=1:2:numphases*2
            y_matrix = [y_matrix line_y(:,j)+1i*line_y(:,j+1)];
        end
        Y_DD = y_matrix % diagonal matrix
        
        [U_Y, Lambda_Y] = eig(Y_DD);%./mpcint.baseMVA);
        [val_Y, ind_Y]=sort(diag(Lambda_Y),'ascend');
        U_Y_arr = U_Y(:,ind_Y');
        
        % change into complex orthogonal matrices:
        r_vect = zeros(length(Y_DD),1);
        for ii=1:length(Y_DD)
            r_vect(ii) = 1./sqrt(transpose(U_Y_arr(:,ii))*U_Y_arr(:,ii));
        end
        
        U_Y_new = U_Y_arr*diag(r_vect);
        disp(U_Y)
        
%         break;
    end
end
% x = [ 1 1 1 ; 1 1*exp(1i*2*pi/3) 1*exp(1i*4*pi/3) ; 1 1*exp(1i*4*pi/3) 1*exp(1i*8*pi/3)];
% Y_DD = transpose(x); % diagonal matrix
%         
% [U_Y, Lambda_Y] = eig(Y_DD);%./mpcint.baseMVA);
% [val_Y, ind_Y]=sort(diag(Lambda_Y),'ascend');
% U_Y_arr = U_Y(:,ind_Y');
% 
% % change into complex orthogonal matrices:
% r_vect = zeros(length(Y_DD),1);
% for ii=1:length(Y_DD)
%     r_vect(ii) = 1./sqrt(transpose(U_Y_arr(:,ii))*U_Y_arr(:,ii));
% end
% 
% U_Y_new = U_Y_arr*diag(r_vect);
% % U_Y_new = 
% % disp(U_Y)



