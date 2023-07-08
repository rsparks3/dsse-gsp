function [voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G)

nb = G.numnodes();
if length(tildeW) ~= nb
    W = U_k * tildeW *(U_k)';
else
    W= tildeW;
end
voltage_angle= zeros(nb,1);
for i = 2:nb
    path = shortestpath(G,1,i);
    angle_sum = 0;
    for j = 1:length(path)-1
        angle_sum = angle_sum - angle(W(path(j),path(j+1)));
    end
    voltage_angle(i,1) = angle_sum;
    
end

voltage_mag = sqrt(abs(diag(W)));
end

