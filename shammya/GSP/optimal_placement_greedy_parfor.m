function [P_S,M_s] = optimal_placement_greedy_parfor(U_K,M,M_s)
%% greedy implementation of optimal placement:
% M_s is a set of already present indices/locations-just need to add onto it
% each time-column vector
M_tilde = M-length(M_s); % total number to be found
[N,~]=size(U_K);
i=1;
j_tilde = 1:N;
j_tilde(M_s)=[];
while i <=M_tilde
  %sigma_min=0;sigma_ind=1;
  sigma_min=zeros(length(j_tilde),1);
  parfor j=1: length(j_tilde)
      A = U_K([M_s;j_tilde(j)],:);
      [~,sigma_min(j),~,flag] = svds(A,1,'smallest');
      if(flag)
        continue;
      end
      % sigma_temp = svds(A,1,‘largest’); %min the largest
%       if(sigma_temp>sigma_min)
%       %if(sigma_temp<sigma_min)
%         sigma_min = sigma_temp;
%         sigma_ind=j;
%       end    
  end
  [~,min_ind]=max(sigma_min);
  M_s = [M_s;j_tilde(min_ind)];
  j_tilde(min_ind)=[];
  i=i+1;
end
P_S = eye(N);
P_S=P_S(:,M_s');
end