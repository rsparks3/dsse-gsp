
clc;
A = randi(10,[5,3])+1i*randi(10,[5,3]);
[m,n] = size(A);
Q =zeros(m,n);
R= zeros(n,n);

for j = 1:n
    v = A(:,j);
    for i = 1:j-1
        R (i,j) = Q(:,i)'*A(:,j);
        v = v - R(i,j)*Q(:,i);
    end
    R(j,j)= norm(v);
    Q(:,j) = v/R(j,j);
end

fprintf('Check whether Q is orthonormal \n');
fprintf('Check Q^H Q == I');
Q'*Q
fprintf('Check Norm for each column \n')
[norm(Q(:,1)) norm(Q(:,2)) norm(Q(:,3))]

fprintf('Check Graham Schmidt is Working by checking A = QR \n')
fprintf('A Matrix \n')
A
fprintf('QR Matrix \n')
Q*R


fprintf('Calculate Q_perp using the definition I-QQ^H \n');
inter_1 = Q*Q';
Q_perp_1 = eye(size(inter_1))-inter_1; 
fprintf('Check for the equations in the paper \n')
fprintf('Check Q_perp^H*Q == 0 \n ')
Q_perp_1'*Q
fprintf('This can be considered as 0 matrix \n')

fprintf('Check Q_perp^H*Q_perp == I \n ')
Q_perp_1'*Q_perp_1
fprintf('Expected to be Identity, but not found \n\n')


fprintf('Calculate Q_perp using the definition I-QQ^{dagger} \n');
inter_2 = Q*pinv(Q);
Q_perp_2 = eye(size(inter_2))-inter_2; 
fprintf('Check for the equations in the paper \n')
fprintf('Check Q_perp^H*Q == 0\n ')
Q_perp_2'*Q
fprintf('This can be considered as 0 matrix \n')

fprintf('Check Q_perp^H*Q_perp == I \n ')
Q_perp_2'*Q_perp_2
fprintf('Expected to be Identity, but not found \n \n')

fprintf('Calculate Q_perp using the definition that orthogonal complement is the basis of the null space of Q^H \n');
Q_perp_3 = null(Q');
fprintf('Check for the equations in the paper \n')
fprintf('Check Q_perp^H*Q == 0\n ')
Q_perp_3'*Q
fprintf('This can be considered as 0 matrix \n')

fprintf('Check Q_perp^H*Q_perp == I \n ')
Q_perp_3'*Q_perp_3
fprintf('Expected to be Identity, and is found \n \n')

