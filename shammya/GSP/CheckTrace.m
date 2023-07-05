clc;
clear all;
close all;
A = randi(10,[4,2]);
nb=4;
syms w11
syms w12
syms w22
W= [1 1+1i; 1-1i 2];
A = A+1i*A;
Y = randi(15,[nb,nb]);
Y = Y+1i*Y;
C=A';
% C = randi(8,[2,4]);
% C = C+1i*C;
test_1 = diag(A*W*C);
test_3 = diag(Y*A*W*C);

Uk_UkH =  zeros(2,2*2);
YUk_UkH =  zeros(2,2*2);
for k = 1:nb
   Uk_UkH(k,:) =  transpose(vec(C(:,k)*A(k,:)));
end

A = Y*A;
for k = 1:nb
   YUk_UkH(k,:) =  transpose(vec(C(:,k)*A(k,:)));
end


for k = 1:nb
    test_2= Uk_UkH(k,:)*vec(transpose(W));
    
end

for k = 1:nb
    test_4= YUk_UkH(k,:)*vec(transpose(W))
    
end

