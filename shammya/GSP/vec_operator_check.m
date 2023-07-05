clc;
clear all;
rng(1)
syms w11;
syms w12;
syms w22;syms w13; syms w33;syms w23;syms w21

A = randi(10,3,2)+1i*randi(10,3,2);
W =[w11 w12 ;w21 w22];
C = A';
% C = randi(10,2,3)+1i*randi(10,2,3);
diag(A*W*C)

transpose(vec(C(:,1).*A(1,:)))*vec(transpose(W))
transpose(vec(C(:,2).*A(2,:)))*vec(transpose(W))
transpose(vec(C(:,3).*A(3,:)))*vec(transpose(W))

% A = randi(10,4,3)+1i*randi(10,4,3);
% W =[w11 w12 w13;conj(w12) w22 w23;conj(w13) conj(w23) w33];
% C=A';
% diag(A*W*C)
% transpose(vec(C(:,1).*A(1,:)))*vec(transpose(W))
x = [transpose(vec(C(:,1).*A(1,:)));
    transpose(vec(C(:,2).*A(2,:)));
    transpose(vec(C(:,3).*A(3,:)))]



% This is written in the paper
a = A(:,1);

