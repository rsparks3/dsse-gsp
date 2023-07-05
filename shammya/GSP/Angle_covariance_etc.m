%% act on the parsed data for ACTIVSG
clear variables;
% close all;
clc;
%% mask of participating buses
fname   = 'case141';
fprintf('Fetching sorted Laplacian...')
tstart = tic;
data = getmpc_laplacian(fname);
fprintf('Complete (%0.3f sec)\n', toc(tstart))
nb = size(data.L,1);
ngbus = sum(data.genmask);
%%
bus_names = cell(nb, 1);
figure(1);
subplot(133)
spy(data.L,10)
hold on;

plot([1, nb], [ngbus, ngbus], '--k');
plot([ngbus, ngbus], [1, nb], '--k');
title('Sorted Laplacian')
%% use k-means/spectral clustering on Laplacian directly:
ncluster = 8;
Q = full(data.B);
[QV, QD]=eigs(Q, ncluster,'sm');
[Qidx, QC] = kmeans(QV(:,1:ncluster),ncluster);
%% plot the adjacency:
Adjacency_clusters=zeros(length(Qidx));
for i=1:length(Qidx)
    edges = find(Qidx==Qidx(i)); %from the same group
    Adjacency_clusters(i,edges')=1;
end
ngbus=0;
figure;
spy(Adjacency_clusters);
hold on;
plot([1, nb], [ngbus, ngbus], '--k');
plot([ngbus, ngbus], [1, nb], '--k');
%% angle covariance
load('ACTIVSg2000Parsed.mat');
%% unwrap angle
va = unwrap(va*pi/180,[],2);

angcov = va*va.'/length(t);
%figure;
%imagesc(angcov);

%% clustering of laplacian and covariance

[LV,LD] = eigs(data.L, ncluster+1 ,'sm');
[Lidx, LC] = kmeans(LV,ncluster);


[AngV, AngD]   = eigs(angcov, ncluster+1, 'lm');
[Angidx, AngC] = kmeans(AngV, ncluster);
%% plot the adjacency:
Adjacency_clusters=zeros(length(Angidx));
for i=1:length(Angidx)
    edges = find(Angidx==Angidx(i)); %from the same group
    Adjacency_clusters(i,edges')=1;
end

figure;
spy(Adjacency_clusters);
hold on;
plot([1, nb], [ngbus, ngbus], '--k');
plot([ngbus, ngbus], [1, nb], '--k');

%% 
phasor = vm.*exp(1i.*va);

covariance_matrix = phasor*transpose(conj(phasor))/length(t);




%% clustering 'cosine' of the angle...
cos_va = sin(va);

cos_ang_cov = cos_va*cos_va.'/length(t);
[CosAngV,CosAngD] = eigs(cos_ang_cov, ncluster+1,'lm');
[CosAngidx, CosAngC]= kmeans(CosAngV, ncluster);

%% plot the adjacency:
Adjacency_clusters_cos=zeros(length(CosAngidx));
for i=1:length(CosAngidx)
    edges = find(CosAngidx==CosAngidx(i)); %from the same group
    Adjacency_clusters_cos(i,edges')=1;
end

figure;
spy(Adjacency_clusters_cos);
hold on;
plot([1, nb], [ngbus, ngbus], '--k');
plot([ngbus, ngbus], [1, nb], '--k');