clc
clear all
warning off
%% Load the data.
% Please uncomment the dataset that you want to test.
load("3Sources.mat")
disp('Dateset：3Sources')

viewnum=length(X);

opt.maxiter=50;

data=X;
k=length(unique(Y));% k is the number of cluster

% Parameter setting
d=[ 3*k ]; 
l=[ 3*k];
num_p=[3]; % num_p is the number of LES
lambda=[0.001 ]; %this is default

for v=1:viewnum
    data{v}=data{v};
    data{v} = zscore(data{v})';
    %data{v}=data{v}';
end


for pp=1:length(num_p)
for dd=1:length(d)
for ll=1:length(l)
for lm=1:length(lambda)

tic;
%Run MSC-DOLES
[L,output]=DOLES(data,num_p(pp),d(dd),l(ll),k,lambda(lm),opt);
toc
%res save the scores of ACC, NMI, PUR and FScore.
res = Clustering8Measure(Y,L);

fprintf(' p:%d\t d:%d\t l:%d\t lambda:%12.5f\t ACC:%12.5f\t nmi:%12.5f\t Purity:%12.5f\t Fscore:%12.5f \t\n',[ num_p(pp) d(dd) l(ll) lambda(lm) res(1) res(2) res(3) res(4) ]);

end

end
end
end

