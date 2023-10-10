function [L,output]=DOLES(X,num_p,d,l,k,lambda,options)
% 1/2*alpha^2*|X(v)-Hp*Zp|+beta*|Zp-PAS|+lambda*|S-G*F|
% X: m*n
% num_p: number of spaces
% d: number of anchor
% l: number of bases
% k: cluster number
% lambda: parameters
% options: maxiter,


if isfield(options,'maxiter') maxiter=options.maxiter;else maxiter=100;end

%% initialize 
num_view = length(X);
num_sample = size(X{1},2);
alpha = ones(num_p,1)/num_p;
beta = ones(num_p,1)*sqrt(1/num_p);
Z = cell(num_p,1);
H = cell(num_p,num_view);
m=k;
cnt = zeros(num_p,num_view);
cnt1 = ones(num_p,1);
P=cell(num_p,1);
for i=1:num_p
    Z{i} = eye(i*m,num_sample);
    P{i} = eye(i*m,d);
end

A=eye(d,l);
S=zeros(l,num_sample);
S(:,1:l) = eye(l);
G = eye(l,k);
F = eye(k,num_sample); 

flag = 1;
iter = 0;
%%
while flag
    iter = iter + 1;
%% Update H
    
    H = update_H(X,Z,H,num_p);
    
%% Update Z
    
    Z = update_Z(X,Z,H,P,A,S,alpha,beta,num_p);
    

%% Update P_i
    
    AS = A*S; 
    parfor iv=1:num_p
        C = Z{iv}*AS';      
        [U,~,V] = svd(C,'econ');
        P{iv} = U*V';
    end
   
%% Update A
    
    sumAlpha = 0;
    part1 = 0;
    for ia = 1:num_p
        al2 = alpha(ia)^2;
        sumAlpha = sumAlpha + al2;
        part1 = part1 + al2 * P{ia}' * Z{ia} * S';
    end
    [Unew,~,Vnew] = svd(part1,'econ');
    A = Unew*Vnew';
    
%% Update S
    
    HS = 2*sumAlpha*eye(l)+2*lambda*eye(l);
    HS = (HS+HS')/2;
    options = optimset( 'Algorithm','interior-point-convex','Display','off'); % interior-point-convex
    parfor ji=1:num_sample
        ff=0;
        e = F(:,ji)'*G';
        for j=1:num_p
            C = P{j} * A;
            ff = ff - 2*Z{j}(:,ji)'*C - 2*lambda*e;
        end
        S(:,ji) = quadprog(HS,ff',[],[],ones(1,l),1,zeros(l,1),ones(l,1),[],options);
    end
   
%% Update alpha beta
    
    for i=1:num_p
        cnt(i)=0.5*norm(Z{i}-P{i}*A*S,'fro');
        for j=1:num_view
            cnt(i,j) = sum(sum((X{j}-H{i,j}*Z{i}).^2));
        end
    end

    alpha = update_alpha(cnt,num_p,num_view);
    cnt_tmp1= cal_obj(cnt,cnt1,alpha,beta,num_p,num_view);
    beta = update_beta(cnt1,num_p);
    
%% Update G
   
    J = S*F';      
    [Ug,~,Vg] = svd(J,'econ');
    G = Ug*Vg';
    
    %% Update F
   
    F=zeros(k,num_sample);
    for iff=1:num_sample
        Dis=zeros(k,1);
        for jf=1:k
            Dis(jf)=(norm(S(:,iff)-G(:,jf)))^2;
        end
        [~,r]=min(Dis);
        F(r(1),iff)=1;
    end
    

%% compute obj
    
    cnt_tmp2= cal_obj(cnt,cnt1,alpha,beta,num_p,num_view);
   cnt_tmp1-cnt_tmp2;
    obj(iter) = cal_obj(cnt,cnt1,alpha,beta,num_p,num_view);
    obj(iter)=obj(iter)+lambda*norm(S - G * F,'fro')^2;
    if (iter>2) && (abs((obj(iter)-obj(iter-1))/(obj(iter)))<1e-5 || iter>maxiter)
        flag =0;
    end
    
end

%% classifier

[~,L]=max(F);
output.S=S;
output.F=F;
output.loss=obj;

end

%%


%% 
function [H] = update_H(fea,Z,H,num_p)

num_view = length(fea);
T = cell(num_p,num_view);
for i=1:num_p
    for j=1:num_view
        T{i,j} = fea{j};
        H{i,j} = T{i,j}*Z{i}';
    end
end
end


%%
function [alpha] = update_alpha(cnt,num_p,num_view)

tmp = zeros(num_p,1);
for i=1:num_p
    for j=1:num_view
        tmp(i) = tmp(i)+cnt(i,j);
    end
end
tmp = ones(num_p,1)./sqrt(tmp);
total = 0;
for i=1:num_p
    total = total+tmp(i);
end
alpha = tmp/total;
end

%%
function [beta] = update_beta(cnt1,num_p)

cnt = 0;
for i=1:num_p
    cnt = cnt+cnt1(i)^2;
end
beta = cnt1/sqrt(cnt);
end
%%
function [Z] = update_Z(fea,Z,H,P,A,S,alpha,beta,num_p)

num_view = length(fea);
U = cell(num_p,1);
V = cell(num_p,1);
T = cell(num_p,num_view);
M = cell(num_p,1);
CNT = cell(num_p,1);
for i=1:num_p
    CNT{i} = zeros(size(Z{1},2),size(Z{i},1));
    for j=1:num_view
        T{i,j} = fea{j};
        CNT{i} = CNT{i}+alpha(i)^2*T{i,j}'*H{i,j};
    end
    CNT{i} = CNT{i}+beta(i)*S'*A'*P{i}';
    [U{i},~,V{i}] = svd(CNT{i},'econ');
    M{i} = U{i}*V{i}';
    Z{i} = M{i}';
end
end

%%
function [obj] = cal_obj(cnt,cnt1,alpha,beta,num_p,num_view)

obj=0;
for i=1:num_p
    tmp=0;
    for j=1:num_view
        tmp = tmp+cnt(i,j);
    end
    obj = obj+1/2*alpha(i)^2*tmp-beta(i)*cnt1(i);
end
end


