function [Cluster_elem,M,EW_history,W,z]=CGFFCM_ablation_v2(X,M,k,p_init,p_max,p_step,t_max,beta_memory,N,fuzzy_degree,d,beta_z,landa,v,G,cfg)
% Ablation-enabled CGFFCM based on your CGFFCM.m
% cfg fields (logical):
%   cfg.learnW, cfg.learnZ, cfg.useMemory, cfg.adaptiveP, cfg.useV, cfg.useKernel
% defaults: all true

if nargin<16 || isempty(cfg), cfg=struct(); end
if ~isfield(cfg,'learnW'),    cfg.learnW=true;    end
if ~isfield(cfg,'learnZ'),    cfg.learnZ=true;    end
if ~isfield(cfg,'useMemory'), cfg.useMemory=true; end
if ~isfield(cfg,'adaptiveP'), cfg.adaptiveP=true; end
if ~isfield(cfg,'useV'),      cfg.useV=true;      end
if ~isfield(cfg,'useKernel'), cfg.useKernel=true; end

% --- handle ablation switches ---
if ~cfg.useV
    v = ones(1,d);
end

if ~cfg.adaptiveP
    % fix p at p_max (or you can choose p_init)
    p_init = p_max;
    p_step = 0;
end

% ---------------- checks (same as your code) ----------------
if p_init<0 || p_init>=1, error('p_init must take a value in [0,1)'); end
if p_max<0 || p_max>=1,  error('p_max must take a value in [0,1)'); end
if p_max<p_init,         error('p_max must be greater or equal to p_init'); end
if p_step<0,             error('p_step must be a non-negative number'); end
if beta_memory<0 || beta_memory>1, error('beta must take a value in [0,1]'); end
if beta_z==0, error('beta_z must be a non-zero number'); end

if p_init==p_max
    p_flag=0; p_step=0;
elseif p_step==0
    p_flag=0; p_max=p_init;
else
    p_flag=1;
end

% weights init
W = ones(1,k)/k;
z(:,1:3)=ones(k,3)/3;
z(:,4:6)=ones(k,3)/3;
z(:,7:8)=ones(k,2)/2;

p = p_init;
p_prev = p-1e-8;
empty=0; Iter=1;
E_w_old=inf;
EW_history=[];
Cluster_elem_history=[]; W_history=[]; z_history=[];
nn=1;

while true
    % -------- Step 1: update memberships --------
    dNK = zeros(N,k);
    for j=1:k
        diff = X - repmat(M(j,:),N,1);
        if cfg.useKernel
            dist_base = 1 - exp((-1.*repmat(landa,N,1)).*(diff.^2));
        else
            dist_base = (repmat(landa,N,1)).*(diff.^2); % non-kernel variant
        end
        dist = dist_base .* repmat(v,N,1);

        WBETA = transpose(z(j,:).^beta_z);
        WBETA(WBETA==inf)=0; WBETA(isnan(WBETA))=0;

        if cfg.learnW
            dNK(:,j) = (W(1,j).^p) * (dist * WBETA);
        else
            % W fixed uniform
            dNK(:,j) = ( (1/k)^p ) * (dist * WBETA);
        end
    end

    tmp1 = zeros(N,k);
    for j=1:k
        tmp2 = (dNK./repmat(dNK(:,j),1,k)).^(1/(fuzzy_degree-1));
        tmp2(tmp2==inf)=0; tmp2(isnan(tmp2))=0;
        tmp1 = tmp1 + tmp2;
    end
    Cluster_elem = transpose(1./tmp1);
    Cluster_elem(isnan(Cluster_elem))=1;
    Cluster_elem(Cluster_elem==inf)=1;

    if nnz(dNK==0)>0
        for i=1:N
            if nnz(dNK(i,:)==0)>0
                Cluster_elem(find(dNK(i,:)==0),i) = 1/nnz(dNK(i,:)==0);
                Cluster_elem(find(dNK(i,:)~=0),i) = 0;
            end
        end
    end

    % -------- objective (use the consistent version) --------
    if cfg.useKernel
        E_w = object_fun1(N,d,k,Cluster_elem,landa,M,fuzzy_degree,W,z,beta_z,p,X,v);
    else
        E_w = object_fun_linear(N,d,k,Cluster_elem,landa,M,fuzzy_degree,W,z,beta_z,p,X,v);
    end
    EW_history(nn)=E_w; nn=nn+1;

    % -------- empty/singleton check (same logic) --------
    for i=1:k
        I=find(Cluster_elem(i,:)<=0.05);
        if length(I)==N-1 || length(I)==N
            E_w=NaN; empty=empty+1;
            if empty>1, p=p-p_step; else, p=p_prev; end
            p_flag=0;

            if p<p_init || p_step==0
                M=NaN(k,size(X,2));
                return;
            end

            a=(k*empty)-(k-1); b=k*empty;
            Cluster_elem=Cluster_elem_history(a:b,:);
            W=W_history(empty,:);
            aa=(k*empty)-(k-1); bb=k*empty;
            z=z_history(aa:bb,:);
            break;
        end
    end

    % -------- convergence --------
    if ~isnan(E_w) && ~isnan(E_w_old) && (abs(1-E_w/E_w_old)<1e-6 || Iter>=t_max)
        break;
    end
    E_w_old=E_w;

    % -------- Step 2: update centers --------
    mf = Cluster_elem.^fuzzy_degree;
    for j=1:k
        diff = X - repmat(M(j,:),N,1);
        if cfg.useKernel
            E = exp((-1.*repmat(landa,N,1)).*(diff.^2));
            numer = (mf(j,:) * (X .* E));
            denom = (mf(j,:) * E);
            M(j,:) = numer ./ denom;
        else
            % fallback to standard weighted mean
            numer = mf(j,:) * X;
            denom = sum(mf(j,:));
            M(j,:) = numer ./ denom;
        end
    end

    % -------- Step 3: increase p --------
    if p_flag==1
        Cluster_elem_history=[Cluster_elem;Cluster_elem_history];
        W_history=[W;W_history];
        z_history=[z;z_history];

        p_prev=p;
        p=p+p_step;
        if p>=p_max
            p=p_max; p_flag=0;
        end
    end

    W_old=W; z_old=z;

    % -------- Step 4: update feature weights z --------
    if cfg.learnZ
        dWkm = zeros(k,d);
        for j=1:k
            diff = X - repmat(M(j,:),N,1);
            if cfg.useKernel
                dist_base = 1 - exp((-1.*repmat(landa,N,1)).*(diff.^2));
            else
                dist_base = (repmat(landa,N,1)).*(diff.^2);
            end
            dist_w = dist_base .* repmat(v,N,1);
            dWkm(j,:) = (Cluster_elem(j,:).^fuzzy_degree) * dist_w;
        end

        z_new = z;

        % group1: 1:3
        tmp1=zeros(k,3);
        for jj=1:3
            tmp2=(dWkm(:,1:3)./repmat(dWkm(:,jj),1,3)).^(1/(beta_z-1));
            tmp2(tmp2==inf)=0; tmp2(isnan(tmp2))=0;
            tmp1=tmp1+tmp2;
        end
        z_new(:,1:3)=1./tmp1;

        % group2: 4:6
        tmp1=zeros(k,3);
        for jj=1:3
            tmp2=(dWkm(:,4:6)./repmat(dWkm(:,3+jj),1,3)).^(1/(beta_z-1));
            tmp2(tmp2==inf)=0; tmp2(isnan(tmp2))=0;
            tmp1=tmp1+tmp2;
        end
        z_new(:,4:6)=1./tmp1;

        % group3: 7:8
        tmp1=zeros(k,2);
        for jj=1:2
            tmp2=(dWkm(:,7:8)./repmat(dWkm(:,6+jj),1,2)).^(1/(beta_z-1));
            tmp2(tmp2==inf)=0; tmp2(isnan(tmp2))=0;
            tmp1=tmp1+tmp2;
        end
        z_new(:,7:8)=1./tmp1;

        z_new(isnan(z_new))=1; z_new(z_new==inf)=1;
        z = z_new;
    end

    % -------- Step 5: update cluster weights W --------
    if cfg.learnW
        Dw=zeros(1,k);
        for j=1:k
            diff = X - repmat(M(j,:),N,1);
            if cfg.useKernel
                dist_base = 1 - exp((-1.*repmat(landa,N,1)).*(diff.^2));
            else
                dist_base = (repmat(landa,N,1)).*(diff.^2);
            end
            dist = dist_base .* repmat(v,N,1);

            WBETA = transpose(z(j,:).^beta_z);
            WBETA(WBETA==inf)=0; WBETA(isnan(WBETA))=0;

            val = dist * WBETA; % N x 1
            Dw(1,j) = sum( val .* (Cluster_elem(j,:)'.^fuzzy_degree) );
        end

        tmp = sum((repmat(Dw,k,1)./transpose(repmat(Dw,k,1))).^(1/(p-1)));
        tmp(tmp==inf)=0; tmp(isnan(tmp))=0;
        W_new = 1./tmp;
        W_new(isnan(W_new))=1; W_new(W_new==inf)=1;

        if nnz(Dw==0)>0
            W_new(find(Dw==0)) = 1/nnz(Dw==0);
            W_new(find(Dw~=0)) = 0;
        end
        W = W_new;
    else
        W = ones(1,k)/k;
    end
    % -------- Step 6: memory effect --------
    if cfg.useMemory
        if cfg.learnW
            % Eq. memory: W^{t+1} = (1-β)W^{t} + βW_{new}
            W = (1-beta_memory)*W_old + beta_memory*W;
        else
            W = ones(1,k)/k;
        end
        if cfg.learnZ
            z = (1-beta_memory)*z_old + beta_memory*z;
        else
            z(:,1:3)=ones(k,3)/3;
            z(:,4:6)=ones(k,3)/3;
            z(:,7:8)=ones(k,2)/2;
        end
    end
    Iter=Iter+1;
end
end
