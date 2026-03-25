function [E_w] = object_fun_linear(N,d,k,Cluster_elem,landa,M,fuzzy_degree,W,z,beta_z,p,X,v)
% Linear (non-kernel) objective consistent with the non-kernel distance:
%   dist = landa .* (X-M).^2
% and the same weighting scheme used in CGFFCM (v and z^beta_z).

E_w = 0;

for i=1:k
    WBETA = transpose(z(i,:).^beta_z);
    WBETA(WBETA==inf)=0;
    WBETA(isnan(WBETA))=0;

    for j=1:N
        diff = (X(j,:)-M(i,:)).^2;
        dist = (landa .* diff) .* v;     % 1 x d
        E_w = E_w + (Cluster_elem(i,j).^fuzzy_degree) * (W(1,i).^p) * (dist * WBETA);
    end
end

end
