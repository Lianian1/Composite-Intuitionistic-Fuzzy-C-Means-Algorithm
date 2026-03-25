function [E_w] = object_fun1(N,d,k,Cluster_elem,landa,M,fuzzy_degree,W,z,beta_z,p,X,v)

%Initialize the objective.
E_w=0;

%Calculate the objective.
for i=1:k
    WBETA = transpose(z(i,:).^beta_z);
    WBETA(WBETA==inf)=0;
    WBETA(isnan(WBETA))=0;
    for j=1:N
        E_w=E_w + Cluster_elem(i,j).^fuzzy_degree * W(1,i).^p * sum(sum((1-exp((-1.*landa).*((X(j,:)-M(i,:)).^2))).*v.*WBETA));
    end
end

end