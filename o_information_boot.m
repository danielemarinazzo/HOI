function o=o_information_boot(X,indsample,indvar)
X=X(indsample,indvar);
% this function takes thw whole X as input, and additionally the indices
% convenient for bootstrap
%X size is N (samples) x M(variables)
[N, M]=size(X);

o= (M-2)*ent_g(X, 1);
for j=1:M
    X1=X;X1(:,j)=[];
    o=o+ent_g(X(:,j), 1)- ent_g(X1, 1);
end