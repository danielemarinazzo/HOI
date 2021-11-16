function o=o_information_lagged_boot(Y,X,m,indstart,chunklength,indvar)
% evaluates the o_information flow
% Y Nx1 target vector .  X NxM drivers
% m order of the model


%R=randi(N-b+1,Nb,1);
if chunklength==0
    indsample=1:length(Y);
else
    nchunks=floor(length(Y)/chunklength);
    indstart=indstart(1:nchunks);
    for istart=1:nchunks
        indsample(1+(istart-1)*chunklength:istart*chunklength)=...
            indstart(istart):indstart(istart)+chunklength-1;
    end
end

Y=Y(indsample,:);
X=X(indsample,indvar);
[N, M]=size(X);
%Y=copnorm(Y);
%X=copnorm(X);
n=N-m;
X0=zeros(n,m,M);
Y0=zeros(n,m);
y=Y(m+1:end);
for i=1:n
    for j=1:m
        Y0(i,j)=Y(m-j+i);
        for k=1:M
            X0(i,j,k)=X(m-j+i,k);
        end
    end
end

o=-(M-1)*gccmi_ccc_nocopnorm(y,reshape(X0,n,m*M),Y0);
for k=1:M
    X=X0;X(:,:,k)=[];
    o=o+gccmi_ccc_nocopnorm(y,reshape(X,n,m*(M-1)),Y0);
end