% demend and covariance matrix
function [outBoot, outJack] = hoi_createBootsData_lagged(X,chunklength,modelorder,nboot,pathOut)
    
    [n,p] = size(X);
    
    covBootst = zeros(p*(modelorder+1),p*(modelorder+1),nboot+1);
        
    parfor i = 1:nboot+1
        if i==1
            idx = 1:n;
        else
            idx = zeros(n,1);
            indstart = randsample(n-chunklength+1,n,1);
            % if the statistics toolbox is not installed use
            % tmp = 1:n-chunklength+1;
            % idx = tmp(randi(n, k))   %select SIZE elements WITH REPETITION
            nchunks = floor(n/chunklength);
            indstart = indstart(1:nchunks);
            for istart = 1:nchunks
                idx(1+(istart-1)*chunklength:istart*chunklength) = ...
                    indstart(istart):indstart(istart)+chunklength-1;
            end
        end
        Xt = X(idx,:);
        Xb = cell(modelorder+1,1);
        for m = 1:modelorder+1
            Xb{m} = Xt(m:end-modelorder-1+m,:);
        end
        Xb = cat(2,Xb{:});
        covBootst(:,:,i) = cov(Xb);
    end
    outBoot = [pathOut filesep 'covBootst.mat']; 
    save(outBoot, 'covBootst');
    
    covJack = zeros(p*(modelorder+1),p*(modelorder+1),n);
    parfor i = 1:n
        Xt = X;
        Xt(i,:) = [];
        Xb = cell(modelorder+1,1);
        for m = 1:modelorder+1
            Xb{m} = Xt(m:end-modelorder-1+m,:);
        end
        Xb = cat(2,Xb{:});        
        covJack(:,:,i) = cov(Xb);
    end
    outJack = [pathOut filesep 'covJack.mat']; 
    save(outJack, 'covJack');
end
