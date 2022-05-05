% demend and covariance matrix
function [outBoot, outJack] = hoi_createBootsData(X,nboot,pathOut)
    
    [n,p] = size(X);
    
    covBootst = zeros(p,p,nboot);
    parfor i = 1:nboot+1
        if i==1
            idx = 1:n;
        else
            idx = randsample(n,n,1); 
        end
        Xb = X(idx,:);
        covBootst(:,:,i) = cov(Xb);
    end
    outBoot = [pathOut filesep 'covBootst.mat']; 
    save(outBoot, 'covBootst');
    
    covJack = zeros(p,p,n);
    parfor i = 1:n
        Xi = X;
        Xi(i,:) = [];
        covJack(:,:,i) = cov(Xi);
    end
    outJack = [pathOut filesep 'covJack.mat']; 
    save(outJack, 'covJack');

end