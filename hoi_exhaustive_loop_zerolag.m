function [Otot, O_val_size_tot] = hoi_exhaustive_loop_zerolag(ts, maxsize, n_best, biascorrection, pathTmp, groups)

% ts= input (observations x variables), time series or static/behavioral data
% maxsize = max number of variables in the multiplet
% n_best = number of most informative multiplets retained for statistical test
% biascorrection: apply or not bias correction for entropy calculation
% pathTmp: directory where to save the covariance matrix for the bootstrap
% groups: if you want to constrain the search to multiplets of variables belonging to different groups, provide a vector of length equal to the number of variables, whose entries are the group assignment of each variable

Xfull = copnorm(ts);

%% Parameters

% matrix size
[N, nvartot] = size(Xfull);
X = Xfull;

nboot = 1000; % number of bootstrap samples
alphaval = .05;
Otot(maxsize) = struct('index_var_red', [], 'sorted_red', [], 'index_red', [], 'bootsig_red', [], 'bootsigCI_red', [],...
    'index_var_syn', [], 'sorted_syn', [], 'index_syn', [], 'bootsig_syn', [], 'bootsigCI_syn', []);
O_val_size_tot(maxsize) = struct('multiplet_value',[]);
if nargin<6
    groups = ones(nvartot,1);
end

%% create data for CI estimation by Bootstrap

[outBoot, outJack] = hoi_createBootsData(X,nboot,pathTmp);

%% this section is for the expansion of redundancy, so maximizing the O
% there's no need to fix the target here

load(outBoot, 'covBootst'); % varible: covBootst

for isize = 3:maxsize
    
    boolPos = false; boolNeg = false; isLoadBoot = false;
    
    for iboot = 1:nboot+1
        if iboot == 1 % compute O, observer data
            C = nchoosek(1:nvartot,isize);
            if max(groups)>1
                C = select_combinations(C,groups);
            end
            
            ncomb = size(C,1);
            O_val_size = zeros(ncomb,1);
            XX = covBootst(:,:,iboot);
            
            %----- for bias correction
            if biascorrection
                psiterms{1} = psi((N - (1:isize))/2) / 2;
                psiterms{2} = psi((N - (1:isize-1))/2) / 2;
                psiterms{3} = psi((N - (1))/2) / 2;
                ln2 = log(2);
                dterm = (ln2 - log(N-1)) / 2;
            else
                psiterms = cell(3,1);
                dterm = [];
            end
            %------
            
            parfor icomb = 1:ncomb
                O_val_size(icomb) = hoi_o_information_boot(XX,C(icomb,:),biascorrection,psiterms,dterm);
            end
            O_val_size_tot(isize).multiplet_val=O_val_size; % here we save all the values
            % select n_best highest values
            ind_pos = find(O_val_size > 0);
            ind_neg = find(O_val_size < 0);
            O_pos = O_val_size(O_val_size>0);
            O_neg = O_val_size(O_val_size<0);
            [Osort_pos, ind_pos_sort] = sort(O_pos,'descend');
            [Osort_neg, ind_neg_sort] = sort(O_neg);
            if ~isempty(Osort_pos)
                n_sel_pos = min(n_best,length(Osort_pos));
            end
            if ~isempty(Osort_neg)
                n_sel_neg = min(n_best,length(Osort_neg));
            end
        else
            if ~isempty(Osort_pos) && boolPos==false
                boot_sig_pos = zeros(n_sel_pos,1); p_pos = boot_sig_pos;
                bstats_pos = zeros(nboot,n_sel_pos);
                boolPos = true;
            end
            if ~isempty(Osort_neg) && boolNeg==false
                boot_sig_neg = zeros(n_sel_neg,1); p_neg = boot_sig_neg;
                bstats_neg = zeros(nboot,n_sel_neg);
                boolNeg = true;
            end
            % compute bootstrapped measure for CI estimation
            if (boolPos || boolNeg) && isLoadBoot==false
                covJack = [];
                load(outJack, 'covJack'); % variable: covJack
                isLoadBoot = true;
            end
            if boolPos
                for isel = 1:n_sel_pos
                    indvar = C(ind_pos(ind_pos_sort(isel)),:);
                    bstats_pos(iboot-1,isel) = hoi_o_information_boot(covBootst(:,:,iboot),indvar,...
                        biascorrection,psiterms,dterm);
                end
            end
            if boolNeg
                for isel = 1:n_sel_neg
                    indvar = C(ind_neg(ind_neg_sort(isel)),:);
                    bstats_neg(iboot-1,isel) = hoi_o_information_boot(covBootst(:,:,iboot),indvar,...
                        biascorrection,psiterms,dterm);
                end
            end
        end
    end
    
    % estimate CI and apply FDR
    
    %----- for bias correction\
    if biascorrection
        psiterms{1} = psi((N - 1 - (1:isize))/2) / 2;
        psiterms{2} = psi((N -1 - (1:isize-1))/2) / 2;
        psiterms{3} = psi((N -1 - (1))/2) / 2;
        dterm = (ln2 - log(N-2)) / 2;
    else
        psiterms = cell(3,1);
        dterm = [];
    end
    %-----
    
    if boolPos
        ci = zeros(n_sel_pos,2);
        for isel = 1:n_sel_pos
            
            indvar = C(ind_pos(ind_pos_sort(isel)),:);
            
            p_pos(isel) = (1+sum(bstats_pos(:,isel)<0)) / (nboot+1);
            
            % jackniffe
            z_0 = fz0(bstats_pos(:,isel),Osort_pos(isel));
            jstat = zeros(N,1);
            parfor i = 1:N
                jstat(i) = hoi_o_information_boot(covJack(:,:,i),indvar,biascorrection,psiterms,dterm);
            end
            ci(isel,:) = bootci_jack(N,jstat,alphaval,z_0,Osort_pos(isel),bstats_pos(:,isel));
            
            boot_sig_pos(isel) = ~((ci(isel,1)<0) && (ci(isel,2)>0));
        end
        
        h = fdr_bh(p_pos);
        Otot(isize).index_var_red = C(ind_pos(ind_pos_sort(1:n_sel_pos)),:);
        Otot(isize).sorted_red = Osort_pos(1:n_sel_pos);
        Otot(isize).index_red = ind_pos(ind_pos_sort(1:n_sel_pos));
        Otot(isize).bootsig_red = h.*boot_sig_pos;
        Otot(isize).bootsigCI_red = ci;
    end
    if boolNeg
        ci = zeros(n_sel_neg,2);
        for isel = 1:n_sel_neg
            
            indvar = C(ind_neg(ind_neg_sort(isel)),:);
            
            p_neg(isel) = (1+sum(bstats_neg(:,isel)>0)) / (nboot+1);
            
            % jackniffe
            z_0 = fz0(bstats_neg(:,isel),Osort_neg(isel));
            jstat = zeros(N,1);
            parfor i = 1:N
                jstat(i) = hoi_o_information_boot(covJack(:,:,i),indvar,biascorrection,psiterms,dterm);
            end
            ci(isel,:) = bootci_jack(N,jstat,alphaval,z_0,Osort_neg(isel),bstats_neg(:,isel));
            
            boot_sig_neg(isel) = ~((ci(isel,1)<0) && (ci(isel,2)>0));
        end
        
        h = fdr_bh(p_neg);
        Otot(isize).index_var_syn = C(ind_neg(ind_neg_sort(1:n_sel_neg)),:);
        Otot(isize).sorted_syn = Osort_neg(1:n_sel_neg);
        Otot(isize).index_syn = ind_neg(ind_neg_sort(1:n_sel_neg));
        Otot(isize).bootsig_syn = h.*boot_sig_neg;
        Otot(isize).bootsigCI_syn = ci;
    end
end

for isize = 1:maxsize
    if ~isempty(Otot(isize).bootsig_red)
        Otot(isize).index_var_red(Otot(isize).bootsig_red==0,:)=[];
        Otot(isize).sorted_red(Otot(isize).bootsig_red==0)=[];
        Otot(isize).index_red(Otot(isize).bootsig_red==0)=[];
        Otot(isize).bootsigCI_red(Otot(isize).bootsig_red==0,:)=[];
        Otot(isize).bootsig_red(Otot(isize).bootsig_red==0)=[];
    end
    if ~isempty(Otot(isize).bootsig_syn)
        Otot(isize).index_var_syn(Otot(isize).bootsig_syn==0,:)=[];
        Otot(isize).sorted_syn(Otot(isize).bootsig_syn==0)=[];
        Otot(isize).index_syn(Otot(isize).bootsig_syn==0)=[];
        Otot(isize).bootsigCI_syn(Otot(isize).bootsig_syn==0,:)=[];
        Otot(isize).bootsig_syn(Otot(isize).bootsig_syn==0)=[];
    end
end
% and now flag the multiplets which don't have a significant increase of
% info with respect to their lower order composants
Otot=find_carryover_significance_zerolag(Otot);

for isize = 1:maxsize
    if ~isempty(Otot(isize).inc_sig_red)
        Otot(isize).index_var_red(Otot(isize).inc_sig_red==0,:)=[];
        Otot(isize).sorted_red(Otot(isize).inc_sig_red==0)=[];
        Otot(isize).index_red(Otot(isize).inc_sig_red==0)=[];
        Otot(isize).bootsigCI_red(Otot(isize).inc_sig_red==0,:)=[];
        Otot(isize).bootsig_red(Otot(isize).inc_sig_red==0)=[];
        Otot(isize).inc_sig_red(Otot(isize).inc_sig_red==0)=[];
    end
    if ~isempty(Otot(isize).inc_sig_syn)
        Otot(isize).index_var_syn(Otot(isize).inc_sig_syn==0,:)=[];
        Otot(isize).sorted_syn(Otot(isize).inc_sig_syn==0)=[];
        Otot(isize).index_syn(Otot(isize).inc_sig_syn==0)=[];
        Otot(isize).bootsigCI_syn(Otot(isize).inc_sig_syn==0,:)=[];
        Otot(isize).bootsig_syn(Otot(isize).inc_sig_syn==0)=[];
        Otot(isize).inc_sig_syn(Otot(isize).inc_sig_syn==0)=[];
    end
end

end

%% internal functions

function ci = bootci_jack(N,jstat,alpha,z_0,stat,bstat) % from bootci function
weights = repmat(1/N,N,1);
% acceleration finding, see DiCiccio and Efron (1996)
mjstat = sum(bsxfun(@times,jstat,weights),1); % mean along 1st dim.
score = bsxfun(@minus,mjstat,jstat); % score function at stat; ignore (N-1) factor because it cancels out in the skew
iszer = all(score==0,1);
skew = sum(bsxfun(@times,score.^3,weights),1) ./ ...
    (sum(bsxfun(@times,score.^2,weights),1).^1.5) /sqrt(N); % skewness of the score function
skew(iszer) = 0;
acc = skew/6;  % acceleration

% transform back with bias corrected and acceleration
z_alpha1 = norminv(alpha/2);
z_alpha2 = -z_alpha1;
pct1 = 100*normcdf(z_0 +(z_0+z_alpha1)./(1-acc.*(z_0+z_alpha1)));
pct1(z_0==Inf) = 100;
pct1(z_0==-Inf) = 0;
pct2 = 100*normcdf(z_0 +(z_0+z_alpha2)./(1-acc.*(z_0+z_alpha2)));
pct2(z_0==Inf) = 100;
pct2(z_0==-Inf) = 0;

% inverse of ECDF
m = numel(stat);
lower = zeros(1,m);
upper = zeros(1,m);
for i=1:m
    lower(i) = prctile(bstat(:,i),pct2(i),1);
    upper(i) = prctile(bstat(:,i),pct1(i),1);
end

% return
ci = sort([lower;upper],1);

end

% -------------------------
function z0 = fz0(bstat,stat) % from bootci function
% Compute bias-correction constant z0
z0 = norminv(mean(bsxfun(@lt,bstat,stat),1) + mean(bsxfun(@eq,bstat,stat),1)/2);
end



