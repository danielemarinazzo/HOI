function [Otot, O_val_size_tot] = goi_exhaustive_loop_zerolag1(ts, maxorder, n_best, biascorrection, groups)

% ts: input (observations x variables), time series or static/behavioral data
% maxorder: maximum order for the gradients  
% n_best: number of most informative multiplets retained for statistical test
% biascorrection: apply or not bias correction for entropy calculation
% groups: if you want to constrain the search to multiplets of variables belonging to different groups, 
%             provide a vector of length equal to the number of variables, whose entries are the group assignment 
%             of each variable

X = copnorm(ts);
[N, nvartot] = size(X); % matrix size

%% Parameters

nboot = 1000; % number of bootstrap samples
alphaval = .05;
Otot(maxorder) = struct('index_var_red', [], 'sorted_red', [], 'index_red', [], 'bootsig_red', [], 'bootsigCI_red', [],...
    'index_var_syn', [], 'sorted_syn', [], 'index_syn', [], 'bootsig_syn', [], 'bootsigCI_syn', []);
O_val_size_tot(maxorder) = struct('multiplet_value',[]);
if nargin<5
    groups = ones(nvartot,1);
end
pathTmp = pwd;

%% create data for CI estimation by Bootstrap

[outBoot, outJack] = hoi_createBootsData(X,nboot,pathTmp);

%% this section is for the expansion of redundancy, so maximizing the O
% there's no need to fix the target here

load(outBoot, 'covBootst'); % varible: covBootst
nvartot_Vec = 1:nvartot;

for isize = 1:maxorder
    
    boolPos = false; boolNeg = false; isLoadBoot = false;
    
    for iboot = 1:nboot+1
        if iboot == 1 % compute O, observer data
            C = nchoosek(1:nvartot,isize);
            if max(groups)>1
                C = select_combinations(C,groups);
            end
            subset_index = dec2bin(0:2^isize-1) - '0';

            ncomb = size(C,1);
            O_val_size = zeros(ncomb,1);
            XX = covBootst(:,:,iboot);
            
            %----- for bias correction
            if biascorrection
               psiterms{3} = psi((N - (1))/2) / 2;
                ln2 = log(2);
                dterm = (ln2 - log(N-1)) / 2;
            else
                psiterms = cell(3,1);
                dterm = [];
            end
            %------
            
           parfor icomb = 1:ncomb
                c = C(icomb,:);
                if biascorrection
                    psitermsP = psiterms;
                end
                    
                for n = 1:length(subset_index) 
                    idx = c(subset_index(n,:)>0);
                    idx = setdiff(nvartot_Vec,idx);
            
                    if biascorrection
                        psitermsP{1} = psi((N - (1:length(idx)))/2) / 2;
                        psitermsP{2} = psi((N - (1:length(idx)-1))/2) / 2;
                    end
                    O_val_size(icomb) = O_val_size(icomb) + ...
                            (-1)^numel(idx)*hoi_o_information_boot(XX, idx, biascorrection, psitermsP, dterm);
                end                
           end
            
            O_val_size_tot(isize).multiplet_val = O_val_size; % here we save all the values
            
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
            
            %%%%%% check here, this loop is differentm until line 133
            if boolPos
                for isel = 1:n_sel_pos
                    indvar = C(ind_pos(ind_pos_sort(isel)),:);
                    for n = 1:length(subset_index)
                        idx = indvar(subset_index(n,:)>0);
                        idx_isel_pos{isel,n} = setdiff(nvartot_Vec,idx);
                    end
                end
            end
            if boolNeg
                for isel = 1:n_sel_neg
                    indvar = C(ind_neg(ind_neg_sort(isel)),:);
                    for n = 1:length(subset_index) 
                        idx = indvar(subset_index(n,:)>0);
                        idx_isel_neg{isel,n} = setdiff(nvartot_Vec,idx);
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            if boolPos
                for isel = 1:n_sel_pos
                    for n = 1:length(subset_index)
                        if biascorrection
                            psiterms{1} = psi((N - (1:length( idx_isel_pos{isel,n})))/2) / 2;
                            psiterms{2} = psi((N - (1:length( idx_isel_pos{isel,n})-1))/2) / 2;
                        end
                        bstats_pos(iboot-1,isel) = bstats_pos(iboot-1,isel) + ...
                                (-1)^numel(idx_isel_pos{isel,n})*hoi_o_information_boot(covBootst(:,:,iboot),  idx_isel_pos{isel,n}, biascorrection, psiterms, dterm);
                    end                                
                end
            end
            if boolNeg
                for isel = 1:n_sel_neg
                    for n = 1:length(subset_index) 
                        if biascorrection
                            psiterms{1} = psi((N - (1:length( idx_isel_neg{isel,n})))/2) / 2;
                            psiterms{2} = psi((N - (1:length( idx_isel_neg{isel,n})-1))/2) / 2;
                        end
                        bstats_neg(iboot-1,isel) = bstats_neg(iboot-1,isel) + ...
                                (-1)^numel(idx_isel_neg{isel,n})*hoi_o_information_boot(covBootst(:,:,iboot), idx_isel_neg{isel,n}, biascorrection, psiterms, dterm);
                    end
                end
            end
        end
    end
    
    % estimate CI and apply FDR
    
    %----- for bias correction\
    if biascorrection
        psiterms{1} = psi((N - 1 - (1:isize))/2) / 2;
        dterm = (ln2 - log(N-2)) / 2;
    else
        psiterms = cell(3,1);
        dterm = [];
    end
    %-----
    
    if boolPos
        ci = zeros(n_sel_pos,2);
        for isel = 1:n_sel_pos
            
            p_pos(isel) = (1+sum(bstats_pos(:,isel)<0)) / (nboot+1);
            
            % jackniffe
            z_0 = fz0(bstats_pos(:,isel),Osort_pos(isel));
            jstat = zeros(N,1);
            for i = 1:N
                if biascorrection
                    psitermsP1 = psiterms;
                end
                    
                for n = 1:length(subset_index) 
                    if biascorrection
                        psitermsP1{1} = psi((N - 1 - (1:length( idx_isel_pos{isel,n})))/2) / 2;
                        psitermsP1{2} = psi((N - 1-  (1:length( idx_isel_pos{isel,n})-1))/2) / 2;
                    end
                    jstat(i) = jstat(i) + ...
                            (-1)^numel(idx)*hoi_o_information_boot(covJack(:,:,i), idx_isel_pos{isel,n}, biascorrection, psitermsP1, dterm);
                end
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
                       
            p_neg(isel) = (1+sum(bstats_neg(:,isel)>0)) / (nboot+1);
            
            % jackniffe
            z_0 = fz0(bstats_neg(:,isel),Osort_neg(isel));
            jstat = zeros(N,1);
            for i = 1:N
                if biascorrection
                    psitermsP2 = psiterms;
                end
                    
                for n = 1:length(subset_index) 
                    if biascorrection
                        psitermsP2{1} = psi((N - 1 - (1:length( idx_isel_neg{isel,n})))/2) / 2;
                        psitermsP2{2} = psi((N - 1-  (1:length( idx_isel_neg{isel,n})-1))/2) / 2;
                    end
                    jstat(i) = jstat(i) + ...
                            (-1)^numel(idx)*hoi_o_information_boot(covJack(:,:,i), idx_isel_neg{isel,n}, biascorrection, psitermsP2, dterm);
                end
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

for isize = 1:maxorder
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

for isize = 4:maxorder
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



