function [Otot_lagged, O_val_size_tot_lagged] = hoi_exhaustive_loop_lagged_all(ts, modelorder, maxsize, biascorrection, groups)

% ts: input (observations x variables), time series or static/behavioral data
% model order: 
% maxsize:  max number of drivers in the multiplet (i.e. size 2 means a multipet of 3 variables, 
%                the target and two drivers). This is different from the zero lag version
% biascorrection: apply or not bias correction for entropy calculation
% groups: if you want to constrain the search to multiplets of variables belonging to different groups, 
%             provide a vector of length equal to the number of variables, whose entries are the group assignment 
%             of each variable

X = copnorm(ts);
[N, nvartot] = size(X); % matrix size

%% Parameters

nboot = 1000; % number of bootstrap samples
alphaval = .05;
chunklength = round(N/5); % can play around with this
Otot_lagged(nvartot,maxsize) = struct('index_var_red', [], 'sorted_red', [], 'index_red', [], 'bootsig_red', [], 'bootsigCI_red', [], ...
    'index_var_syn', [], 'sorted_syn', [], 'index_syn', [], 'bootsig_syn', [], 'bootsigCI_syn', []);
O_val_size_tot_lagged(nvartot,maxsize) = struct('multiplet_value',[]);
if nargin<5
    groups = ones(nvartot,1);
end
pathTmp = pwd;

%% create data for CI estimation by Bootstrap

[outBoot, outJack] = hoi_createBootsData_lagged(X,chunklength,modelorder,nboot,pathTmp);

%% this section is for the expansion of redundancy, so maximizing the O
% there's no need to fix the target here

load(outBoot, 'covBootst'); % varible: covBootst
isLoadBoot = false;

for itarget = 1:nvartot
    
    for isize = 2:maxsize
        
        boolPos = false; boolNeg =  false; 
        
        C = nchoosek(setdiff(1:nvartot,itarget),isize);
        if max(groups)>1
            C = select_combinations(C,groups);
        end

        ncomb = size(C,1);
        O_val_size_lagged = zeros(ncomb,1);
        XX = covBootst(:,:,1);

        %----- for bias correction
        if biascorrection
            Ntrl = N-modelorder;
            Nvarxyz = size(C,2)*modelorder+modelorder+1;
            psiterms = psi((Ntrl - (1:Nvarxyz))/2) / 2;
            ln2 = log(2);
            dterm = (ln2 - log(Ntrl-1)) / 2;
        else
            psiterms = [];
            dterm = [];
        end
        %------
                
        parfor icomb = 1:ncomb
            O_val_size_lagged(icomb) = hoi_o_information_lagged_boot(XX,itarget,C(icomb,:),modelorder, ...
                biascorrection,psiterms,dterm);
        end
        O_val_size_tot_lagged(itarget,isize).multiplet_val=O_val_size_lagged; % here we save all the values

        % select positive and negative values
        ind_pos = find(O_val_size_lagged > 0);
        ind_neg = find(O_val_size_lagged < 0);
        O_pos = O_val_size_lagged(O_val_size_lagged>0);
        O_neg = O_val_size_lagged(O_val_size_lagged<0);
        [Osort_pos, ind_pos_sort] = sort(O_pos,'descend');
        [Osort_neg, ind_neg_sort] = sort(O_neg);
         if ~isempty(Osort_pos)
            n_pos = length(Osort_pos);
         end
        if ~isempty(Osort_neg)
            n_neg = length(Osort_neg);
        end
        
        if ~isempty(Osort_pos) && boolPos==false
            boot_sig_pos = zeros(n_pos,1); p_pos = boot_sig_pos;
            boolPos = true;
        end
        if ~isempty(Osort_neg) && boolNeg==false
            boot_sig_neg = zeros(n_neg,1); p_neg = boot_sig_neg;
            boolNeg = true;
        end
        
        % compute bootstrapped measure for CI estimation
        if (boolPos || boolNeg) && isLoadBoot==false
            covJack = [];
            load(outJack, 'covJack'); % variable: covJack
            isLoadBoot = true;
        end
        
        %----- for bias correction\
        if biascorrection
            Ntrl = N-1-modelorder;
            Nvarxyz = size(C,2)*modelorder+modelorder+1;
            psiterms = psi((Ntrl - (1:Nvarxyz))/2) / 2;
            ln2 = log(2);
            dterm = (ln2 - log(Ntrl-1)) / 2;
        else
            psiterms = [];
            dterm = [];
        end
        %-----
        
        if boolPos
            ci = zeros(n_pos,2);
            icount_pos=0;
            for isel = 1:n_pos
                bstats_pos = zeros(nboot,1);
                indvar = C(ind_pos(ind_pos_sort(isel)),:);
                
                for iboot=1:nboot
                    bstats_pos(iboot) = hoi_o_information_lagged_boot(covBootst(:,:,iboot+1),itarget, indvar,...
                        modelorder, biascorrection,psiterms,dterm);
                end
                
                p_pos(isel) = (1+sum(bstats_pos(:,1)<0)) / (nboot+1);

                % jackniffe
                z_0 = fz0(bstats_pos,Osort_pos(isel));
                jstat = zeros(N,1);
                parfor i = 1:N
                    jstat(i) = hoi_o_information_lagged_boot(covJack(:,:,i),itarget,indvar,modelorder,biascorrection,psiterms,dterm);
                end
                ci(isel,:) = bootci_jack(N,jstat,alphaval,z_0,Osort_pos(isel),bstats_pos);

                boot_sig_pos(isel) = ~((ci(isel,1)<0) && (ci(isel,2)>0));
                h = fdr_bh(p_pos(1:isel));
                if h(end)==0
                    icount_pos=icount_pos+1;
                else
                    icount_pos=0;
                end
                if icount_pos>10
                    break;
                end            
            end
            
            htot = zeros(n_pos,1);
            htot(1:isel) = h;
            Otot_lagged(itarget,isize).index_var_red = C(ind_pos(ind_pos_sort(1:n_pos)),:);
            Otot_lagged(itarget,isize).sorted_red = Osort_pos(1:n_pos);
            Otot_lagged(itarget,isize).index_red = ind_pos(ind_pos_sort(1:n_pos));
            Otot_lagged(itarget,isize).bootsig_red = htot.*boot_sig_pos;
            Otot_lagged(itarget,isize).bootsigCI_red = ci(1:n_pos,:);
        end
                
        if boolNeg
            ci = zeros(n_neg,2);
            icount_neg=0;
        
            for isel = 1:n_neg
                bstats_neg = zeros(nboot,1);
                indvar = C(ind_neg(ind_neg_sort(isel)),:);
                
                 for iboot=1:nboot
                        bstats_neg(iboot) = hoi_o_information_lagged_boot(covBootst(:,:,iboot+1),itarget, indvar,...
                            modelorder, biascorrection,psiterms,dterm);
                 end
                 
                 p_neg(isel) = (1+sum(bstats_neg>0)) / (nboot+1);
                
                % jackniffe
                z_0 = fz0(bstats_neg,Osort_neg(isel));
                jstat = zeros(N,1);
                parfor i = 1:N
                    jstat(i) = hoi_o_information_lagged_boot(covJack(:,:,i),itarget,indvar,modelorder,biascorrection,psiterms,dterm);
                end
                ci(isel,:) = bootci_jack(N,jstat,alphaval,z_0,Osort_neg(isel),bstats_neg);
                
                boot_sig_neg(isel) = ~((ci(isel,1)<0) && (ci(isel,2)>0));
                
                 h = fdr_bh(p_neg(1:isel));
                if h(end)==0
                    icount_neg=icount_neg+1;
                else
                    icount_neg=0;
                end
                if icount_neg>10
                    break
                end
            end
            
            htot = zeros(n_neg,1);
            htot(1:isel) = h;
            Otot_lagged(itarget,isize).index_var_syn = C(ind_neg(ind_neg_sort(1:n_neg)),:);
            Otot_lagged(itarget, isize).sorted_syn = Osort_neg(1:n_neg);
            Otot_lagged(itarget,isize).index_syn = ind_neg(ind_neg_sort(1:n_neg));
            Otot_lagged(itarget,isize).bootsig_syn = htot.*boot_sig_neg(1:n_neg);
            Otot_lagged(itarget,isize).bootsigCI_syn = ci(1:n_neg,:);
        end
    end
end                
                   
for itarget = 1:nvartot
    for isize = 1:maxsize
        if ~isempty(Otot_lagged(itarget,isize).bootsig_red)
            Otot_lagged(itarget,isize).index_var_red(Otot_lagged(itarget,isize).bootsig_red==0,:)=[];
            Otot_lagged(itarget,isize).sorted_red(Otot_lagged(itarget,isize).bootsig_red==0)=[];
            Otot_lagged(itarget,isize).index_red(Otot_lagged(itarget,isize).bootsig_red==0)=[];
            Otot_lagged(itarget,isize).bootsigCI_red(Otot_lagged(itarget,isize).bootsig_red==0,:)=[];
            Otot_lagged(itarget,isize).bootsig_red(Otot_lagged(itarget,isize).bootsig_red==0)=[];
        end
        if ~isempty(Otot_lagged(itarget,isize).bootsig_syn)
            Otot_lagged(itarget,isize).index_var_syn(Otot_lagged(itarget,isize).bootsig_syn==0,:)=[];
            Otot_lagged(itarget,isize).sorted_syn(Otot_lagged(itarget,isize).bootsig_syn==0)=[];
            Otot_lagged(itarget,isize).index_syn(Otot_lagged(itarget,isize).bootsig_syn==0)=[];
            Otot_lagged(itarget,isize).bootsigCI_syn(Otot_lagged(itarget,isize).bootsig_syn==0,:)=[];
            Otot_lagged(itarget,isize).bootsig_syn(Otot_lagged(itarget,isize).bootsig_syn==0)=[];
        end
    end
end

% and now flag the multiplets which don't have a significant increase of
% info with respect to their lower order composants
Otot_lagged = find_carryover_significance_lagged(Otot_lagged);
for itarget = 1:nvartot
    for isize = 4:maxsize
        if ~isempty(Otot_lagged(isize).inc_sig_red)
            Otot_lagged(itarget,isize).index_var_red(Otot_lagged(isize).inc_sig_red==0,:)=[];
            Otot_lagged(itarget,isize).sorted_red(Otot_lagged(isize).inc_sig_red==0)=[];
            Otot_lagged(itarget,isize).index_red(Otot_lagged(isize).inc_sig_red==0)=[];
            Otot_lagged(itarget,isize).bootsigCI_red(Otot_lagged(isize).inc_sig_red==0,:)=[];
            Otot_lagged(itarget,isize).bootsig_red(Otot_lagged(isize).inc_sig_red==0)=[];
            Otot_lagged(itarget,isize).inc_sig_red(Otot_lagged(isize).inc_sig_red==0)=[];
        end
        if ~isempty(Otot_lagged(isize).inc_sig_syn)
            Otot_lagged(itarget,isize).index_var_syn(Otot_lagged(isize).inc_sig_syn==0,:)=[];
            Otot_lagged(itarget,isize).sorted_syn(Otot_lagged(isize).inc_sig_syn==0)=[];
            Otot_lagged(itarget,isize).index_syn(Otot_lagged(isize).inc_sig_syn==0)=[];
            Otot_lagged(itarget,isize).bootsigCI_syn(Otot_lagged(isize).inc_sig_syn==0,:)=[];
            Otot_lagged(itarget,isize).bootsig_syn(Otot_lagged(isize).inc_sig_syn==0)=[];
            Otot_lagged(itarget,isize).inc_sig_syn(Otot_lagged(isize).inc_sig_syn==0)=[];
        end
    end
end

end


%% internal functions

function ci = bootci_jack(N,jstat,alpha,z_0,stat,bstat) % from matlab bootci function
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
for i = 1:m
    lower(i) = prctile(bstat(:,i),pct2(i),1);
    upper(i) = prctile(bstat(:,i),pct1(i),1);
end

% return
ci = sort([lower;upper],1);
end

% -------------------------
function z0 = fz0(bstat,stat) % from matlab bootci function
% Compute bias-correction constant z0
z0 = norminv(mean(bsxfun(@lt,bstat,stat),1) + mean(bsxfun(@eq,bstat,stat),1)/2);
end
