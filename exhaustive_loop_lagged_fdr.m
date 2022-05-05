function [Otot_lagged, O_val_size_tot_lagged] = exhaustive_loop_lagged_fdr(ts,modelorder,maxsize,n_best)

Xfull=copnorm(ts);
%Xfull=X;
%% Parameters
% matrix size
[N, nvartot]=size(Xfull);
X=Xfull;
nboot = 1000; % number of bootstrap samples
chunklength = round(N/5); %can play around with this
alphaval=.05;
Otot_lagged(nvartot,maxsize) = struct('index_var_red', [], 'sorted_red', [], 'index_red', [], 'bootsig_red', [], 'bootsigCI_red', [],...
    'index_var_syn', [], 'sorted_syn', [], 'index_syn', [], 'bootsig_syn', [], 'bootsigCI_syn', []);
O_val_size_tot_lagged(nvartot,maxsize) = struct('multiplet_value',[]);

%% this section is for the expansion of redundancy, so maximizing the O
for itarget=1:nvartot
    %tic
    t=X(:,itarget);
    for isize = 2:maxsize
        C = nchoosek(setdiff(1:nvartot,itarget),isize);
        ncomb = size(C,1);
        O_val_size_lagged = zeros(ncomb,1);
        for icomb = 1:ncomb
            %             if mod(icomb,floor(ncomb/10))==1
            %                 fprintf('target %d, size %d, fraction completed %4.2f.\n', itarget, isize, icomb/ncomb);
            %             end
            %disp([itarget isize icomb/ncomb])%Y,X,m,indstart,chunklength,indvar
            O_val_size_lagged(icomb) = o_information_lagged_boot(t,X,modelorder,1:N,0,C(icomb,:));
        end
        O_val_size_tot_lagged(itarget,isize).multiplet_val=O_val_size_lagged; % here we save all the values
        %disp('ranking multiplets and bootstrapping')
        ind_pos = find(O_val_size_lagged > 0);
        ind_neg = find(O_val_size_lagged < 0);
        O_pos = O_val_size_lagged(O_val_size_lagged>0);
        O_neg = O_val_size_lagged(O_val_size_lagged<0);
        [Osort_pos, ind_pos_sort] = sort(O_pos,'descend');
        [Osort_neg, ind_neg_sort] = sort(O_neg);
        if ~isempty(Osort_pos)
            n_sel = min(n_best,length(Osort_pos));
            boot_sig=zeros(n_sel,1);p=boot_sig;
            ci = zeros(n_sel,2);
            for isel=1:n_sel
                indvar=C(ind_pos(ind_pos_sort(isel)),:);
                f  = @(xsamp) o_information_lagged_boot(t,X,modelorder,xsamp,chunklength,indvar);
                [ci(isel,:),bstats] = bootci(nboot,{f,1:1:N-chunklength+1},'alpha',alphaval,'Options',statset('UseParallel',true));
                p(isel) = (1+sum(bstats < 0)) / (nboot+1);
                boot_sig(isel) = ~((ci(isel,1)<0) && (ci(isel,2)>0));
            end
            h=fdr_bh(p);
            Otot_lagged(itarget,isize).sorted_red = Osort_pos(1:n_sel);
            Otot_lagged(itarget,isize).index_red = ind_pos(ind_pos_sort(1:n_sel));
            Otot_lagged(itarget,isize).bootsig_red = h.*boot_sig;
            Otot_lagged(itarget,isize).bootsigCI_red = ci;
        end
        if ~isempty(Osort_neg)
            n_sel = min(n_best,length(Osort_neg));
            boot_sig=zeros(n_sel,1);p=boot_sig;
            ci = zeros(n_sel,2);
            for isel=1:n_sel
                indvar=C(ind_neg(ind_neg_sort(isel)),:);
                f  = @(xsamp) o_information_lagged_boot(t,X,modelorder,xsamp,chunklength,indvar);
                [ci(isel,:),bstats] = bootci(nboot,{f,1:1:N-chunklength+1},'alpha',alphaval,'Options',statset('UseParallel',true));
                p(isel) = (1+sum(bstats < 0)) / (nboot+1);
                boot_sig(isel) = ~((ci(isel,1)<0) && (ci(isel,2)>0));
            end
            h=fdr_bh(p);
            Otot_lagged(itarget,isize).sorted_syn = Osort_neg(1:n_sel);
            Otot_lagged(itarget,isize).index_syn = ind_neg(ind_neg_sort(1:n_sel));
            Otot_lagged(itarget,isize).bootsig_syn = h.*boot_sig;
            Otot_lagged(itarget,isize).bootsigCI_syn = ci;
        end
    end
    %toc
end

%%
for itarget=1:nvartot
    for isize=3:maxsize;if ~isempty(Otot_lagged(itarget,isize).bootsig_red); Otot_lagged(itarget,isize).sorted_red(Otot_lagged(itarget,isize).bootsig_red==0)=[];end;end
    for isize=3:maxsize;if ~isempty(Otot_lagged(itarget,isize).bootsig_red); Otot_lagged(itarget,isize).index_red(Otot_lagged(itarget,isize).bootsig_red==0)=[];end;end
    for isize=3:maxsize;if ~isempty(Otot_lagged(itarget,isize).bootsig_red); Otot_lagged(itarget,isize).bootsig_red(Otot_lagged(itarget,isize).bootsig_red==0)=[];end;end
    for isize=3:maxsize;if ~isempty(Otot_lagged(itarget,isize).bootsig_syn); Otot_lagged(itarget,isize).sorted_syn(Otot_lagged(itarget,isize).bootsig_syn==0)=[];end;end
    for isize=3:maxsize;if ~isempty(Otot_lagged(itarget,isize).bootsig_syn); Otot_lagged(itarget,isize).index_syn(Otot_lagged(itarget,isize).bootsig_syn==0)=[];end;end
    for isize=3:maxsize;if ~isempty(Otot_lagged(itarget,isize).bootsig_syn); Otot_lagged(itarget,isize).bootsig_syn(Otot_lagged(itarget,isize).bootsig_syn==0)=[];end;end
end


%%

Otot_lagged=find_carryover_significance_lagged(Otot_lagged);

for itarget = 1:nvartot
    for isize = 1:maxsize
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
