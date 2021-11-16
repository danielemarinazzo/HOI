function Otot=exhaustive_loop_lagged_fdr(ts,modelorder,maxsize,n_best,isig)
if nargin<5
    isig=0;
end
Xfull=copnorm(ts);
%Xfull=X;
%% Parameters
% matrix size
[N, nvartot]=size(Xfull);
X=Xfull;
nboot = 1000; % number of bootstrap samples
chunklength = round(N/5); %can play around with this
alphaval=.05;
Otot(nvartot,maxsize) = struct();


%% this section is for the expansion of redundancy, so maximizing the O
for itarget=1:nvartot
    %tic
    t=X(:,itarget);
    for isize = 2:maxsize
        C = nchoosek(setdiff(1:nvartot,itarget),isize);
        ncomb = size(C,1);
        Osize = zeros(ncomb,1);
        for icomb = 1:ncomb
%             if mod(icomb,floor(ncomb/10))==1
%                 fprintf('target %d, size %d, fraction completed %4.2f.\n', itarget, isize, icomb/ncomb);
%             end
            %disp([itarget isize icomb/ncomb])%Y,X,m,indstart,chunklength,indvar
            Osize(icomb) = o_information_lagged_boot(t,X,modelorder,1:N,0,C(icomb,:));
        end
        %disp('ranking multiplets and bootstrapping')
        ind_pos = find(Osize > 0);
        ind_neg = find(Osize < 0);
        O_pos = Osize(Osize>0);
        O_neg = Osize(Osize<0);
        [Osort_pos, ind_pos_sort] = sort(O_pos,'descend');
        [Osort_neg, ind_neg_sort] = sort(O_neg);
        if ~isempty(Osort_pos)
            n_sel = min(n_best,length(Osort_pos));
            boot_sig=zeros(n_sel,1);p=boot_sig;
            for isel=1:n_sel
                indvar=C(ind_pos(ind_pos_sort(isel)),:);
                f  = @(xsamp) o_information_lagged_boot(t,X,modelorder,xsamp,chunklength,indvar);
                [ci,bstats] = bootci(nboot,{f,1:1:N-chunklength+1},'alpha',alphaval,'Options',statset('UseParallel',true));
                p(isel) = (1+sum(bstats < 0)) / (nboot+1);
                boot_sig(isel) = ~((ci(1)<0) && (ci(2)>0));
            end
            h=fdr_bh(p);
            Otot(itarget,isize).sorted_red = Osort_pos(1:n_sel);
            Otot(itarget,isize).index_red = ind_pos(ind_pos_sort(1:n_sel));
            Otot(itarget,isize).bootsig_red = h.*boot_sig;
        end
        if ~isempty(Osort_neg)
            n_sel = min(n_best,length(Osort_neg));
            boot_sig=zeros(n_sel,1);p=boot_sig;
            for isel=1:n_sel
                indvar=C(ind_neg(ind_neg_sort(isel)),:);
                f  = @(xsamp) o_information_lagged_boot(t,X,modelorder,xsamp,chunklength,indvar);
                [ci,bstats] = bootci(nboot,{f,1:1:N-chunklength+1},'alpha',alphaval,'Options',statset('UseParallel',true));
                p(isel) = (1+sum(bstats < 0)) / (nboot+1);
                boot_sig(isel) = ~((ci(1)<0) && (ci(2)>0));
            end
            h=fdr_bh(p);
            Otot(itarget,isize).sorted_syn = Osort_neg(1:n_sel);
            Otot(itarget,isize).index_syn = ind_neg(ind_neg_sort(1:n_sel));
            Otot(itarget,isize).bootsig_syn = h.*boot_sig;
        end
    end
    %toc
end
%%
if isig
    for itarget=1:nvartot
        for isize=3:maxsize;if ~isempty(Otot(itarget,isize).bootsig_red); Otot(itarget,isize).sorted_red(Otot(itarget,isize).bootsig_red==0)=[];end;end
        for isize=3:maxsize;if ~isempty(Otot(itarget,isize).bootsig_red); Otot(itarget,isize).index_red(Otot(itarget,isize).bootsig_red==0)=[];end;end
        for isize=3:maxsize;if ~isempty(Otot(itarget,isize).bootsig_red); Otot(itarget,isize).bootsig_red(Otot(itarget,isize).bootsig_red==0)=[];end;end
        for isize=3:maxsize;if ~isempty(Otot(itarget,isize).bootsig_syn); Otot(itarget,isize).sorted_syn(Otot(itarget,isize).bootsig_syn==0)=[];end;end
        for isize=3:maxsize;if ~isempty(Otot(itarget,isize).bootsig_syn); Otot(itarget,isize).index_syn(Otot(itarget,isize).bootsig_syn==0)=[];end;end
        for isize=3:maxsize;if ~isempty(Otot(itarget,isize).bootsig_syn); Otot(itarget,isize).bootsig_syn(Otot(itarget,isize).bootsig_syn==0)=[];end;end
    end
end