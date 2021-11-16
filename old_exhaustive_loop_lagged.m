function Otot=exhaustive_loop_lagged(ts)
Xfull=copnorm(ts);
%Xfull=X;
%% Parameters
% matrix size
[N, nvartot]=size(Xfull);
X=Xfull;
modelorder = 3; % check this
maxsize = 6; % max number of variables in the multiplet
n_best = 10;  % number of most informative multiplets retained
nboot = 100; % number of bootstrap samples
chunklength = round(N/5); %can play around with this
alphaval=.05;
o_b=zeros(nboot,1);
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
            boot_sig=zeros(n_sel,1);
            for isel=1:n_sel
                indvar=C(ind_pos(ind_pos_sort(isel)),:);
                f  = @(xsamp) o_information_lagged_boot(t,X,modelorder,xsamp,chunklength,indvar);
                ci = bootci(nboot,{f,1:N-chunklength+1},'alpha',alphaval,'Options',statset('UseParallel',true));
                boot_sig(isel) = ~((ci(1)<0) && (ci(2)>0));
            end
            Otot(itarget,isize).sorted_red = Osort_pos(1:n_sel);
            Otot(itarget,isize).index_red = ind_pos(ind_pos_sort(1:n_sel));
            Otot(itarget,isize).bootsig_red = boot_sig;
        end
        if ~isempty(Osort_neg)
            n_sel = min(n_best,length(Osort_neg));
            boot_sig=zeros(n_sel,1);
            for isel=1:n_sel
                indvar=C(ind_neg(ind_neg_sort(isel)),:);
                f  = @(xsamp) o_information_lagged_boot(t,X,modelorder,xsamp,chunklength,indvar);
                ci = bootci(nboot,{f,1:1:N-chunklength+1},'alpha',alphaval);
                boot_sig(isel) = ~((ci(1)<0) && (ci(2)>0));
            end
            Otot(itarget,isize).sorted_syn = Osort_neg(1:n_sel);
            Otot(itarget,isize).index_syn = ind_neg(ind_neg_sort(1:n_sel));
            Otot(itarget,isize).bootsig_syn = boot_sig;
        end
    end
    %toc
end