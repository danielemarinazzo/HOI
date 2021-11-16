function Otot=exhaustive_loop_zerolag(ts,isig)
if nargin<2
    isig=0;
end
Xfull=copnorm(ts);
%Xfull=X;
%% Parameters
% matrix size
[N, nvartot]=size(Xfull);
X=Xfull;
maxsize = 6; % max number of variables in the multiplet
n_best = 20;  % number of most informative multiplets retained
nboot = 100; % number of bootstrap samples
alphaval=.005;
%o_b=zeros(nboot,1);
Otot(maxsize) = struct();


%% this section is for the expansion of redundancy, so maximizing the O
% there's no need to fix the target here
for isize = 3:maxsize
    C = nchoosek(1:nvartot,isize);
    ncomb = size(C,1);
    Osize = zeros(ncomb,1);
    for icomb = 1:ncomb
        %         if mod(icomb,floor(ncomb/10))==1
        %             fprintf('size %d, fraction completed %4.2f.\n',  isize, icomb/ncomb);
        %         end
        Osize(icomb) = o_information_boot(X,1:N,C(icomb,:));
    end
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
            f  = @(xsamp) o_information_boot(X,xsamp,indvar);
            ci = bootci(nboot,{f,1:N},'alpha',alphaval,'Options',statset('UseParallel',true));
            boot_sig(isel) = ~((ci(1)<0) && (ci(2)>0));
        end
        Otot(isize).sorted_red = Osort_pos(1:n_sel);
        Otot(isize).index_red = ind_pos(ind_pos_sort(1:n_sel));
        Otot(isize).bootsig_red = boot_sig;
    end
    if ~isempty(Osort_neg)
        n_sel = min(n_best,length(Osort_neg));
        boot_sig=zeros(n_sel,1);
        for isel=1:n_sel
            indvar=C(ind_neg(ind_neg_sort(isel)),:);
            f  = @(xsamp) o_information_boot(X,xsamp,indvar);
            ci = bootci(nboot,{f,1:N},'alpha',alphaval);
            boot_sig(isel) = ~((ci(1)<0) && (ci(2)>0));
        end
        Otot(isize).sorted_syn = Osort_neg(1:n_sel);
        Otot(isize).index_syn = ind_neg(ind_neg_sort(1:n_sel));
        Otot(isize).bootsig_syn = boot_sig;
    end
end
if isig
    for i=3:maxsize;if ~isempty(Otot(i).bootsig_red); Otot(i).sorted_red(Otot(i).bootsig_red==0)=[];end;end
    for i=3:maxsize;if ~isempty(Otot(i).bootsig_red); Otot(i).index_red(Otot(i).bootsig_red==0)=[];end;end
    for i=3:maxsize;if ~isempty(Otot(i).bootsig_red); Otot(i).bootsig_red(Otot(i).bootsig_red==0)=[];end;end
    for i=3:maxsize;if ~isempty(Otot(i).bootsig_syn); Otot(i).sorted_syn(Otot(i).bootsig_syn==0)=[];end;end
    for i=3:maxsize;if ~isempty(Otot(i).bootsig_syn); Otot(i).index_syn(Otot(i).bootsig_syn==0)=[];end;end
    for i=3:maxsize;if ~isempty(Otot(i).bootsig_syn); Otot(i).bootsig_syn(Otot(i).bootsig_syn==0)=[];end;end
end