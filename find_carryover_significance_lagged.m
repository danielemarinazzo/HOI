function Otot=find_carryover_significance_lagged(Otot)
% this checks whether the increment of the O information going one order up
% is significant with respect to the value at the lower order

[nvartot, maxsize]= size(Otot);
for itarget = 1:nvartot
    for isize=maxsize:-1:3
        higher=isize;
        nlower_red=size(Otot(itarget,higher-1).index_var_red,1);
        nlower_syn=size(Otot(itarget,higher-1).index_var_syn,1);
        nhigher_red=size(Otot(itarget,higher).index_var_red,1);
        nhigher_syn=size(Otot(itarget,higher).index_var_syn,1);
        sig_red=zeros(nhigher_red,1);
        sig_syn=zeros(nhigher_syn,1);
        for ilower=1:nlower_red
            for ihigher=1:nhigher_red
                if sum(ismember(Otot(itarget,higher-1).index_var_red(ilower,:),Otot(itarget,higher).index_var_red(ihigher,:)))==higher-1
                    v1=Otot(itarget,higher-1).bootsigCI_red(ilower,:);
                    v2=Otot(itarget,higher).bootsigCI_red(ihigher,:);
                    l1=v1(1);l2=v2(1);r1=v1(2);r2=v2(2);
                    r = bsxfun(@min, r1, r2);
                    l = bsxfun(@max, l1, l2);
                    sig_red(ihigher)= l >= r;
                end
            end
        end
        Otot(itarget,isize).inc_sig_red=sig_red;
        for ilower=1:nlower_syn
            for ihigher=1:nhigher_syn
                if sum(ismember(Otot(itarget,higher-1).index_var_syn(ilower,:),Otot(itarget,higher).index_var_syn(ihigher,:)))==higher-1
                    v1=Otot(itarget,higher-1).bootsigCI_syn(ilower,:);
                    v2=Otot(itarget,higher).bootsigCI_syn(ihigher,:);
                    l1=v1(1);l2=v2(1);r1=v1(2);r2=v2(2);
                    r = bsxfun(@min, r1, r2);
                    l = bsxfun(@max, l1, l2);
                    sig_syn(ihigher)= l >= r;
                end
            end
        end
        Otot(itarget,isize).inc_sig_syn=sig_syn;
    end
end