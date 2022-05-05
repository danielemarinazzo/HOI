maxsize=length(O_val_size_tot);
nbest_tot=[20 50 100 300 500];
mycol=parula(length(nbest_tot));
for isize=3:maxsize
    figure(isize);
    O=O_val_size_tot(isize).multiplet_val;
    Opos=O(O>0);
    Oneg=O(O<0);
    [Npos,edgespos]=histcounts(Opos);
    subplot(2,1,1);hold on
    pred(1)=bar(edgespos(2:end)-mean(diff(edgespos))/2,Npos);
    Cpos=cumsum(Npos(end:-1:1));
    [Nneg,edgesneg]=histcounts(Oneg);
    Cneg=cumsum(Nneg);
    subplot(2,1,2);hold on
    psyn(1)=bar(edgesneg(2:end)-mean(diff(edgesneg))/2,Nneg);
    nbest_label=string(nbest_tot);
    for ibest=1:length(nbest_tot)
        nbest=nbest_tot(ibest);
        if(max(Cpos>0))
            xcut_pos=edgespos(max(1,(length(Cpos)-find(Cpos<=nbest, 1, 'last' ))));
            subplot(2,1,1);hold on;
            pred(ibest+1)=line([xcut_pos xcut_pos]+mean(diff(edgespos))/2,[0 max(Npos)],'color',mycol(ibest,:),'LineWidth',2);
        end
        if(max(abs(Cneg))>0)
            xcut_neg=edgesneg(find(Cneg<=nbest, 1, 'last' )+1);
            subplot(2,1,2);hold on;
            psyn(ibest+1)=line([xcut_neg xcut_neg]+mean(diff(edgesneg))/2,[0 max(Nneg)],'color',mycol(ibest,:),'LineWidth',2);
        end
    end
    subplot(2,1,1);title('redundancy');legend(pred(2:length(nbest_tot)+1),nbest_label)
    subplot(2,1,2);title('synergy');legend(psyn(2:length(nbest_tot)+1),nbest_label)
end