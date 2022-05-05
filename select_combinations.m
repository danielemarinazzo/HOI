function Cg = select_combinations(C,groups)
    ngr = max(groups);
    dumm_gr = false([size(C), ngr]);
    for i = 1:ngr
        dumm_gr(:,:,i) = ismember(C,find(groups==i));
    end
    Cg = squeeze(sum(dumm_gr,2));
    Cg = sum(Cg>0,2) == ngr;
    Cg = C(Cg,:);
end