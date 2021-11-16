function o = hoi_o_information_boot(X, indvar, biascorrection, psiterms, dterm)
    M = length(indvar);
    X = X(indvar,indvar);
    o = (M-2)*hoi_ent_g_COV(X,biascorrection,psiterms{1},dterm);
    for j = 1:M
        X1 = X; X1(:,j) = []; X1(j,:) = [];
        o = o + hoi_ent_g_COV(X(j,j),biascorrection,psiterms{3},dterm) - ...
                hoi_ent_g_COV(X1,biascorrection,psiterms{2},dterm);
    end
end