function o = hoi_o_information_lagged_boot(X, itarget, indvar, modelorder, biascorrection, psiterms, dterm)
% evaluates the o_information flow
% Y Nx1 target vector .  X NxM drivers
% m order of the model
M = length(indvar);
nvar = size(X,1)/(modelorder+1);

indvarAll = [];
itargetAll = [];
for i = 1:modelorder
    indvarAll = [indvarAll indvar+(nvar*(i-1))];
    itargetAll = [itargetAll itarget+(nvar*(i-1))];
end

itarget = itarget+(nvar*(modelorder));
X1 = X([itarget, itargetAll indvarAll],[itarget, itargetAll indvarAll]);
o = -(M-1)*hoi_ent_ggg_COV(X1,modelorder,biascorrection, psiterms, dterm);
for k = 1:M
    indvar1 = indvar;
    indvar1(k) = [];
    
    indvarAll = [];
    for i = 1:modelorder
        indvarAll = [indvarAll indvar1+(nvar*(i-1))];
    end
    
    X1 = X([itarget, itargetAll indvarAll],[itarget, itargetAll indvarAll]);
    o = o+hoi_ent_ggg_COV(X1,modelorder,biascorrection, psiterms, dterm);
end