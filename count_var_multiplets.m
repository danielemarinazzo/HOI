nvartot=7;
maxorder=size(Otot,2);
for iorder=3:maxorder
    eval(['count_red_order_' num2str(iorder) '=zeros(nvartot,1);']);
    eval(['count_syn_order_' num2str(iorder) '=zeros(nvartot,1);']);
    for ivar=1:nvartot
        eval(['count_red_order_' num2str(iorder) '(ivar,1)=length(find(Otot(' num2str(iorder) ...
            ').index_var_red==ivar));'])
        eval(['count_syn_order_' num2str(iorder) '(ivar,1)=length(find(Otot(' num2str(iorder) ...
            ').index_var_syn==ivar));'])
        eval(['save count_red_order_' num2str(iorder) '.txt count_red_order_' num2str(iorder) ...
            ' /ascii'])
        eval(['save count_syn_order_' num2str(iorder) '.txt count_syn_order_' num2str(iorder) ...
            ' /ascii'])
    end
end