function gamma=CoeffOperGraph(op,alpha,beta,a,b)
% added weighted l1 norm: xiaoqun zhang 04/26/2016
[nD Level]=size(alpha);
% add weighted l1norm
if strcmp(op,'norm1') || strcmp(op,'norm2')||strcmp(op,'wnorm1')
    gamma=0;
end
for j=1:Level
    for k=1:nD
        if op=='-'
            gamma{k,j}=alpha{k,j}-beta{k,j};
        elseif op=='+'
            gamma{k,j}=alpha{k,j}+beta{k,j};
        elseif strcmp(op,'*c')
            gamma{k,j}=alpha{k,j}*beta;
        elseif strcmp(op,'*+')
            gamma{k,j}=alpha{k,j}*a+beta{k,j}*b;
        elseif op=='p' % projection onto l-infity ball,d{k}=CoeffOperGraph('p',d{k},Thresh); 
            if k==1
                % gamma{k,j} = beta{k,j};
                % modified by jiayong Liu,2017.10.23.Reference:Wavele frame
                % based multiphase segmentation,page 2531
                gamma{k,j} = beta{k,j}.*(beta{k,j} < abs(alpha{k,j})) + alpha{k,j} .* (beta{k,j} >= abs(alpha{k,j}));
            else
                % gamma{k,j} = beta{k,j}.*( abs(alpha{k,j}) >= beta{k,j} );
                % modified by jiayong Liu,2017.10.23.Reference:Wavele frame
                % based multiphase segmentation,page 2531
                gamma{k,j} = beta{k,j} .* sign(alpha{k,j});
                gamma{k,j} = gamma{k,j} + alpha{k,j}.*( abs(alpha{k,j}) < beta{k,j} );                
            end
        elseif op=='h'
            if k==1
                gamma{k,j}=alpha{k,j};
            else
                gamma{k,j}=alpha{k,j}.*(abs(alpha{k,j})>=beta{k,j});
            end
        elseif op=='s'
            if k==1
                gamma{k,j}=alpha{k,j};
            else
                gamma{k,j}=(alpha{k,j}-beta{k,j}.*sign(alpha{k,j})).*(abs(alpha{k,j})>=beta{k,j});
            end
        elseif strcmp(op,'s_band')
            for l=1:length(a)
                if k==a(l)
                    gamma{k,j}=(alpha{k,j}-beta{k,j}.*sign(alpha{k,j})).*(abs(alpha{k,j})>=beta{k,j});
                else
                    gamma{k,j}=alpha{k,j};
                end
            end
        elseif strcmp(op,'norm1')
            gamma=gamma+norm(alpha{k,j},1);
        elseif strcmp(op,'norm2')
            gamma=gamma+norm(alpha{k,j});
        elseif strcmp(op,'wnorm1')
            gamma=gamma+norm(alpha{k,j}.*beta{k,j},1);            
        end
    end
end