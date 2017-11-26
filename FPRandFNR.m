function [fpr,fnr] = FPRandFNR(ytrue,ypre)
    indexTrue = find(ytrue);
    indexPre = find(ypre);
    fpr = length(setdiff(indexPre,indexTrue))/length(indexPre);

    indexTrue = find(ytrue==0);
    indexPre = find(ypre==0);
    fnr = length(setdiff(indexPre,indexTrue))/length(indexPre);
end