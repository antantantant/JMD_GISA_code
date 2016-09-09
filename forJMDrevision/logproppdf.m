function p = logproppdf(x,y,An,dn)
%     e = An*x';

global eee
    if isempty(eee)
        eee = ones(size(An,1),1);
    end
    p = sum(log(chi2pdf(eee,ones(size(An,1),1))));
%     y = log(mvnpdf(x, w0, Sigma));