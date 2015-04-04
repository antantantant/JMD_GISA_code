% conditional probabilities if all queries known

% s = 1e1;
% scale = 1e-1;
% keep = true(nx,1);
% A = bsxfun(@times,dX,2*(dX*w'>0)-1);
% A(abs(dX*w')<1e-3,:) = [];
% rid = [1:4,6:9,11:14,16:19,21:24,26:29];
% A = unique(A,'rows');
% AA = A(:,rid);
% bb = -A(:,setdiff(1:30,rid))*w(setdiff(1:30,rid))'*scale;
% W = zeros(s,d);
% count = 0;

% while count<s && iter <10
%     ww = mvtnconstrained(s,length(rid),[],[],AA,bb,...
%         w(rid)*scale,scale);
%     WC = zeros(s,d);
%     WC(:,5) = sum(w(1:5))-sum(ww(:,1:4),2);
%     WC(:,10) = sum(w(6:10))-sum(ww(:,5:8),2);
%     WC(:,15) = sum(w(11:15))-sum(ww(:,9:12),2);
%     WC(:,20) = sum(w(16:20))-sum(ww(:,13:16),2);
%     WC(:,25) = sum(w(21:25))-sum(ww(:,17:20),2);
%     WC(:,30) = sum(w(26:30))-sum(ww(:,21:24),2);
%     WC(:,rid) = ww;
%     WC = WC(sum(WC*A'>0,2)==size(A,1),:);
%     W(count+(1:size(WC,1)),:) = WC;
%     count = count + size(WC,1);
%     iter = iter + 1;
% end
% W = W(1:min(count,s),:);

W = truncnormrnd([1e4,1],0,1,0,100)*w;
% ww = randn(1e5,30);
% rr = sqrt(sum(ww.^2,2));
% rr = chi2rnd(30,[1e4,1]);
% W = rr*w/norm(w);
pref = W*Xf';
area = (sum(bsxfun(@eq, pref, max(pref,[],2)))/s)';
const = 1/sum(area>0)./area;
const(area==0) = 0;
h = sparse(bsxfun(@eq, pref, max(pref,[],2))*const);

obj = bsxfun(@times,1./(1+exp(-W*Xf(keep,:)')),c(keep)');
f(keep) = sum(bsxfun(@times, bsxfun(@eq, obj, max(obj,[],2))/s, h))';
probability_obj = f/sum(f);
