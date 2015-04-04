function W = sampling(s,d,w,Xf,dX,dXID,pairs)
    
    if isempty(pairs)
        W =  randn(s,d); % generate candidates
    else
        W = zeros(s,d);
        count = 0;
        iter = 0;
        A = dX(dXID(sub2ind(size(dXID),pairs(:,1),pairs(:,2))),:);
        A = bsxfun(@times,A,sign(pairs(:,2)-pairs(:,1)));
        b = zeros(size(pairs,1),1);
%         rid = [1:4,6:9,11:14,16:19,21:24,26:29];
        rid = 1:30;
        if size(w,2)>1
%             w0 = feasible(w,A,b);
%             if ~isempty(w0)
%                 w0 = w0(1,:);
%             else
%                 model = train([ones(size(A,1),1);-ones(size(A,1),1)],...
%                     sparse([A;-A]),'-s 0 -c 1e6 -e 1e-6 -q');
%                 w0 = model.w;
%             end
            W = feasible(w,A,b);
            if isempty(W)
                model = train([ones(size(A,1),1);-ones(size(A,1),1)],...
                sparse([A;-A]),'-s 0 -c 1e6 -e 1e-6 -q');
                w0 = model.w;
                w0 = w0/norm(w0);
%                 while count<s && iter <10
%                     try
% %                         WC = mvtnconstrained(s,length(rid),[],[],A(:,rid),b,...
% %                             w0);
                        WC = mvtrcnml(s,length(rid),A(:,rid)',b,w0);
                        W = WC(:,s+1:end)';
%                     catch
%                         W = [];
%                         break;
%                     end
%         %             WC = generate(WC,s,d);
%                     WC = feasible(WC,A,b);
%                     W(count+(1:size(WC,1)),:) = WC;
%                     count = count + size(WC,1);
%                     iter = iter + 1;
%                 end
%                 if ~isempty(W)
%                     W = W(1:min(count,s),:);
%                 end
            end
        else
            model = train([ones(size(A,1),1);-ones(size(A,1),1)],...
                sparse([A;-A]),'-s 0 -c 1e6 -e 1e-6 -q');
            w0 = model.w;
            w0 = w0/norm(w0);
%             while count<s && iter <10
%                 try
%                     WC = mvtnconstrained(s,length(rid),[],[],A(:,rid),b,...
%                         w0);
                    WC = mvtrcnml(s,length(rid),A(:,rid)',b,w0);
                    W = WC(:,s+1:end)';
%                 catch
%                     W = [];
%                     break;
%                 end
    %             WC = generate(WC,s,d);
%                 WC = feasible(WC,A,b);
%                 W(count+(1:size(WC,1)),:) = WC;
%                 count = count + size(WC,1);
%                 iter = iter + 1;
%             end
%             if ~isempty(W)
%                 W = W(1:min(count,s),:);
%             end
        end
    end
        
%     % deal with large matrices
%     if s>1e6
%         segl = 1e3;
%         segn = ceil(s/segl);
%         area = zeros(size(Xf,1),1);
%         pref_id = zeros(s,1);
%         for i = 1:segn
%             id = ((i-1)*segl+1):min(i*segl,s);
%             pref = W(id,:)*Xf';
%             temp = bsxfun(@eq, pref, max(pref,[],2));
%             pref_id(id) = temp*(1:size(Xf,1))';
%             area = area + (sum(temp)/s)';
%         end
%         const = 1/sum(area>0)./area;
%         const(area==0) = 0;
%         h = sparse(area(pref_id));
%     else
%         pref = W*Xf';
%         area = (sum(bsxfun(@eq, pref, max(pref,[],2)))/s)';
%         const = 1/sum(area>0)./area;
%         const(area==0) = 0;
%         h = sparse(bsxfun(@eq, pref, max(pref,[],2))*const);
%     end
    
    function WC = generate(W,s,d)
%         WC = randn(s,d); % generate candidates
%         WC = zeros(s,d);
%         WC(:,5) = -0.02-sum(W(:,1:4),2);
%         WC(:,10) = -0.02-sum(W(:,5:8),2);
%         WC(:,15) = -0.02-sum(W(:,9:12),2);
%         WC(:,20) = -0.02-sum(W(:,13:16),2);
%         WC(:,25) = -0.02-sum(W(:,17:20),2);
%         WC(:,30) = -0.02-sum(WC(:,21:24),2);
        WC = W;
    
%     function WC = feasible(WC,dX,dXID,pairs)
    function WC = feasible(WC,A,b)
%         segn = toobig(WC,size(pairs,1));
%         if (segn>1)
%             feasible = true(s,1);
%             for i = 1:segn
%                 segl = round(size(dX,1)/segn);
%                 pref = zeros(size(WC,1),segl);
%                 id = ((i-1)*segl+1):min(i*segl,size(dX,1));
%                 pref(feasible,:) = WC(feasible,:)*dX(dXID(sub2ind(size(dXID),pairs(id,1),pairs(id,2))),:)';
%                 ff = bsxfun(@times,pref,(pairs(id,2)-pairs(id,1))')>0;
%                 feasible = feasible&(sum(ff,2)==segl);
%             end
%         else
%             pref = WC*dX(dXID(sub2ind(size(dXID),pairs(:,1),pairs(:,2))),:)';
%             feasible = bsxfun(@times,pref,(pairs(:,2)-pairs(:,1))')>0;
%             feasible = sum(feasible,2)==size(pairs,1);
%         end
        feasible = sum(bsxfun(@lt, A*WC', b),1)==0;
        WC = WC(feasible,:);

    function n = toobig(A,B)
        l1 = size(A,1);
        l2 = size(B,1);
        n = 1;
        if l1*l2>1e6
            n = round(l2/100);
        end