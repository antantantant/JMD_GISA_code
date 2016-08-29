function f = appDistribution(info,pairs,p)
    nx = info.nx;
    X = info.X;
    d = info.d;
    f = zeros(nx,1);
    
    if ~isempty(pairs)
        best = pairs(:,1);
        for i = 1:length(best)
            if find(pairs(:,2)==best(i))
                best(i)=0;
            end
        end
        best = unique(best);
        best(best==0)=[];
        worse = unique([find(p<1e-3/info.nx);pairs(:,2)]);
        rest = 1:nx;
        rest([best',worse'])=[];

        for i = 1:length(rest)
            Z = zeros(size(pairs,1)+nx-1,d);
            count = 1;
            for j = [1:rest(i)-1,rest(i)+1:nx]
                Z(count,:) = X(rest(i),:)-X(j,:);
                count = count + 1;
            end
            for j = 1:size(pairs,1)
                Z(count,:) = X(pairs(j,1),:)-X(pairs(j,2),:);
                count = count + 1;
            end
            Z = bsxfun(@rdivide,Z,sqrt(dot(Z,Z,2)));
            Z = sparse(Z);
            model = train(ones(size(Z,1),1),Z,'-s 0 -c 1e6 -e 1e-6 -q');
            wh = model.w';
            wh = wh/norm(wh);
            f(rest(i))=max(min(Z*wh),0)^(info.d-1);
        end

        for i = 1:length(best)
            pairid = find(pairs(:,1)~=best(i));
            Z = zeros(length(pairid)+nx-1,d);
            count = 1;
            for j = [1:best(i)-1,best(i)+1:nx]
                Z(count,:) = X(best(i),:)-X(j,:);
                count = count + 1;
            end
            for j = 1:length(pairid)
                Z(count,:) = X(pairs(pairid(j),1),:)-X(pairs(pairid(j),2),:);
                count = count + 1;
            end
            Z = bsxfun(@rdivide,Z,sqrt(dot(Z,Z,2)));
            Z = sparse(Z);
            model = train(ones(size(Z,1),1),Z,'-s 0 -c 1e6 -e 1e-6 -q');
            wh = model.w';
            wh = wh/norm(wh);
            f(best(i))=max(min(Z*wh),0)^(info.d-1);
        end
    else
        for i = 1:nx
            Z = zeros(nx-1,d);
            count = 1;
            for j = [1:(i-1),(i+1):nx]
                Z(count,:) = X(i,:)-X(j,:);
                count = count + 1;
            end
            Z = bsxfun(@rdivide,Z,sqrt(dot(Z,Z,2)));
            Z = sparse(Z);
            model = train(ones(size(Z,1),1),Z,'-s 0 -c 1e6 -e 1e-6 -q');
            wh = model.w';
            wh = wh/norm(wh);
            f(i)=max(min(Z*wh),0)^(info.d-1);
        end        
    end