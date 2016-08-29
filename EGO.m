function [E,wh,pairs] = EGO(info,pairs,w,bestID,nq)

    % calculate A
    if ~isempty(pairs)
        count_set = appDistribution(info,pairs,[]);
        A = info.theta'*count_set;
        probability_obj = info.theta.*count_set/A;
    else
        probability_obj = info.theta;
        A = 1;
    end
    
    % calculate current entropy
    E = -sum(probability_obj.*log2(probability_obj+1e-99));
    disp('*');
    
    % find the next query
    if  nq>0
        if sum(probability_obj<(1/info.nx*1e-3))==info.nx-1&&sum(pairs(:,1)==bestID)>0
            Z = zeros(size(pairs,1),info.d);
            for i = 1:size(pairs,1)
                Z(i,:) = info.X(pairs(i,1),:)-info.X(pairs(i,2),:);
            end
            Z = bsxfun(@rdivide,Z,sqrt(dot(Z,Z,2)));
            Z = sparse(Z);
            model = train(ones(size(Z,1),1),Z,'-s 0 -c 1e6 -e 1e-6 -q');
            wh = model.w';
        else
            % pick the next query that minimizes the expectation of entropy
            % current wh
            Z = zeros(size(pairs,1),info.d);
            unique_pairs = unique(pairs);
            XX = zeros(length(unique_pairs),info.d);
            for i = 1:size(pairs,1)
                Z(i,:) = info.X(pairs(i,1),:)-info.X(pairs(i,2),:);
            end
            for i = 1:length(unique_pairs)
                XX(i,:) = info.X(unique_pairs(i),:);
            end
            Z = bsxfun(@rdivide,Z,sqrt(dot(Z,Z,2)));
            Z = sparse(Z);
            model = train(ones(size(Z,1),1),Z,'-s 0 -c 1e6 -e 1e-6 -q');
            wh = model.w';
            C1 = info.X*wh;
            C2 = zeros(info.nx,1);
%             for i = 1:info.nx
%                 temp = info.X(i,:)*ones(nq+1,1)-XX;
%                 C2(i) = min(diag(temp*temp'));
%             end
            C = C1 + C2;
            a = 1:info.nx;
            C(unique(pairs))=[];a(unique(pairs))=[];
            
            options = [pairs(end,1),a(C==max(C))];

            u = (info.X(options(1),:)-info.X(options(2),:))*w;
            if u>0
                pairs = [pairs;options(1),options(2)];
            elseif u<0
                pairs = [pairs;options(2),options(1)];
            else
                fuck = 1;
            end

            if exist('E_','var')
                [E__,wh_] = EGO(info,pairs,w,bestID,nq+1);
                wh = wh + wh_;
                if length(E__) == length(E_)
                    E_ = E_ + E__;
                    weight = weight + ones(1,length(E_));
                else
                    if length(E_)<length(E__)
                        weight = [weight,zeros(1,length(E__)-length(E_))]+ones(1,length(E__));
                        E_ = [E_,zeros(1,length(E__)-length(E_))]+E__;
                    else
                        weight = [ones(1,length(E__)),zeros(1,length(E_)-length(E__))]+weight;
                        E_ = [E__,zeros(1,length(E_)-length(E__))]+E_;
                    end
                end
            else
                [E_,wh,pairs] = EGO(info,pairs,w,bestID,nq+1);
                weight = ones(1,length(E_));
            end
            E_ = E_./weight;
            E = [E,E_];
            wh = wh/length(options);
        end
    elseif nq == 0
        [whatever, id] = sort(probability_obj,'descend');
        u = (info.X(id(1),:)-info.X(id(2),:))*w;
        if u>0
            pairs = [id(1),id(2)];
        elseif u<0
            pairs = [id(2),id(1)];
        else
            fuck = 1;
        end
        [E_,wh,pairs] = EGO(info,pairs,w,bestID,nq+1);
        E = [E,E_];
    end 