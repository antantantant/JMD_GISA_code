function E = pairwiseFullClayton_v2(info,p_id,q_id,q,nx,np)
    [Pr,Pq] = updatePr(info,q_id,q,nx); % calculate current probability for each object
    Pr(Pr==0)=1e-18;
    E = sum(-Pr.*log2(Pr));
    if ~sum(Pr==1)>0
        % pick the next query that minimizes the expectation of entropy
        a = 1:np;a(q_id)=[];
        Pr_pos = zeros(nx,length(a));
        Pr_q = zeros(length(a),1);
        rho_k = zeros(nx,length(a));
        rho = zeros(length(a),1);
        C = rho;
        
        for i = 1:length(a)
            % calculate the probability of a new query returning 1
            if ~isempty(q_id)
                q_q = info.Q;q_q = q_q(q_q(:,a(i))==1,q_id);
                if isempty(q_q)
                    Pr_q(i) = 0;
                else
                    t = q_q - ones(size(q_q,1),1)*q;
                    Pr_q(i) = sum(diag(t*t')==0)/size(q_q,1)*...
                        sum(info.Q(:,a(i))==1)/size(info.Q,1)/sum(Pq);
                end
            else
                Pr_q(i) = sum(info.Q(:,a(i))==1)/size(info.Q,1);
            end
            
            Pr_pos(:,i) = updatePr(info,[q_id,a(i)],[q,1],nx);
            Pr_pos(Pr_pos(:,i)==0,i)=1e-18;
            if Pr_q(i)>0.5
                rho(i) = Pr_q(i);
            else
                rho(i) = 1-Pr_q(i);
            end
            for j = 1:nx
                if Pr_pos(j,i)*Pr_q(i)/Pr(j)>0.5
                    rho_k(j,i) = Pr_pos(j,i)*Pr_q(i)/Pr(j);
                else
                    rho_k(j,i) = 1 - Pr_pos(j,i)*Pr_q(i)/Pr(j);
                end
            end
            C(i) = rho(i)*log2(rho(i))-sum(Pr.*rho_k(:,i).*log2(rho_k(:,i)));
        end
        
        options = a(C==min(C));
        for i = 1:length(options)
            q_id_ = [q_id,options(i)];
            q_ = [q,info.Q(p_id,options(i))];
            if exist('E_','var')
                E__ = pairwiseFullClayton_v2(info,p_id,q_id_,q_,nx,np);
                if length(E__) == length(E_)
                    E_ = E_ + E__;
                    w = w + ones(1,length(E_));
                else
                    if length(E_)<length(E__)
                        w = [w,zeros(1,length(E__)-length(E_))]+ones(1,length(E__));
                        E_ = [E_,zeros(1,length(E__)-length(E_))]+E__;
                    else
                        w = [ones(1,length(E__)),zeros(1,length(E_)-length(E__))]+w;
                        E_ = [E__,zeros(1,length(E_)-length(E__))]+E_;
                    end
                end
            else
                E_ = pairwiseFullClayton_v2(info,p_id,q_id_,q_,nx,np);
                w = ones(1,length(E_));
            end
        end
        try
        E_ = E_./w;
        catch
            fuck = 1;
        end
        E = [E,E_];
    end