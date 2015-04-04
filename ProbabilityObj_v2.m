% conditional probability of obj i under given query
function p = ProbabilityObj_v2(info,A_set,ii)
    o = info.unique_label(ii,:);
    q = zeros(info.r-1,1);
    q_id = zeros(info.r-1,1);
    for i = 1:info.r-1
        if o(i)>o(i+1)
            q_id(i) = info.Z_id(o(i+1),o(i));
            q(i) = -1;
        else
            q_id(i) = info.Z_id(o(i),o(i+1));
            q(i) = 1;
        end
    end

    a = 1:info.nx;
    a(o)=[];
    last = o(end);
%         pair = [ones(info.nx-info.r,1)*last,a'];
    query = zeros(info.nx-info.r,1);
    query_id = zeros(info.nx-info.r,1);
    for i = 1:info.nx-info.r
        if a(i)>last
            query_id(i) = info.Z_id(last,a(i));
            query(i) = 1;
        else
            query_id(i) = info.Z_id(a(i),last);
            query(i) = -1;
        end
    end
    temp_set = orderingID(info,[q;query],[q_id;query_id]);
    try
    temp_set = (temp_set*ones(1,size(A_set,2))) & A_set;
    p = sum((info.beta*ones(1,size(A_set,2))).*temp_set);
    catch
        fuck = 1;
    end
    
