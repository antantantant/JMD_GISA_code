% conditional probability of obj i under given query
function p = ProbabilityObj(info,Z,Z_id,A_set,ii)
    if sum(A_set)>0
        o = info.unique_label(ii,:);
        Z_o = info.Zobject{ii};
        a = 1:info.nx;
        a(o)=[];
        last = o(end);
        Z_obj = zeros(info.nx-info.r,info.d);
        Z_obj(a<last,:) = -Z(Z_id(a(a<last),last),:);
        Z_obj(a>last,:) = Z(Z_id(last,a(a>last)),:);
        temp_set = sum(reshape(info.sample*[Z_obj;Z_o]'>0,size(info.sample,1),info.nx-1),2)==info.nx-1;
        temp_set = temp_set & A_set;
        p = sum(probability_w(info.sample_label_id(logical(temp_set),:)',info.theta))*info.sample_c;
    else
        p = 0;
    end