rng(0);
cor_step = 1:-0.1:-1;
wn = w/norm(w);
price = z(26:30);
cv = 3;
c = Xf(:,26:30)*price'-cv;
bar = zeros(length(cor_step)-1,1);
obj = 1./(1+exp(-wn*Xf')).*c';
true_best = find(obj==max(obj));

size = 1e4;
bin = zeros(size,length(cor_step)-1);
while sum(bar<size) > 0
    disp('generating random partworth...')
    
%     W = bsxfun(@plus, randn(min(size,sum(sum(bin==0))),30)*0.2, wn);
%     W = bsxfun(@rdivide, W, sqrt(dot(W,W,2)));
    cor = W*wn';
    obj = bsxfun(@times,1./(1+exp(-W*Xf')),c');
    opt_id = bsxfun(@eq, obj, max(obj,[],2));
    opt_id(:,true_best) = 0;
    s = sum(opt_id,2);
    for i = 1:(length(cor_step)-1)
        sprintf('*');
        if bar(i)<size
            id = (cor>cor_step(i+1))&(cor<cor_step(i));
            id = find(id==1);
            l = min((bar(i)+length(id)),size)-bar(i);
            bin(bar(i)+(1:l),i) = s(id(1:l));
            bar(i) = bar(i) + l;
        end
    end
end