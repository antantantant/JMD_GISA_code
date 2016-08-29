rng(0);
level = 5;
step = -2:0.01:2;
% wn = w/norm(w);
wn = w;
price = z(26:30);
cv = 3;
c = Xf(:,26:30)*price'-cv;
obj = 1./(1+exp(-wn*Xf')).*c';
true_best = find(obj==max(obj));

bin = zeros(30,length(step));
for j = 1:30
    groupid = floor((j-1)/level);
    groupsum = sum(wn(groupid*level+(1:level)))-wn(j);
    for i = 1:length(step)
%         if j == 9 && i == 198
%             a = 1;
%         end
        wb = wn;
%         wb(groupid*level+(1:level)) =...
%             wn(groupid*level+(1:level))...
%             /groupsum * (sum(wn(groupid*level+(1:level)))-step(i));
        wb(j) = step(i);
        obj = 1./(1+exp(-wb*Xf')).*c';
        best = find(obj==max(obj));
        if length(best)==1 && best == true_best
            bin(j,i) = 1;
        else
            bin(j,i) = 0;
        end
    end
end

figure;
imagesc(bin*0.6);
colormap bone;
hold on;

id = (round(wn*100)+200)+1;
for j = 1:30
    plot([id(j),id(j)],[j-0.5,j+0.5],'r','linewidth',3);
    plot([-0.5,length(step)+.5],[j-0.5,j-0.5],'Color',[0.5,0.5,0.5]);
end
