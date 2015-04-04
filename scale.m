% clear;
% 
% w = [-0.534,0.129,0.228,0.104,0.052,...
%     -0.058,0.253,0.278,-0.025,-0.467,...
%     0.015,-0.098,0.049,0.047,-0.033,...
%     -0.366,-0.164,0.215,0.194,0.100,...
%     -0.744,-0.198,0.235,0.291,0.396,...
%     0.719,0.482,0.054,-0.368,-0.908];
% % attributes and levels
% z1 = 200:50:400;
% z2 = [0.75,0.88,1.00,1.14,1.33];
% z3 = 100:10:140;
% z4 = 1/16:1/32:6/32;
% z5 = 0.75:0.25:1.75;
% z6 = 10:5:30;
% z = [z1,z2,z3,z4,z5,z6];
% 
% % design parameters
% y1 = 0.3;
% y2 = 0.5;
% y3 = 1.9;
% y4 = 0.25;
% y5 = 3;
% y6 = 2;
% y7 = 1.13;
% y8 = 1;
% y9 = 1.1;
% y10 = 0.31;
% y11 = 16;
% y12 = 1.29;
% y13 = 4;
% y = [y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13];
% % linear constraints
% A = zeros(6,15);
% b = zeros(6,1);
% A(1,[12,14]) = [1,-1];
% b(1) = -2*y1;
% A(2,[7,12,13]) = [1,1,-1];
% b(2) = -2*y1 - y9;
% A(3,[4,5,13]) = [1,1,-1];
% b(3) = -2*y1;
% A(4,[2,5]) = [-1,1];
% A(5,[7,8,11,13]) = [1,1,1,-1];
% b(5) = -2*y1 - y9;
% A(6,[7,8,10,12,13]) = [1,1,1,0.5,-1];
% b(6) = -2*y1 - y7 - y9;
% 
% lb = [0.125*ones(5,1);1;0.5;1;0.25;0.5;0.5;1;1;1;10];
% ub = [36;36;24;24;36;200;12;36;24;1.9;1.9;36;36;36;30];
% 
% % get feasible designs
% X = zeros(5^6,30);
% count = 1;
% for i1 = 1:5
%     for i2 = 1:5
%         for i3 = 1:5
%             for i4 = 1:5
%                 for i5 = 1:5
%                     for i6 = 1:5
%                         X(count,[i1,5+i2,10+i3,15+i4,20+i5,25+i6]) = 1;
%                         count = count + 1;
%                     end
%                 end
%             end
%         end
%     end
% end
% [design_set, feasible] = checkFeasibility(X,z);
% feasible = feasible==1;
% Xf = X(feasible,:);
% 
% % find the true best design
% true_obj = zeros(size(Xf,1),1);
% for i = 1:size(Xf,1)
%     attribute = z(Xf(i,:)==1);
%     true_obj(i) = scaleObj([],y,w,attribute);
% end
% true_best = find(true_obj==min(true_obj));
% 
% % parameter study
% test_range = -1:0.1:1;
% best_set = zeros(length(w),length(test_range));
% df2true = zeros(length(w),length(test_range));
load trueSolution.mat;
for i = 1:length(w);
    disp(['\n parameter study on w ',num2str(i)]);
    x_c = zeros(length(test_range),15);
    w_c = repmat(w,length(test_range),1);
    w_c(:,i) = test_range';
    for j = 1:length(test_range)
        disp('.');
        obj = zeros(size(Xf,1),1);
        for k = 1:size(Xf,1)
            attribute = z(Xf(k,:)==1);
            obj(k) = scaleObj([],y,w_c(j,:),attribute);
        end
        best_set(i,j) = find(obj==min(obj));
        attribute = z(Xf(best_set(i,j),:)==1);
        df2true(i,j) = abs(scaleObj([],y,w,attribute)-fval);
    end
end


% test_range = -1:0.1:1;
% dist2true = zeros(length(w),length(test_range));
% f2true = zeros(length(w),length(test_range));
% for i = 1:length(w);
%     disp(['\n parameter study on w ',num2str(i)]);
%     x_c = zeros(length(test_range),15);
%     w_c = repmat(w,length(test_range),1);
%     w_c(:,i) = test_range';
%     if exist('x0','var')
%         clear x0;
%     end
%     for j = 1:length(test_range)
%         disp('.');
%         if exist('x0','var')
%             x0 = x_c(j-1,:);
%         else
%             x0 = x;
%         end
% %         x_c(j,:) = fmincon(@(x)scaleObj(x,y,w_c(j,:)),x0,A,b,[],[],lb,ub,@(x)scaleCon(x,y));
%         options = gaoptimset('InitialPopulation',x0);
%         [xx,ff,exitflag] = ga(@(x)scaleObj(x',y,w_c(j,:)),15,A,b,[],[],lb,ub,@(x)scaleCon(x',y),options);
%         x_c(j,:) = xx;
%         
%         dist2true(i,j) = norm(x_c(j,:)-x);
%         f2true(i,j) = abs(ff-fval);
%     end
% end
    