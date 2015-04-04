load basedata.mat;

%% preprocessing
dXunique = unique(dX, 'rows'); %keep unique rows
A = bsxfun(@times, dXunique, (abs(dXunique*w')>1e-6).*sign(dXunique*w')); %flip the pair if the second has higher utility
A(sum(A.^2,2)<1e-6,:) = []; %remove all-zero rows

%% get w distribution
% A = bsxfun(@times,dX,2*(dX*w'>0)-1);
% A(abs(dX*w')<1e-3,:) = [];
addpath('C:\Users\Max Yi Ren\Documents\Research\Code\Tools\liblinear\matlab');
model = train([ones(size(A,1),1);-ones(size(A,1),1)],...
                    sparse([A;-A]),'-s 0 -c 1e6 -e 1e-6 -q');
w_all = model.w;



%% sample according to the distribution of w
%sample with c implementation
sample_target_c = mvtrcnml(1e5,30,A',zeros(size(A,1),1),w);

% %sample with matlab implementation
% sample_target_matlab = mvtnconstrained(2,30,[],[],A,zeros(size(A,1),1),w);

% % find duplicates in [A;-A]
% AA = [A;-A];
% [AAu,IC,IA] = unique(AA,'rows');
% duplicate_ID = setdiff(1:size(AA,1),IC); %should be empty
% % --debugged, sign(dXunique*w') can be zero

s = 1e5;
W = sample_target_c(:,(1e5+1):end)';
obj = bsxfun(@times,1./(1+exp(-W*Xf')),c');
f = sum(bsxfun(@eq, obj, max(obj,[],2))/s)';

