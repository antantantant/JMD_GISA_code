A = rand(10,2);
[w0, C] = crossvalidation(sparse(A),ones(size(A,1),1)); 
As = bsxfun(@times, A, sqrt(exp(A*w0'))./(1+exp(A*w0')));
As(isnan(As))=0;
Hessian = (eye(length(w0))/C+(As'*As));
Sigma = inv(Hessian);
[X,Y] = meshgrid(-2:0.1:2, -2:0.1:2);
Z1 = zeros(size(X));
Z2 = zeros(size(X));
for i = 1:size(X,1)
    for j = 1:size(X,2)
        w = [X(i,j)+w0(1),Y(i,j)+w0(2)];
        Z1(i,j) = mvnpdf(w, w0, Sigma);
        dW = norm(w)^2;
        Z2(i,j) = exp(-dW/2/C).*prod(1./(1+exp(-w*A')),2);
    end
end
figure;surf(X,Y,Z1);xlabel('x');ylabel('y');
figure;surf(X,Y,Z2);xlabel('x');ylabel('y');