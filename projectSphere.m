function x = projectSphere(X)
    x = zeros(size(X));
    for i = 1:size(x,1)
        x(i,:) = X(i,:)/norm(X(i,:));
    end