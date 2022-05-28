function [x] = gradientDescent(f, X, t, threshold)
% gradient descent.
prex = X;
x = prex - t.*f(prex);
while norm(prex - x, "fro") / norm(prex, "fro") >= threshold
    prex = x;
    x = prex - t.*f(prex);
    % disp(norm(prex - x, "fro") / norm(prex, "fro"));
end
end