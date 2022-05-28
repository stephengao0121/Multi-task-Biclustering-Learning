function [result] = mtcost(Y, X, W)
% calculates sum(|y - xw|_2^2)
[~, l] = size(W);
result = 0;
for i = 1: l
    result = result + norm(Y{i} - X{i}*W(:, i), 2)^2;
end
end