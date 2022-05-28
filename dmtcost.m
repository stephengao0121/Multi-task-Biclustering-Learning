function [dw] = dmtcost(Y, X, W)
% gradient of mtcost.
[p, l] = size(W);
dw = zeros(p, l);
for i = 1: l
    dw(:, i) = X{i}.' * (X{i}*W(:, i) - Y{i});
end
end