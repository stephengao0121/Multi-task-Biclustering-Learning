function [U] = COBRA(X)
% DLPA in matrix scope.
[r, c] = size(X);
U = X;
P = zeros(r, c);
Q = zeros(r, c);
preU = zeros(r, c);
while norm(U - preU, "fro") < 10
    preU = U;
    preP = P;
    preQ = Q;
    Y = prox(U.T + P.T);
    P = preU + preP - Y.T;
    U = prox(Y.T + preQ.T);
    Q = Y + preQ - U.T;
end
end