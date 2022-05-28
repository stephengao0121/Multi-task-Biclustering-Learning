function [W] = admm(Y, X, rho1, rho2, threshold, alpha, beta)
% Realization of admm algorithm.
%% Initilization
[~, p] = size(X);
[~, l] = size(Y);
W = zeros(p, l);
U = zeros(p*(p-1)/2, l);
V = zeros(p, l*(l-1)/2);
lambda1 = zeros(p*(p-1)/2, l);
lambda2 = zeros(p, l*(l-1)/2);
R = zeros(p*(p-1)/2, p);
C = zeros(l, l*(l-1)/2);
t = 1e-4;
gdthreshold = 1e-4;
tau1 = alpha / rho1;
tau2 = beta / rho2;
text1 = 'initilization done';
disp(text1);
fvalue = [];
%% Contruct R & C
counter = 1;
for i = 1: p
    for j = i+1: p
        R(counter, i) = 1;
        R(counter, j) = -1;
        counter = counter + 1;
    end
end
counter = 1;
for i = 1: l
    for j = i+1: l
        C(i, counter) = 1;
        C(j, counter) = -1;
        counter = counter + 1;
    end
end
text2 = 'r&c done';
disp(text2);
%% Iteration
counter = 1;
while true
    preW = W;
    preU = U;
    preV = V;
    prel1 = lambda1;
    prel2 = lambda2;
    df = @(W) X.'*(X*W-Y) + rho1.*R.'*(R*W-preU-prel1) + rho2.*(W*C-preV-prel2)*C.';
    W = gradientDescent(df, preW, t, gdthreshold);
    U = proximalL12norm(W.'*R.' - prel1.', tau1).';
    V = proximalL12norm(W*C - prel2, tau2);
    lambda1 = prel1 + rho1.*(U - R*W);
    lambda2 = prel2 + rho2.*(V - W*C);
    f = @(w) 0.5*norm(Y-X*w, "fro")^2 + (rho1/2)*norm(U-R*w+lambda1, "fro")^2 + (rho2/2)*norm(V-w*C+lambda2, "fro")^2;
    fprev = @(w) 0.5*norm(Y-X*w, "fro")^2 + (rho1/2)*norm(preU-R*w+prel1, "fro")^2 + (rho2/2)*norm(preV-w*C+prel2, "fro")^2;
    fvalue(counter) = f(W);
    counter = counter + 1;
    if abs(f(W) - fprev(preW)) / abs(fprev(preW)) < threshold
        break;
    end
    if counter > 500
       break;
    end
end
%% Plotting
plot(fvalue);
end