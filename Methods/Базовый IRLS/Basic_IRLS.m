function [x, iter_n] = Basic_IRLS(A, b, p, KK)
% KK - количество итераций, A - начальная матрица, b - значения
% приближающей функции.
if nargin < 4, KK =10; end;
x = pinv(A)*b;
E = [];
iter_n = 0;

for k = 1:KK
    iter_n = iter_n + 1;
    e = A*x - b;
    w = abs(e).^((p-2)/2);
    W = diag(w/sum(w));
    WA = W*A;
    x = (WA'*WA)\(WA'*W)*b;
    ee = norm(e, p); E = [E; ee];

end