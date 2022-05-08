function [P,sgP] = LinApproximator(y,r,funcs)
N = size(y, 2);
% K = size(r, 1);
M  = size(funcs, 1);
% y is array of N experimantal points: [y1, y2, ..., yN]
% r is array of N K-dimensional vectors: [(x_11; x_21; x_31; ...; x_k1), ...
% , (x_1N; x_2N; x_3N; ...; x_kN)]
% funcs is array of M K-dimensional functions: [(f_11; f_21; f_31; ...; f_k1), ...
% , (f_1M; f_2M; f_3M; ...; f_kM)]

%intermediate matrix with N rows and M columns
g = zeros(N, M);
for ii = 1 : N
    for jj = 1 : M
        f = cell2mat(funcs(jj));
        vec = num2cell(r(:, ii));
        g(ii, jj) =  f(vec{:});
    end
end

Y = y*g;
Y = Y';

G = g'*g;

a = G\Y;
P = a;

Err_0 = sqrsum(y - g*a)/N;

Gkk = trace(G);
Delta_a = sqrt(Err_0/(2*Gkk));
sgP = Delta_a;

end

function[sum] = sqrsum(x)
    len = length(x);
    sum = 0;
    for kk = 1 : len
        sum = sum + x(kk)^2;
    end
end
