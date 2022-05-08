function [P, sgP] = NonLinApproximator (y,r,fun, P_0)
% P_0 is a line-array with starting guess for values of parameters
N = size(y, 2);

% p_number = find_arg_number(fun, N); % number of parameters P_i
p_number = length(P_0);

for iterations = 1 : 1000

F = zeros(1, N);
for ii = 1 : N
    F(ii) = fun(r(:, ii), P_0); 
end

% find jacobian of fun
J = jacobi(fun, r, p_number, P_0);

% equation is G*(P - P_0) = J*(y - fun(r, P_0))
Y = J'*(y' - F');
G = J'*J;
P = P_0' + G\Y;

if isnan(P)
    break
end

if iterations > 1
    if abs(P_0 - P_previous) < 1e-10
        break
    end
end

P_previous = P;
P_0 = P';
end

sgP = 'bla bla';
end

function [number] = find_arg_number(fun, N)
ii = 1;
exception_identifier = 'MATLAB:badsubscript';
while strmatch(exception_identifier, 'MATLAB:badsubscript')
    exception_identifier = ' ';
    try fun(zeros(1, N), zeros(1, ii))
    catch E
        exception_identifier = E.identifier;
        ii = ii + 1;
    end
end

number = ii;
end

function [J] = jacobi(fun, r, p_number, P_0)
% numerical method of finding derivative
h = 1e-6; %step along axis
N = length(r(1, :));

J = zeros(N, p_number);

for ii = 1 : N
    for jj = 1 : p_number
        delta = zeros(1, p_number);
            delta(jj) = h*abs(P_0(jj));
            J(ii, jj) = (fun(r(:, ii),P_0 + delta) - fun(r(:, ii),P_0 - delta))/(2*h);
    end   
end

end
