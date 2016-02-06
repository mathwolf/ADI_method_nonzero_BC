function z = u(x,y,t)
% Exact solution u(x,y,t) used for testing

% Constants used to select different test functions
EXPONENT_0 = 0;
EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLY = 4;
global test_solution

if test_solution == EXPONENT_0
    z = x * (1-x) * y * (1-y) * exp(x+y+t);
elseif test_solution == TRIG
    z = 10 * cos(16*x^2 + 4*y^2 + t);
elseif test_solution == EXPONENT_1
    z = exp(x + y + t);
elseif test_solution == EXPONENT_2
    z = exp(3*x + 2*y + t);
elseif test_solution == POLY
    z = 40*x^3 - 60*x^2*y + 20*y^3 + 10*t;
end

end

