function f = u(x,y,z,t)
% Exact solution u(x,y,t) used for testing

% Constants used to select different test functions
EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLY = 4;
EXPONENT_0 = 5;
global test_solution

if test_solution == TRIG
    f = 10 * cos(16*x^2 + 4*y^2 + z^2 + t);
elseif test_solution == EXPONENT_1
    f = exp(x + y + z + t);
elseif test_solution == EXPONENT_2
    f = exp(3*x + 2*y + 4*z + t);
elseif test_solution == POLY
    f = 40*x^3 - 60*x^2*y + 20*y^3 + - x*y*z + 10*t;
elseif test_solution == EXPONENT_0
    f = x * (x-1) * y * (y-1) * z * (z-1) * exp(x + y + z + t);
end

end

