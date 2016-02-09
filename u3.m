function w = u3(x,y,z,t)
% Exact solution u(x,y,t) used for testing

% Constants used to select different test functions
EXPONENT_0 = 0;
EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLY = 4;
global test_solution

if test_solution == EXPONENT_0
    w = x * (1-x) * y * (1-y) * z * (1-z) * exp(x+y+z+t);
elseif test_solution == TRIG
    w = 10 * cos(16*x^2 + 4*y^2 + z^2 + t);
elseif test_solution == EXPONENT_1
    w = exp(x + y + z + t);
elseif test_solution == EXPONENT_2
    w = exp(4*x + 3*y + 2*z + t);
elseif test_solution == POLY
    w = 40*x^3 - 60*x^2*y + 20*y^3 + 10*t;
end

end

