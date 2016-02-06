function z = f(x,y,t)
% Forcing term for test function

% Constants used to select different test functions
EXPONENT_0 = 0;
EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLY = 4;
global test_solution

if test_solution == TRIG
    z = (10240*x^2 + 640*y^2) * cos(16*x^2 + 4*y^2 + t) + ...
        390 * sin(16*x^2 + 4*y^2 + t);
elseif test_solution == EXPONENT_0
    z = - x * y * (-7 + 3*y + 3*x + x*y) * exp(x + y + t);
elseif test_solution == EXPONENT_1
    z = - exp(x+y+t);
elseif test_solution == EXPONENT_2
    z = - 12 * exp(3*x + 2*y + t);
elseif test_solution == POLY
    z = 10 - 240*x;
end

end

