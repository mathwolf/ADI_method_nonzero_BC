function w = f3(x,y,z,t)
% Forcing term for test function

% Constants used to select different test functions
EXPONENT_0 = 0;
EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLY = 4;
global test_solution

if test_solution == TRIG
    w = (10240*x^2 + 640*y^2) * cos(16*x^2 + 4*y^2 + t) + ...
        390 * sin(16*x^2 + 4*y^2 + t);
elseif test_solution == EXPONENT_0
    w = 2 * x * y * z * (5 + y*(z-3) - 3*z + x*(y*z+y+z-3)) * ...
        exp(x + y + z + t);
elseif test_solution == EXPONENT_1
    w = - exp(x+y+t);
elseif test_solution == EXPONENT_2
    w = - 12 * exp(3*x + 2*y + t);
elseif test_solution == POLY
    w = 10 - 240*x;
end

end

