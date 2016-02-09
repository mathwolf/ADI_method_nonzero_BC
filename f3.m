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
    w = 40 * (256*x^2 + 16*y^2 + z^2) * cos(16*x^2 + 4*y^2 + z^2 + t) + ...
        410 * sin(16*x^2 + 4*y^2 + z^2 + t);
elseif test_solution == EXPONENT_0
    w = 2 * x * y * z * (5 + y*(z-3) - 3*z + x*(y*z+y+z-3)) * ...
        exp(x + y + z + t);
elseif test_solution == EXPONENT_1
    w = - 2 * exp(x+y+z+t);
elseif test_solution == EXPONENT_2
    w = - 28 * exp(4*x + 3*y + 2*z + t);
elseif test_solution == POLY
    w = 10 - 240*x;
end

end

