function w = f3(x,y,z,t)
% Forcing term for test function

% Constants used to select different test functions
EXPONENT_0 = 0;
EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLY = 4;
EXPONENT_3 = 5;
global test_solution

switch test_solution 
    case TRIG
        w = 40 * (256*x^2 + 16*y^2 + z^2) * cos(16*x^2 + 4*y^2 + z^2 + t) + ...
            410 * sin(16*x^2 + 4*y^2 + z^2 + t);
    case EXPONENT_0
        w = 2 * x * y * z * (5 + y*(z-3) - 3*z + x*(y*z+y+z-3)) * ...
            exp(x + y + z + t);
    case EXPONENT_1
        w = - 2 * exp(x+y+z+t);
    case EXPONENT_2
        w = - 28 * exp(4*x + 3*y + 2*z + t);
    case POLY
        w = 10 - 240*x - 40*y^3 + 120*y*z - 120*y*z^2;
    case EXPONENT_3
        w = exp(x*y*z*t)*(x*y*z - t^2 * ( x^2*y^2 + x^2*z^2 + y^2*z^2));
end

end

