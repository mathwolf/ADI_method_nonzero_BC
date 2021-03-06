function Dyakonov_2D_table()

% Matlab R2013a

% ADI method for solution of parabolic PDE with Dirichlet boundary
% conditions. This program tests the method using the function
% e^(2x+3y+t).  It uses a range of grid sizes and displays
% the associated error measurements in a table.

% Constants used to select different testing parameters
EXPONENT_0 = 0;
EXPONENT_1 = 1;
EXPONENT_2 = 2;
EXPONENT_3 = 3;
TRIG = 4;

global test_solution 
test_solution = TRIG;

OFF = 0;
ON = 1;
PERTURBATION = OFF;

% Table for storing error data
table_data = zeros(5,6);

% Check five different grid sizes. Each step will decrease the size by 
% half.

for p = 1:5
    
    % Use a spatial grid of 0.1 times 2 to the power p-1
    N = 5 * 2^(p-1);
    h = 1./N;
    
    % Use a temporal spacing of the same size, go to 1 second
    M = N;
    tau = 1./M;

    % Set up matrix of gridpoints
    % Use g1 for initial conditions
    U = zeros(N+1,N+1);
    for i = 0:N
        for j = 0:N
            U(i+1,j+1) = g1(h*i, h*j); 
        end
    end

    % Tridiagonal matrices used on LHS of ADI method
    tridiagonal = spdiags([-ones(N-1,1), 2*ones(N-1,1), -ones(N-1,1)], ...
                        [-1,0,1], N-1, N-1);
    A = speye(N-1) + (tau/(2*h^2))*tridiagonal;

    % Matrix V used on RHS of matrix formulas, interior points only
    V = zeros(N-1,N-1);

    % Array for forcing term, used at interior points only
    load_vector = zeros(N-1,N-1);

    % Step through the scheme
    for m = 1:M

        % Create an array of values for V at each interior gridpoint
        % Used on RHS of step 1 of the Dyakonov method
        for i = 1:N-1
            for j = 1:N-1
                % Take derivatives
                V(i,j) = U(i+1,j+1) + ...
                    (tau / (2*h^2)) * (U(i+1,j) - 2*U(i+1,j+1) + U(i+1,j+2)) + ...
                    (tau / (2*h^2)) * (U(i,j+1) - 2*U(i+1,j+1) + U(i+2,j+1)) + ...
                    (tau^2 / (4*h^4)) * (U(i,j) - 2*U(i+1,j) + U(i+2,j) ...
                        - 2*U(i,j+1) + 4*U(i+1,j+1) - 2*U(i+2,j+1) ...
                        + U(i,j+2) - 2*U(i+1,j+2) + U(i+2,j+2));
                % Add the effect of the forcing term
                V(i,j) = V(i,j) + (tau/2)*(f(h*i, h*j, tau*(m-1)) + ...
                    f(h*i, h*j, tau*m));
            end
        end
        
        % Find the approximate values at the half timestep
        % Need to solve system once for each j = 1,...,N-1
        for i = 1:N-1  
            if PERTURBATION == ON
                % Calculate the boundary values for the first and last
                % points with perturbation terms
                btf = g2(h*i, 0, tau*m) ...
                    - (tau/(2*h^2)) * (g2(h*(i-1), 0, tau*m) ...
                    - 2*g2(h*i, 0, tau*m) + g2(h*(i+1), 0, tau*m));
                btl = g2(h*i, 1, tau*m) ...
                    - (tau/(2*h^2)) * (g2(h*(i-1), 1, tau*m) ...
                    - 2*g2(h*i, 1, tau*m) + g2(h*(i+1), 1, tau*m));
            else
                % Calculate the boundary values for the first and last 
                % points, no perturbation terms
                btf = g2(h*i, 0, tau*m);
                btl = g2(h*i, 1, tau*m);
            end
            % Find rhs of vector equation
            b = V(i,1:N-1)';
            % Add effect of boundary values of rhs
            b(1) = b(1) + (tau/(2*h^2)) * btf;
            b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
            % Solve vector equation to get points for half timestep
            U(i+1,2:N) = (A\b)';
        end
    
        % Now find the interior points for next whole timestep
        % Solve system one for each j = 1, ... , N-1
        for j = 1:N-1 
            % Find boundary values for first and last points
            btf = g2(0, h*j, tau*m);
            btl = g2(1, h*j, tau*m);
            % Calculate rhs of next vector equation
            b = U(2:N,j+1);
            % Add effect of boundary values
            b(1) = b(1) + (tau/(2*h^2)) * btf;
            b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
            % Solve vector equation to get points for next full timestep
            U(2:N,j+1) = (A\b)';       
        end
    
        % Update boundary points
        for i = 0:N
            U(1,i+1) = g2(0, h*i, tau*m); 
            U(N+1,i+1) = g2(h*N, h*i, tau*m);
            U(i+1,1) = g2(h*i, 0, tau*m);
            U(i+1,N+1) = g2(h*i, h*N, tau*m);
        end

    end
    
    % Evaluate error at interior gridpoints at final timestep
    max_nodal_error = 0;
    nodal_error_squared = 0;    
    root_mean_sq_error = 0;

    for i = 1:N-1
       for j = 1:N-1
           current_error = abs(U(i+1,j+1) - u(h*i, h*j, 1));
           if current_error > max_nodal_error
              max_nodal_error = current_error; 
           end
           nodal_error_squared = nodal_error_squared + current_error^2;
       end
    end
    root_mean_sq_error = sqrt(h^2 * nodal_error_squared);

    
    % Write data to table
    table_data(p,1) = h;
    table_data(p,2) = max_nodal_error;
    if p ~= 1
       % estimate order of convergence for max err
       table_data(p,3) = log(table_data(p-1,2) / table_data(p,2)) / ...
           log(table_data(p-1,1) / table_data(p,1)); 
    end   
    table_data(p,4) = h;
    table_data(p,5) = root_mean_sq_error;
    if p ~= 1
       % estimate order of convergence for rms err
       table_data(p,6) = log(table_data(p-1,5) / table_data(p,5)) / ...
           log(table_data(p-1,4) / table_data(p,4));
    end
    

end

% Display table of error information
disp('    h                  max error          order of conv');
disp(table_data(:,1:3));
disp('    h                  rms error          order of conv');
disp(table_data(:,4:6));
end

