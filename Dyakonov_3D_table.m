function Dyakonov_3D_table()

% Matlab R2013a

% ADI method for solution of parabolic PDE with Dirichlet boundary
% conditions in 3D.

FALSE = 0;
TRUE = 1;

% Constants used to select different test functions
EXPONENT_0 = 0;
EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLY = 4;
global test_solution
test_solution = EXPONENT_1; % exp test function with 0 BC

% Turn on or off perturbation on half step boundary conditions
perturbation = TRUE;

% Table for storing error data
table_data = zeros(5,6);

% Check five different grid sizes. Each step will decrease the size by 
% half.

for p = 1:5
    
    % Use a spatial grid of 0.1 times 2 to the power p-1
    N = 10 * 2^(p-1);
    h = 1./N;
    
    % Use a temporal spacing of the same size, go to 1 second
    M = N;
    tau = 1./M;

    % Set up matrix of gridpoints
    % Use g1 for initial conditions
    U = zeros(N+1,N+1,N+1);
    for i = 0:N
       for j = 0:N
           for k = 0:N
                U(i+1, j+1, k+1) = g3(h*i, h*j, h*k); 
           end
       end
    end

    % Tridiagonal matrices used for ADI method
    tridiagonal = spdiags([-ones(N-1,1), 2*ones(N-1,1), -ones(N-1,1)], ...
                        [-1,0,1], N-1, N-1);
    A = speye(N-1) + (tau/(2*h^2)) * tridiagonal;
    D = speye(N-1) - (tau/(2*h^2)) * tridiagonal;

    % Array V used on RHS of matrix formulas, interior points only
    V = zeros(N+1,N+1,N+1);

    % Array W used on RHS of matrix formulas, interior points only
    W = zeros(N+1,N+1,N+1);

    % Array for forcing term, used at interior points only
    load_vector = zeros(N+1,N+1,N+1);

    % Step through the scheme
    for m = 1:M
        % Create an array for the forcing term at each interior
        % gridpoint at one half of a timestep
        for i = 1:N+1
           for j = 1:N+1
               for k = 1:N+1
                load_array(i,j,k) = ...
                    (1/2)*(f3(h*(i-1), h*(j-1), h*(k-1), tau*(m-1)) + ...
                    f3(h*(i-1), h*(j-1), h*(k-1), tau*m)); 
               end
           end
        end

        % Create an array of values for W at each interior gridpoint
        % and also at the extreme values of i=0, i=N, j=0, j=N.
        % Fix the row and column and get the z-derivative
        for i = 1:N+1
            for j = 1:N+1
                W(i,j,2:N) = permute(D*squeeze(U(i,j,2:N)), [3 2 1]);
                btf = U(i,j,1);
                btl = U(i,j,N+1);
                W(i,j,2) = W(i,j,2) + (tau/(2*h^2)) * btf;
                W(i,j,N) = W(i,j,N) + (tau/(2*h^2)) * btl;      
            end
        end

        % Create array of values for V at each interior gridpoint
        % and also the extreme values of i=0, i=N.
        % Fix the row and page and get the y-derivative
        for i = 1:N+1
            for k = 2:N
                V(i,2:N,k) = (D*W(i,2:N,k)')';
                btf = W(i,1,k);
                btl = W(i,N+1,k);
                V(i,2,k) = V(i,2,k) + (tau/(2*h^2)) * btf;
                V(i,N,k) = V(i,N,k) + (tau/(2*h^2)) * btl;      
            end
        end

        % Find the approximation at the one-third timestep
        % Need to solve system once for each j = 2,...,N
        % and k = 2,...,N
        for j = 2:N
            for k = 2:N
                % Find boundary values for first and last point, no perturbation
                % term
                if perturbation == TRUE
                    btf = g4(0, h*(j-1), h*(k-1), tau*m) ...
                        - (tau/(2*h^2)) * ...
                        (g4(0, h*(j-1), h*(k-2), tau*m) ...
                        - 2 * g4(0, h*(j-1), h*(k-1), tau*m) ...
                        + g4(0, h*(j-1), h*k, tau*m) ...
                        + g4(0,h*(j-2), h*(k-1), tau*m) ...
                        - 2 * g4(0, h*(j-1), h*(k-1), tau*m) ...
                        + g4(0, h*j, h*(k-1), tau*m))...
                        + (tau^2/(4*h^4)) ...
                        * (g4(0, h*(j-2), h*(k-2), tau*m) ...
                        + g4(0, h*(j-2), h*k, tau*m) ...
                        + g4(0, h*j, h*(k-2), tau*m) ...
                        + g4(0, h*j, h*k, tau*m) ...
                        -2 * g4(0, h*(j-1), h*(k-2), tau*m) ...
                        -2 * g4(0, h*(j-1), h*k, tau*m) ...
                        -2 * g4(0, h*(j-2), h*(k-1), tau*m) ...
                        -2 * g4(0, h*j, h*(k-1), tau*m) ...
                        + 4 * g4(0, h*(j-1), h*(k-1), tau*m)) ...
                        + V(1,j,k);
                    btl = g4(1, h*(j-1), h*(k-1), tau*m) ...
                        - (tau/(2*h^2)) * ...
                        (g4(1, h*(j-1), h*(k-2), tau*m) ...
                        - 2 * g4(1, h*(j-1), h*(k-1), tau*m) ...
                        + g4(1, h*(j-1), h*k, tau*m) ...
                        + g4(1,h*(j-2), h*(k-1), tau*m) ...
                        - 2 * g4(1, h*(j-1), h*(k-1), tau*m) ...
                        + g4(1, h*j, h*(k-1), tau*m))...
                        + (tau^2/(4*h^4)) ...
                        * (g4(1, h*(j-2), h*(k-2), tau*m) ...
                        + g4(1, h*(j-2), h*k, tau*m) ...
                        + g4(1, h*j, h*(k-2), tau*m) ...
                        + g4(1, h*j, h*k, tau*m) ...
                        -2 * g4(1, h*(j-1), h*(k-2), tau*m) ...
                        -2 * g4(1, h*(j-1), h*k, tau*m) ...
                        -2 * g4(1, h*(j-2), h*(k-1), tau*m) ...
                        -2 * g4(1, h*j, h*(k-1), tau*m) ...
                        + 4 * g4(1, h*(j-1), h*(k-1), tau*m)) ...
                        + V(N+1, j, k);

                else
                    btf = g4(0, h*(j-1), h*(k-1), tau*m) ...
                        + V(1,j,k);
                    btl = g4(1, h*(j-1), h*(k-1), tau*m) ...
                        + V(N+1, j, k);
                end

                % Find rhs of vector equation
                b = D*V(2:N,j,k) + tau * load_array(2:N,j,k);
                % Add effect of the boundary points
                b(1) = b(1) + (tau/(2*h^2)) * btf;
                b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
                % Solve for approximate values at one-third timestep
                U(2:N,j,k) = A\b;
            end
        end

        % Find the approximation at the two-thirds timestep
        % Need to solve system once for each i = 2,...,N
        % and k = 2,...,N
        for i = 2:N
            for k = 2:N
                % Find boundary values for first and last point, no perturbation
                % term
                if perturbation == TRUE            
                    btf = g4(h*(i-1), 0, h*(k-1), tau*m) ...
                        - (tau/(2*h^2)) * ...
                        (g4(h*(i-1), 0, h*(k-2), tau*m) ...
                        - 2 * g4(h*(i-1), 0, h*(k-1), tau*m) ...
                        + g4(h*(i-1), 0, h*k, tau*m));
                    btl = g4(h*(i-1), h*N, h*(k-1), tau*m) ...
                        - (tau/(2*h^2)) * ...
                        (g4(h*(i-1), h*N, h*(k-2), tau*m) ...
                        - 2 * g4(h*(i-1), h*N, h*(k-1), tau*m) ...
                        + g4(h*(i-1), h*N, h*k, tau*m));
                else
                    btf = g4(h*(i-1), 0, h*(k-1), tau*m);
                    btl = g4(h*(i-1), 1, h*(k-1), tau*m);
                end

                % Find rhs of vector equation
                b = U(i,2:N,k)';
                % Add effect of the boundary points
                b(1) = b(1) + (tau/(2*h^2)) * btf;
                b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
                % Solve for approximate values at two-thirds timestep
                U(i,2:N,k) = (A\b)';
            end
        end

        % Now find the interior points for next whole timestep
        % Solve system one for each fixed i and j
        for i = 2:N
            for j = 2:N

                btf = g4(h*(i-1), h*(j-1), 0, tau*m);
                btl = g4(h*(i-1), h*(j-1), h*N, tau*m);
                % Get the rhs of the vector equation
                b = squeeze(U(i,j,2:N));
                % Add effect of the boundary points
                b(1) = b(1) + (tau/(2*h^2)) * btf;
                b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
                % Solve vector equation for the next full timestep
                U(i,j,2:N) = permute(A\b, [3 2 1]);       
            end
        end

        % Update boundary points
        for i = 0:N
           for j = 0:N
              U(1,i+1,j+1) = g4(0, h*i, h*j, tau*m);
              U(N+1, i+1, j+1) = g4(h*N, h*i, h*j, tau*m);
              U(i+1, 1, j+1) = g4(h*i, 0, h*j, tau*m);
              U(i+1, N+1, j+1) = g4(h*i, h*N, h*j, tau*m);
              U(i+1,j+1,1) = g4(h*i, h*j, 0, tau*m);
              U(i+1, j+1, N+1) = g4(h*i, h*j, h*N, tau*m);
           end
        end
    end
    
    % Evaluate error at interior gridpoints at final timestep
    max_nodal_error = 0;
    nodal_error_squared = 0;    
    root_mean_sq_error = 0;

    for i = 1:N-1
       for j = 1:N-1
           for k = 1:N-1
               current_error = abs(U(i+1,j+1,k+1) - u3(h*i, h*j, h*k, 1));
               if current_error > max_nodal_error
                  max_nodal_error = current_error; 
               end
               nodal_error_squared = nodal_error_squared + current_error^2;
           end
       end
    end
    root_mean_sq_error = sqrt(h^3 * nodal_error_squared);

    
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

