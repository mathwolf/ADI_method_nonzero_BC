function Douglas_Rachford_3D_table()

% Matlab R2013a

% ADI method for solution of parabolic PDE with Dirichlet boundary
% conditions in 3D.

% Constants used to select different test functions
EXPONENT_0 = 0;
EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLY = 4;
EXPONENT_3 = 5;

global test_solution
test_solution = EXPONENT_3; 

% Turn on or off perturbation on half step boundary conditions
OFF = 0;
PARTIAL = 1;
FULL = 2;
perturbation = OFF;

% Table for storing error data
table_data = zeros(5,6);

% Check five different grid sizes. Each step will decrease the size by 
% half.

for p = 1:5
    
    % Use a spatial grid of 0.2 times 2 to the power p-1
    N = 5 * 2^(p-1);
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

    % Array used to hold boundary values for the next step
    % Only the boundary points are used, but I am too lazy to figure
    % out how to store this efficiently right now
    U_next = zeros(N+1,N+1,N+1);

    % Create tridiagonal matrix used on left hand side of ADI method
    tridiagonal = spdiags([-ones(N-1,1), 2*ones(N-1,1), -ones(N-1,1)], ...
                            [-1,0,1], N-1, N-1);
    A = speye(N-1) + (tau/(2*h^2)) * tridiagonal;

    % Create tridiagonal matrix used on right hand side of method
    D = speye(N-1) - (tau/(2*h^2)) * tridiagonal;

    % Array for forcing term, used at interior points only
    load_array = zeros(N-1,N-1,N-1);

    % Step through the scheme
    for m = 1:M
        % Create an array for the forcing term at each interior
        % gridpoint at one half of a timestep
        for i = 1:N-1
           for j = 1:N-1
               for k = 1:N-1
                load_array(i,j,k) = (1/2)*(f3(h*i, h*j, h*k, tau*(m-1)) + ...
                  f3(h*i, h*j, h*k, tau*m)); 
               end
           end
        end
        
        % Find only the boundary values of U at
        % the next full timestep.  We will use these values to get the
        % BCs on the subsequent steps of the method.
        for i = 0:N
           for j = 0:N
              U_next(1,i+1,j+1) = g4(0, h*i, h*j, tau*m);
              U_next(N+1, i+1, j+1) =  g4(h*N, h*i, h*j, tau*m);
              U_next(i+1, 1, j+1) = g4(h*i, 0, h*j, tau*m);
              U_next(i+1, N+1, j+1) = g4(h*i, h*N, h*j, tau*m);
              U_next(i+1,j+1,1) = g4(h*i, h*j, 0, tau*m);
              U_next(i+1, j+1, N+1) = g4(h*i, h*j, h*N, tau*m);
           end
        end    

        % Find the approximation for V at the one-third timestep
        % Need to solve system once for each j = 2,...,N
        % and k = 2,...,N
        for j = 2:N
            for k = 2:N
                % Find boundary values for first and last point, no perturbation
                % term
                if perturbation == FULL            
    %                 btf = V(1,j,k) - (tau/(2*h^2)) * ...
    %                     (V(1,j,k-1) - 2 * V(1,j,k) + V(1,j,k+1)) ...
    %                     - (tau/(2*h^2)) * ...
    %                     (V(1,j-1,k) - 2 * V(1,j,k) + V(1,j+1,k)) ...
    %                     +(tau^2/(4*h^4)) * ...
    %                     (V(1,j-1,k-1) - 2 * V(1,j,k-1) + V(1,j+1,k-1) ...
    %                     - 2 * V(1,j-1,k) + 4 * V(1,j,k) - 2 * V(1,j+1,k) ...
    %                     + V(1,j-1,k+1) - 2 * V(1,j,k+1) + V(1,j+1,k+1));
    %                 btl = V(N+1,j,k) - (tau/(2*h^2)) * ...
    %                     (V(N+1,j,k-1) - 2 * V(N+1,j,k) + V(N+1,j,k+1)) ...
    %                     - (tau/(2*h^2)) * ...
    %                     (V(N+1,j-1,k) - 2 * V(N+1,j,k) + V(N+1,j+1,k)) ...
    %                     +(tau^2/(4*h^4)) * ...
    %                     (V(N+1,j-1,k-1) - 2 * V(N+1,j,k-1) + V(N+1,j+1,k-1) ...
    %                     - 2 * V(N+1,j-1,k) + 4 * V(N+1,j,k) - 2 * V(N+1,j+1,k) ...
    %                     + V(N+1,j-1,k+1) - 2 * V(N+1,j,k+1) + V(N+1,j+1,k+1));
    %             elseif perturbation == PARTIAL
    %                 btf = V(1,j,k) - (tau/(2*h^2)) * ...
    %                     (V(1,j,k-1) - 2 * V(1,j,k) + V(1,j,k+1));
    %                 btl = V(N+1,j,k) - (tau/(2*h^2)) * ...
    %                     (V(N+1,j,k-1) - 2 * V(N+1,j,k) + V(N+1,j,k+1));
                else
                    btf = U_next(1,j,k) + U(1,j,k);
                    btl = U_next(N+1,j,k) + U(N+1,j,k);

                end

                % Find rhs of vector equation
                % First, create a vector describing the three terms that 
                % contain a difference quotients.
                b = D*U(2:N,j,k);
                for i = 2:N
                    % Add effects of y derivative
                    b(i-1) = b(i-1) + (tau / (h^2)) * ...
                        (U(i,j-1,k) - 2*U(i,j,k) + U(i,j+1,k));
                    % Add effects of z derviative
                    b(i-1) = b(i-1) + (tau / (h^2)) * ...
                        (U(i,j,k-1) - 2*U(i,j,k) + U(i,j,k+1));
                end
                % Add the effect of the forcing term
                b = b + tau * load_array(1:N-1,j-1,k-1);
                % Add effect of the boundary points
                b(1) = b(1) + (tau/(2*h^2)) * btf;
                b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
                % Solve for approximate values of V at one-third timestep
                U_next(2:N,j,k) = A\b;
            end
        end

        % Store plot data for debugging
    %     for i = 1:N+1
    %         for j = 1:N+1
    %             plot_data(i,j,1) = U(i,j,3);
    %             plot_data(i,j,2) = u3(h*(i-1),h*(j-1),0.5,tau*(m-(1./2.)));
    %             plot_data(i,j,3) = plot_data(i,j,1) - plot_data(i,j,2);
    %         end
    %     end

        % Find the approximation for V at the two-thirds timestep
        % Need to solve system once for each i = 2,...,N
        % and k = 2,...,N
        for i = 2:N
            for k = 2:N
                % Find boundary values for first and last point, no perturbation
                % term
                if perturbation == FULL            
    %                 btf = V(i,1,k) - (tau/(2*h^2)) * ...
    %                     (V(i,1,k-1) - 2 * V(i,1,k) + V(i,1,k+1));
    %                 btl = V(i,N+1,k) - (tau/(2*h^2)) * ...
    %                     (V(i,N+1,k-1) - 2 * V(i,N+1,k) + V(i,N+1,k+1));
    %             elseif perturbation == PARTIAL
    %                 btf = V(i,1,k) - (tau/(2*h^2)) * ...
    %                     (V(i,1,k-1) - 2 * V(i,1,k) + V(i,1,k+1));
    %                 btl = V(i,N+1,k) - (tau/(2*h^2)) * ...
    %                     (V(i,N+1,k-1) - 2 * V(i,N+1,k) + V(i,N+1,k+1));
                else
                    btf = U_next(i,1,k) - U(i,1,k);
                    btl = U_next(i,N+1,k) - U(i,N+1,k);
                end

                % Find rhs of vector equation
                b = U_next(i,2:N,k)';
                % Add effect of y derivative
                b = b + (tau/(2*h^2)) * tridiagonal * U(i,2:N,k)';
                % Add effect of the boundary points
                b(1) = b(1) + (tau/(2*h^2)) * btf;
                b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
                % Solve for approximate values of V at two-thirds timestep
                U_next(i,2:N,k) = (A\b)';
            end
        end

        % Store plot data for debugging
    %     for i = 1:N+1
    %         for j = 1:N+1
    %             plot_data(i,j,4) = U(i,j,3);
    %             plot_data(i,j,5) = u3(h*(i-1),h*(j-1),0.5,tau*(m-(1./3.)));
    %             plot_data(i,j,6) = plot_data(i,j,4) - plot_data(i,j,5);
    %         end
    %     end

        % Now find the approximation for U at next whole timestep
        % Solve system one for each fixed i and j
        for i = 2:N
            for j = 2:N
                btf = U_next(i,j,1) - U(i,j,1);
                btl = U_next(i,j,N+1) - U(i,j,N+1);

                % Get the rhs of the vector equation
                b = squeeze(U_next(i,j,2:N));
                b = b + (tau/(2*h^2)) * tridiagonal * squeeze(U(i,j,2:N));
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
              U(1,i+1,j+1) = U_next(1,i+1,j+1);
              U(N+1, i+1, j+1) = U_next(N+1,i+1,j+1);
              U(i+1, 1, j+1) = U_next(i+1,1,j+1);
              U(i+1, N+1, j+1) = U_next(i+1,N+1,j+1);
              U(i+1,j+1,1) = U_next(i+1,j+1,1);
              U(i+1, j+1, N+1) = U_next(i+1,j+1,N+1);
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

