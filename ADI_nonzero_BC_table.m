function ADI_nonzero_BC_table()
% Summary of this function goes here
%   Detailed explanation goes here

% Table for storing error data
table_data = zeros(5,5);

% Check five different grid sizes. Each step will decrease the size by 
% half.

for p = 1:5
    
    % Use a spatial grid of 0.1 times 2 to the power p-1
    N = 10 * 2^(p-1) + 1;
    h = 1./(N-1);
    
    % Use a temporal spacing of the same size, go to 1 second
    tau = h;
    M = floor(1./tau);

    % Permittivity parameter for material
    a = 1.;

    % Set up matrix of gridpoints
    % Use g1 for initial conditions
    U = sparse(N-2,N-2);
    for i = 1:N-2
        for j = 1:N-2
            U(i,j) = g1(h*i, h*j); 
        end
    end


    % Matrices used in ADI method
    Bh = (1/h^2) * spdiags([-ones(N-2,1), 2*ones(N-2,1), -ones(N-2,1)], ...
                        [-1,0,1], N-2, N-2);
    R = chol(speye(N-2) + (a*tau/2)*Bh);

    % Matrix used to store alternating step
    U_twiddles = sparse(N-2,N-2);

    % Array for forcing term
    load_vector = sparse(N-2,1);

    % Step through the scheme
    for m = 1:M
        % Create an array holding value of the forcing term, used on 
        % next step
        for i = 1:N-2
            for j = 1:N-2
                load_vector(i,j) = (1/2)*(f(h*i, h*j, tau*(m-1)) + ...
                    f(h*i, h*j, tau*m)); 
            end
        end
    
        % Create the next U_twiddles
        for i = 1:N-2  
            if i == 1
                b = (1 - 2*(a*tau)/(2*h^2)) * U(:,i) + ...
                    (a*tau)/(2*h^2) * U(:,i+1) + ...
                    (tau/2)*load_vector(:,i);
            elseif i == N-2
                b = (a*tau)/(2*h^2) * U(:,i-1) + ...
                    (1 - 2*(a*tau)/(2*h^2)) * U(:,i) + ...
                    (tau/2)*load_vector(:,i);
            else
                b = (a*tau)/(2*h^2) * U(:,i-1) + ...
                    (1 - 2*(a*tau)/(2*h^2)) * U(:,i) + ...
                    (a*tau)/(2*h^2) * U(:,i+1) + ...
                    (tau/2)*load_vector(:,i);
            end
            b = R'\b;
            U_twiddles(:,i) = R\b;
        end
        
        % Now find next U
        for i = 1:N-2
            b = 2*U_twiddles(i,:)' - (speye(N-2) - ...
                (a*tau/2)*Bh) *  U(i,:)';
        b = R'\b;
        U(i,:) = (R\b)';       
        end    
    end
    
    % Evaluate error at final timestep
    max_nodal_error = 0;
    mean_sq_error = 0;
    nodal_error_squared = 0;
    
    for i = 1:N-2
       for j = 1:N-2
           current_error = abs(U(i,j) - u(h*i, h*j, m*tau));
           if current_error > max_nodal_error
              max_nodal_error = current_error; 
           end
           nodal_error_squared = nodal_error_squared + current_error^2;
       end
    end
    
    mean_sq_error = sqrt(nodal_error_squared / (N-2)^2);
    
    % Write data to table
    table_data(p,1) = h;
    table_data(p,2) = tau;
    table_data(p,3) = max_nodal_error;
    table_data(p,4) = mean_sq_error;
    if p ~= 1
       table_data(p,5) = log2(table_data(p-1,4) / table_data(p,4)); 
    end
end
disp('    h               tau                     max error                       rms errror                            order of conv');
disp(table_data);

end

