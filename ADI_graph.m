function ADI_graph()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N= 17;
h = 0.125;
x_min = -1.;
y_min = -1.;

stencil = make_stencil(h, N, x_min, y_min);


% Set up matrix of gridpoints to represent domain of problem
for i = 1:N
   % Go along the columns starting at the left
   for j = 1:N
       if stencil(i,j) == 0
           % Not an interior point
           interior_pts(i,j).on = 0; % FALSE
       else
           % Interior point
           interior_pts(i,j).on = 1; % TRUE
           interior_pts(i,j).x = x_min + h * (i-1);
           interior_pts(i,j).y = y_min + h * (j-1);
           interior_pts(i,j).U = 0.;
           interior_pts(i,j).V = 0.;
       end
   end
end

% Set up row structure for matrix equations
n_rows = 0;
for j = 1:N
    % Go along row until we reach the end or find an interior point
    i = 1;
    while (i <= N) && (interior_pts(i,j).on == 0) % FALSE
       i = i + 1; 
    end

    if i > N
       % We are at the end of an empty row, go to the next row
       continue;
    end
    
    % We are at the first interior point in this row
    n_rows = n_rows + 1;
    row(n_rows).j = j;
    row(n_rows).starti = i;
    while (i <= N) && (interior_pts(i,j).on == 1) % TRUE
       i = i + 1; 
    end
    row(n_rows).endi = i - 1;
    
    % Add boundary value terms to each nonempty row
    % Use 0 BC for now
    row(n_rows).btf = 0.;
    row(n_rows).btl = 0.;
end

% Set up col structure for matrix equations
n_cols = 0;
for i = 1:N
    % Go through col until we reach the end or find an interior point
    j = 1;
    while (j <= N) && (interior_pts(i,j).on == 0) % FALSE
       j = j + 1; 
    end

    if j > N
       % We are at the end of an empty col, go to the next col
       continue;
    end
    
    % We are at the first interior point in this col
    n_cols = n_cols + 1;
    col(n_cols).i = i;
    col(n_cols).startj = j;
    while (j <= N) && (interior_pts(i,j).on == 1) % TRUE
       j = j + 1; 
    end
    col(n_cols).endj = j - 1;
    
    % Add boundary value terms to each nonempty col
    % Use 0 BC for now
    col(n_cols).btf = 0.;
    col(n_cols).btl = 0.;
end

disp(n_rows);
disp(n_cols);
%{
for i = 1:N
    for j = 1:N
        if interior_pts(i,j).on == 1
            disp(interior_pts(i,j));
        end
    end
end
%}

for j=1:n_rows
   disp(row(j)); 
end

for i=1:n_cols
   disp(col(i)); 
end

end
