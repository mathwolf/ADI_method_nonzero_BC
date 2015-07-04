function ADI_graph()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N= 5;
h = 0.5;
x_min = -1.;
y_min = -1.;

[n_pts, stencil] = make_stencil(h, N, x_min, y_min);

% Vector holding gridpoints
% Entries in each cell: i, j, U, V, next_pt

% Array with one entry for each gridpoints
interior_pts = zeros(n_pts,6);
current_pt = 0;
n_rows = 0;
n_cols = 0;
row_ptr = zeros(N,1);
col_ptr = zeros(N,1);

% Look at each row of stencil, start with row 1
for i = 1:N
   % Go along the columns until we get the first 1
   j = 1;
   while (j <= N) && (stencil(i,j) == 0)
      j = j + 1;
   end
   
   if j > N
      % We have reached the end of the current row without finding any 1s
      % If we are at the beginning of the stencil, then go to the next row
      if current_pt == 0
          continue;
      % If we are not at the beginning, an empty row means we have set up 
      % all interior points
      else
          break; 
      end
   end
   
   % We found the first 1 in the current row 
   % Write the next row to the data structure
   % Each 1 in the current row maatches up to one entry in the data vector
   n_rows = n_rows + 1;
   row_ptr(i) = current_pt + 1;   % grid pt containing beginning pt of row
   while stencil(i,j) == 1
       current_pt = current_pt + 1;
       interior_pts(current_pt,1) = i;
       interior_pts(current_pt,2) = j;
       interior_pts(current_pt,3) = x_min + h * (i-1);
       interior_pts(current_pt,4) = y_min + h * (j-1);
       interior_pts(current_pt,5) = 0;  % initial condition for PDE
       j = j + 1;
       if j > N
           % We reached the end of the row without finding a 0
           break;
       end
   end
      
   % Go back through the current row and create the column structure
   % for our grid
   if n_rows == 1
       % This is the first row - every column is a new column
       current_row_end = current_pt;
       current_pt = row_ptr(i) - 1;     % go back to beginning of row
       while current_pt < current_row_end
          current_pt = current_pt + 1;
          n_cols = n_cols + 1;
          current_col = interior_pts(current_pt,2);
          col_ptr(current_col) = current_pt;
       end
   else
      % This is not the first row, so we need to align the columns of the
      % current row with the columns of the previous row
      current_row_end = current_pt;
      current_pt = row_ptr(i);     % go back to beginning of row
      prev_row_pt = row_ptr(i-1);
      prev_row_end = current_pt;
      
      current_col = interior_pts(current_pt, 2);
      % Take care of all points in the current row that have 0 above them
      while (current_col < interior_pts(prev_row_pt, 2))
            % Make a new column
            n_cols = n_cols + 1;
            col_ptr(current_col) = current_pt;
            current_pt = current_pt + 1;
            current_col = current_col + 1;
      end
      
      % Next we have all points with a 1 in the previous row
      while (prev_row_pt <= prev_row_end) && ...
              (current_pt <= current_row_end)
          % Update the pointer for this col in the previous row
          interior_pts(prev_row_pt, 6) = current_pt;
          prev_row_pt = prev_row_pt + 1;
          current_pt = current_pt + 1;
      end
      
      % See if we have any new columns at the end of the current row
      while current_pt <= current_row_end
          % Make a new column
          current_col = interior_pts(current_pt, 2);
          n_cols = n_cols + 1;
          col_ptr(current_col) = current_pt;
          current_pt = current_pt + 1;
      end
      current_pt = current_pt - 1;
   end
   
end

disp(n_rows);
disp(n_cols);
for i = 1:current_pt
    disp([i, interior_pts(i,:)]);
end
end

