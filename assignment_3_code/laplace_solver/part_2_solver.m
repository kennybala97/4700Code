function [V,E_x,E_y,J_x,J_y,current_end,current_start] = part_2_solver(nx,ny,delta,a_cond,b_cond,box_w,box_gap,V0)
    x_c = nx/2;
    y_c = ny/2;

    [X,Y] = meshgrid(1:1:nx,1:1:ny);
    X = X';
    Y = Y';

    c_map = ones(nx,ny)*b_cond;

    box_top_idx = (X > x_c - box_w/2) & (X < x_c + box_w/2) & (Y > y_c + box_gap/2);
    box_bottom_idx = (X > x_c - box_w/2) & (X < x_c + box_w/2) & (Y < y_c - box_gap/2);

    c_map(box_top_idx | box_bottom_idx) = a_cond;

    [V,E_x,E_y,J_x,J_y] = laplace_solver_2d_part_2(nx,ny,delta,c_map,V0);
    
    J_mag = J_x.^2 + J_y.^2;
    current_end = sum(J_mag(end,:));
    current_start = sum(J_mag(1,:));
end