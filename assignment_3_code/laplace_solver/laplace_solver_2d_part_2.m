
function [Vmap,E_x,E_y,J_x,J_y] = laplace_solver_2d_part_2(nx,ny,delta,c_map,V0)
    G = sparse(nx*ny);
    B = zeros(1,nx*ny);
    
    [co_ords.i,co_ords.j] = get_co_ord_list(nx,ny);    
    co_ords.n = ij_to_n(ny,co_ords.i,co_ords.j);

    for u = 1:1:length(co_ords.n)
        i = co_ords.i(u);
        j = co_ords.j(u);

        [G,B] = set_values(nx,ny,i,j,c_map,G,B,V0);
    end

    V = G\B';
    
    Vmap = reshape(V,ny,nx)';
    
    [E_x,E_y] = gradient(Vmap,delta);
    E_x = -E_x;
    E_y = -E_y;
    
    J_x = c_map.*E_x;
    J_y = c_map.*E_y;
    
end

function [G,B] = set_values(nx,ny,i,j,c_map,G,B,V0)
    l = ij_to_n(ny,i-1,j);
    r = ij_to_n(ny,i+1,j);
    t = ij_to_n(ny,i,j+1);
    b = ij_to_n(ny,i,j-1);
    c = ij_to_n(ny,i,j);
    
    if i == 1
        G(c,:) = 0;
        G(c,c) = 1/V0;
        B(c) = 1;
    elseif i == nx
        G(c,:) = 0;
        G(c,c) = 1;
    elseif j == 1
        s_l = (c_map(i,j) + c_map(i-1,j))/2;
        s_r = (c_map(i,j) + c_map(i+1,j))/2;
        s_t = (c_map(i,j) + c_map(i,j+1))/2;

        G(c,c) = -(s_l+s_r+s_t);
        G(c,l) = s_l;
        G(c,r) = s_r;
        G(c,t) = s_t;
    elseif j == ny
        s_l = (c_map(i,j) + c_map(i-1,j))/2;
        s_r = (c_map(i,j) + c_map(i+1,j))/2;
        s_b = (c_map(i,j) + c_map(i,j-1))/2;

        G(c, c) = -(s_l + s_r + s_b);
        G(c, l) = s_l;
        G(c, r) = s_r;
        G(c, b) = s_b;        
    else
        s_l = (c_map(i,j) + c_map(i-1,j))/2;
        s_r = (c_map(i,j) + c_map(i+1,j))/2;
        s_b = (c_map(i,j) + c_map(i,j-1))/2;
        s_t = (c_map(i,j) + c_map(i,j+1))/2;

        G(c,c) = -(s_l+s_r+s_b+s_t);
        G(c,l) = s_l;
        G(c,r) = s_r;
        G(c,b) = s_b;
        G(c,t) = s_t;
    end

end


function [i,j] = get_co_ord_list(nx,ny)
    i = repmat([1:1:nx]',1,ny);
    j = repmat([1:1:ny],nx,1);

    i = reshape(i',1,[]);
    j = reshape(j',1,[]);
end

function [n] = ij_to_n(ny,i,j)
    n = ny*(i-1) + j;
end