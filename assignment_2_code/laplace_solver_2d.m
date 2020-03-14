
function [Vmap] = laplace_solver_2d(nx,ny,delta,initial_mesh)
    G = sparse(nx*ny);
    B = zeros(1,nx*ny);
    
    [co_ords.i,co_ords.j] = get_co_ord_list(nx,ny);    
    co_ords.n = ij_to_n(ny,co_ords.i,co_ords.j);
    
    
    for u = 1:1:length(co_ords.n)
        i = co_ords.i(u);
        j = co_ords.j(u);
        n = co_ords.n(u);
        
        G(n, :) = 0;

        if ~isnan(initial_mesh(i,j))
            G(n, n) = 1;
            B(n) = initial_mesh(i,j);
        else
            i_vals = ij_to_n(ny,[i, i-1, i+1, i, i,],[j, j, j, j-1, j+1]);
            vals = [-4, 1, 1, 1, 1]; 
                        
            if (j == 1) || (j == ny)
                i_vals(4:5) = [];
                vals(4:5) = [];         
            end
            
            vals = vals./(delta*delta);
            G(n,i_vals) = vals;
        end
    end

    V = G\B';
    
    Vmap = zeros(nx,ny);
    
    for i = 1:nx
        for j = 1:ny
            n = j + (i - 1) * ny;

            Vmap(i, j) = V(n);
        end
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