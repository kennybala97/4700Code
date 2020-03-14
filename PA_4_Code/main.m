clear all;
close all;
clc;

% Mesh size

N = 100;
M = 100;

% Number of iterations

n_max = 1000; 

% Set up the mesh
[X,Y] = meshgrid(1:1:M,1:1:N);

finite_diff_mesh = zeros(M,N);

% Set up boundary conditions

finite_diff_mesh(:,1) = 100;
finite_diff_mesh(:,end) = 100;

wrap_n = @(x, N) (1 + mod(x-1, N));

n = 0;

while n < n_max

    for i = 2:1:M-1
        for j = 2:1:N-1            
            i_vec = wrap_n([i + 1, i - 1],N);
            finite_diff_mesh(i,j) = (finite_diff_mesh(i_vec(1), j + 1) + finite_diff_mesh(i_vec(1), j - 1) + finite_diff_mesh(i_vec(2), j + 1) + finite_diff_mesh(i_vec(2), j - 1))/4;

        end

    end

    n = n+1;
end
figure;
h = surf(X,Y,finite_diff_mesh);
set(h,'LineStyle','none');
colormap 'jet';

[dx, dy] = gradient(finite_diff_mesh);

figure;
quiver(X,Y,dx,dy);

