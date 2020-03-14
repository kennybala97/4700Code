
function [num_out] = uniform_random_gen(a,b,N)
    num_out = a + (b-a).*rand(N,1);
end
