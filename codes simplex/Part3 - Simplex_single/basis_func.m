function [eval] = basis_func(m,gamma,bary_coordinate_vertex)

b_k   = 1;
k_fac = 1;
for i_k=1:length(bary_coordinate_vertex)
    b_k=(bary_coordinate_vertex(i_k)^gamma(i_k))*b_k;
    k_fac=(factorial(gamma(i_k)))*k_fac;
end
eval=factorial(m)*b_k/k_fac;
end

