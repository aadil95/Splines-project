function [eval] = multinomial_coeff(d,k,bary_coordinate)

b_k   = 1;
k_fac = 1;
for i_k=1:length(bary_coordinate)
    b_k=(bary_coordinate(i_k)^k(i_k))*b_k;
    k_fac=(factorial(k(i_k)))*k_fac;
end
eval=factorial(d)*b_k/k_fac;
end

