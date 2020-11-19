function [est_coeff] = est_coeff_simplex(d,r,,total_bcoeffs_one_simp,Ntr,Ntst,kappa_mat,b_tr,b_tst)
for tj=1:n_simplices
for i_index=1:Ntr{tj}
    B_tr{tj}(i_index,:)=(factorial(d) ./ (factorial(kappa_mat(:,1)') ...
                    .*factorial(kappa_mat(:,2)').*factorial(kappa_mat(:,3))')) ...
                    .*((b_tr{tj}(i_index,1).^(kappa_mat(:,1)')) .* (b_tr{tj}(i_index,2).^(kappa_mat(:,2)')) ...
                    .* (b_tr{tj}(i_index,3).^(kappa_mat(:,3)')));


end
clear i_index
for i_index=1:Ntst{tj}
    B_tst{tj}(i_index,:)=(factorial(d) ./ (factorial(kappa_mat(:,1)') ...
                .*factorial(kappa_mat(:,2)').*factorial(kappa_mat(:,3))')) ...
                .*((b_tst{tj}(i_index,1).^(kappa_mat(:,1)')) .* (b_tst{tj}(i_index,2).^(kappa_mat(:,2)')) ...
                .* (b_tst{tj}(i_index,3).^(kappa_mat(:,3)')));;

end
clear i_index
end
Btr = blkdiag(B_tr{:});
Btst = blkdiag(B_tst{:});

%creating the smoothness matrix
H=[];
i_eqn=1;
% r=order of continuity
for m=0:r
   for tj=1:n_simplices
   for i_pair_row=1:length(nb(:,1)) %equations for each pair of simplices
   if (~isnan(nb(tj,i_pair_row)) && tj < nb(tj,i_pair_row))
   t1=nb(tj,i_pair_row);
   t2=tj;
   [non_edge_vertex_t1, loc_non_edge_t1]=setdiff(TRI(t1,:), TRI(t2,:));
   [non_edge_vertex_t2, loc_non_edge_t2]=setdiff(TRI(t2,:), TRI(t1,:));
   %create multi-indices for Left hand side
   [kappa_m_mat,k0,k1] = creat_multi_index_n_2_m(d,m,loc_non_edge_t1);
   %create multi-indices for Right hand side
   [kappa_zero_mat,gamma_mat] = creat_multi_index_n_2_zero(m,k0,k1,loc_non_edge_t2);  
   %coloumn locations of t1 in the kappa_mat matrix
   [t1_col_locs] = get_smoothness_equation_locations(kappa_mat,kappa_m_mat);

   bary_out_edge = cartesianToBarycentric(TRI, t2, V(non_edge_vertex_t1,:));%CHECK
   
   for i_row=1:length(t1_col_locs)
       H=[H; zeros(1,total_bcoeffs_one_simp*n_simplices)];
       H(i_eqn,total_bcoeffs_one_simp*(t1-1)+t1_col_locs(i_row))=-1;

       [t2_col_locs] = get_smoothness_equation_locations(kappa_mat,kappa_zero_mat{i_row});
        for rhs_row=1:length(t2_col_locs)
            H(i_eqn,total_bcoeffs_one_simp*(t2-1)+t2_col_locs(rhs_row))=+basis_func(m,gamma_mat(rhs_row,:),bary_out_edge);
        end
    i_eqn=i_eqn+1;
   end
   end
   end
   end
end
%arrange the smoothness matrix in the order of simplices
[rows_H, cols_H]=size(H);
F=Btr'*Btr;
M=[F H';H zeros(length(H(:,1)))];
if rank (H)< rows_H
    inv_M=pinv(M);
else
    inv_M=inv(M);
end
C1=inv_M(1:size(F),1:size(F));
est_coeff=C1*Btr'*Cmtr;% estimate the coeffients
end

