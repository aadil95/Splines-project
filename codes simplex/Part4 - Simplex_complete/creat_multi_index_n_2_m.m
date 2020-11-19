function [kappa_m_mat,k0,k1] = creat_multi_index_n_2_m(d,m,loc_non_edge_t1)

index=1;
for i=d:-1:0
for j=d:-1:0
if(i+m+j)==d
    if loc_non_edge_t1==1
        kappa_m_mat(index,1)=m;
        kappa_m_mat(index,2)=i;
        kappa_m_mat(index,3)=j;
    elseif loc_non_edge_t1==2
        kappa_m_mat(index,1)=i;
        kappa_m_mat(index,2)=m;
        kappa_m_mat(index,3)=j;
    elseif loc_non_edge_t1==3
        kappa_m_mat(index,1)=i;
        kappa_m_mat(index,2)=j;
        kappa_m_mat(index,3)=m;
    end        
    
    k0(index,1)=i;
    k1(index,1)=j;
    index=index+1;
end
end
end
clear i j index

end

