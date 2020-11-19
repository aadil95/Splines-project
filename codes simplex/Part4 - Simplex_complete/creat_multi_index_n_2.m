function [kappa_mat] = creat_multi_index_n_2(d)

index=1;
for i=d:-1:0
for j=d:-1:0
for k=d:-1:0
if(i+j+k)==d
kappa_mat(index,:)=[i j k];
index=index+1;
end
end
end
end
clear i j k index

end

