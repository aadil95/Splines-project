function Hx = kf_calc_Hx(t, x, u)
u=x(1);
v=x(2);
w=x(3);
C_alfa_up=x(4);

% shortening squares
sum_sq=u^2+v^2+w^2;
u_w_sq=u^2+w^2;
u_w_sqrt=sqrt(u_w_sq);

Hx = zeros(3,4);

Hx(1,:)=[-w*(1+C_alfa_up)/u_w_sq 0 u*(1+C_alfa_up)/u_w_sq atan(w/u)];
Hx(2,:)=[-3*u*v/(sum_sq*u_w_sqrt) u_w_sqrt/sum_sq -3*w*v/(sum_sq*u_w_sqrt) 0];
Hx(3,:)=[u v w 0]/sqrt(sum_sq);
               
end
