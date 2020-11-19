function zpred = kf_calcHx(t, x, u)
    u=x(1);
    v=x(2);
    w=x(3);
    C_alfa_up=x(4);

    alfa_true=atan(w/u);
    beta_true=atan(v/sqrt(u^2+w^2));
    V_true=sqrt(u^2+v^2+w^2);
    
    zpred(1) = alfa_true*(1+C_alfa_up);
    zpred(2) = beta_true;
    zpred(3) = V_true;
    end
    