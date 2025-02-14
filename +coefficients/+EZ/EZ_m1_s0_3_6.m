function EZ = EZ_m1_s0_3_6(dt,Dm0,Dm1,Ds0)
%EZ_m1_s0_3_6
%    EZ = EZ_m1_s0_3_6(DT,Dm0,Dm1,Ds0)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    08-Mar-2024 12:53:36

%E(Z) function of refined Hermite, J=3  K=6 ndiff_mu = 1, ndiff_sigma = 0
t2 = Dm0.^2;
t3 = Dm0.^3;
t4 = Dm1.^2;
t5 = Dm1.^3;
t7 = Dm1.^5;
t8 = Ds0.^2;
t9 = dt.^2;
t10 = dt.^3;
t12 = dt.^5;
t13 = sqrt(dt);
t6 = t4.^2;
t11 = t9.^2;
t14 = t13.^3;
t15 = t13.^5;
t16 = t13.^7;
t17 = t14.^3;
mt1 = [(Dm0.*t13.*7.2e+2+Dm0.*Dm1.*t14.*3.6e+2+Dm0.*t4.*t15.*1.2e+2+Dm0.*t5.*t16.*3.0e+1+Dm0.*t6.*t17.*6.0+Dm0.*t7.*t13.^11)./(Ds0.*7.2e+2);(dt.*t2.*3.6e+2+t8.*(Dm1.*dt.*3.6e+2+t4.*t9.*2.4e+2+t5.*t10.*1.2e+2+t6.*t11.*4.8e+1+t7.*t12.*1.6e+1+3.6e+2)+Dm1.*t2.*t9.*3.6e+2+t2.*t4.*t10.*2.1e+2+t2.*t5.*t11.*9.0e+1+t2.*t6.*t12.*3.1e+1)./(t8.*3.6e+2)];
mt2 = [(1.0./Ds0.^3.*(t3.*t14.*1.2e+2+t8.*(Dm0.*t13.*3.6e+2+Dm0.*Dm1.*t14.*5.4e+2+Dm0.*t4.*t15.*4.8e+2+Dm0.*t5.*t16.*3.15e+2+Dm0.*t6.*t17.*1.66e+2)+Dm1.*t3.*t15.*1.8e+2+t3.*t4.*t16.*1.5e+2+t3.*t5.*t17.*9.0e+1))./1.2e+2];
EZ = [mt1;mt2];
end
