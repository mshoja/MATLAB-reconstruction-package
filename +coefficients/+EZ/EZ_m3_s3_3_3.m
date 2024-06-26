function EZ = EZ_m3_s3_3_3(dt,Dm0,Dm1,Dm2,Dm3,Ds0,Ds1,Ds2,Ds3)
%EZ_M3_S3_3_3
%    EZ = EZ_M3_S3_3_3(DT,DM0,DM1,DM2,DM3,DS0,DS1,DS2,DS3)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    12-May-2023 16:47:41

%E(Z) function of refined Hermite, J=3  K=3 ndiff_mu = 3, ndiff_sigma = 3
t2 = Dm0.^2;
t3 = Dm0.^3;
t4 = Dm1.^2;
t5 = Ds1.^2;
t6 = Ds1.^3;
t7 = dt.^2;
t8 = Dm0.*2.4e+1;
t9 = Ds1.*1.2e+1;
t10 = 1.0./Ds0.^3;
t11 = sqrt(dt);
t12 = t11.^3;
t13 = t11.^5;
EZ = [(t10.*(Ds0.*(Ds0.*(t11.*(t8+Dm0.*Dm1.*dt.*1.2e+1+Dm0.*dt.*t5.*1.2e+1+Dm0.*t4.*t7.*4.0+Dm2.*t2.*t7.*4.0+Dm0.*t5.^2.*t7.*4.0+Ds1.*Ds2.*t2.*t7.*2.0e+1+Dm1.*t5.*t7.*t8)-Ds0.*(t11.*(t9+Dm0.*Ds2.*dt.*1.2e+1+Dm1.*dt.*t9+Dm1.*t6.*t7.*4.0+Ds1.*t4.*t7.*8.0+Ds3.*t2.*t7.*6.0+Dm0.*Dm1.*Ds2.*t7.*1.8e+1+Dm0.*Dm2.*Ds1.*t7.*1.0e+1+Dm0.*Ds2.*t5.*t7.*1.0e+1)+Ds0.*(Ds0.*(t11.*(Ds3.*dt.*3.0+Dm1.*Ds3.*t7.*6.0+Dm2.*Ds2.*t7.*5.0+Ds3.*t5.*t7)+Ds0.*Ds2.*Ds3.*t13)-t11.*(Dm2.*dt.*6.0+Dm0.*Dm3.*t7.*4.0+Dm1.*Dm2.*t7.*6.0+Dm2.*t5.*t7.*2.0+Dm0.*Ds2.^2.*t7.*4.0+Dm0.*Ds1.*Ds3.*t7.*2.0+Dm1.*Ds1.*Ds2.*t7.*8.0))))-t11.*(Ds2.*t3.*t7.*4.0+dt.*t2.*t9+t2.*t6.*t7.*1.6e+1+Dm1.*t2.*t7.*t9))+t3.*t5.*t13.*8.0))./2.4e+1;(t10.*(Ds0.*(dt.*t2.*2.4e+1-Ds0.*(-Ds0.*(Dm1.*dt.*2.4e+1+dt.*t5.*6.0+t4.*t7.*1.6e+1-Ds0.*(Ds2.*dt.*1.2e+1-Ds0.*(Dm3.*t7.*8.0-Ds1.*Ds3.*t7.*5.0)+Dm0.*Ds3.*t7.*2.2e+1+Dm1.*Ds2.*t7.*3.2e+1+Dm2.*Ds1.*t7.*1.4e+1)+Dm0.*Dm2.*t7.*2.8e+1+Dm1.*t5.*t7.*2.8e+1+Dm0.*Ds1.*Ds2.*t7.*4.4e+1+2.4e+1)+Dm0.*Ds1.*dt.*4.8e+1+Dm0.*t6.*t7.*2.8e+1+Ds2.*t2.*t7.*4.0e+1+Dm0.*Dm1.*Ds1.*t7.*8.4e+1)+Dm1.*t2.*t7.*2.4e+1+t2.*t5.*t7.*6.8e+1)-Ds1.*t3.*t7.*2.4e+1))./2.4e+1;(t10.*(t3.*t12.*8.0-Ds0.*(Ds0.*(Ds0.*(t11.*(t9+dt.*t6+Dm0.*Ds2.*dt.*3.2e+1+Dm1.*Ds1.*dt.*3.2e+1)-Ds0.*(t11.*(Dm2.*dt.*1.4e+1+Ds1.*Ds2.*dt.*2.0)-Ds0.*Ds3.*t12.*7.0))-t11.*(t8+Dm0.*Dm1.*dt.*3.6e+1+Dm0.*dt.*t5.*3.8e+1))+Ds1.*t2.*t12.*4.8e+1)))./8.0];
