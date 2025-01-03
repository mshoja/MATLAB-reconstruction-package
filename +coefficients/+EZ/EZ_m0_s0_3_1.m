function EZ = EZ_m0_s0_3_1(dt,Dm0,Ds0)
%EZ_m0_s0_3_1
%    EZ = EZ_m0_s0_3_1(DT,Dm0,Ds0)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    08-Mar-2024 12:52:11

%E(Z) function of refined Hermite, J=3  K=1 ndiff_mu = 0, ndiff_sigma = 0
t2 = Ds0./1.0e+100;
EZ = [(Dm0.*sqrt(dt))./Ds0;t2+1.0;t2];
end
