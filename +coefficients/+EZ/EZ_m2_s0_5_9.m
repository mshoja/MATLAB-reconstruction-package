function EZ = EZ_m2_s0_5_9(dt,Dm0,Dm1,Dm2,Ds0)
%EZ_m2_s0_5_9
%    EZ = EZ_m2_s0_5_9(DT,Dm0,Dm1,Dm2,Ds0)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    25-Mar-2024 19:40:02

%E(Z) function of refined Hermite, J=5  K=9 ndiff_mu = 2, ndiff_sigma = 0
t2 = Dm0.^2;
t3 = Dm0.^3;
t4 = Dm1.^2;
t6 = Dm1.^3;
t7 = Dm2.^2;
t8 = Dm0.^5;
t10 = Dm2.^3;
t12 = Dm1.^5;
t15 = Dm2.^5;
t16 = Dm1.^7;
t17 = Ds0.^2;
t18 = dt.^2;
t19 = dt.^3;
t21 = dt.^5;
t23 = dt.^7;
t26 = sqrt(dt);
t5 = t2.^2;
t9 = t4.^2;
t11 = t2.^3;
t13 = t7.^2;
t14 = t4.^3;
t20 = t18.^2;
t22 = t18.^3;
t27 = t26.^3;
t28 = t26.^5;
t29 = t26.^7;
t31 = t26.^11;
t32 = t26.^13;
t34 = t26.^17;
t35 = t6.*t19.*4.8384e+5;
t36 = Dm0.*t26.*2.90304e+6;
t24 = t9.^2;
t25 = t20.^2;
t30 = t27.^3;
t33 = t27.^5;
et1 = Dm2.*t27.*7.2576e+5+t17.*(t17.*(t15.*t33.*3.9609e+4+Dm1.*t15.*t34.*1.11483e+5)+t10.*t30.*1.27008e+5+Dm1.*t10.*t31.*2.44944e+5+Dm0.*t13.*t32.*1.29168e+5+t4.*t10.*t32.*2.52288e+5+t6.*t10.*t33.*1.8306e+5+t2.*t15.*t34.*7.8858e+4+t9.*t10.*t34.*1.04454e+5+Dm0.*Dm1.*t13.*t33.*3.02778e+5+Dm0.*t4.*t13.*t34.*3.73446e+5)+Dm1.*Dm2.*t28.*7.2576e+5+Dm2.*t4.*t29.*4.2336e+5+Dm0.*t7.*t29.*4.2336e+5+Dm2.*t6.*t30.*1.8144e+5+Dm2.*t9.*t31.*6.2496e+4+Dm2.*t12.*t32.*1.8144e+4+Dm2.*t14.*t33.*4.572e+3+Dm2.*t16.*t34.*1.02e+3;
et2 = t2.*t10.*t31.*1.7136e+5+t3.*t13.*t33.*5.8032e+4+t2.*t4.*t10.*t33.*3.19896e+5+t2.*t6.*t10.*t34.*2.23488e+5+Dm0.*Dm1.*t7.*t30.*6.16896e+5+Dm0.*t4.*t7.*t31.*4.87872e+5+Dm0.*t6.*t7.*t32.*2.74752e+5+Dm1.*t2.*t10.*t32.*3.21408e+5+Dm0.*t7.*t9.*t33.*1.2258e+5+Dm1.*t3.*t13.*t34.*1.32672e+5+Dm0.*t7.*t12.*t34.*4.5828e+4;
et3 = t36+t17.*(et1+et2)+Dm0.*Dm1.*t27.*1.45152e+6+Dm0.*t4.*t28.*4.8384e+5+Dm2.*t2.*t28.*4.8384e+5+Dm0.*t6.*t29.*1.2096e+5+Dm0.*t9.*t30.*2.4192e+4+Dm0.*t12.*t31.*4.032e+3+Dm0.*t14.*t32.*5.76e+2+Dm0.*t16.*t33.*7.2e+1+Dm0.*t24.*t34.*8.0+t3.*t7.*t30.*9.6768e+4+t5.*t10.*t32.*1.9584e+4+t8.*t13.*t34.*3.968e+3+t3.*t4.*t7.*t32.*1.0368e+5+t3.*t6.*t7.*t33.*5.5296e+4+t3.*t7.*t9.*t34.*2.3232e+4+t4.*t5.*t10.*t34.*3.4304e+4+Dm1.*Dm2.*t2.*t29.*4.8384e+5+Dm2.*t2.*t4.*t30.*2.66112e+5+Dm2.*t2.*t6.*t31.*1.04832e+5+Dm1.*t3.*t7.*t31.*1.37088e+5;
et4 = Dm2.*t2.*t9.*t32.*3.2832e+4+Dm1.*t5.*t10.*t33.*3.5712e+4+Dm2.*t2.*t12.*t33.*8.64e+3+Dm2.*t2.*t14.*t34.*1.976e+3;
et5 = t35+Dm1.*dt.*1.45152e+6+t4.*t18.*9.6768e+5+t9.*t20.*1.93536e+5+t12.*t21.*6.4512e+4+t14.*t22.*1.8432e+4+t16.*t23.*4.608e+3+t24.*t25.*1.024e+3;
et6 = t17.*(t7.*t19.*6.3504e+5+t17.*(t13.*t22.*3.16872e+5+Dm1.*t13.*t23.*9.63738e+5+Dm0.*t15.*t25.*5.63787e+5+t4.*t13.*t25.*1.527714e+6)+Dm1.*t7.*t20.*1.342656e+6+Dm0.*t10.*t21.*9.04176e+5+t4.*t7.*t21.*1.521072e+6+t6.*t7.*t22.*1.212192e+6+t2.*t13.*t23.*7.09722e+5+t7.*t9.*t23.*7.57026e+5+t7.*t12.*t25.*3.92508e+5+Dm0.*Dm1.*t10.*t22.*2.293056e+6+Dm0.*t4.*t10.*t23.*3.058236e+6+Dm0.*t6.*t10.*t25.*2.838528e+6+Dm1.*t2.*t13.*t25.*2.093238e+6)+Dm0.*Dm2.*t18.*1.69344e+6;
et7 = t2.*t7.*t20.*1.02816e+6+t3.*t10.*t22.*4.64256e+5+t5.*t13.*t25.*1.76896e+5+t2.*t4.*t7.*t22.*2.23776e+6+t2.*t6.*t7.*t23.*1.691496e+6+t3.*t4.*t10.*t25.*1.45792e+6+t2.*t7.*t9.*t25.*1.001724e+6+Dm0.*Dm1.*Dm2.*t19.*2.66112e+6+Dm0.*Dm2.*t4.*t20.*2.310336e+6+Dm0.*Dm2.*t6.*t21.*1.435392e+6+Dm0.*Dm2.*t9.*t22.*7.05888e+5+Dm0.*Dm2.*t12.*t23.*2.89872e+5+Dm0.*Dm2.*t14.*t25.*1.02748e+5+Dm1.*t2.*t7.*t21.*2.078496e+6+Dm1.*t3.*t10.*t23.*1.136016e+6+1.45152e+6;
et8 = dt.*t2.*1.45152e+6+t17.*(et5+et6+et7)+Dm1.*t2.*t18.*1.45152e+6+Dm2.*t3.*t19.*4.8384e+5+t2.*t4.*t19.*8.4672e+5+t2.*t6.*t20.*3.6288e+5+t2.*t9.*t21.*1.24992e+5+t5.*t7.*t21.*1.37088e+5+t2.*t12.*t22.*3.6288e+4+t2.*t14.*t23.*9.144e+3+t8.*t10.*t23.*3.5712e+4+t2.*t16.*t25.*2.04e+3+t4.*t5.*t7.*t23.*2.73024e+5+t5.*t6.*t7.*t25.*1.9584e+5+Dm1.*Dm2.*t3.*t20.*7.2576e+5+Dm2.*t3.*t4.*t21.*5.88672e+5+Dm2.*t3.*t6.*t22.*3.38688e+5+Dm1.*t5.*t7.*t22.*2.66112e+5+Dm2.*t3.*t9.*t23.*1.53792e+5;
et9 = Dm2.*t3.*t12.*t25.*5.832e+4+Dm1.*t8.*t10.*t25.*8.448e+4;
et10 = Dm2.*t27.*1.69344e+6+t17.*(t10.*t30.*1.478736e+6+Dm1.*t10.*t31.*4.71744e+6+Dm0.*t13.*t32.*3.356316e+6+t4.*t10.*t32.*7.8813e+6+t6.*t10.*t33.*9.126e+6+t15.*t17.*t33.*1.280727e+6+Dm0.*Dm1.*t13.*t33.*1.2127698e+7)+Dm1.*Dm2.*t28.*3.6288e+6+Dm2.*t4.*t29.*4.29408e+6+Dm0.*t7.*t29.*3.532032e+6+Dm2.*t6.*t30.*3.6288e+6+Dm2.*t9.*t31.*2.421216e+6+Dm2.*t12.*t32.*1.34568e+6+Dm2.*t14.*t33.*6.4414e+5+t2.*t10.*t31.*3.699072e+6+t3.*t13.*t33.*2.691748e+6;
et11 = t2.*t4.*t10.*t33.*1.8393456e+7+Dm0.*Dm1.*t7.*t30.*9.332064e+6+Dm0.*t4.*t7.*t31.*1.3084992e+7+Dm0.*t6.*t7.*t32.*1.283148e+7+Dm1.*t2.*t10.*t32.*1.141128e+7+Dm0.*t7.*t9.*t33.*9.820416e+6;
et12 = t36+t17.*(et10+et11)+Dm0.*Dm1.*t27.*4.35456e+6+Dm0.*t4.*t28.*3.87072e+6+Dm2.*t2.*t28.*3.14496e+6+Dm0.*t6.*t29.*2.54016e+6+Dm0.*t9.*t30.*1.338624e+6+Dm0.*t12.*t31.*5.92704e+5+Dm0.*t14.*t32.*2.26944e+5+Dm0.*t16.*t33.*7.668e+4+t3.*t7.*t30.*2.052288e+6+t5.*t10.*t32.*1.022688e+6+t3.*t4.*t7.*t32.*7.019568e+6+t3.*t6.*t7.*t33.*6.57648e+6+Dm1.*Dm2.*t2.*t29.*6.53184e+6+Dm2.*t2.*t4.*t30.*7.346304e+6+Dm2.*t2.*t6.*t31.*5.854464e+6+Dm1.*t3.*t7.*t31.*5.225472e+6+Dm2.*t2.*t9.*t32.*3.673296e+6;
et13 = Dm1.*t5.*t10.*t33.*3.0528e+6+Dm2.*t2.*t12.*t33.*1.9188e+6;
et14 = t35+Dm1.*dt.*4.8384e+5+t4.*t18.*5.6448e+5+t9.*t20.*3.33312e+5+t12.*t21.*1.93536e+5+t14.*t22.*9.7536e+4+t16.*t23.*4.352e+4+t17.*(t17.*(t13.*t22.*6.38163e+5+Dm1.*t13.*t23.*2.720322e+6)+t7.*t19.*4.85856e+5+Dm1.*t7.*t20.*1.574496e+6+Dm0.*t10.*t21.*1.449312e+6+t4.*t7.*t21.*2.702544e+6+t6.*t7.*t22.*3.23928e+6+t2.*t13.*t23.*2.152514e+6+t7.*t9.*t23.*3.026198e+6+Dm0.*Dm1.*t10.*t22.*5.38056e+6+Dm0.*t4.*t10.*t23.*1.0416208e+7)+Dm0.*Dm2.*t18.*8.8704e+5;
et15 = t2.*t7.*t20.*1.328544e+6+t3.*t10.*t22.*1.263792e+6+t2.*t4.*t7.*t22.*6.937272e+6+t2.*t6.*t7.*t23.*8.01152e+6+Dm0.*Dm1.*Dm2.*t19.*2.33856e+6+Dm0.*Dm2.*t4.*t20.*3.337152e+6+Dm0.*Dm2.*t6.*t21.*3.372096e+6+Dm0.*Dm2.*t9.*t22.*2.680704e+6+Dm0.*Dm2.*t12.*t23.*1.7728e+6+Dm1.*t2.*t7.*t21.*4.182528e+6+Dm1.*t3.*t10.*t23.*4.56e+6+2.4192e+5;
et16 = t17.*(dt.*t2.*4.8384e+5+t17.*(et14+et15)+Dm1.*t2.*t18.*9.6768e+5+Dm2.*t3.*t19.*5.6448e+5+t2.*t4.*t19.*1.08864e+6+t2.*t6.*t20.*8.8704e+5+t2.*t9.*t21.*5.76576e+5+t5.*t7.*t21.*4.01856e+5+t2.*t12.*t22.*3.14496e+5+t2.*t14.*t23.*1.4852e+5+t8.*t10.*t23.*2.1952e+5+t4.*t5.*t7.*t23.*1.966e+6+Dm1.*Dm2.*t3.*t20.*1.45152e+6+Dm2.*t3.*t4.*t21.*1.994496e+6+Dm2.*t3.*t6.*t22.*1.927296e+6+Dm1.*t5.*t7.*t22.*1.227744e+6+Dm2.*t3.*t9.*t23.*1.45976e+6)+t5.*t18.*8.064e+4;
et17 = Dm1.*t5.*t19.*1.6128e+5+Dm2.*t8.*t20.*5.376e+4+t4.*t5.*t20.*1.7472e+5+t5.*t6.*t21.*1.344e+5+t5.*t9.*t22.*8.1648e+4+t5.*t12.*t23.*4.144e+4+t7.*t11.*t22.*2.4192e+4+Dm1.*Dm2.*t8.*t21.*1.344e+5+Dm2.*t4.*t8.*t22.*1.77408e+5+Dm2.*t6.*t8.*t23.*1.6352e+5+Dm1.*t7.*t11.*t23.*7.168e+4;
et18 = Dm0.*t26.*2.4192e+5+t17.*(Dm2.*t27.*2.2176e+5+t17.*(t10.*t30.*4.99968e+5+Dm1.*t10.*t31.*2.1609e+6+Dm0.*t13.*t32.*1.933161e+6+t4.*t10.*t32.*4.86375e+6)+Dm1.*Dm2.*t28.*7.056e+5+Dm2.*t4.*t29.*1.21128e+6+Dm0.*t7.*t29.*9.35088e+5+Dm2.*t6.*t30.*1.47e+6+Dm2.*t9.*t31.*1.402632e+6+Dm2.*t12.*t32.*1.11328e+6+t2.*t10.*t31.*1.801728e+6+Dm0.*Dm1.*t7.*t30.*3.49104e+6+Dm0.*t4.*t7.*t31.*6.857712e+6+Dm0.*t6.*t7.*t32.*9.3701e+6+Dm1.*t2.*t10.*t32.*7.6036e+6);
et19 = Dm0.*Dm1.*t27.*6.048e+5+Dm0.*t4.*t28.*8.4672e+5+Dm2.*t2.*t28.*6.4512e+5+Dm0.*t6.*t29.*8.568e+5+Dm0.*t9.*t30.*6.91488e+5+Dm0.*t12.*t31.*4.6872e+5+Dm0.*t14.*t32.*2.7544e+5+t3.*t7.*t30.*8.72256e+5+t5.*t10.*t32.*8.0924e+5+t3.*t4.*t7.*t32.*6.07352e+6+Dm1.*Dm2.*t2.*t29.*2.016e+6+Dm2.*t2.*t4.*t30.*3.366048e+6+Dm2.*t2.*t6.*t31.*3.95136e+6+Dm1.*t3.*t7.*t31.*3.18024e+6+Dm2.*t2.*t9.*t32.*3.6346e+6;
et20 = t8.*t28.*1.6128e+4+t17.*(t3.*t27.*1.6128e+5+t17.*(et18+et19)+Dm1.*t3.*t28.*4.032e+5+Dm2.*t5.*t29.*2.0832e+5+t3.*t4.*t29.*5.5104e+5+t3.*t6.*t30.*5.376e+5+t3.*t9.*t31.*4.15296e+5+t7.*t8.*t31.*1.62288e+5+t3.*t12.*t32.*2.6824e+5+Dm1.*Dm2.*t5.*t30.*6.384e+5+Dm2.*t4.*t5.*t31.*1.035216e+6+Dm2.*t5.*t6.*t32.*1.1732e+6+Dm1.*t7.*t8.*t32.*5.7736e+5)+Dm1.*t8.*t29.*4.032e+4+Dm2.*t11.*t30.*1.344e+4+t4.*t8.*t30.*5.376e+4+t6.*t8.*t31.*5.04e+4+t8.*t9.*t32.*3.7072e+4+Dm0.^7.*t7.*t32.*7.168e+3;
et21 = Dm1.*Dm2.*t11.*t31.*4.032e+4+Dm2.*t4.*t11.*t32.*6.3392e+4;
mt1 = [(et3+et4)./(Ds0.*2.90304e+6);(et8+et9)./(t17.*1.45152e+6)];
mt2 = [(1.0./Ds0.^3.*(t3.*t27.*9.6768e+5+t17.*(et12+et13)+Dm1.*t3.*t28.*1.45152e+6+Dm2.*t5.*t29.*4.8384e+5+t3.*t4.*t29.*1.2096e+6+t3.*t6.*t30.*7.2576e+5+t3.*t9.*t31.*3.46752e+5+t7.*t8.*t31.*1.77408e+5+t3.*t12.*t32.*1.39104e+5+t3.*t14.*t33.*4.84e+4+t10.*t11.*t33.*5.632e+4+t4.*t7.*t8.*t33.*5.6064e+5+Dm1.*Dm2.*t5.*t30.*9.6768e+5+Dm2.*t4.*t5.*t31.*1.032192e+6+Dm2.*t5.*t6.*t32.*7.74144e+5+Dm1.*t7.*t8.*t32.*4.35456e+5+Dm2.*t5.*t9.*t33.*4.5552e+5))./9.6768e+5];
mt3 = [(1.0./t17.^2.*(et16+et17))./8.064e+4;(1.0./Ds0.^5.*(et20+et21))./1.6128e+4];
EZ = [mt1;mt2;mt3];