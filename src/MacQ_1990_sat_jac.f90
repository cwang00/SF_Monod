      MacQ_1990_sat_jac(1,1) = mpars_kmax*A*S*X/((mpars_Ka+A)*mpars_Rs*(
     1   S+mpars_Ks)**2)-mpars_kmax*A*X/((mpars_Ka+A)*mpars_Rs*(S+mpars_
     2   Ks))
      MacQ_1990_sat_jac(1,2) = mpars_kmax*A*S*X/((mpars_Ka+A)**2*mpars_R
     1   s*(S+mpars_Ks))-mpars_kmax*S*X/((mpars_Ka+A)*mpars_Rs*(S+mpars_
     2   Ks))
      MacQ_1990_sat_jac(1,3) = -mpars_kmax*A*S/((mpars_Ka+A)*mpars_Rs*(S
     1   +mpars_Ks))
      MacQ_1990_sat_jac(2,1) = mpars_kmax*mpars_r*A*S*X/((mpars_Ka+A)*mp
     1   ars_Rm*(S+mpars_Ks)**2)-mpars_kmax*mpars_r*A*X/((mpars_Ka+A)*mp
     2   ars_Rm*(S+mpars_Ks))
      MacQ_1990_sat_jac(2,2) = mpars_kmax*mpars_r*A*S*X/((mpars_Ka+A)**2
     1   *mpars_Rm*(S+mpars_Ks))-mpars_kmax*mpars_r*S*X/((mpars_Ka+A)*mp
     2   ars_Rm*(S+mpars_Ks))
      MacQ_1990_sat_jac(2,3) = -mpars_kmax*mpars_r*A*S/((mpars_Ka+A)*mpa
     1   rs_Rm*(S+mpars_Ks))
      MacQ_1990_sat_jac(3,1) = mpars_kmax*A*mpars_Y*X/((mpars_Ka+A)*(S+m
     1   pars_Ks))-mpars_kmax*A*mpars_Y*S*X/((mpars_Ka+A)*(S+mpars_Ks)**
     2   2)
      MacQ_1990_sat_jac(3,2) = mpars_kmax*mpars_Y*S*X/((mpars_Ka+A)*(S+m
     1   pars_Ks))-mpars_kmax*A*mpars_Y*S*X/((mpars_Ka+A)**2*(S+mpars_Ks
     2   ))
      MacQ_1990_sat_jac(3,3) = mpars_kmax*A*mpars_Y*S/((mpars_Ka+A)*(S+m
     1   pars_Ks))-mpars_decay
