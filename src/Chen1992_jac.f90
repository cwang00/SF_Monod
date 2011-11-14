      Chen1992_jac(1,1) = mpars_kmax_tol*tol*A*X1/((mpars_Kdo+A)*(mpars_
     1   Ktol+tol)**2*mpars_Rtol)-mpars_kmax_tol*A*X1/((mpars_Kdo+A)*(mp
     2   ars_Ktol+tol)*mpars_Rtol)
      Chen1992_jac(1,2) = 0
      Chen1992_jac(1,3) = mpars_kmax_tol*tol*A*X1/((mpars_Kdo+A)**2*(mpa
     1   rs_Ktol+tol)*mpars_Rtol)-mpars_kmax_tol*tol*X1/((mpars_Kdo+A)*(
     2   mpars_Ktol+tol)*mpars_Rtol)
      Chen1992_jac(1,4) = -mpars_kmax_tol*tol*A/((mpars_Kdo+A)*(mpars_Kt
     1   ol+tol)*mpars_Rtol)
      Chen1992_jac(1,5) = 0
      Chen1992_jac(2,1) = 0
      Chen1992_jac(2,2) = ben*mpars_kmax_ben*A*X2/((mpars_Kben+ben)**2*(
     1   mpars_Kdo+A)*mpars_Rben)-mpars_kmax_ben*A*X2/((mpars_Kben+ben)*
     2   (mpars_Kdo+A)*mpars_Rben)
      Chen1992_jac(2,3) = ben*mpars_kmax_ben*A*X2/((mpars_Kben+ben)*(mpa
     1   rs_Kdo+A)**2*mpars_Rben)-ben*mpars_kmax_ben*X2/((mpars_Kben+ben
     2   )*(mpars_Kdo+A)*mpars_Rben)
      Chen1992_jac(2,4) = 0
      Chen1992_jac(2,5) = -ben*mpars_kmax_ben*A/((mpars_Kben+ben)*(mpars
     1   _Kdo+A)*mpars_Rben)
      Chen1992_jac(3,1) = mpars_kmax_tol*mpars_r_tol_do*tol*A*X1/((mpars
     1   _Kdo+A)*(mpars_Ktol+tol)**2*mpars_Rtol)-mpars_kmax_tol*mpars_r_
     2   tol_do*A*X1/((mpars_Kdo+A)*(mpars_Ktol+tol)*mpars_Rtol)
      Chen1992_jac(3,2) = ben*mpars_kmax_ben*mpars_r_ben_do*A*X2/((mpars
     1   _Kben+ben)**2*(mpars_Kdo+A)*mpars_Rben)-mpars_kmax_ben*mpars_r_
     2   ben_do*A*X2/((mpars_Kben+ben)*(mpars_Kdo+A)*mpars_Rben)
      Chen1992_jac(3,3) = -ben*mpars_kmax_ben*mpars_r_ben_do*X2/((mpars_
     1   Kben+ben)*(mpars_Kdo+A)*mpars_Rben)+ben*mpars_kmax_ben*mpars_r_
     2   ben_do*A*X2/((mpars_Kben+ben)*(mpars_Kdo+A)**2*mpars_Rben)-mpar
     3   s_kmax_tol*mpars_r_tol_do*tol*X1/((mpars_Kdo+A)*(mpars_Ktol+tol
     4   )*mpars_Rtol)+mpars_kmax_tol*mpars_r_tol_do*tol*A*X1/((mpars_Kd
     5   o+A)**2*(mpars_Ktol+tol)*mpars_Rtol)
      Chen1992_jac(3,4) = -mpars_kmax_tol*mpars_r_tol_do*tol*A/((mpars_K
     1   do+A)*(mpars_Ktol+tol)*mpars_Rtol)
      Chen1992_jac(3,5) = -ben*mpars_kmax_ben*mpars_r_ben_do*A/((mpars_K
     1   ben+ben)*(mpars_Kdo+A)*mpars_Rben)
      Chen1992_jac(4,1) = mpars_kmax_tol*A*mpars_Y*X1/((mpars_Kdo+A)*(mp
     1   ars_Ktol+tol))-mpars_kmax_tol*tol*A*mpars_Y*X1/((mpars_Kdo+A)*(
     2   mpars_Ktol+tol)**2)
      Chen1992_jac(4,2) = 0
      Chen1992_jac(4,3) = mpars_kmax_tol*tol*mpars_Y*X1/((mpars_Kdo+A)*(
     1   mpars_Ktol+tol))-mpars_kmax_tol*tol*A*mpars_Y*X1/((mpars_Kdo+A)
     2   **2*(mpars_Ktol+tol))
      Chen1992_jac(4,4) = mpars_kmax_tol*tol*A*mpars_Y/((mpars_Kdo+A)*(m
     1   pars_Ktol+tol))-mpars_decay
      Chen1992_jac(4,5) = 0
      Chen1992_jac(5,1) = 0
      Chen1992_jac(5,2) = mpars_kmax_ben*A*mpars_Y*X2/((mpars_Kben+ben)*
     1   (mpars_Kdo+A))-ben*mpars_kmax_ben*A*mpars_Y*X2/((mpars_Kben+ben
     2   )**2*(mpars_Kdo+A))
      Chen1992_jac(5,3) = ben*mpars_kmax_ben*mpars_Y*X2/((mpars_Kben+ben
     1   )*(mpars_Kdo+A))-ben*mpars_kmax_ben*A*mpars_Y*X2/((mpars_Kben+b
     2   en)*(mpars_Kdo+A)**2)
      Chen1992_jac(5,4) = 0
      Chen1992_jac(5,5) = ben*mpars_kmax_ben*A*mpars_Y/((mpars_Kben+ben)
     1   *(mpars_Kdo+A))-mpars_decay
