      real*8 function get_kappa(temp_kap, rho_kap)

      include "../lib/flux_cal.h"
      
      real*8 Rkap, logt, kapl, kapr, kapm, temp_kap, rho_kap
      integer ir, it, i,j

      Rkap=log10(rho_kap/((temp_kap/1.e6)**3))
      logt=log10(temp_kap)

      ir=1
      do i=1, numr
         if(Rkap.ge.gridR(i).and.Rkap.le.gridR(i+1)) then
            ir=i
            exit
         end if
      end do

      it=1
      do i=1, numo-1
         if(logt.le.gridT(i).and.logt.ge.gridT(i+1)) then
            it=i
            exit
         end if
      end do

      kapl=kap(ir,it)+(kap(ir,it+1)-kap(ir,it))*
     &     (logt-gridT(it))/(gridT(it+1)-gridT(it))
      kapr=kap(ir+1,it)+(kap(ir+1,it+1)-kap(ir+1,it))*
     &     (logt-gridT(it))/(gridT(it+1)-gridT(it))

      kapm=kapl+(kapr-kapl)*(Rkap-gridR(ir))/(gridR(ir+1)-gridR(ir))
      get_kappa=10.**kapm
      
      return
      end
