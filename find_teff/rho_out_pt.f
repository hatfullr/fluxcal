c     this function creates a table to get back density if T and P are known
      subroutine create_density_array
      include "../lib/flux_cal.h"
      
      real*8 tmin_rho, tmax_rho, pmin_rho, pmax_rho,dt_tabrho, dp_tabrho
      integer nt_tabrho, np_tabrho
      real*8 rho_pt(1000,1000), derr(1000,1000)
      common /rho_lookup/ rho_pt, tmin_rho, tmax_rho,pmin_rho,
     &     pmax_rho, dt_tabrho, dp_tabrho, nt_tabrho, np_tabrho
      real*8 umin, umax, rhomin, rhomax
      real*8 drho, du, dt, dp, tloc, ploc, rholoc, uloc, err_tp
      integer nrho, nu, ip, it, np, nt
      integer  jrho, ju, jp, jt
c     real*8 useeostable

      write(o,*) "Creating envfit density lookup table"

      nt=150
      np=340
      nt_tabrho=nt
      np_tabrho=np
      tmin_rho=3.
      tmax_rho=4.5
      pmin_rho=-4.
      pmax_rho=13.
      dt=(tmax_rho-tmin_rho)/nt
      dp=(pmax_rho-pmin_rho)/np
      dt_tabrho=dt
      dp_tabrho=dp

      rhomin=-15.
      rhomax=0.
      umin=5.
      umax=15.
      nrho=3000
      nu=4000
      drho=(rhomax-rhomin)/nrho
      du=(umax-umin)/nu


      do jt=1,nt
         do jp=1,np
            rho_pt(jp, jt)=-100.
            derr(jp, jt)=100.
         end do
      end do

      do jrho=1,nrho
         do ju=1,nu
            rholoc=10.**(rhomin+drho*(jrho-1))
            uloc=10.**(umin+du*(ju-1))
            tloc=useeostable(uloc,rholoc,1)
            ploc=useeostable(uloc,rholoc,3)
            tloc=log10(tloc)
            ploc=log10(ploc)
            jt=ceiling(1.d0*(tloc-tmin_rho)/dt)
            jp=ceiling(1.d0*(ploc-pmin_rho)/dp)

            if(jt.ge.1.and.jt.le.nt.and.jp.ge.1.and.jp.le.np) then
               err_tp=sqrt( (abs(tloc-tmin_rho-dt*jt)/dt)**2+
     &              (abs(ploc-pmin_rho-dp*jp)/dp)**2)
               if(err_tp.lt.derr(jp,jt)) then
                  derr(jp,jt)=err_tp
                  rho_pt(jp,jt)=rhomin+drho*(jrho-1)
               end if
            end if
         end do
      end do

c      write(o,*) "look up table for density is created"

      return
      end

      

c     this function used useeos to find closest density that is work for this T and P
      real*8 function get_density(p, t)
      
      real*8 tmin_rho, tmax_rho, pmin_rho, pmax_rho,dt_tabrho, dp_tabrho
      integer nt_tabrho, np_tabrho
      real*8 rho_pt(1000,1000)
      common /rho_lookup/ rho_pt, tmin_rho, tmax_rho,pmin_rho,
     &     pmax_rho, dt_tabrho, dp_tabrho, nt_tabrho, np_tabrho

      real*8 rho, p, t, Rgas, a_rad
      integer ip,it, found
      Rgas=8.31e7
      a_rad=7.565e-15

      ip=ceiling((log10(p)-pmin_rho)/dp_tabrho)
      it=ceiling((log10(t)-tmin_rho)/dt_tabrho)

      rho=-40.
      
      if(ip.ge.1.and.ip.le.np_tabrho.and.
     &     it.ge.1.and.it.le.nt_tabrho) then
         rho=rho_pt(ip,it)
      end if

      if(rho.ge.-40) then
         rho=10.**rho
      else
         if(t.le.5000) then
            rho=(p-a_rad*(t**4)/3.)/Rgas*1.293/t
         else
            rho=(p-a_rad*(t**4)/3.)/Rgas*0.6/t            
         end if
      end if
      
      get_density=rho
      
      return
      end
