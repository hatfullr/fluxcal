      subroutine output(filename,pid)
c     This routine is for writing output to file unit fid for particle
c     pid.
      include 'flux_cal.h'
      character*255 filename
      integer pid
 200  format(13ES22.14,i22)

      if(len(trim(adjustl(filename))).gt.0) then
         open(10,file=trim(adjustl(filename)),action='write',
     $        position='append')
         write(10,200)x(pid)/runit_out,y(pid)/runit_out,z(pid)/runit_out,
     $        am(pid)/munit_out,hp(pid)/runit_out,rho(pid)/rhounit_out,
     $        a(pid)/Eunit_out*munit_out,wmeanmolecular(pid)/muunit_out,
     $        localg(pid)/gunit_out,tempp(pid)/tempunit_out,
     $        Teff(pid)/tempunit_out,pp(pid)/punit_out,tauA(pid),pid
         close(10)
      else
         write(*,200)x(pid)/runit_out,y(pid)/runit_out,z(pid)/runit_out,
     $        am(pid)/munit_out,hp(pid)/runit_out,rho(pid)/rhounit_out,
     $        a(pid)/Eunit_out*munit_out,wmeanmolecular(pid)/muunit_out,
     $        localg(pid)/gunit_out,tempp(pid)/tempunit_out,
     $        Teff(pid)/tempunit_out,pp(pid)/punit_out,tauA(pid),pid
      end if
      
      end subroutine
