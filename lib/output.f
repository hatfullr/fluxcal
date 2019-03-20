      subroutine output(filename,pid,binary)
c     This routine is for writing output to file unit fid for particle
c     pid.
      include 'flux_cal.h'
      character*255 filename
      integer pid
      logical binary
      
 200  format(15ES22.14,i22)

      if(len(trim(adjustl(filename))).gt.0) then
         if(binary)then
            open(10,file=trim(adjustl(filename)),action='write',
     $        position='append',form='unformatted')
            write(10)
     $           x(pid)/runit_out,
     $           y(pid)/runit_out,
     $           z(pid)/runit_out,
     $           am(pid)/munit_out,
     $           hp(pid)/runit_out,
     $           rho(pid)/rhounit_out,
     $           a(pid)/Eunit_out*munit_out,
     $           wmeanmolecular(pid)/muunit_out,
     $           localg(pid)/gunit_out,
     $           tempp(pid)/tempunit_out,
     $           Teff(pid)/tempunit_out,
     $           pp(pid)/punit_out,
     $           opac_sph(pid)/kunit_out,
     $           tauA(pid),
     $           entropy(pid)/sunit_out,
     $           pid
            close(10)
         else
            open(10,file=trim(adjustl(filename)),action='write',
     $        position='append')
            write(10,200)
     $           x(pid)/runit_out,
     $           y(pid)/runit_out,
     $           z(pid)/runit_out,
     $           am(pid)/munit_out,
     $           hp(pid)/runit_out,
     $           rho(pid)/rhounit_out,
     $           a(pid)/Eunit_out*munit_out,
     $           wmeanmolecular(pid)/muunit_out,
     $           localg(pid)/gunit_out,
     $           tempp(pid)/tempunit_out,
     $           Teff(pid)/tempunit_out,
     $           pp(pid)/punit_out,
     $           opac_sph(pid)/kunit_out,
     $           tauA(pid),
     $           entropy(pid)/sunit_out,
     $           pid
            close(10)
         endif
      else
         write(*,200)
     $        x(pid)/runit_out,
     $        y(pid)/runit_out,
     $        z(pid)/runit_out,
     $        am(pid)/munit_out,
     $        hp(pid)/runit_out,
     $        rho(pid)/rhounit_out,
     $        a(pid)/Eunit_out*munit_out,
     $        wmeanmolecular(pid)/muunit_out,
     $        localg(pid)/gunit_out,
     $        tempp(pid)/tempunit_out,
     $        Teff(pid)/tempunit_out,
     $        pp(pid)/punit_out,
     $        opac_sph(pid)/kunit_out,
     $        tauA(pid),
     $        entropy(pid)/sunit_out,
     $        pid
      end if
      
      end subroutine
