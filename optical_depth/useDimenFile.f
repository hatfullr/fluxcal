      subroutine useDimenFile
c     Check to see if the dimen file for this viewing angle and iteration
c     already exists. If so, read from it. If not, throw an error.
      include 'optical_depth.h'
      character*255 dimenfname
      logical fileexists

 193  format ('dimen',i5.5,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
 194  format ('dimen',i6.6,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
      if(innit.le.99999) then
         write(dimenfname,193) innit,nint(anglez*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
      else
         write(dimenfname,194) innit,nint(anglez*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
      end if

      inquire(file=trim(adjustl(dimenfname)),exist=fileexists)
      if(.not.fileexists) then
         write(o,*) "Could not find dimen file '",
     $        trim(adjustl(dimenfname)),"'."
         write(o,*) "Try setting get_fluxes=.true."
         error stop "useDimenFile.f"
      end if
      
      write(o,*) "Reading pre-existing dimen file '" //
     $     trim(adjustl(dimenfname)) // "'"
      open (73,file=trim(adjustl(dimenfname)))
      read(73,*) xminmap,hxmap,nxmap,yminmap,hymap,nymap
      close(73)
      xmaxmap=xminmap+hxmap*(nxmap-1)
      ymaxmap=yminmap+hymap*(nymap-1)
      if(nxmap.gt.max_Nx) then
         nxmap=max_Nx
         hxmap=(xmaxmap-xminmap)/(nxmap-1)
         print *,'Reset nxmap=',nxmap,' hxmap=',hxmap
      endif
      if(nymap.gt.max_Ny) then
         nymap=max_Ny
         hymap=(ymaxmap-yminmap)/(nymap-1)
         print *,'Reset nymap=',nymap,' hymap=',hymap
      endif
      
      end subroutine
