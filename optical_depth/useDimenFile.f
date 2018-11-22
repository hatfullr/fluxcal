      subroutine useDimenFile
      include 'optical_depth.h'

c     Check to see if appropriate dimen*.dat already exists.  If so, use
c     it to set the grid size!
c      inquire(file=trim(adjustl(fname2)),exist=dimenFileAlreadyExists)

c      if(dimenFileAlreadyExists) then
      write(*,*) "Reading pre-existing dimen file '" //
     $     trim(adjustl(fname2)) // "'"
      open (73,file=trim(adjustl(fname2)))
      read(73,*)xminmap,hxmap,nxmap,yminmap,hymap,nymap
      close(73)
      xmaxmap=xminmap+hxmap*(nxmap-1)
      ymaxmap=yminmap+hymap*(nymap-1)
      if(nxmap.gt.nxmapmax) then
         nxmap=nxmapmax
         hxmap=(xmaxmap-xminmap)/(nxmap-1)
         print *,'Reset nxmap=',nxmap,' hxmap=',hxmap
      endif
      if(nymap.gt.nymapmax) then
         nymap=nymapmax
         hymap=(ymaxmap-yminmap)/(nymap-1)
         print *,'Reset nymap=',nymap,' hymap=',hymap
      endif
      
      end subroutine
