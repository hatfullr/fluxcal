      subroutine useDimenFile
c     L296 to L344
      include 'optical_depth.h'

c     Check to see if appropriate dimen*.sph already exists.  If so, use
c     it to set the grid size!
      inquire(file=FNAME2,exist=dimenFileAlreadyExists)

      if(dimenFileAlreadyExists) then
         write(*,*) "Reading pre-existing dimen file '" //
     $        trim(adjustl(fname2)) // "'"
         OPEN (73,FILE=FNAME2)
         read(73,*)xminmap,hxmap,nxmap,yminmap,hymap,nymap
         close(73)
         xmaxmap=xminmap+hxmap*(nxmap-1)
         ymaxmap=yminmap+hymap*(nymap-1)
c         print *,'Set grid size according to pre-existing file ',FNAME2
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

      else
       
c     NXMAP and NYMAP are how many grid cells to use in the X and Y
c     directions, respectively:
         NXMAP=NXMAPnext
         NYMAP=NYMAPnext

c     xminmap and xmaxmap are the left and rightmost coordinates of grid
c     yminmap and ymaxmap are the bottom and topmost coordinates of grid
         xminmap=xminmapnext
         xmaxmap=xmaxmapnext
         yminmap=yminmapnext
         ymaxmap=ymaxmapnext

C     Compute cell widths:             
         HXMAP=(XMAXMAP-XMINMAP)/DBLE(NXMAP-1)                           
         HYMAP=(YMAXMAP-YMINMAP)/DBLE(NYMAP-1)
      endif
      end subroutine
