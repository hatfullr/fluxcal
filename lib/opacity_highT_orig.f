      subroutine init_highT_opac(verbose)
c     Here, we are initializing logopacity_highT_tables by,
c     1) Reading in the opacityfiles
c     2) Interpolating across the opacityfiles to assign values to
c        logopacity_highT_tables (this array is larger than the table)
      use Grid_Interpolation
      use bilinear_interpolation
      include 'flux_cal.h'
      integer i,j,k

      integer form, version
      integer nrows(numopacityfiles), ncols(numopacityfiles)
      real*8 dummy,logR_max,logR_min,logT_max,logT_min

      real*8 logopacitytables(numopacityfiles,nrows_opac_table,
     $     ncols_opac_table)
      real*8 logTs(numopacityfiles,nrows_opac_table),
     $     logRs(numopacityfiles,ncols_opac_table)
      real*8 index_array_T(nrows_opac_table),
     $     index_array_T_int(nrows_opac_table)
      real*8 index_array_R(ncols_opac_table),
     $     index_array_R_int(ncols_opac_table)
      integer IERR
      integer numtabused
      real*8 logTminima(numopacityfiles),
     $     logTmaxima(numopacityfiles),
     $     logRminima(numopacityfiles),
     $     logRmaxima(numopacityfiles)
      integer testsizeT,testsizeR
      parameter (testsizeT=138,testsizeR=37)
      real*8 testlogT(testsizeT), testlogR(testsizeR)
      real*8 testlogopacity(testsizeT,testsizeR)
      real*8 testlogrho
      real*8 logTmin, logTmax, logRmin, logRmax
      real*8 testlogTmin,testlogTmax,testlogRmin,testlogRmax
      character*255 fmt1,fmt2
      real*8 dT, dR
      real*8 minT, minR, maxT, maxR
      real*8 weights(numopacityfiles)
      real*8 Tboundary, Rboundary, Twindow, Rwindow
      real*8 weightT, weightR, weight, sumweights
      real*8 logTarrays(numopacityfiles,nrows_opac_table),
     $     logRarrays(numopacityfiles,ncols_opac_table)
      real*8 opacity,tempopacity,rhocgs

      external cop
      real*8 eD_R, eG_R, eD_P, eG_P
      DIMENSION eD_R(5,6), eD_P(5,6)
      DIMENSION eG_R(71,71), eG_P(71,71)
      common/opac/eD_R,eG_R,eD_P,eG_P

      integer verbose
      

      do j=1,nrows_opac_table
         logT_opac_array(j) = -1.d30
         index_array_T_int(j) = float(j)/float(nrows_opac_table)
      end do
      do k=1,ncols_opac_table
         logrho_opac_array(k) = -1.d30
         index_array_R_int(k) = float(k)/float(ncols_opac_table)
      end do

c     Initialize T and R arrays
      do i=1,numopacityfiles
         do j=1,nrows_opac_table
            logTs(i,j) = -1.d30 ! Zero
         end do
         do k=1, ncols_opac_table
            logRs(i,k) = -1.d30 ! Zero
         end do
      end do
      
c     Initialize the array to hold all the opacity tables
      do i=1,numopacityfiles
         nrows(i) = 0
         ncols(i) = 0
         do j=1,nrows_opac_table
            do k=1,ncols_opac_table
               logopacitytables(i,j,k) = -1.d30 ! Zero
               logopacity_highT_tables(i,j,k) = -1.d30 ! Zero
            end do
         end do
      end do

c     Loop through the opacityfiles, which should all have
c     the format of MESA opacity files (I think OPAL?)

 100  format("logTmins(",I1,"), logTmaxs(",I1,") = ",E10.4,", ",E10.4)
 101  format("logRmins(",I1,"), logRmaxs(",I1,") = ",E10.4,", ",E10.4)
 102  format("    logTmin,     logTmax = ",E10.4,", ",E10.4)
 103  format("    logRmin,     logRmax = ",E10.4,", ",E10.4)


      numtabused = 0
      if(len(trim(adjustl(opacityfiles(1)))).gt.0)
     $     numtabused=numtabused+1
      if(len(trim(adjustl(opacityfiles(2)))).gt.0)
     $     numtabused=numtabused+1
      
      do i=3,numopacityfiles
         if(len(trim(adjustl(opacityfiles(i)))).le.0) cycle
         open(42,file=trim(adjustl(opacityfiles(i))),status='old')
         read(42,*)             ! Descriptor
         read(42,*)             ! Descriptor
         read(42,*)form,version,dummy,dummy,ncols(i),
     $        logR_min,logR_max,nrows(i),
     $        logT_min,logT_max
         if ( logTmins(i).le.-1.d30 ) logTmins(i) = logT_min
         if ( logTmaxs(i).le.-1.d30 ) logTmaxs(i) = logT_max
         if ( logRmins(i).le.-1.d30 ) logRmins(i) = logR_min
         if ( logRmaxs(i).le.-1.d30 ) logRmaxs(i) = logR_max
         read(42,*)             ! Newline
         read(42,*)             ! Descriptor
         read(42,*) (logRs(i,k),k=1,ncols(i)) ! This is still R, not rho
         read(42,*)             ! Newline

         if ( ncols(i).gt.ncols_opac_table .or.
     $        nrows(i) .gt. nrows_opac_table ) then
            write(*,*) "i = ",i
            write(*,*) "opacityfiles(i) = '"//
     $           trim(adjustl(opacityfiles(i)))//"'"
            write(*,*) "nrows in opacityfiles(i) = ",nrows(i)
            write(*,*) "ncols in opacityfiles(i) = ",ncols(i)
            write(*,*) "nrows_opac_table = ",nrows_opac_table
            write(*,*) "ncols_opac_table = ",ncols_opac_table
            write(*,*) "The size of the high temperature"//
     $           " opacity table must be greater than the maximum"//
     $           " possible size of all input opacity tables."
            error stop "opacity_highT.f"
         end if
         
         do j=1, nrows(i)
            read(42,*) logTs(i,j),
     $           (logopacitytables(i,j,k),k=1,ncols(i))
         end do

         close(42)
         
         if ( logTmins(i).lt.logT_min .or.
     $        logTmaxs(i).gt.logT_max .or.
     $        logRmins(i).lt.logR_min .or.
     $        logRmaxs(i).gt.logR_max ) then
            write(*,100) i,i,logTmins(i),logTmaxs(i)
            write(*,101) i,i,logRmins(i),logRmaxs(i)
            write(*,102) logT_min,logT_max
            write(*,103) logR_min,logR_max
            write(*,*) "The specified logT,logR domain is out "//
     $           "of bounds for '",trim(adjustl(opacityfiles(i))),"'."
            error stop "opacity_highT.f"
         end if







c     Interpolate each opacityfile up to a universal resolution defined in the code
         do j=1,nrows(i)
            index_array_T(j) = float(j)/float(nrows(i))
         end do
         do k=1,ncols(i)
            index_array_R(k) = float(k)/float(ncols(i))
         end do


! We need to interpolate across T to get new Ts for interpolating opacity
! This is because the original Ts do not have to be linear
         call interp_linear(1,
     $        nrows(i),
     $        index_array_T(1:nrows(i)),
     $        logTs(i,1:nrows(i)),
     $        nrows_opac_table,
     $        index_array_T_int,
     $        logTarrays(i,:))
         call interp_linear(1,
     $        ncols(i),
     $        index_array_R(1:ncols(i)),
     $        logRs(i,1:ncols(i)),
     $        ncols_opac_table,
     $        index_array_R_int,
     $        logRarrays(i,:))
         call rgsf3p(0,
     $        nrows(i),
     $        ncols(i),
     $        logTs(i,1:nrows(i)),
     $        logRs(i,1:ncols(i)),
     $        logopacitytables(i,1:nrows(i),1:ncols(i)),
     $        nrows_opac_table,
     $        logTarrays(i,:),
     $        ncols_opac_table,
     $        logRarrays(i,:),
     $        logopacity_highT_tables(i,:,:),
     $        IERR)
         if(IERR.ne.0) then
            write(*,*) "IERR = ", IERR
            error stop "opacity_highT.f: Interpolation with rgsf3p "//
     $           "encountered an error."
         end if
         
         numtabused = numtabused + 1
      end do
      
      if(numtabused.eq.0) then
         write(*,*) "No opacity tables are being used."
         error stop "opacity_highT.f"
      end if

      minT = 1.d30
      maxT = -1.d30
      minR = 1.d30
      maxR = -1.d30
      do i=1,numopacityfiles
         if ( logTmins(i).gt.-1.d30.and.logTmins(i).lt.minT )
     $        minT = logTmins(i)
         if ( logTmaxs(i).gt.-1.d30.and.logTmaxs(i).gt.maxT )
     $        maxT = logTmaxs(i)
         if ( logRmins(i).gt.-1.d30.and.logRmins(i).lt.minR )
     $        minR = logRmins(i)
         if ( logRmaxs(i).gt.-1.d30.and.logRmaxs(i).gt.maxR )
     $        maxR = logRmaxs(i)
      end do

      if ( minT.eq.1.d30 .or.
     $     maxT.eq.-1.d30 .or.
     $     minR.eq.1.d30 .or.
     $     maxR.eq.-1.d30 ) then
         write(*,*) "minT,maxT = ",minT,maxT
         write(*,*) "minR,maxR = ",minR,maxR
         write(*,*) "Could not find proper R,T bounds. Check "//
     $        "your domain for your opacityfiles in your input file."
         error stop "opacity_highT.f"
      end if


      ! Automatically set the blending regions that are -1.d30 to the
      ! edges of the file
      do i=1,numopacityfiles
         if(logT_blend1(i).eq.-1.d30) logT_blend1(i)=logTmins(i)
         if(logT_blend2(i).eq.-1.d30) logT_blend2(i)=logTmaxs(i)
         if(logR_blend1(i).eq.-1.d30) logR_blend1(i)=logRmins(i)
         if(logR_blend2(i).eq.-1.d30) logR_blend2(i)=logRmaxs(i)
      end do

      

      do j=1,nrows_opac_table
c         logT_opac_array(j) = minT+
c     $        (maxT-minT)/float(nrows_opac_table-1)*(j-1)
         logT_opac_array(j) = log10(opacity_low_T +
     $        (opacity_high_T-opacity_low_T)/
     $        float(nrows_opac_table-1)*(j-1))
         do k=1,ncols_opac_table
c            logR_opac_array(k) = minR+
c     $           (maxR-minR)/float(ncols_opac_table-1)*(k-1)
            logrho_opac_array(k) = log10(opacity_low_rho +
     $           (opacity_high_rho-opacity_low_rho)/
     $           float(ncols_opac_table-1)*(k-1))
            logR_opac_array(k) = logrho_opac_array(k) +
     $           3.d0*logT_opac_array(j) - 18.d0

            opacity = 0.d0 ! Calculating this
            sumweights = 0.d0
            do i=1,numopacityfiles
               if(len(trim(adjustl(opacityfiles(i)))).le.0) cycle ! Skip non--existent files
               if ( .not. (logT_opac_array(j).ge.logTmins(i) .and. ! Skip when R,T isn't within domain
     $              logT_opac_array(j).le.logTmaxs(i) .and.
     $              logR_opac_array(k).ge.logRmins(i) .and.
     $              logR_opac_array(k).le.logRmaxs(i) )) cycle

               ! Catch input file errors
               if ((logT_opac_array(j).lt.logT_blend1(i) .and.
     $              logT_opac_array(j).gt.logT_blend2(i)) .or.
     $             (logR_opac_array(k).lt.logR_blend1(i) .and.
     $              logR_opac_array(k).gt.logR_blend2(i)))then
                  write(*,*) "logT, logR = ",logT_opac_array(j),
     $                 logR_opac_array(k)
                  write(*,*) "logT_blend1(i) = ",logT_blend1(i)
                  write(*,*) "logT_blend2(i) = ",logT_blend2(i)
                  write(*,*) "logR_blend1(i) = ",logR_blend1(i)
                  write(*,*) "logR_blend2(i) = ",logR_blend2(i)
                  write(*,*) "Blending regions are not set" //
     $                 " properly in your input file."
                  error stop "opacity_highT.f"
               end if
               
               weightT = 1.d0
               weightR = 1.d0

               ! Runtime errors have been caught, so this should be safe
               if ( logT_opac_array(j).lt.logT_blend1(i) ) then
                  weightT = (logT_opac_array(j)-logTmins(i))/
     $                 (logT_blend1(i)-logTmins(i))
               else if ( logT_opac_array(j).gt.logT_blend2(i) ) then
                  weightT = (logTmaxs(i)-logT_opac_array(j))/
     $                 (logTmaxs(i)-logT_blend2(i))
               end if
               
               if ( logR_opac_array(k).lt.logR_blend1(i) ) then
                  weightR = (logR_opac_array(k)-logRmins(i))/
     $                 (logR_blend1(i)-logRmins(i))
               else if ( logR_opac_array(k).gt.logR_blend2(i) ) then
                  weightR = (logRmaxs(i)-logR_opac_array(k))/
     $                 (logRmaxs(i)-logR_blend2(i))
               end if
               
               weight = weightT*weightR
                  
               if ( i.gt.2 ) then ! High T opacities
c                  write(*,*) "weight = ",weight
c                  opacity = opacity+weight
                  opacity = opacity+weight*10.d0**bilinear_interpolate(
     $                 nrows_opac_table,
     $                 logTarrays(i,:),
     $                 ncols_opac_table,
     $                 logRarrays(i,:),
     $                 logopacity_highT_tables(i,:,:),
     $                 logT_opac_array(j), logR_opac_array(k))
               else ! Low T opacities
                  if ( i.eq.1 ) then ! Planck opacity
                     rhocgs = 10.d0**(logR_opac_array(k)+
     $                    3.d0*logT_opac_array(j)-18.d0)
                     call cop(eD_P,eG_P,rhocgs,10.d0**logT_opac_array(j),
     $                    tempopacity,verbose)
                  else if ( i.eq.2 ) then ! Rosseland opacity
                     rhocgs = 10.d0**(logR_opac_array(k)+
     $                    3.d0*logT_opac_array(j)-18.d0)
                     call cop(eD_R,eG_R,rhocgs,10.d0**logT_opac_array(j),
     $                    tempopacity,verbose)
                  end if
                  
                  opacity = opacity+weight*tempopacity
                  
               end if
               
               sumweights = sumweights + weight
               
            end do ! opacityfiles loop


            logopacity_highT_table(j,k) = log10(opacity/sumweights)
            
         end do ! R loop
      end do ! T loop


      write(fmt1,'("("A","I3"E15.4)")') '"               ", ',
     $     ncols_opac_table
      write(fmt2,'("("I3"E15.4)")') ncols_opac_table+1

 400  format("               ",138E15.4)
 401  format(138E15.4)
      open(45,file='test_opacity_file.txt',status='unknown')
c      write(45,fmt1) (testlogR(i),i=1,size(testlogR))
      write(45,fmt1) (logR_opac_array(i),i=1,ncols_opac_table)
      
      
      do j=1,nrows_opac_table
         write(45,fmt2) logT_opac_array(j),
     $        (logopacity_highT_table(j,k),k=1,ncols_opac_table)
      end do
      close(45)


      write(*,'(A,I3," x ",I3)') " Created master opacity table with "//
     $     "resolution (rows, cols) = ",nrows_opac_table,
     $     ncols_opac_table
      
      end subroutine
