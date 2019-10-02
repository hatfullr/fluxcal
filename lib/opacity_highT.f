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
      real*8 dummy,logR_max,logR_min,logT_max,logT_min

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
      parameter (testsizeT=nrows_opac_table,testsizeR=ncols_opac_table)
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
         nrows_tables(i) = 0
         ncols_tables(i) = 0
         do j=1,nrows_opac_table
            do k=1,ncols_opac_table
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
         read(42,*)form,version,dummy,dummy,ncols_tables(i),
     $        logR_min,logR_max,nrows_tables(i),
     $        logT_min,logT_max
         if ( logTmins(i).le.-1.d30 ) logTmins(i) = logT_min
         if ( logTmaxs(i).le.-1.d30 ) logTmaxs(i) = logT_max
         if ( logRmins(i).le.-1.d30 ) logRmins(i) = logR_min
         if ( logRmaxs(i).le.-1.d30 ) logRmaxs(i) = logR_max
         read(42,*)             ! Newline
         read(42,*)             ! Descriptor
         read(42,*) (logRs(i,k),k=1,ncols_tables(i)) ! This is still R, not rho
         read(42,*)             ! Newline

         if ( ncols_tables(i).gt.ncols_opac_table .or.
     $        nrows_tables(i) .gt. nrows_opac_table ) then
            write(*,*) "i = ",i
            write(*,*) "opacityfiles(i) = '"//
     $           trim(adjustl(opacityfiles(i)))//"'"
            write(*,*) "nrows in opacityfiles(i) = ",nrows_tables(i)
            write(*,*) "ncols in opacityfiles(i) = ",ncols_tables(i)
            write(*,*) "nrows_opac_table = ",nrows_opac_table
            write(*,*) "ncols_opac_table = ",ncols_opac_table
            write(*,*) "The size of the high temperature"//
     $           " opacity table must be greater than the maximum"//
     $           " possible size of all input opacity tables."
            error stop "opacity_highT.f"
         end if
         
         do j=1, nrows_tables(i)
            read(42,*) logTs(i,j),
     $           (logopacity_highT_tables(i,j,k),k=1,ncols_tables(i))
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





c     Don't interpolate up to a higher resolution.

cc     Interpolate each opacityfile up to a universal resolution defined in the code
c         do j=1,nrows_tables(i)
c            index_array_T(j) = float(j)/float(nrows_tables(i))
c         end do
c         do k=1,ncols_tables(i)
c            index_array_R(k) = float(k)/float(ncols_tables(i))
c         end do
c
c
c! We need to interpolate across T to get new Ts for interpolating opacity
c! This is because the original Ts do not have to be linear
c         call interp_linear(1,
c     $        nrows_tables(i),
c     $        index_array_T(1:nrows_tables(i)),
c     $        logTs(i,1:nrows_tables(i)),
c     $        nrows_opac_table,
c     $        index_array_T_int,
c     $        logTarrays(i,:))
c         call interp_linear(1,
c     $        ncols_tables(i),
c     $        index_array_R(1:ncols_tables(i)),
c     $        logRs(i,1:ncols_tables(i)),
c     $        ncols_opac_table,
c     $        index_array_R_int,
c     $        logRarrays(i,:))
c         call rgsf3p(0,
c     $        nrows_tables(i),
c     $        ncols_tables(i),
c     $        logTs(i,1:nrows_tables(i)),
c     $        logRs(i,1:ncols_tables(i)),
c     $        logopacitytables(i,1:nrows_tables(i),1:ncols_tables(i)),
c     $        nrows_opac_table,
c     $        logTarrays(i,:),
c     $        ncols_opac_table,
c     $        logRarrays(i,:),
c     $        logopacity_highT_tables(i,:,:),
c     $        IERR)
c         if(IERR.ne.0) then
c            write(*,*) "IERR = ", IERR
c            error stop "opacity_highT.f: Interpolation with rgsf3p "//
c     $           "encountered an error."
c         end if
         
         numtabused = numtabused + 1
      end do


c      do i=1,numopacityfiles
c         write(*,'("logTmins(",I1,") = ",ES11.4)') i,logTmins(i)
c         write(*,'("logTmaxs(",I1,") = ",ES11.4)') i,logTmaxs(i)
c         write(*,'("logRmins(",I1,") = ",ES11.4)') i,logRmins(i)
c         write(*,'("logRmaxs(",I1,") = ",ES11.4)') i,logRmaxs(i)
c         write(*,*)
c      end do

      
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
c
c      
c
c      do j=1,nrows_opac_table
cc         logT_opac_array(j) = minT+
cc     $        (maxT-minT)/float(nrows_opac_table-1)*(j-1)
c         logT_opac_array(j) = log10(opacity_low_T +
c     $        (opacity_high_T-opacity_low_T)/
c     $        float(nrows_opac_table-1)*(j-1))
c         do k=1,ncols_opac_table
cc            logR_opac_array(k) = minR+
cc     $           (maxR-minR)/float(ncols_opac_table-1)*(k-1)
c            logrho_opac_array(k) = log10(opacity_low_rho +
c     $           (opacity_high_rho-opacity_low_rho)/
c     $           float(ncols_opac_table-1)*(k-1))
c            logR_opac_array(k) = logrho_opac_array(k) +
c     $           3.d0*logT_opac_array(j) - 18.d0
c
c            opacity = 0.d0 ! Calculating this
c            sumweights = 0.d0
c            do i=1,numopacityfiles
c               if(len(trim(adjustl(opacityfiles(i)))).le.0) cycle ! Skip non--existent files
c               if ( .not. (logT_opac_array(j).ge.logTmins(i) .and. ! Skip when R,T isn't within domain
c     $              logT_opac_array(j).le.logTmaxs(i) .and.
c     $              logR_opac_array(k).ge.logRmins(i) .and.
c     $              logR_opac_array(k).le.logRmaxs(i) )) cycle
c
c               ! Catch input file errors
c               if ((logT_opac_array(j).lt.logT_blend1(i) .and.
c     $              logT_opac_array(j).gt.logT_blend2(i)) .or.
c     $             (logR_opac_array(k).lt.logR_blend1(i) .and.
c     $              logR_opac_array(k).gt.logR_blend2(i)))then
c                  write(*,*) "logT, logR = ",logT_opac_array(j),
c     $                 logR_opac_array(k)
c                  write(*,*) "logT_blend1(i) = ",logT_blend1(i)
c                  write(*,*) "logT_blend2(i) = ",logT_blend2(i)
c                  write(*,*) "logR_blend1(i) = ",logR_blend1(i)
c                  write(*,*) "logR_blend2(i) = ",logR_blend2(i)
c                  write(*,*) "Blending regions are not set" //
c     $                 " properly in your input file."
c                  error stop "opacity_highT.f"
c               end if
c               
c               weightT = 1.d0
c               weightR = 1.d0
c
c               ! Runtime errors have been caught, so this should be safe
c               if ( logT_opac_array(j).lt.logT_blend1(i) ) then
c                  weightT = (logT_opac_array(j)-logTmins(i))/
c     $                 (logT_blend1(i)-logTmins(i))
c               else if ( logT_opac_array(j).gt.logT_blend2(i) ) then
c                  weightT = (logTmaxs(i)-logT_opac_array(j))/
c     $                 (logTmaxs(i)-logT_blend2(i))
c               end if
c               
c               if ( logR_opac_array(k).lt.logR_blend1(i) ) then
c                  weightR = (logR_opac_array(k)-logRmins(i))/
c     $                 (logR_blend1(i)-logRmins(i))
c               else if ( logR_opac_array(k).gt.logR_blend2(i) ) then
c                  weightR = (logRmaxs(i)-logR_opac_array(k))/
c     $                 (logRmaxs(i)-logR_blend2(i))
c               end if
c               
c               weight = weightT*weightR
c                  
c               if ( i.gt.2 ) then ! High T opacities
cc                  write(*,*) "weight = ",weight
cc                  opacity = opacity+weight
c                  opacity = opacity+weight*10.d0**bilinear_interpolate(
c     $                 nrows_opac_table,
c     $                 logTarrays(i,:),
c     $                 ncols_opac_table,
c     $                 logRarrays(i,:),
c     $                 logopacity_highT_tables(i,:,:),
c     $                 logT_opac_array(j), logR_opac_array(k))
c               else ! Low T opacities
c                  if ( i.eq.1 ) then ! Planck opacity
c                     rhocgs = 10.d0**(logR_opac_array(k)+
c     $                    3.d0*logT_opac_array(j)-18.d0)
c                     call cop(eD_P,eG_P,rhocgs,10.d0**logT_opac_array(j),
c     $                    tempopacity,verbose)
c                  else if ( i.eq.2 ) then ! Rosseland opacity
c                     rhocgs = 10.d0**(logR_opac_array(k)+
c     $                    3.d0*logT_opac_array(j)-18.d0)
c                     call cop(eD_R,eG_R,rhocgs,10.d0**logT_opac_array(j),
c     $                    tempopacity,verbose)
c                  end if
c                  
c                  opacity = opacity+weight*tempopacity
c                  
c               end if
c               
c               sumweights = sumweights + weight
c               
c            end do ! opacityfiles loop
c
c            logopacity_highT_table(j,k) = log10(opacity/sumweights)
c            
c         end do ! R loop
c      end do ! T loop
c
c


c      do j=1,nrows_opac_table
c         testlogT(j) = log10(500.d0) +
c     $        (8.d0-log10(500.d0))/float(nrows_opac_table-1)*(j-1)
c      end do
c      do k=1,ncols_opac_table
c         testlogR(k) = -13.d0 +
c     $        (3.d0+13.d0)/float(ncols_opac_table-1)*(k-1)
c      end do
c
c      do j=1,nrows_opac_table
c         do k=1,ncols_opac_table
c            testlogrho = testlogR(k) + 3.d0*testlogT(j) - 18.d0
c
c            testlogopacity(j,k)=log10(
c     $           getOpacity(10.d0**testlogT(j),10.d0**testlogrho,0.7d0))
c         end do
c      end do
c      
c      write(fmt1,'("("A","I3"E15.4)")') '"               ", ',
c     $     ncols_opac_table
c      write(fmt2,'("("I3"E15.4)")') ncols_opac_table+1
c
c 400  format("               ",138E15.4)
c 401  format(138E15.4)
c      open(45,file='test_opacity_file.txt',status='unknown')
c      write(45,fmt1) (testlogR(k),k=1,ncols_opac_table)
c      
c      do j=1,nrows_opac_table
c         write(45,fmt2) testlogT(j),
c     $        (testlogopacity(j,k),k=1,ncols_opac_table)
c      end do
c      close(45)

c      stop
      
      end subroutine
