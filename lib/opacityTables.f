      subroutine init_highT_opac(verbose)
c     Here, we are initializing logopacitytables by,
c     1) Reading in the opacityfiles
c     2) Assign values to logopacitytables from opacityfiles
      include 'flux_cal.h'
      integer i,j,k

      integer form, version
      real*8 dummy,logR_max,logR_min,logT_max,logT_min
      integer numtabused
      integer testsizeT,testsizeR
      real*8 testlogTmin,testlogTmax,testlogRmin,testlogRmax
      parameter (testsizeT=nrows_opac_table,testsizeR=ncols_opac_table)
      real*8 testlogT(testsizeT), testlogR(testsizeR)
      real*8 testlogopacity(testsizeT,testsizeR)
      real*8 testlogrho
      real*8 logTmin, logTmax, logRmin, logRmax
      character*255 fmt1,fmt2
      real*8 minT, minR, maxT, maxR
      integer verbose
      
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
               logopacitytables(i,j,k) = -1.d30 ! Zero
            end do
         end do
      end do



c     Need to set the bounds for the first two opacityfiles if they
c     haven't been set yet
      logTmins(1) = log10(500.d0+1.d-8)
      logRmins(1) = (-19.d0+1.d-8) - 3.d0*logTmins(1) + 18.d0
      logTmaxs(2) = 4.d0-1.d-8
      logRmaxs(2) = (-7.d0-1.d-8) - 3.d0*(logTmaxs(2)+1.d-8) + 18.d0

      
      if ( logT_blend2(2).eq.-1.d30 ) logT_blend2(2)=logTmaxs(2)
      if ( logT_blend1(1).eq.-1.d30 ) logT_blend1(1)=logTmins(1)
      if ( logR_blend1(1).eq.-1.d30 ) logR_blend1(1)=logRmins(1)
      if ( logR_blend2(2).eq.-1.d30 ) logR_blend2(2)=logRmaxs(2)

c     Now logT_blend2(2), logT_blend1(1), logR_blend1(1), and
c     logR_blend2(2) have valid values

c     Here, the user has set no preference on transitioning from
c     Planck to Rosseland opacities, so we set a blending by default
      if ( logTmaxs(1).eq.-1.d30 .and. logTmins(2).eq.-1.d30 ) then
         logTmaxs(1) = 3.041393d0
         logTmins(2) = 2.995635d0
         
c     The user has set a minimum for Rosseland opacity, so the
c     maximum for Planck opacity should equal the minimum for
c     Rosseland opacity
      else if ( logTmaxs(1).eq.-1.d30 .and. logTmins(2).ne.-1.d30 ) then
         logTmaxs(1) = logTmins(2)
         
c     The user has set a maximum for Planck opacity, so the
c     minimum for Rosseland opacity should equal the minimum for
c     Planck opacity
      else if ( logTmaxs(1).ne.-1.d30 .and. logTmins(2).eq.-1.d30 ) then
         logTmins(2) = logTmaxs(1)
      end if

c     Please do not use R blending for Planck and Rosseland opacities...
      logRmaxs(1) = (-7.d0-1.d-8) - 3.d0*logT_blend2(2) + 18.d0
      logRmins(2) = (-19.d0+1.d-8) - 3.d0*logT_blend1(1) + 18.d0
      logR_blend1(1) = logRmins(1)
      logR_blend2(1) = logRmaxs(1)
      logR_blend1(2) = logRmins(2)
      logR_blend2(2) = logRmaxs(2)

c     Turning off blending if blending not specified
      if ( logT_blend2(1).eq.-1.d30 ) logT_blend2(1)=logTmaxs(1)
      if ( logT_blend1(2).eq.-1.d30 ) logT_blend1(2)=logTmins(2)
c      if ( logR_blend2(1).eq.-1.d30 ) logR_blend2(1)=logRmaxs(1)
c      if ( logR_blend1(2).eq.-1.d30 ) logR_blend1(2)=logRmins(2)


      logR_blend1(1) = -1.d30
      logR_blend2(1) = -1.d30
      logR_blend1(2) = -1.d30
      logR_blend2(2) = -1.d30
      
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

         if ( nrows_tables(i).gt.nrows_opac_table .or.
     $        ncols_tables(i).gt.ncols_opac_table ) then
            write(*,*) "opacityfiles(i) = ",
     $           trim(adjustl(opacityfiles(i)))
            write(*,*) "nrows_tables(i), ncols_tables(i) = ",
     $           nrows_tables(i),ncols_tables(i)
            write(*,*) "nrows_opac_table, ncols_opac_table = ",
     $           nrows_opac_table, ncols_opac_table
            write(*,*) "Must have nrows_tables(i) <= nrows_opac_table"//
     $           " and ncols_tables(i) <= ncols_opac_table"
            write(*,*) "Try manually increasing the nrows_opac_table "//
     $           "and ncols_opac_table parameters in lib/flux_cal.h "//
     $           "to fit your opacityfiles, or use smaller "//
     $           "opacityfiles."
            error stop "opacityTables.f"
         end if
         
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
            error stop "opacityTables.f"
         end if
         
         do j=1, nrows_tables(i)
            read(42,*) logTs(i,j),
     $           (logopacitytables(i,j,k),k=1,ncols_tables(i))
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
            error stop "opacityTables.f"
         end if
         
         numtabused = numtabused + 1
      end do


      if(numtabused.eq.0) then
         write(*,*) "No opacity tables are being used."
         error stop "opacityTables.f"
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
         error stop "opacityTables.f"
      end if


      ! Automatically set the blending regions that are -1.d30 to the
      ! edges of the file
      do i=1,numopacityfiles
         if(logT_blend1(i).eq.-1.d30) logT_blend1(i)=logTmins(i)
         if(logT_blend2(i).eq.-1.d30) logT_blend2(i)=logTmaxs(i)
         if(logR_blend1(i).eq.-1.d30) logR_blend1(i)=logRmins(i)
         if(logR_blend2(i).eq.-1.d30) logR_blend2(i)=logRmaxs(i)
      end do

      
c$$$c     This is for printing out an opacity file for debugging purposes
c$$$c
c$$$      testlogTmin = log10(200.d0)
c$$$      testlogTmax = 10.d0
c$$$      testlogRmin = -20.d0
c$$$      testlogRmax = 10.d0
c$$$      
c$$$      do j=1,nrows_opac_table
c$$$         testlogT(j) = testlogTmin +
c$$$     $        (testlogTmax-testlogTmin)/float(nrows_opac_table-1)*(j-1)
c$$$      end do
c$$$      do k=1,ncols_opac_table
c$$$         testlogR(k) = testlogRmin +
c$$$     $        (testlogRmax-testlogRmin)/float(ncols_opac_table-1)*(k-1)
c$$$      end do
c$$$
c$$$      do j=1,nrows_opac_table
c$$$         do k=1,ncols_opac_table
c$$$            testlogrho = testlogR(k) + 3.d0*testlogT(j) - 18.d0
c$$$
c$$$            testlogopacity(j,k)=log10(
c$$$     $           getOpacity(10.d0**testlogT(j),10.d0**testlogrho,0.7d0))
c$$$         end do
c$$$      end do
c$$$      
c$$$      write(fmt1,'("("A","I3"E15.4)")') '"               ", ',
c$$$     $     ncols_opac_table
c$$$      write(fmt2,'("("I3"E15.4)")') ncols_opac_table+1
c$$$
c$$$ 400  format("               ",138E15.4)
c$$$ 401  format(138E15.4)
c$$$      open(45,file='test_opacity_file.txt',status='unknown')
c$$$      write(45,fmt1) (testlogR(k),k=1,ncols_opac_table)
c$$$      
c$$$      do j=1,nrows_opac_table
c$$$         write(45,fmt2) testlogT(j),
c$$$     $        (testlogopacity(j,k),k=1,ncols_opac_table)
c$$$      end do
c$$$      close(45)
c$$$
c$$$      write(*,*) "BREAK"
c$$$      write(*,*) ""
c$$$      write(*,*) getOpacity(10.d0**(3.753),10.d0**(4.849+3.*3.753-18),0.7d0)
c$$$      write(*,*) ""
c$$$      write(*,*) "WE ARE OUT"
c$$$      
c$$$      stop
      
      end subroutine
