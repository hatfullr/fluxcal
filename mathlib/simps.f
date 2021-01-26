      subroutine simpson(xpos,ypos,z1,z2,Tthin4,tau_thin,kount1)
      include '../lib/flux_cal.h'
      integer kount1,ncycle,flow,detect_flow
      real*8 z1,z2,dz,Tthin4,tau_thin,opacit
      real*8 dtaudz1,dtaudz2,t1,t2,dtau,dTthin4,frac_Tthin4
      real*8 xpos,ypos,zpos
      real*8 rhocgs,xhp,ucgs,gcgs,pcgs,tcgs
      real xh,t6
      common/localQuantities/ rhocgs,xh,t6, xhp,ucgs,gcgs,pcgs,tcgs
      real*8 z1temp,z2temp,tempdz1,tempdz2
      real start_time,finish_time
      real time_get_dtaudz_T,time_increase_dz,time_decrease_dz
      common/simps_timing/time_get_dtaudz_T,time_increase_dz,
     $     time_decrease_dz
      integer lastpart
      common/lastparticle/ lastpart
      integer intout
      common/intoutput/ intout
      
      external getLocalQuantities
      external get_dtaudz_T
      external increase_dz,decrease_dz

      if ( debug ) then
         write(o,*) "simps.f: Using Simpson's Rule integration "//
     $        "(integrator=1)"
      end if
      
      call cpu_time(start_time)
      time_get_dtaudz_T = 0.
      time_increase_dz = 0.
      time_decrease_dz = 0.
      
      ncycle = 0
      ndtau_increase = 0
      ndF_increase = 0
      ndz_increase = 0
      ndtau_decrease = 0
      ndF_decrease = 0
      ndz_decrease = 0

c     Because we are calculating the integrated flux, we always must
c     ensure that z1 > z2. Otherwise, we get a negative tau, a negative
c     flux, and the magnitude of the flux is not equal to what it would
c     have been if z1 > z2. Negative taus and fluxes make no sense, so
c     we strictly enforce z1 > z2.
      z1temp = z1
      z2temp = z2
      z1 = max(z1temp,z2temp)
      z2 = min(z1temp,z2temp)

      kount1 = 0
      tau_thin = 0.d0
      Tthin4 = 0.d0
      dtau = 0.d0
      dtaudz1 = 0.d0
      dtaudz2 = 0.d0
      t1 = 0.d0
      t2 = 0.d0
      zpos = z1

      call get_dtaudz_T(xpos,ypos,zpos,dtaudz1,t1)
      
      ! Determine the initial stepsize
      dz = (z2-z1)/float(abs(MAXSTP))
      if (simps_min_dz.ne.0 .and. simps_max_dz.ne.0) then
         dz = min(max(dz,sign(dabs(simps_min_dz),dz)),simps_max_dz)
      else if (simps_min_dz.ne.0 .and. simps_max_dz.eq.0) then
         dz = max(dz,sign(dabs(simps_min_dz),dz))
      else if (simps_min_dz.eq.0 .and. simps_max_dz.ne.0) then
         dz = min(dz,sign(dabs(simps_max_dz),dz))
      end if

      if(printIntegrationSteps) then
         call getLocalQuantities(xpos,ypos,zpos)
         opacit = getOpacity(dble(t6*1d6),rhocgs,xh)
         write(intout,'(12ES22.14,I22)')xpos/runit_out,
     $        ypos/runit_out,
     $        zpos/runit_out,xhp/runit_out,rhocgs/rhounit_out,
     $        ucgs/Eunit_out*munit_out,xh/muunit_out,
     $        gcgs/gunit_out,tcgs/tempunit_out,pcgs/punit_out,
     $        opacit,tau_thin,lastpart
      end if
      
      flow=0
      nstp=0                    ! Start at the zero'th step
      do while (ncycle.le.1)
         ! If we are moving out of the integration bounds for z,
         ! make this step the last step
         call detect_last_step(dz,zpos,z2,ndz_increase,ndz_decrease)
c         if ( zpos+dz.lt.z2 ) then
c            if (dabs(dz).lt.dabs(z2-zpos)) then
c               ndz_increase = ndz_increase + 1
c            else if (dabs(dz).gt.dabs(z2-zpos)) then
c               ndz_decrease = ndz_decrease + 1
c            end if
c            ! This signals to the code below that this is the last step
c            dz = z2-zpos
c         end if

         flow = detect_flow(dz,zpos)
         if ( flow.ne.0 ) goto 100 ! Exit main loop and throw error
         
         ! Calculate this step
         call get_dtaudz_T(xpos,ypos,zpos+dz,dtaudz2,t2)
c         write(o,*) "dtaudz2 = ",dtaudz2
         do while (dtaudz2.ge.1.d30 .and. flow.eq.0)
            call decrease_dz(dz,simps_min_dz,simps_alpha,0)
            flow = detect_flow(dz,zpos)
            call get_dtaudz_T(xpos,ypos,zpos+dz,dtaudz2,t2)
c            write(o,*) "We are in here"
         end do
         if ( flow.eq.2 ) then
            write(e,*) "WARNING (simps.f): Stepsize underflow "//
     $           "encountered in a region for which the opacity, "//
     $           "density, or temperature could not be calculated."//
     $           " This means there exists regions that are "//
     $           "outside the provided opacity tables or EOS "//
     $           "tables. The integrator will now stop at "//
     $           "whatever depth it is currently at and assume "//
     $           "all is well. Your final result may not be "//
     $           "physical!"
            flow=0
            goto 100            ! Exit main loop without throwing an error
         end if
         
         dtau = 0.5*(dtaudz1+dtaudz2)*dz
         dTthin4 = 0.5*(t1**4.d0*exp(-tau_thin)+
     $        t2**4.d0*exp(-(tau_thin+dtau)))*dtau

         ! Only enter here if we have taken a step that is allowed by
         ! the dtau limits
         if ( dtau .le. simps_max_dtau .and.
     $        dtau .ge. simps_min_dtau ) then
         
            ! Check against the taulimit
            if ( taulimit.ne.1.d30 .and.
     $           tau_thin+dtau.gt.taulimit-taulimit_threshold ) then
               if ( debug ) then
                  write(o,*) "tau,dtau = ",tau_thin,dtau
                  write(o,*) "tau+dtau > taulimit, refining dz"
               end if
               ! Refine dz if tau is not within the taulimit_threshold
               do while(tau_thin+dtau.gt.taulimit+taulimit_threshold.or.
     $              tau_thin+dtau.lt.taulimit-taulimit_threshold)
                  
                  call get_dtaudz_T(xpos,ypos,zpos+dz,dtaudz2,t2)
                  
                  ! We can only get here if we have already found a 
                  ! dtaudz2 value that is not 1.d30 (error). Our only
                  ! goal in this while loop is to get close to the
                  ! taulimit. Hopefully, changing dz does not move us
                  ! into a region where dtaudz2=1.d30. However, if it
                  ! does, then there is not much that we can do. The
                  ! best option in this case is to just stop the 
                  ! integration where it is, as this will be the closest
                  ! to the taulimit.
                  if ( dtaudz2.ge.1.d30 ) then
                     write(e,*) "WARNING (simps.f): Integration "//
     $                    "stopped in a region where the opacity, "//
     $                    "density, and/or temperature could not "//
     $                    "be calculated when trying to refine the "//
     $                    "stepsize to meet the taulimit requirement."//
     $                    " As a result, the integrator was not able "//
     $                    "to meet the taulimit requirement to within"//
     $                    " the provided taulimit_threshold."
                     goto 100   ! Exit main loop without throwing an error
                  end if
                  
                  dtau = 0.5*(dtaudz1+dtaudz2)*dz
               
                  if ( tau_thin+dtau.gt.taulimit+taulimit_threshold)then
                     ! Need to use 0 as the lower limit here in case we
                     ! need to take a stepsize smaller than that limit.
                     ! The stepsize will depend on taulimit_threshold.
                     ndtaulimit_decrease = ndtaulimit_decrease + 1
                     call decrease_dz(dz,0.d0,simps_alpha,0)
                  else if(tau_thin+dtau.lt.
     $                    taulimit-taulimit_threshold)then
                     ndtaulimit_increase = ndtaulimit_increase + 1
                     call increase_dz(dz,simps_min_dz,simps_alpha,0)
                  end if
                  
                  flow = detect_flow(dz,zpos)
                  if(flow.ne.0) goto 100 ! Exit main loop and throw error
                  
               end do
               
               ! The integration step has been accepted
               ! (we have reached the taulimit)
               dTthin4 = 0.5*(t1**4.d0*exp(-tau_thin)+
     $              t2**4.d0*exp(-(tau_thin+dtau)))*dtau
               ncycle = 0
               Tthin4 = Tthin4 + dTthin4
               tau_thin = tau_thin + dtau
               t2 = t1
               zpos = zpos+dz
               nstp = nstp + 1
               
               if(min_step_size.eq.0.d0)then
                  min_step_size = dz
               else
                  min_step_size = min(min_step_size,dz)
               end if
               if(max_step_size.eq.0.d0)then
                  max_step_size = dz
               else
                  max_step_size = max(max_step_size,dz)
               end if

               if ( debug ) then
                  write(o,*) "Integration exited at the taulimit "//
     $                 "tau,dtau = ",tau_thin,dtau
               end if

               if(printIntegrationSteps) then
                  call getLocalQuantities(xpos,ypos,zpos)
                  opacit = getOpacity(dble(t6*1d6),rhocgs,xh)
                  write(intout,'(12ES22.14,I22)')xpos/runit_out,
     $                 ypos/runit_out,
     $                 zpos/runit_out,xhp/runit_out,rhocgs/rhounit_out,
     $                 ucgs/Eunit_out*munit_out,xh/muunit_out,
     $                 gcgs/gunit_out,tcgs/tempunit_out,pcgs/punit_out,
     $                 opacit,tau_thin,lastpart
               end if
               
               exit             ! Stop integration
            end if
         end if
         
         ! Refine the step
         frac_Tthin4 = dabs(dTthin4/(Tthin4+simps_F_cutoff))
         
         if ( zpos+dz.ne.z2 ) then ! Not the last step
            ! Conditions for decreasing dz
            if ( simps_max_dtau.ne.0 .and.
     $           dtau.gt.simps_max_dtau ) then
               ndtau_decrease = ndtau_decrease+1
               call decrease_dz(dz,simps_min_dz,simps_alpha,ncycle)
               cycle
            else if(simps_max_frac_dF.ne.0 .and.
     $              frac_Tthin4.gt.simps_max_frac_dF) then
               ndF_decrease = ndF_decrease+1
               call decrease_dz(dz,simps_min_dz,simps_alpha,ncycle)
               cycle
            else if(simps_max_dz.ne.0 .and.
     $              dabs(dz).gt.dabs(simps_max_dz)) then
               ndz_decrease = ndz_decrease+1
               call decrease_dz(dz,simps_min_dz,simps_alpha,ncycle)
               cycle
            ! If there's no max limits on anything, never decrease dz
            ! Conditions for increasing dz
            else if((simps_min_dtau.ne.0 .and. ! Min dtau
     $              simps_min_frac_dF.ne.0) .and. ! Min dF/F
     $              (dtau.lt.simps_min_dtau .and.
     $              frac_Tthin4.lt.simps_min_frac_dF)) then
               ndtau_increase = ndtau_increase+1
               ndF_increase = ndF_increase+1
               call increase_dz(dz,simps_max_dz,simps_alpha,ncycle)
               cycle
            else if((simps_min_dtau.ne.0 .and. ! Min dtau
     $              simps_min_frac_dF.eq.0) .and. ! No min dF/F
     $              dtau.lt.simps_min_dtau) then
               ndtau_increase = ndtau_increase+1
               call increase_dz(dz,simps_max_dz,simps_alpha,ncycle)
               cycle
            else if((simps_min_dtau.eq.0 .and. ! No min dtau
     $              simps_min_frac_dF.ne.0) .and. ! Min dF/F
     $              frac_Tthin4.lt.simps_min_frac_dF) then
               ndF_increase = ndF_increase+1
               call increase_dz(dz,simps_max_dz,simps_alpha,ncycle)
               cycle
            else if((simps_min_dtau.eq.0 .and. ! No min dtau
     $              simps_min_frac_dF.eq.0 .and. ! No min dF/F
     $              simps_min_dz.ne.0) .and. ! Min dz
     $              dabs(dz).lt.dabs(simps_min_dz)) then
               ndz_increase = ndz_increase+1
               call increase_dz(dz,simps_max_dz,simps_alpha,ncycle)
               cycle
            end if
            ! If there's no min limits on anything, never increase dz
         end if

         ! The integration step has been accepted
         ncycle = 0
         Tthin4 = Tthin4 + dTthin4
         tau_thin = tau_thin + dtau

         if ( tau_thin.lt.0 ) then
            write(e,*) "dtaudz1 = ",dtaudz1
            write(e,*) "dtaudz2 = ",dtaudz2
            write(e,*) "dtau = ",dtau
            write(e,*) "Somehow we got a negative tau: tau = ",tau_thin
            error stop "simps.f"
         end if

         
         
         ! Take the step forward
         dtaudz1=dtaudz2
         t2=t1
         zpos = zpos+dz
         nstp=nstp+1

         if(printIntegrationSteps) then
            call getLocalQuantities(xpos,ypos,zpos)
            opacit = getOpacity(dble(t6*1d6),rhocgs,xh)
            write(intout,'(12ES22.14,I22)')xpos/runit_out,
     $           ypos/runit_out,
     $           zpos/runit_out,xhp/runit_out,rhocgs/rhounit_out,
     $           ucgs/Eunit_out*munit_out,xh/muunit_out,
     $           gcgs/gunit_out,tcgs/tempunit_out,pcgs/punit_out,
     $           opacit,tau_thin,lastpart
         end if
         
         if(min_step_size.eq.0.d0)then
            min_step_size = dz
         else
            min_step_size = min(min_step_size,dz)
         end if
         if(max_step_size.eq.0.d0)then
            max_step_size = dz
         else
            max_step_size = max(max_step_size,dz)
         end if
c         min_step_size = min(min_step_size,dz)
c         max_step_size = max(max_step_size,dz)

         if ( debug ) then
            write(o,*) "Took integration step, zpos,tau = ",
     $           zpos/runit_out,tau_thin
         end if
         
         ! Check to see if we have detected optical thickness
         if (tau_thin .ge. tau_thick_integrator) kount1 = 1

         
         
         ! Most of the stopping conditions:

         ! If we are outside integration bounds
         ! 1.d-8 is just some small number so that our numerical check
         ! is free from floating point number shenanigans
         if (dabs(zpos-z2)/runit.le.1.d-8) exit
         if ( (z2-z1.gt.0 .and. zpos.gt.z2) .or.
     $        (z2-z1.lt.0 .and. zpos.lt.z2) ) then
            write(e,*)dabs(zpos-z2)/runit
            write(e,*) "This should never happen, but somehow the"//
     $           " integrator managed to get outside the integration"//
     $           " bounds... Good luck!"
            write(e,*) "xpos,ypos,zpos,tau,z1,z2=",
     $        xpos/runit_out,ypos/runit_out,zpos/runit_out,
     $        tau_thin,z1/runit_out,z2/runit_out
            error stop "simps.f"
         end if
         
         ! If we have exceeded the maximum allowed iterations
         if (MAXSTP.gt.0 .and. nstp.ge.abs(MAXSTP)) then
            if (simps_max_step_error) then
               write(o,*) "Reached maximum number of steps "//
     $              "MAXSTP=",MAXSTP,"at xpos,ypos,zpos,tau,z1,z2=",
     $              xpos/runit_out,ypos/runit_out,zpos/runit_out,
     $              tau_thin,z1/runit_out,z2/runit_out
               write(o,*) "You should do one of the following "//
     $              "(in order of safety):"
               write(o,*) "   (1) Increase MAXSTP"
               write(o,*) "   (2) Set MAXSTP to a negative value"
               write(o,*) "   (3) Adjust your integration limits"
               write(o,*) "   (4) Set simps_max_step_error=.false. "//
     $              "in flux_cal.input"
               error stop "simps.f"
            else if (simps_max_step_warning) then
               write(o,*) "WARNING (simps.f): Reached maximum number "//
     $              "integration steps MAXSTP=",MAXSTP,"at "//
     $              "xpos,ypos,zpos,tau,z1,z2=",
     $              xpos/runit_out,ypos/runit_out,zpos/runit_out,
     $              tau_thin,z1/runit_out,z2/runit_out
            end if
            exit
         end if
      end do
 100  continue
      call cpu_time(finish_time)
      
      if ( debug ) then
         write(o,*) ""
         write(o,*) "simps.f: Finished integration"
         write(o,*) "Total time        = ",finish_time-start_time
         write(o,*) "time_get_dtaudz_T = ",time_get_dtaudz_T
         write(o,*) "time_increase_dz  = ",time_increase_dz
         write(o,*) "time_decrease_dz  = ",time_decrease_dz
         write(o,*) "xpos,ypos,zpos = ",
     $        xpos/runit_out,ypos/runit_out,zpos/runit_out
         write(o,*) "z1,z2 = ",z1/runit_out,z2/runit_out
         write(o,*) "nstp = ",nstp
         write(o,*) "ndz_decrease   = ",ndz_decrease
         write(o,*) "ndtau_decrease = ",ndtau_decrease
         write(o,*) "ndF_decrease   = ",ndF_decrease
         write(o,*) "ndz_increase   = ",ndz_increase
         write(o,*) "ndtau_increase = ",ndtau_increase
         write(o,*) "ndF_increase   = ",ndF_increase
         write(o,*) "tau_thin = ",tau_thin
      end if

c      flow = detect_flow(dz,zpos)
      if (flow.eq.1) then
         write(e,*) "xpos,ypos,zpos,dz,tau,z1,z2=",
     $        xpos/runit_out,ypos/runit_out,zpos/runit_out,
     $        dz/runit_out,tau_thin,z1/runit_out,z2/runit_out
         write(e,*) "Stepsize overflow (1/dz = 0)"
         error stop "simps.f"
      else if (flow.eq.2) then
         write(e,*) "xpos,ypos,zpos,dz,tau,z1,z2=",
     $        xpos/runit_out,ypos/runit_out,zpos/runit_out,
     $        dz/runit_out,tau_thin,z1/runit_out,z2/runit_out
         write(e,*) "Stepsize underflow (zpos+dz = zpos)"
         error stop "simps.f"
      end if

      if (ncycle.gt.1) then ! Detect falling out of limits
         write(e,*) "xpos,ypos,zpos,dz,tau,z1,z2=",
     $        xpos/runit_out,ypos/runit_out,zpos/runit_out,
     $        dz/runit_out,tau_thin,z1/runit_out,z2/runit_out
         write(e,*) "Could not increase/decrease the step "//
     $        "size above/below limits simps_max_dz/simps_min_dz."
         write(e,*) "Current dz=",dz
         tempdz1 = dz
         tempdz2 = dz
         call decrease_dz(tempdz1,simps_min_dz,simps_alpha,-1)
         call increase_dz(tempdz2,simps_max_dz,simps_alpha,-1)
         if (tempdz1.le.simps_min_dz) then
            write(e,*) "Desired dz=",tempdz1
         else if (tempdz2.ge.simps_max_dz) then
            write(e,*) "Desired dz=",tempdz2
         end if
         error stop "simps.f"
      end if
      end subroutine

      subroutine get_dtaudz_T(xpos,ypos,zpos,dtaudz,temp)
      implicit none
      real*8 xpos,ypos,zpos,temp,dtaudz,
     $     opacit,getOpacity
      real*8 rhocgs,xhp,ucgs,gcgs,pcgs,tcgs
      real xh,t6
      common/localQuantities/ rhocgs,xh,t6, xhp,ucgs,gcgs,pcgs,tcgs
      external getLocalQuantities
      real start_time,finish_time
      real time_get_dtaudz_T,time_increase_dz,time_decrease_dz
      common/simps_timing/time_get_dtaudz_T,time_increase_dz,
     $     time_decrease_dz
      call cpu_time(start_time)
      call getLocalQuantities(xpos,ypos,zpos)
      temp = t6*1d6
      opacit = getOpacity(temp,rhocgs,xh)
      if ( opacit.eq.-1.d30 ) then
         dtaudz = 1.d30         ! Signal to main code badness is happening
      else
         dtaudz = -opacit*rhocgs
      end if
      call cpu_time(finish_time)
      time_get_dtaudz_T = time_get_dtaudz_T + finish_time-start_time
      end subroutine
      
      subroutine increase_dz(dz,simps_max_dz,simps_alpha,ncycle)
      implicit none
      real*8 dz,simps_max_dz,simps_alpha
      integer ncycle
      real start_time,finish_time
      real time_get_dtaudz_T,time_increase_dz,time_decrease_dz
      common/simps_timing/time_get_dtaudz_T,time_increase_dz,
     $     time_decrease_dz
      logical debug
      common/dbg/debug
      if ( debug ) then
         write(*,*) "Increasing dz"
      end if
      call cpu_time(start_time)
      dz = (1.d0+simps_alpha)*dz
      if (ncycle.lt.0) return   ! -1 for debugging
      if (simps_max_dz.ne.0 .and. dz.ge.simps_max_dz) then
         dz = simps_max_dz
         ncycle = ncycle + 1
      end if
      call cpu_time(finish_time)
      time_increase_dz = time_increase_dz + finish_time-start_time
      end subroutine
      
      subroutine decrease_dz(dz,simps_min_dz,simps_alpha,ncycle)
      implicit none
      real*8 dz,simps_min_dz,simps_alpha
      integer ncycle
      real start_time,finish_time
      real time_get_dtaudz_T,time_increase_dz,time_decrease_dz
      common/simps_timing/time_get_dtaudz_T,time_increase_dz,
     $     time_decrease_dz
      logical debug
      common/dbg/debug
      if ( debug ) then
         write(*,*) "Decreasing dz"
      end if
      call cpu_time(start_time)
      dz = simps_alpha*dz
      if (ncycle.lt.0) return   ! -1 for debugging
      if (simps_min_dz.ne.0 .and. dz.le.simps_min_dz) then
         dz = simps_min_dz
         ncycle = ncycle + 1
      end if
      call cpu_time(finish_time)
      time_decrease_dz = time_decrease_dz + finish_time-start_time
      end subroutine

      integer function detect_flow(dz,zpos)
      include '../lib/flux_cal.h'
      real*8 dz,xpos,ypos,zpos,tau_thin,z1,z2
      detect_flow=0
      if (1.d0/dz.eq.0) then    ! Detect stepsize overflow
         detect_flow=1
      else if (zpos+dz.eq.zpos) then ! Detect stepsize underflow
         detect_flow=2
      end if
      end function

      subroutine detect_last_step(dz,zpos,z2,ndz_increase,ndz_decrease)
      implicit none
      real*8 dz,zpos,z2
      integer ndz_increase,ndz_decrease
      if ( zpos+dz.lt.z2 ) then
         if (dabs(dz).lt.dabs(z2-zpos)) then
            ndz_increase = ndz_increase + 1
         else if (dabs(dz).gt.dabs(z2-zpos)) then
            ndz_decrease = ndz_decrease + 1
         end if
         dz = z2-zpos
      end if
      end subroutine



c     You can use the following commented out code to test thsi integrator
c     outside of FluxCal. Just copy every line with "c$$$" in front and
c     paste into a new document. Then, compile the code with gfortran and
c     run a.out.
c$$$      program main
c$$$      implicit none
c$$$      real*8 xpos,ypos,h
c$$$      common/share/xpos,ypos
c$$$
c$$$      real*8 pi
c$$$      parameter(pi=3.14159265359d0)
c$$$      
c$$$      real*8 z1,z2,zp,zg,hole,dummy,tau_thin
c$$$      integer kount1,nstp,MAXSTP
c$$$      integer ndtau_increase,ndF_increase,ndz_increase,
c$$$     $     ndtau_decrease,ndF_decrease,ndz_decrease,ndtaulimit_decrease,
c$$$     $     ndtaulimit_increase
c$$$      
c$$$      external simpson,getLocalQuantities
c$$$
c$$$      real*8 tau_thick_integrator,taulimit,
c$$$     $     simps_max_dtau,simps_max_frac_dF,taulimit_threshold,
c$$$     $     simps_max_dz,simps_min_dz,simps_min_frac_dF,
c$$$     $     simps_min_dtau,simps_F_cutoff,zpos
c$$$      logical simps_max_step_error,simps_max_step_warning
c$$$      common/share1/MAXSTP,nstp,ndtau_increase,ndF_increase,
c$$$     $     ndz_increase,ndtau_decrease,ndF_decrease,ndz_decrease,
c$$$     $     ndtaulimit_decrease,ndtaulimit_increase
c$$$      common/share2/tau_thick_integrator,taulimit,
c$$$     $     simps_max_step_error,simps_max_step_warning,
c$$$     $     simps_max_dtau,simps_max_frac_dF,taulimit_threshold,
c$$$     $     simps_max_dz,simps_min_dz,simps_min_frac_dF,
c$$$     $     simps_min_dtau,simps_F_cutoff,zpos
c$$$      real*8 sphoto
c$$$
c$$$      real*8 time1,time2
c$$$      integer i
c$$$      real*8 Tthin4,my_tau_thin
c$$$
c$$$      common/intlimits/z1,z2,zp,zg,hole
c$$$      integer trial
c$$$      common/trialno/trial
c$$$
c$$$c     It does not matter what is chosen for z1 and z2, because we
c$$$c     strictly enforce z1 > z2 in the integrator. This is to ensure that
c$$$c     we always get a positive optical zpos and integrated flux.
c$$$c      z1=pi/2.d0                ! Integration start
c$$$c      z2=0.d0                   ! Integration end
c$$$      z1=log(11.d0)             ! For use with trial=2
c$$$      z2=0.d0                   ! For use with trial=2
c$$$      
c$$$      MAXSTP = -2
c$$$      simps_F_cutoff = 1.d-4
c$$$      simps_min_frac_dF = 0
c$$$      simps_max_frac_dF = 0
c$$$      simps_min_dtau = 0
c$$$      simps_max_dtau = 0
c$$$      simps_min_dz = 1.d-2
c$$$      simps_max_dz = 1.d-1
c$$$      taulimit = 0
c$$$      taulimit_threshold = 1.d-4
c$$$      tau_thick_integrator = 1.d1
c$$$      simps_max_step_error = .true.
c$$$      simps_max_step_warning = .false.
c$$$      
c$$$c     both trial=0 and trial=1 use tau=-sin(z)dz, but 1 has a "hole"
c$$$c     centered on (z2+z1)/2 where density is zero. The hole size can be
c$$$c     set using the "hole" variable
c$$$c     trial=2 is for integrating the function rho=e^z. Analytically, if
c$$$c     we want tau=10 to be the result, we can set z2=0 and z1=ln(11).
c$$$      trial = 2
c$$$      hole = dabs(z2-z1)*0.2d0
c$$$      
c$$$c     If trial = 0:
c$$$c     With the settings for kappa and rho below in getLocalQuantities,
c$$$c     we will be doing tau = -int_z1^z2 sin(z) dz and 
c$$$c     Tthin4 = int e^(-tau) dtau = -e^(-tau) + C. We know Tthin4 = 0
c$$$c     when tau=0, so C=1. Thus, Tthin4 = 1 - e^(-tau). Thus, we expect,
c$$$c     tau = cos(z2) - cos(z1),
c$$$c     Tthin4 = 1 - e^[cos(z1)-cos(z2)].
c$$$c      
c$$$c     If trial = 1:
c$$$c     We get the same thing as in trial = 0 except with a hole of size
c$$$c     "hole" centered on (z2-z1)/2. Density is zero in the hole.
c$$$c      
c$$$c     This calculation is done to make it easier to calculate the
c$$$c     expected results later:
c$$$      if (trial.eq.1) then
c$$$         if (z2.lt.z1) then
c$$$            zg = 0.5d0*(z1-z2+hole)
c$$$            zp = 0.5d0*(z1-z2-hole)
c$$$         else
c$$$            zg = 0.5d0*(z2-z1+hole)
c$$$            zp = 0.5d0*(z2-z1-hole)
c$$$         end if
c$$$      end if
c$$$
c$$$      xpos=0.d0                 ! These don't matter
c$$$      ypos=0.d0                 ! These don't matter
c$$$
c$$$
c$$$c     These lines will go in init.f
c$$$c     ----------------------------------------------------
c$$$      if (MAXSTP.eq.0) then
c$$$         write(o,*) "Cannot have MAXSTP = 0"
c$$$         error stop "init.f"
c$$$      end if
c$$$      if (abs(MAXSTP).eq.1) then
c$$$         write(o,*) "Cannot have |MAXSTP| = 1"
c$$$         error stop "init.f"
c$$$      end if
c$$$      if ( (simps_min_dz.ne.0 .and. simps_max_dz.ne.0) .and.
c$$$     $     (dabs(simps_min_dz).gt.dabs(simps_max_dz)) ) then
c$$$         write(o,*) "|simps_min_dz| > |simps_max_dz|"
c$$$         write(o,*) "simps_min_dz = ",simps_min_dz
c$$$         write(o,*) "simps_max_dz = ",simps_max_dz
c$$$         error stop "init.f"
c$$$      end if
c$$$      if ( (simps_min_dtau.ne.0 .and. simps_max_dtau.ne.0) .and.
c$$$     $     (dabs(simps_min_dtau).gt.dabs(simps_max_dtau)) ) then
c$$$         write(o,*) "|simps_min_dtau| > |simps_max_dtau|"
c$$$         write(o,*) "simps_min_dtau = ",simps_min_dtau
c$$$         write(o,*) "simps_max_dtau = ",simps_max_dtau
c$$$         error stop "init.f"
c$$$      end if
c$$$      if ( (simps_min_frac_dF.ne.0 .and. simps_max_frac_dF.ne.0) .and.
c$$$     $     (dabs(simps_min_frac_dF).gt.dabs(simps_max_frac_dF)) ) then
c$$$         write(o,*) "|simps_min_frac_dF| > |simps_max_frac_dF|"
c$$$         write(o,*) "simps_min_frac_dF = ",simps_min_frac_dF
c$$$         write(o,*) "simps_max_frac_dF = ",simps_max_frac_dF
c$$$         error stop "init.f"
c$$$      end if
c$$$      if ( simps_min_frac_dF.ne.0 .and. simps_max_frac_dF.ne.0 .and.
c$$$     $     simps_F_cutoff.eq.0 ) then
c$$$         write(o,*) "Must give a non-zero value for simps_F_cutoff "//
c$$$     $        "when simps_min_frac_dF != 0 and simps_max_frac_dF != 0"
c$$$         error stop "init.f"
c$$$      end if
c$$$      if ( taulimit.ne.1.d30 .and. taulimit_threshold.le.0 ) then
c$$$         write(o,*) "When taulimit=1.d30, you must give a positive, "//
c$$$     $        "non-zero value for taulimit_threshold"
c$$$         error stop "init.f"
c$$$      end if
c$$$c     ----------------------------------------------------
c$$$
c$$$      
c$$$      call cpu_time(time1)
c$$$      call simpson(xpos,ypos,z1,z2,Tthin4,tau_thin,kount1)
c$$$      call cpu_time(time2)
c$$$
c$$$c     Print the reuslts
c$$$      write(o,*) "Integrator results:"
c$$$      write(o,*) "That took ",time2-time1,"seconds"
c$$$      write(o,*) "Final zpos=",zpos
c$$$      write(o,*) "z1,z2 = ",z1,z2
c$$$      write(o,*) "nstp = ",nstp
c$$$      if (kount1.eq.1) then
c$$$         write(o,*) "photosphere detected"
c$$$      else
c$$$         write(o,*) "No photosphere detected"
c$$$      end if
c$$$      write(o,*) "tau_thin = ",tau_thin
c$$$      write(o,*) "Tthin4 = ",Tthin4
c$$$      write(o,*) ""
c$$$      write(o,*) "ndtaulimit_increase = ",ndtaulimit_increase
c$$$      write(o,*) "ndtaulimit_decrease = ",ndtaulimit_decrease
c$$$      write(o,*) ""
c$$$      write(o,*) "ndtau_increase      = ",ndtau_increase
c$$$      write(o,*) "ndtau_decrease      = ",ndtau_decrease
c$$$      write(o,*) ""
c$$$      write(o,*) "ndF_increase        = ",ndF_increase
c$$$      write(o,*) "ndF_decrease        = ",ndF_decrease
c$$$      write(o,*) ""
c$$$      write(o,*) "ndz_increase        = ",ndz_increase
c$$$      write(o,*) "ndz_decrease        = ",ndz_decrease
c$$$      write(o,*) ""
c$$$      write(o,*) "----------------------------------------"
c$$$      write(o,*) ""
c$$$      write(o,*) "Expected analytic results:"
c$$$      if (trial.eq.0) then
c$$$         ! Just need to be careful to make sure we are integrating
c$$$         ! "downward" into the "fluid". That is, so that the integration
c$$$         ! ends at a smaller z than when it started.
c$$$         if (z2.lt.z1) then     ! Start tau=0 at z1 instead of z2
c$$$            my_tau_thin = cos(z2)-cos(z1)
c$$$         else                   ! Start tau=0 at z2
c$$$            my_tau_thin = cos(z1)-cos(z2)
c$$$         end if
c$$$      else if (trial.eq.1) then
c$$$         if (z2.lt.z1) then
c$$$            my_tau_thin = cos(zg)-cos(z1)+cos(z2)-cos(zp)
c$$$         else
c$$$            my_tau_thin = cos(zg)-cos(z2)+cos(z1)-cos(zp)
c$$$         end if
c$$$      else if (trial.eq.2) then
c$$$         if (z2.lt.z1) then
c$$$            my_tau_thin = exp(z1)-exp(z2)
c$$$         else
c$$$            my_tau_thin = exp(z2)-exp(z1)
c$$$         end if
c$$$      end if
c$$$
c$$$      write(o,*) "tau_thin = ",my_tau_thin
c$$$      write(o,*) "Tthin4 = ",1.d0-exp(-my_tau_thin)
c$$$      write(o,*) "% difference in tau_thin = ",
c$$$     $     dabs(my_tau_thin-tau_thin)/my_tau_thin * 100.d0
c$$$      write(o,*) "% difference in Tthin4 = ",
c$$$     $     dabs((1.d0-exp(-my_tau_thin))-Tthin4)/
c$$$     $     (1.d0-exp(-my_tau_thin)) * 100.d0
c$$$      
c$$$
c$$$      end program
c$$$
c$$$      subroutine getLocalQuantities(xpos,ypos,zpos)
c$$$      implicit none
c$$$      real*8 xpos,ypos,zpos
c$$$      real*8 rhocgs,xhp,ucgs,gcgs,pcgs,tcgs
c$$$      real xh,t6
c$$$      common/localQuantities/ rhocgs,xh,t6, xhp,ucgs,gcgs,pcgs,tcgs
c$$$      real*8 z1,z2,zp,zg,hole
c$$$      common/intlimits/z1,z2,zp,zg,hole
c$$$      integer trial
c$$$      common/trialno/trial
c$$$
c$$$      if (trial.eq.0) then
c$$$         rhocgs = sin(zpos)
c$$$      else if (trial.eq.1) then
c$$$         if ( zp.le.zpos .and. zpos.le.zg ) then
c$$$            rhocgs = 0.d0
c$$$         else
c$$$            rhocgs = sin(zpos)
c$$$         end if
c$$$      else if (trial.eq.2) then
c$$$         rhocgs = exp(zpos)
c$$$      end if
c$$$      t6 = 1.d-6
c$$$      
c$$$      end subroutine
c$$$
c$$$
c$$$      real*8 function getOpacity(T,rhocgs,xh)
c$$$      implicit none
c$$$      real*8 T,rhocgs
c$$$      real xh
c$$$      getOpacity = 1.d0
c$$$      end function
c$$$      
c$$$      subroutine simpson(xpos,ypos,z1,z2,Tthin4,tau_thin,kount1)
c$$$      implicit none
c$$$      integer kount1,MAXSTP,nstp,ndtau_increase,ndF_increase,
c$$$     $     ndz_increase,ndtau_decrease,ndF_decrease,ndz_decrease,ncycle,
c$$$     $     ndtaulimit_decrease,ndtaulimit_increase
c$$$      real*8 z1,z2,dz,Tthin4,tau_thin
c$$$      real*8 dtaudz1,dtaudz2,t1,t2,dtau,dTthin4,frac_Tthin4
c$$$      real*8 xpos,ypos,zpos
c$$$      real*8 rhocgs,xhp,ucgs,gcgs,pcgs,tcgs
c$$$      real xh,t6
c$$$      common/localQuantities/ rhocgs,xh,t6, xhp,ucgs,gcgs,pcgs,tcgs
c$$$      real*8 tau_thick_integrator,taulimit,
c$$$     $     simps_max_frac_dF,simps_max_dtau,taulimit_threshold,
c$$$     $     simps_max_dz,simps_min_dz,simps_min_frac_dF,
c$$$     $     simps_min_dtau,simps_F_cutoff
c$$$      logical simps_max_step_error,simps_max_step_warning
c$$$      common/share1/MAXSTP,nstp,ndtau_increase,ndF_increase,
c$$$     $     ndz_increase,ndtau_decrease,ndF_decrease,ndz_decrease,
c$$$     $     ndtaulimit_decrease,ndtaulimit_increase
c$$$      common/share2/tau_thick_integrator,taulimit,
c$$$     $     simps_max_step_error,simps_max_step_warning,
c$$$     $     simps_max_dtau,simps_max_frac_dF,taulimit_threshold,
c$$$     $     simps_max_dz,simps_min_dz,simps_min_frac_dF,
c$$$     $     simps_min_dtau,simps_F_cutoff,zpos
c$$$      real*8 z1temp,z2temp,tempdz1,tempdz2
c$$$      
c$$$      real*8 getOpacity
c$$$      external getLocalQuantities
c$$$      external get_dtaudz_T
c$$$      external increase_dz,decrease_dz
c$$$      
c$$$      ncycle = 0
c$$$      ndtau_increase = 0
c$$$      ndF_increase = 0
c$$$      ndz_increase = 0
c$$$      ndtau_decrease = 0
c$$$      ndF_decrease = 0
c$$$      ndz_decrease = 0
c$$$
c$$$c     Because we are calculating the integrated flux, we always must
c$$$c     ensure that z1 > z2. Otherwise, we get a negative tau, a negative
c$$$c     flux, and the magnitude of the flux is not equal to what it would
c$$$c     have been if z1 > z2. Negative taus and fluxes make no sense, so
c$$$c     we strictly enforce z1 > z2.
c$$$      z1temp = z1
c$$$      z2temp = z2
c$$$      z1 = max(z1temp,z2temp)
c$$$      z2 = min(z1temp,z2temp)
c$$$
c$$$      kount1 = 0
c$$$      tau_thin = 0.d0
c$$$      Tthin4 = 0.d0
c$$$      dtau = 0.d0
c$$$      dtaudz1 = 0.d0
c$$$      dtaudz2 = 0.d0
c$$$      t1 = 0.d0
c$$$      t2 = 0.d0
c$$$      zpos = z1
c$$$
c$$$      call get_dtaudz_T(xpos,ypos,zpos,dtaudz1,t1)
c$$$
c$$$      ! Determine the initial stepsize
c$$$      dz = (z2-z1)/float(abs(MAXSTP))
c$$$      if (simps_min_dz.ne.0 .and. simps_max_dz.ne.0) then
c$$$         dz = min(max(dz,sign(dabs(simps_min_dz),dz)),simps_max_dz)
c$$$      else if (simps_min_dz.ne.0 .and. simps_max_dz.eq.0) then
c$$$         dz = max(dz,sign(dabs(simps_min_dz),dz))
c$$$      else if (simps_min_dz.eq.0 .and. simps_max_dz.ne.0) then
c$$$         dz = min(dz,sign(dabs(simps_max_dz),dz))
c$$$      end if
c$$$      
c$$$      nstp=0                    ! Start at the zero'th step
c$$$c     Prevent stepsize underflow and overflow
c$$$      do while (dz.ne.0 .and. 1.d0/dz.ne.0 .and. ncycle.le.1) 
c$$$         ! If we are moving out of the integration bounds for z,
c$$$         ! make this step the last step
c$$$         if ( zpos+dz.lt.z2 ) then
c$$$            if (dabs(dz).lt.dabs(z2-zpos)) then
c$$$               ndz_increase = ndz_increase + 1
c$$$            else if (dabs(dz).gt.dabs(z2-zpos)) then
c$$$               ndz_decrease = ndz_decrease + 1
c$$$            end if
c$$$            ! This signals to the code below that this is the last step
c$$$            dz = z2-zpos
c$$$         end if
c$$$
c$$$         ! Calculate this step
c$$$         call get_dtaudz_T(xpos,ypos,zpos+dz,dtaudz2,t2)
c$$$         dtau = 0.5*(dtaudz1+dtaudz2)*dz
c$$$         dTthin4 = 0.5*(t1**4.d0*exp(-tau_thin)+
c$$$     $        t2**4.d0*exp(-(tau_thin+dtau)))*dtau
c$$$         
c$$$         ! Check against the taulimit
c$$$         if ( taulimit.ne.1.d30 .and.
c$$$     $        tau_thin+dtau.gt.taulimit-taulimit_threshold ) then
c$$$            ! Refine dz if tau is not within the taulimit_threshold
c$$$            do while(tau_thin+dtau.gt.taulimit+taulimit_threshold .or.
c$$$     $           tau_thin+dtau.lt.taulimit-taulimit_threshold)
c$$$               
c$$$               call get_dtaudz_T(xpos,ypos,zpos+dz,dtaudz2,t2)
c$$$               dtau = 0.5*(dtaudz1+dtaudz2)*dz
c$$$               
c$$$               if (tau_thin+dtau.gt.taulimit+taulimit_threshold) then
c$$$                  ndtaulimit_decrease = ndtaulimit_decrease + 1
c$$$                  ! Need to use 0 as the lower limit here in case we
c$$$                  ! need to take a stepsize smaller than that limit.
c$$$                  ! The stepsize will depend on taulimit_threshold.
c$$$                  call decrease_dz(dz,0.d0,0)
c$$$                  cycle
c$$$               else if(tau_thin+dtau.lt.taulimit-taulimit_threshold)then
c$$$                  ndtaulimit_increase = ndtaulimit_increase + 1
c$$$                  call increase_dz(dz,simps_min_dz,0)
c$$$                  cycle
c$$$               end if
c$$$            end do
c$$$            
c$$$            ! The integration step has been accepted
c$$$            ! (we have reached the taulimit)
c$$$            dTthin4 = 0.5*(t1**4.d0*exp(-tau_thin)+
c$$$     $           t2**4.d0*exp(-(tau_thin+dtau)))*dtau
c$$$            ncycle = 0
c$$$            Tthin4 = Tthin4 + dTthin4
c$$$            tau_thin = tau_thin + dtau
c$$$            t2 = t1
c$$$            zpos = zpos+dz
c$$$            exit                ! Stop integration
c$$$         end if
c$$$         
c$$$         ! Refine the step
c$$$         frac_Tthin4 = dabs(dTthin4/(Tthin4+simps_F_cutoff))
c$$$         
c$$$c         write(o,'(A,7E15.7)')
c$$$c     $        " zpos,tau,Tthin4,dz,dtau,dTthin4,frac_Tthin4=",
c$$$c     $        zpos,tau_thin,Tthin4,dz,dtau,dTthin4,frac_Tthin4
c$$$         
c$$$         if ( zpos+dz.ne.z2 ) then ! Not the last step
c$$$            ! Conditions for decreasing dz
c$$$            if ( simps_max_dtau.ne.0 .and.
c$$$     $           dtau.gt.simps_max_dtau ) then
c$$$               ndtau_decrease = ndtau_decrease+1
c$$$               call decrease_dz(dz,simps_min_dz,ncycle)
c$$$               cycle
c$$$            else if(simps_max_frac_dF.ne.0 .and.
c$$$     $              frac_Tthin4.gt.simps_max_frac_dF) then
c$$$               ndF_decrease = ndF_decrease+1
c$$$               call decrease_dz(dz,simps_min_dz,ncycle)
c$$$               cycle
c$$$            else if(simps_max_dz.ne.0 .and.
c$$$     $              dabs(dz).gt.dabs(simps_max_dz)) then
c$$$               ndz_decrease = ndz_decrease+1
c$$$               call decrease_dz(dz,simps_min_dz,ncycle)
c$$$               cycle
c$$$            ! If there's no max limits on anything, never decrease dz
c$$$            ! Conditions for increasing dz
c$$$            else if((simps_min_dtau.ne.0 .and. ! Min dtau
c$$$     $              simps_min_frac_dF.ne.0) .and. ! Min dF/F
c$$$     $              (dtau.lt.simps_min_dtau .and.
c$$$     $              frac_Tthin4.lt.simps_min_frac_dF)) then
c$$$               ndtau_increase = ndtau_increase+1
c$$$               ndF_increase = ndF_increase+1
c$$$               call increase_dz(dz,simps_max_dz,ncycle)
c$$$               cycle
c$$$            else if((simps_min_dtau.ne.0 .and. ! Min dtau
c$$$     $              simps_min_frac_dF.eq.0) .and. ! No min dF/F
c$$$     $              dtau.lt.simps_min_dtau) then
c$$$               ndtau_increase = ndtau_increase+1
c$$$               call increase_dz(dz,simps_max_dz,ncycle)
c$$$               cycle
c$$$            else if((simps_min_dtau.eq.0 .and. ! No min dtau
c$$$     $              simps_min_frac_dF.ne.0) .and. ! Min dF/F
c$$$     $              frac_Tthin4.lt.simps_min_frac_dF) then
c$$$               ndF_increase = ndF_increase+1
c$$$               call increase_dz(dz,simps_max_dz,ncycle)
c$$$               cycle
c$$$            else if((simps_min_dtau.eq.0 .and. ! No min dtau
c$$$     $              simps_min_frac_dF.eq.0 .and. ! No min dF/F
c$$$     $              simps_min_dz.ne.0) .and. ! Min dz
c$$$     $              dabs(dz).lt.dabs(simps_min_dz)) then
c$$$               ndz_increase = ndz_increase+1
c$$$               call increase_dz(dz,simps_max_dz,ncycle)
c$$$               cycle
c$$$            end if
c$$$            ! If there's no min limits on anything, never increase dz
c$$$         end if
c$$$
c$$$         ! The integration step has been accepted
c$$$         ncycle = 0
c$$$         Tthin4 = Tthin4 + dTthin4
c$$$         tau_thin = tau_thin + dtau
c$$$         
c$$$         ! Take the step forward
c$$$         dtaudz1=dtaudz2
c$$$         t2=t1
c$$$         zpos = zpos+dz
c$$$         nstp=nstp+1
c$$$         
c$$$         ! Check to see if we have detected optical thickness
c$$$         if (tau_thin .ge. tau_thick_integrator) kount1 = 1
c$$$
c$$$         
c$$$         
c$$$         ! Most of the stopping conditions:
c$$$
c$$$         ! If we are outside integration bounds
c$$$         if (dabs(zpos-z2).le.1.d-30) exit
c$$$         if ( (z2-z1.gt.0 .and. zpos.gt.z2) .or.
c$$$     $        (z2-z1.lt.0 .and. zpos.lt.z2) ) then
c$$$            write(o,*) "This should never happen, but somehow the"//
c$$$     $           " integrator managed to get outside the integration"//
c$$$     $           " bounds... Good luck!"
c$$$            error stop "simps.f"
c$$$         end if
c$$$         
c$$$         ! If we have exceeded the maximum allowed iterations
c$$$         if (MAXSTP.gt.0 .and. nstp.ge.abs(MAXSTP)) then
c$$$            if (simps_max_step_error) then
c$$$               write(o,*) "Reached maximum number of steps "//
c$$$     $              "MAXSTP=",MAXSTP,"at xpos,ypos,zpos,tau,z1,z2=",
c$$$     $              xpos,ypos,zpos,tau_thin,z1,z2
c$$$               write(o,*) "You should do one of the following "//
c$$$     $              "(in order of safety):"
c$$$               write(o,*) "   (1) Increase MAXSTP"
c$$$               write(o,*) "   (2) Set MAXSTP to a negative value"
c$$$               write(o,*) "   (3) Adjust your integration limits"
c$$$               write(o,*) "   (4) Set simps_max_step_error=.false. "//
c$$$     $              "in flux_cal.input"
c$$$               error stop "simps.f"
c$$$            else if (simps_max_step_warning) then
c$$$               write(o,*) "WARNING: Reached maximum number "//
c$$$     $              "integration steps MAXSTP=",MAXSTP,"at "//
c$$$     $              "xpos,ypos,zpos,tau,z1,z2=",xpos,ypos,zpos,
c$$$     $              tau_thin,z1,z2
c$$$            end if
c$$$            exit
c$$$         end if
c$$$
c$$$         
c$$$
c$$$      end do
c$$$      
c$$$      if (1.d0/dz.eq.0) then    ! Detect stepsize overflow
c$$$         write(o,*) "xpos,ypos,zpos,tau,z1,z2=",xpos,ypos,zpos,
c$$$     $        tau_thin,z1,z2
c$$$         write(o,*) "Stepsize overflow (1/dz = 0)"
c$$$         error stop "simps.f"
c$$$      else if (dz.eq.0) then    ! Detect stepsize underflow
c$$$         write(o,*) "xpos,ypos,zpos,tau,z1,z2=",xpos,ypos,zpos,
c$$$     $        tau_thin,z1,z2
c$$$         write(o,*) "Stepsize underflow (dz = 0)"
c$$$         error stop "simps.f"
c$$$      else if (ncycle.gt.1) then ! Detect falling out of limits
c$$$         write(o,*) "xpos,ypos,zpos,tau,z1,z2=",xpos,ypos,zpos,
c$$$     $        tau_thin,z1,z2
c$$$         write(o,*) "Could not increase/decrease the step size "//
c$$$     $        "above/below limits simps_max_dz/simps_min_dz."
c$$$         write(o,*) "Current dz=",dz
c$$$         tempdz1 = dz
c$$$         tempdz2 = dz
c$$$         call decrease_dz(tempdz1,simps_min_dz,-1)
c$$$         call increase_dz(tempdz2,simps_max_dz,-1)
c$$$         if (tempdz1.le.simps_min_dz) then
c$$$            write(o,*) "Desired dz=",tempdz1
c$$$         else if (tempdz2.ge.simps_max_dz) then
c$$$            write(o,*) "Desired dz=",tempdz2
c$$$         end if
c$$$         error stop "simps.f"
c$$$      end if
c$$$      end subroutine
c$$$
c$$$      subroutine get_dtaudz_T(xpos,ypos,zpos,dtaudz,T)
c$$$      implicit none
c$$$      real*8 xpos,ypos,zpos,T,dtaudz,
c$$$     $     opacit,getOpacity
c$$$      real*8 rhocgs,xhp,ucgs,gcgs,pcgs,tcgs
c$$$      real xh,t6
c$$$      common/localQuantities/ rhocgs,xh,t6, xhp,ucgs,gcgs,pcgs,tcgs
c$$$      external getLocalQuantities
c$$$      call getLocalQuantities(xpos,ypos,zpos)
c$$$      T = t6*1d6
c$$$      opacit = getOpacity(T,rhocgs,xh)
c$$$      dtaudz = -opacit*rhocgs
c$$$      end subroutine
c$$$      
c$$$      subroutine increase_dz(dz,simps_max_dz,ncycle)
c$$$      implicit none
c$$$      real*8 dz,simps_max_dz
c$$$      integer ncycle
c$$$      write(o,*) "Increase"
c$$$      dz = 1.5d0*dz
c$$$      if (ncycle.lt.0) return   ! -1 for debugging
c$$$      if (simps_max_dz.ne.0 .and. dz.ge.simps_max_dz) then
c$$$         dz = simps_max_dz
c$$$         ncycle = ncycle + 1
c$$$      end if
c$$$      return
c$$$      end subroutine
c$$$
c$$$      subroutine decrease_dz(dz,simps_min_dz,ncycle)
c$$$      implicit none
c$$$      real*8 dz,simps_min_dz
c$$$      integer ncycle
c$$$      write(o,*) "Decrease"
c$$$      dz = 0.5d0*dz
c$$$      if (ncycle.lt.0) return   ! -1 for debugging
c$$$      if (simps_min_dz.ne.0 .and. dz.le.simps_min_dz) then
c$$$         dz = simps_min_dz
c$$$         ncycle = ncycle + 1
c$$$      end if
c$$$      return
c$$$      end subroutine
