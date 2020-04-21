      subroutine simpson(xpos,ypos,z1,z2,Tthin4,tau_thin,kount1)
c     In this Simpson's rule integrator, we want to be able to stop where
c     tau_thin ~ tau_limit. To do so, we will use a modified version of the
c     Composite Simpson's rule, which is normally written as
c
c     S = \int_a^b f(z) dz = h/3 \sum_{i=0}^m (f_{2i+2} + 4f_{2i+1} + f_{2i})
c
c     We modify it such that if we detect that tau_thin > tau_limit after
c     adding the i'th contribution, we undo that addition and instead reduce
c     the stepsize by some factor alpha, and then continue integration. That is,
c
c     S =   \sum_{i=0}^m_1   h/3         (f_{2i+2} + 4f_{2i+1} + f_{2i})
c         + \sum_{i=m_1}^m_2 h/3 alpha^1 (f_{2i+2} + 4f_{2i+1} + f_{2i})
c         + \sum_{i=m_2}^m_3 h/3 alpha^2 (f_{2i+2} + 4f_{2i+1} + f_{2i})
c         + ...
c         + \sum_{i=m_{k-1}}^m_k h/3 alpha^{k-1} (f_{2i+2} + 4f_{2i+1} + f_{2i})
c
c     S = h/3 \sum_{k=1} [ \sum_{i=m_{k-1}}^m_k alpha^{k-1} (f_{2i+2} + 4f_{2i+1} + f_{2i}) ]
c       = h/3 \sum_{k=1} [ alpha^{k-1} \sum_{i=m_{k-1}}^m_k (f_{2i+2} + 4f_{2i+1} + f_{2i}) ]
c     
c     where m_k is the number of iterations required to reach tau_thin > tau_limit
c     using the stepsize h*alpha^{k-1}, m_0 = 0, and alpha is any constant with
c     alpha < 1 giving convergent results and alpha = 1 giving exactly the same
c     result as normal Simpson's rule integration.
c
c     The error in the result is numerically unattainable.
c
c     (c) Roger Hatfull, Apr. 2020. University of Alberta

      include '../lib/flux_cal.h'
      integer k,kount1
      real*8 z1,z2,h,new_h,Tthin4,tau_thin,mystep,new_step,
     $     tau_max_acceptable,tau_min_acceptable,ho3alphak,t6d6
      real*8 opacit, rhocgs, simps_alphak, ho3,ho3simps_alphak
      real*8 dtaudz1, dtaudz2, dtaudz3, Tthin4_1, Tthin4_2, Tthin4_3
      real*8 xpos,ypos,zpos,xhp
      real xh,t6
      common/localQuantities/ rhocgs,xh,t6,xhp
      
      external getLocalQuantities
c      write(*,*) ""
      if ( simps_h .eq. 0.d0 ) then
         if ( MAXSTP .lt. 2 ) then
            write(*,*) "ERROR: MAXSTP = ",MAXSTP," but "//
     $           "must be >= 2."
            error stop "simps.f"
         end if
         h = 0.5d0*abs(z2-z1)/float(MAXSTP-2)
      else
         h = simps_h
         ! Prevent stepsizes that give a total number of steps less than 2
         if ( abs(h) .gt. abs(0.5d0*(z2-z1)) ) then
            write(*,*) "WARNING (simps.f): Integration stepsize |h| >"//
     $           " |(z2-z1)/2|, so only one step would be taken. |h| "//
     $           "will be set to |(z2-z1)/2|, but this will result in"//
     $           " inaccurate results."
            h = sign(h,abs(0.5d0*(z2-z1)))
         end if
      end if
      
      if ( simps_alpha.ge.1.d0 .or. simps_alpha.le.0.d0 ) then
         write(*,*) "ERROR (simps.f): Must have stepsize modifier "//
     $        "0 < simps_alpha < 1."
         error stop "simps.f"
      end if

c     Make sure sign is correct on stepsizes
      if ( (z2.lt.z1 .and. h.gt.0) .or.
     $     (z2.gt.z1 .and. h.lt.0) ) then
         h = -h
      end if

      ho3 = h/3.d0

c     Find what the maximum and minimum acceptable tau values are
      tau_max_acceptable = taulimit*(1. + simps_fracacc)
      tau_min_acceptable = taulimit*(1. - simps_fracacc)

      mystep = h
      zpos = z1+mystep            ! Start integration at z1+h

      nstp = 1                 ! We are "on" the first iteration
      tau_thin = 0.d0
      Tthin4 = 0.d0
      kount1 = 0
      
c     Calculate the first integration step
      call getLocalQuantities(xpos,ypos,zpos-mystep)
      t6d6 = t6*1d6
      opacit = getOpacity(t6d6,rhocgs,xh)
      dtaudz1 = -opacit*rhocgs*ho3
      Tthin4_1 = (t6d6)**4.d0*exp(-dtaudz1)*dtaudz1

c     Prevent min_step_size from being the default 0
      if ( min_step_size.eq.0.d0 ) then
         min_step_size = abs(mystep)
      end if
      
      k=1
      do while (nstp .lt. MAXSTP)
         simps_alphak = simps_alpha**(k-1)
         ho3simps_alphak = ho3*simps_alphak
         
         new_step = mystep * simps_alphak
         ! Each time we reduce the step size, we need to move zpos into place
         zpos = zpos - mystep + new_step
         mystep = new_step
         min_step_size = min(min_step_size,abs(mystep))
         max_step_size = max(max_step_size,abs(mystep))
         
         ! Catch an edge case where we could infinitely loop if the stepsize is zero.
         if ( mystep .eq. 0 ) then
            write(*,*) "Failed to find a small enough stepsize"
            error stop "simps.f"
         end if

         do while (nstp .lt. MAXSTP)
            
            ! If we have reached the upper bound on integration,
            if ( (zpos+mystep .lt. z2) .or.
     $           (zpos .le. z2) ) then
c               write(*,*) "Final iteration."
               ! This will be the final iteration
               zpos = 0.5d0*(zpos-mystep+z2)
               mystep = z2-zpos
               min_step_size = min(min_step_size,abs(mystep))
               max_step_size = max(max_step_size,abs(mystep))
            end if

            ! Must add to tau_thin here, as Tthin4 depends on it
            tau_thin = tau_thin + dtaudz1
            
            call getLocalQuantities(xpos,ypos,zpos)
            t6d6 = t6*1d6
            opacit = getOpacity(t6d6,rhocgs,xh)
            dtaudz2 = -4.d0*opacit*rhocgs*ho3simps_alphak
            tau_thin = tau_thin + dtaudz2
            Tthin4_2 = (t6d6)**4.d0*exp(-tau_thin)*dtaudz2
            
            call getLocalQuantities(xpos,ypos,zpos+mystep)
            t6d6 = t6*1d6
            opacit = getOpacity(t6d6,rhocgs,xh)
            dtaudz3 = -opacit*rhocgs*ho3simps_alphak
            tau_thin = tau_thin + dtaudz3
            Tthin4_3 = (t6d6)**4.d0*exp(-tau_thin)*dtaudz3

c     For debugging:
c            write(*,'(2I5,4E15.7,"     ",4E15.7,"     ",3E15.7)')
c     $           nstp,k,zpos,mystep,zpos+mystep,z2,
c     $           dtaudz1,dtaudz2,dtaudz3,tau_thin,
c     $           taulimit,tau_min_acceptable,tau_max_acceptable

            ! Reject this integration step if tau is going to exceed the
            ! taulimit. We add dtaudz3 here to check for if tau will 
            ! exceed taulimit at the beginning of the next iteration. 
            ! If we don't check this way, it is possible for the tau to
            ! exceed taulimit such that even a stepsize of 0 cannot 
            ! reduce tau below taulimit.
            if ( tau_thin+dtaudz3 .gt. tau_max_acceptable ) then
               ! Integration step rejected. Undo this step.
               tau_thin = tau_thin - (dtaudz1+dtaudz2+dtaudz3)
               exit             ! Return to the main loop
            end if
            
            ! Integration step accepted (completed an iteration)
            nstp = nstp + 1
            
            Tthin4 = Tthin4 + (Tthin4_1+Tthin4_2+Tthin4_3)
            
            ! Store these values for the next iteration
            dtaudz1 = dtaudz3
            Tthin4_1 = Tthin4_3

            ! Terminate integration if...
            if ( (zpos+mystep .eq. z2) .or. ! Integration bound reached
     $           (tau_thin .ge. tau_min_acceptable) ) then ! tau_thin ~ taulimit
               goto 100         ! Exit the loops
            end if
            
            zpos = zpos + 2.d0*mystep

            ! If we are about to enter a zero-density region at zpos+mystep,
            ! increase the stepsize to search for the next non-zero density
            ! region.
            if ( rhocgs .eq. 0.d0 ) then
               ! k=1 gives the user's original stepsize. We don't want to take steps
               ! larger than that so as to not compromise the user's control of
               ! minimum integration accuracy.
               k = 0            ! 1 gets added right after
               exit
            end if
            
         end do
         k = k+1
      end do
      
 100  continue
      ! Indicate if we have integrated beyond the definition of optical thickness
      if ( tau_thin .gt. tau_thick_integrator ) kount1 = 1

      if ( nstp.gt.MAXSTP ) then
         write(*,*) "Exceeded MAXSTP during integration for " //
     $        "xpos, ypos = ",xpos,ypos
         write(*,*) "Maximum depth achieved = ",zpos
         write(*,*) "Optical depth = ",tau_thin
         write(*,*) "nstp = ",nstp
         write(*,*) "Try increasing simps_h, MAXSTP, "//
     $        "simps_alpha, or simps_fracacc. You can also try using "//
     $        "a smaller taulimit."
         error stop "simps.f"
      end if

      return
      end subroutine
