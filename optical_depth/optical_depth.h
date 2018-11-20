      include '../lib/flux_cal.h'
      real*8 XMINMAP,XMAXMAP,YMINMAP,YMAXMAP
      real*8 XPOS, YPOS, TMAX, TMIN, anglex, angley, anglez
      COMMON/TAUGRID/XPOS,YPOS
      integer mlog
      PARAMETER (MLOG=1)

C NXMAP and NYMAP ARE NOW SET IN THE CODE
      integer NXMAPMAX,NYMAPMAX,NZMAP
c      PARAMETER (NXMAPMAX=999,NYMAPMAX=999,NZMAP=202)
c      real*8 rhoxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      PARAMETER (NXMAPMAX=499,NYMAPMAX=499,NZMAP=416)
      real rhoxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      real uxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      real hpxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      real xhxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      real gxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      real pxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      real txyz(NXMAPMAX,NYMAPMAX,NZMAP)
      integer last_part(NXMAPMAX,NYMAPMAX,NZMAP)
      integer lastpart
      common/lastparticle/ lastpart,last_part

      real*8 tauthin(NXMAPMAX,NYMAPMAX)
      common/thintau/ tauthin
      
      integer thick_part(nxmapmax,nymapmax)
      common/thickparticle/ thick_part
	
      REAL*8 TXY(NXMAPMAX,NYMAPMAX),TphotoXY(NXMAPMAX,NYMAPMAX),
     $     TthermXY(NXMAPMAX,NYMAPMAX),
     $     TpracticalXYthin(NXMAPMAX,NYMAPMAX),
     $	   TpracticalXYthick(NXMAPMAX,NYMAPMAX)
      real*8 zmin(NXMAPMAX,NYMAPMAX),zmax(NXMAPMAX,NYMAPMAX),
     $     h1(NXMAPMAX,NYMAPMAX),zmax_thick(NXMAPMAX,NYMAPMAX)
      integer closest(NXMAPMAX,NYMAPMAX)
      common/close/ closest
      CHARACTER*33 FNAME,FNAME2,fname3,fname4
      real*8 TOTALLUM,TOTALpracticalLUM
      INTEGER kmax,kount,nbad,nok,nrhs,k
      COMMON /countrhs/nrhs
      integer KMAXX,NVARMAX,NVAR
      PARAMETER (KMAXX=200,NVARMAX=50)
      integer ifilter !,numfiltersmax
      real*8 dssav,eps,s,tau,taustart(NVARMAX)
      COMMON /path/ kmax,kount,dssav,s(KMAXX),tau(NVARMAX,KMAXX)
      external derivs2,rkqs,derivs3
      integer NXMAP, NYMAP
      real*8 HXMAP, HYMAP
      integer J, IP
      real*8 rhocgs,xhp,gcgs,ucgs,pcgs,tcgs
      real xh,t6
      common/localQuantities/ rhocgs,xh,t6,xhp,ucgs,gcgs,pcgs,tcgs
      real*8 sa,sb,smid
      real*8 xmin,xmax,ymin,ymax!,xymax
      real*8 xold,yold,zold
      integer iminglow,imaxglow,jminglow,jmaxglow
      real*8 lastTOTALLUM
      real*8 xminmapnext,xmaxmapnext,yminmapnext,ymaxmapnext
      real*8 hxtmp,hytmp
      integer nxmapnext,nymapnext
      integer deltai,deltaj
      integer cellcount
      real*8 rmid,t4avg,tavg,t2avg,vxold,vyold,vzold,vzavg,r2mid
      real*8 costheta,costhetamax,vPCygni,
     $     xcosthetamax,ycosthetamax,zcosthetamax
      real*8 area, sigmat, sigmar,vz2avg,sigmavz,sigmavr
      real*8 localvr,vravg,vr2avg,localvx,localvy,localvz
      integer kount1 !,kount2
      common /kounts/kount1 !,kount2
      real*8 producthigh,productlow,sphoto,stherm,areaphoto,areatherm
      integer ccphoto, cctherm
      real*8 tphoto4avg,ttherm4avg,tpractical4avg
      real*8 tphotoavg,tthermavg
      real*8 konstant,maxdz
      real*8 na,nb,nphoto,nphotoavg
      integer warning,nsl
      parameter(nsl=10)
      real*8 goodtphoto4avg(nsl),Tfit(0:nsl)
      real*8 peak,lambdamax
      integer goodccphoto(nsl),insl
      integer ilogwavelength
      real*8 logwavelength,wavelengthcm,exponent
      integer ilogwavelengthmin,ilogwavelengthmax
      parameter (ilogwavelengthmin=-150,ilogwavelengthmax=150)
      real*8 spectrum(ilogwavelengthmin:ilogwavelengthmax,0:nsl)
      real*8 fluxdensity(ilogwavelengthmin:ilogwavelengthmax,0:nsl)
      integer ilambdamax
      real*8 loglambdamax,ratio

      integer numrho,numtem,numrhodust,numtemdust
      real*8 temtable1,rhotable1,temtabledust1,rhotabledust1
      integer maxtablesize
      parameter(maxtablesize=500)
      real*8 fluxdensityXY(NXMAPMAX,NYMAPMAX,numfiltersmax),
     $     totalfluxdensity(numfiltersmax),mag(numfiltersmax)
      logical dimenFileAlreadyExists
      real*8 rold,phi
      real*8 eexponent,denominator !,coeff
      logical resolved
      common /resolvedboolean/ resolved
      real*8 xposmin,xposmax,yposmin,yposmax,maxdx
      integer imin,imax,jmin,jmax,izmin,izmax,iz,index
      real*8 r2,wpc,zpos,zposmin,zposmax,hzmap
      common/densitygrid/xminmap,yminmap,zmin,zmax,
     $     hxmap,hymap,rhoxyz,uxyz,hpxyz,xhxyz,gxyz,pxyz,txyz,zmax_thick
      real*8 xhi,zz,zzz
      common/metals/zzz

c     Things Roger Hatfull added
      common/subroutine_setviewingAngle/anglex,angley,anglez
      common/subroutine_createGrid/lastTOTALLUM, NXMAPnext, NYMAPnext,
     $     fname, fname2
      common/subroutine_useDimenFile/nxmap,nymap

      common/universal1/ xminmapnext, xmaxmapnext, yminmapnext,
     $     ymaxmapnext,xmaxmap,ymaxmap,
     $     eps,IMINGLOW,JMINGLOW,IMAXGLOW,JMAXGLOW,rmid,r2mid,
     $     costhetamax,vzavg,vz2avg,nphotoavg,warning
      
      common/universal2/ goodtphoto4avg,goodccphoto,ccphoto
      common/universal3/ h1,nok,nbad,TpracticalXYthick,
     $     TpracticalXYthin,fluxdensityXY,
     $     TphotoXY,ilogwavelength
      common/universal4/ logwavelength
      common/universal5/ spectrum, localvx,localvy,localvz,
     $     xcosthetamax,ycosthetamax,zcosthetamax,vravg,vr2avg,TthermXY,
     $     TXY, TOTALLUM,tavg,t2avg,t4avg,tphoto4avg,
     $     tpractical4avg,ttherm4avg,tphotoavg,tthermavg,
     $     totalfluxdensity,TMAX,TMIN,mag,VPCygni,area,areaphoto,
     $     areatherm,sigmat,sigmar,sigmavz,sigmavr,konstant,
     $     dimenFileAlreadyExists
      common/universal6/ xmax,xmin,ymax,ymin,Tfit,ilambdamax
      common/foroutput/ TOTALpracticalLUM
	
      real*8 taus(nmax)
      common/tauparts/ taus

      real*8 nataT(NXMAPMAX,NYMAPMAX),sphT(NXMAPMAX,NYMAPMAX),
     $     nataThi(NXMAPMAX,NYMAPMAX),nataTlo(NXMAPMAX,NYMAPMAX)
      common/splitTs/ nataT, sphT,nataThi,nataTlo

      logical codeerror
      common/errors/ codeerror

      integer intout
      common/intoutput/ intout
