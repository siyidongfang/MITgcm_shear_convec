#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

C--  File external_forcing.F:
C--   Contents
C--   o EXTERNAL_FORCING_U
C--   o EXTERNAL_FORCING_V
C--   o EXTERNAL_FORCING_T
C--   o EXTERNAL_FORCING_S

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: EXTERNAL_FORCING_U
C     !INTERFACE:
      SUBROUTINE EXTERNAL_FORCING_U(
     I           iMin,iMax, jMin,jMax, bi,bj, kLev,
     I           myTime, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R EXTERNAL_FORCING_U
C     | o Contains problem specific forcing for zonal velocity.
C     *==========================================================*
C     | Adds terms to gU for forcing by external sources
C     | e.g. wind stress, bottom friction etc ...
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     bi,bj     :: Current tile indices
C     kLev      :: Current vertical level index
C     myTime    :: Current time in simulation
C     myThid    :: Thread Id number
      INTEGER iMin, iMax, jMin, jMax, kLev, bi, bj
      _RL myTime
      INTEGER myThid

#ifdef USE_OLD_EXTERNAL_FORCING
C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j       :: Loop counters
C     kSurface  :: index of surface level
      INTEGER i, j
      INTEGER kSurface
CEOP

C=== Yidongfang Si =====
      INTEGER Nshear_smooth_half
      INTEGER NsmoothStart
      INTEGER NsmoothEnd
      REAL*8 shearProfile(Nr)
      REAL*8 shearRatio
      REAL*8 Ndivide

      REAL*8 tauTop
      REAL*8 tauProfile
      REAL*8 dz_const

      REAL*8 amp_tide
      REAL*8 amp_force
      REAL*8 omega_tide

      REAL*8 u_shear
      REAL*8 h_shear
      REAL*8 Nsquare_const
      REAL*8 gravity_const
      REAL*8 temp_Alpha
      REAL*8 Tref_const

      REAL*8 cos_slope
      REAL*8 sin_slope
      REAL*8 topog_slope
      REAL*8 pi_const
C==========================

      IF ( fluidIsAir ) THEN
       kSurface = 0
      ELSEIF ( usingPCoords ) THEN
       kSurface = Nr
      ELSE
       kSurface = 1
      ENDIF


C--   Forcing term
#ifdef ALLOW_AIM
      IF ( useAIM ) CALL AIM_TENDENCY_APPLY_U(
     U                       gU(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_AIM */

#ifdef ALLOW_ATM_PHYS
      IF ( useAtm_Phys ) CALL ATM_PHYS_TENDENCY_APPLY_U(
     U                       gU(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_ATM_PHYS */

#ifdef ALLOW_FIZHI
      IF ( useFIZHI ) CALL FIZHI_TENDENCY_APPLY_U(
     U                       gU(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_FIZHI */

C     Ocean: Add momentum surface forcing (e.g., wind-stress) in surface level
      IF ( kLev .EQ. kSurface ) THEN
c      DO j=1,sNy
C-jmc: Without CD-scheme, this is OK ; but with CD-scheme, needs to cover [0:sNy+1]
       DO j=0,sNy+1
        DO i=1,sNx+1
          gU(i,j,kLev,bi,bj) = gU(i,j,kLev,bi,bj)
     &      +foFacMom*surfaceForcingU(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacW(i,j,kLev,bi,bj)
        ENDDO
       ENDDO
      ELSEIF ( kSurface.EQ.-1 ) THEN
       DO j=0,sNy+1
        DO i=1,sNx+1
         IF ( kSurfW(i,j,bi,bj).EQ.kLev ) THEN
          gU(i,j,kLev,bi,bj) = gU(i,j,kLev,bi,bj)
     &      +foFacMom*surfaceForcingU(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacW(i,j,kLev,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF

#ifdef ALLOW_EDDYPSI
         CALL TAUEDDY_TENDENCY_APPLY_U(
     U                 gU(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
#endif

#ifdef ALLOW_RBCS
      IF (useRBCS) THEN
        CALL RBCS_ADD_TENDENCY(
     U                 gU(1-OLx,1-OLy,kLev,bi,bj),
     I                 kLev, bi, bj, -1,
     I                 myTime, 0, myThid )

      ENDIF
#endif /* ALLOW_RBCS */

#ifdef ALLOW_OBCS
      IF (useOBCS) THEN
        CALL OBCS_SPONGE_U(
     U                   gU(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_OBCS */

#ifdef ALLOW_MYPACKAGE
      IF ( useMYPACKAGE ) THEN
        CALL MYPACKAGE_TENDENCY_APPLY_U(
     U                 gU(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_MYPACKAGE */

C=== Yidongfang Si =====
      omega_tide = 1.454441043328608e-04
      u_shear = 2.1e-03
      h_shear = 250.5
      Nsquare_const = 1.0e-06

C      topog_slope = 4.912113091736391
      topog_slope = 4.d0
      pi_const = 3.141592653589793
      cos_slope = cos(pi_const*topog_slope/180.)
      sin_slope = sin(pi_const*topog_slope/180.)
      gravity_const = 9.81e+00
      temp_Alpha = 2.0e-4
      Tref_const = 0.d0

      amp_tide    = u_shear * h_shear
      amp_force   = amp_tide *
     &     (omega_tide*omega_tide
     &      -f0*f0*cos_slope*cos_slope
     &      -Nsquare_const*sin_slope*sin_slope)
     &                       / (omega_tide*omega_tide)

C      tauTop = 1./3600./100.
C      tauProfile = 0.
C      if((rC(kLev)-rC(Nr)).gt.900.) then
C      tauProfile=tauTop*((rC(kLev)-rC(Nr))-900.)/((rC(1)-rC(Nr))-900.)
C      end if

      dz_const = 3.

      DO i=1,Nr
       if((rC(i)-rC(Nr)).lt.h_shear) then
       shearProfile(i)=(rC(i)-rC(Nr))/h_shear
       else
       shearProfile(i)=1.
       end if
      END DO

      Nshear_smooth_half = 15

      if(kLev.gt.Nshear_smooth_half) then
      if((Nr-kLev).ge.Nshear_smooth_half) then
      NsmoothStart = kLev-Nshear_smooth_half
      NsmoothEnd = kLev+Nshear_smooth_half
      end if
      end if

      if((Nr-kLev).lt.Nshear_smooth_half) then
      NsmoothStart = kLev-(Nr-kLev)
      NsmoothEnd = Nr
      end if

      if(kLev.le.Nshear_smooth_half) then
      NsmoothStart = 1
      NsmoothEnd = kLev+(kLev-1)
      end if

      shearRatio = 0.
      Ndivide = 0.

      if(NsmoothEnd.gt.NsmoothStart) then
      DO i=NsmoothStart,NsmoothEnd
       shearRatio = shearRatio + shearProfile(i)
       Ndivide = Ndivide + 1.
      END DO
      shearRatio = shearRatio/Ndivide
      else
      shearRatio = shearProfile(kLev)
      end if

      DO j=0,sNy+1
       DO i=2,sNx+1
        gU(i,j,kLev,bi,bj) = gU(i,j,kLev,bi,bj)
C     Add tidal forcing
     &  + amp_force*omega_tide*sin(omega_tide*myTime)
     &  *shearRatio ! *cos_slope
C     Add projection of buoyancy production term to u-tendency
     &  + (0.5*(theta(i,j,kLev,bi,bj)+theta(i-1,j,kLev,bi,bj))
     &  -Tref_const)*gravity_const*temp_Alpha*sin_slope
C     Add projection of Coriolis term difference to u-tendency
     &  + (vVel(i,j,kLev,bi,bj)+vVel(i-1,j,kLev,bi,bj))
     &  *f0*0.5*(cos_slope-1.)
C     Add tidal relaxation in sponge layer
C     &  + (amp_tide*cos(omega_tide*myTime)-uVel(i,j,kLev,bi,bj))
C     &  * tauProfile
       ENDDO
C     For zonal grid at the re-entry boundaries
        gU(1,j,kLev,bi,bj) = gU(1,j,kLev,bi,bj)
C     Add tidal forcing
     &  + amp_force*omega_tide*sin(omega_tide*myTime)
     &  *shearRatio ! *cos_slope
C     Add projection of buoyancy production term to u-tendency
     &  + (0.5*(theta(1,j,kLev,bi,bj)+theta(sNx+1,j,kLev,bi,bj))
     &  -Tref_const)*gravity_const*temp_Alpha*sin_slope
C     Add projection of Coriolis term difference to u-tendency
     &  + (vVel(1,j,kLev,bi,bj)+vVel(sNx+1,j,kLev,bi,bj))
     &  *f0*0.5*(cos_slope-1.)
C     Add tidal relaxation in sponge layer
C     &  + (amp_tide*cos(omega_tide*myTime)-uVel(1,j,kLev,bi,bj))
C     &  * tauProfile
      ENDDO
C==========================

#endif /* USE_OLD_EXTERNAL_FORCING */
      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: EXTERNAL_FORCING_V
C     !INTERFACE:
      SUBROUTINE EXTERNAL_FORCING_V(
     I           iMin,iMax, jMin,jMax, bi,bj, kLev,
     I           myTime, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R EXTERNAL_FORCING_V
C     | o Contains problem specific forcing for merid velocity.
C     *==========================================================*
C     | Adds terms to gV for forcing by external sources
C     | e.g. wind stress, bottom friction etc ...
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     bi,bj     :: Current tile indices
C     kLev      :: Current vertical level index
C     myTime    :: Current time in simulation
C     myThid    :: Thread Id number
      INTEGER iMin, iMax, jMin, jMax, kLev, bi, bj
      _RL myTime
      INTEGER myThid

#ifdef USE_OLD_EXTERNAL_FORCING
C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j       :: Loop counters
C     kSurface  :: index of surface level
      INTEGER i, j
      INTEGER kSurface
CEOP

C=== Yidongfang Si =====
      REAL*8 tauTop
      REAL*8 tauProfile

      REAL*8 amp_tide
      REAL*8 omega_tide

      REAL*8 u_shear
      REAL*8 h_shear

      REAL*8 cos_slope
      REAL*8 sin_slope
      REAL*8 topog_slope
      REAL*8 pi_const
C==========================

      IF ( fluidIsAir ) THEN
       kSurface = 0
      ELSEIF ( usingPCoords ) THEN
       kSurface = Nr
      ELSE
       kSurface = 1
      ENDIF

C--   Forcing term
#ifdef ALLOW_AIM
      IF ( useAIM ) CALL AIM_TENDENCY_APPLY_V(
     U                       gV(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_AIM */

#ifdef ALLOW_ATM_PHYS
      IF ( useAtm_Phys ) CALL ATM_PHYS_TENDENCY_APPLY_V(
     U                       gV(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_ATM_PHYS */

#ifdef ALLOW_FIZHI
      IF ( useFIZHI ) CALL FIZHI_TENDENCY_APPLY_V(
     U                       gV(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_FIZHI */

C     Ocean: Add momentum surface forcing (e.g., wind-stress) in surface level
      IF ( kLev .EQ. kSurface ) THEN
       DO j=1,sNy+1
c       DO i=1,sNx
C-jmc: Without CD-scheme, this is OK ; but with CD-scheme, needs to cover [0:sNx+1]
        DO i=0,sNx+1
          gV(i,j,kLev,bi,bj) = gV(i,j,kLev,bi,bj)
     &      +foFacMom*surfaceForcingV(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacS(i,j,kLev,bi,bj)
        ENDDO
       ENDDO
      ELSEIF ( kSurface.EQ.-1 ) THEN
       DO j=1,sNy+1
        DO i=0,sNx+1
         IF ( kSurfS(i,j,bi,bj).EQ.kLev ) THEN
          gV(i,j,kLev,bi,bj) = gV(i,j,kLev,bi,bj)
     &      +foFacMom*surfaceForcingV(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacS(i,j,kLev,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF

#ifdef ALLOW_EDDYPSI
         CALL TAUEDDY_TENDENCY_APPLY_V(
     U                 gV(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
#endif

#ifdef ALLOW_RBCS
      IF (useRBCS) THEN
        CALL RBCS_ADD_TENDENCY(
     U                 gV(1-OLx,1-OLy,kLev,bi,bj),
     I                 kLev, bi, bj, -2,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_RBCS */

#ifdef ALLOW_OBCS
      IF (useOBCS) THEN
        CALL OBCS_SPONGE_V(
     U                   gV(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_OBCS */

#ifdef ALLOW_MYPACKAGE
      IF ( useMYPACKAGE ) THEN
        CALL MYPACKAGE_TENDENCY_APPLY_V(
     U                 gV(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_MYPACKAGE */

C=== Yidongfang Si =====
C      omega_tide = 1.454441043328608e-04
C      u_shear = 1.5e-03
C      h_shear = 250.5

C      topog_slope = 4.912113091736391
      topog_slope = 4.d0
      pi_const = 3.141592653589793e+00
      cos_slope = cos(pi_const*topog_slope/180.)
      sin_slope = sin(pi_const*topog_slope/180.)

C      amp_tide    = u_shear * h_shear * f0/omega_tide

C      tauTop = 1./3600./100.
C      tauProfile = 0.
C      if((rC(kLev)-rC(Nr)).gt.900.) then
C      tauProfile=tauTop*((rC(kLev)-rC(Nr))-900.)/((rC(1)-rC(Nr))-900.)
C      end if

      DO j=0,sNy+1
       DO i=1,sNx
        gV(i,j,kLev,bi,bj) = gV(i,j,kLev,bi,bj)
C     Add projection of Coriolis term difference to v-tendency (along-slope component)
     &   +f0*0.5*(uVel(i,j,kLev,bi,bj)+uVel(i+1,j,kLev,bi,bj))
     &   *(1.-cos_slope)
C     Add tidal relaxation in sponge layer
C     &  + (-amp_tide*sin(omega_tide*myTime)-vVel(i,j,kLev,bi,bj))
C     &  * tauProfile
       ENDDO
C     For zonal grid at the re-entry boundaries
        gV(sNx+1,j,kLev,bi,bj) = gV(sNx+1,j,kLev,bi,bj)
C     Add projection of Coriolis term difference to v-tendency (along-slope component)
     &   +f0*0.5*(uVel(sNx+1,j,kLev,bi,bj)+uVel(1,j,kLev,bi,bj))
     &   *(1.-cos_slope)
C     Add tidal relaxation in sponge layer
C     &  + (-amp_tide*sin(omega_tide*myTime)-vVel(sNx+1,j,kLev,bi,bj))
C     &  * tauProfile
      ENDDO

C     Add projection of Coriolis term to v-tendency (slope-normal component)
C     2D, we do not need to convert v grid to mass grid
      IF ( kLev.GT.1 ) THEN
       DO j=0,sNy+1
        DO i=1,sNx+1
         gV(i,j,kLev,bi,bj) = gV(i,j,kLev,bi,bj)
     &  + f0*0.5*(wVel(i,j,kLev,bi,bj)+wVel(i,j,kLev-1,bi,bj))
     &  *sin_slope
        ENDDO
       ENDDO
      ENDIF
C     Upper boundary grid cell
      IF ( kLev.EQ.1 ) THEN
       DO j=0,sNy+1
        DO i=1,sNx+1
         gV(i,j,kLev,bi,bj) = gV(i,j,kLev,bi,bj)
     &  + f0*0.5*(0.+wVel(i,j,kLev,bi,bj))
     &  *sin_slope
        ENDDO
       ENDDO
      ENDIF


C==========================

#endif /* USE_OLD_EXTERNAL_FORCING */
      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: EXTERNAL_FORCING_T
C     !INTERFACE:
      SUBROUTINE EXTERNAL_FORCING_T(
     I           iMin,iMax, jMin,jMax, bi,bj, kLev,
     I           myTime, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R EXTERNAL_FORCING_T
C     | o Contains problem specific forcing for temperature.
C     *==========================================================*
C     | Adds terms to gT for forcing by external sources
C     | e.g. heat flux, climatalogical relaxation, etc ...
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "SURFACE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     bi,bj     :: Current tile indices
C     kLev      :: Current vertical level index
C     myTime    :: Current time in simulation
C     myThid    :: Thread Id number
      INTEGER iMin, iMax, jMin, jMax, kLev, bi, bj
      _RL myTime
      INTEGER myThid

#ifdef USE_OLD_EXTERNAL_FORCING
C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j       :: Loop counters
C     kSurface  :: index of surface level
      INTEGER i, j
      INTEGER kSurface
      INTEGER km, kc, kp
      _RL tmpVar(1:sNx,1:sNy)
      _RL tmpFac, delPI
      _RL recip_Cp
CEOP

C=== Yidongfang Si ========
      REAL*8 Nsquare_const
      REAL*8 gravity_const
      REAL*8 temp_Alpha

      REAL*8 cos_slope
      REAL*8 sin_slope
      REAL*8 topog_slope
      REAL*8 pi_const
C==========================

#ifdef SHORTWAVE_HEATING
      _RL minusone
      PARAMETER (minusOne=-1.)
      _RL swfracb(2)
      INTEGER kp1
#endif

      IF ( fluidIsAir ) THEN
       kSurface = 0
      ELSEIF ( usingZCoords .AND. useShelfIce ) THEN
       kSurface = -1
      ELSEIF ( usingPCoords ) THEN
       kSurface = Nr
      ELSE
       kSurface = 1
      ENDIF
      recip_Cp = 1. _d 0 / HeatCapacity_Cp

C--   Forcing term
#ifdef ALLOW_AIM
      IF ( useAIM ) CALL AIM_TENDENCY_APPLY_T(
     U                       gT(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_AIM */

#ifdef ALLOW_ATM_PHYS
      IF ( useAtm_Phys ) CALL ATM_PHYS_TENDENCY_APPLY_T(
     U                       gT(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_ATM_PHYS */

#ifdef ALLOW_FIZHI
      IF ( useFIZHI ) CALL FIZHI_TENDENCY_APPLY_T(
     U                       gT(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_FIZHI */

#ifdef ALLOW_ADDFLUID
      IF ( selectAddFluid.NE.0 .AND. temp_addMass.NE.UNSET_RL ) THEN
       IF ( ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &      .OR. convertFW2Salt.EQ.-1. _d 0 ) THEN
         DO j=1,sNy
          DO i=1,sNx
            gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj)
     &        + addMass(i,j,kLev,bi,bj)*mass2rUnit
     &          *( temp_addMass - theta(i,j,kLev,bi,bj) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
C    &          *recip_deepFac2C(kLev)*recip_rhoFacC(kLev)
          ENDDO
         ENDDO
       ELSE
         DO j=1,sNy
          DO i=1,sNx
            gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj)
     &        + addMass(i,j,kLev,bi,bj)*mass2rUnit
     &          *( temp_addMass - tRef(kLev) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
C    &          *recip_deepFac2C(kLev)*recip_rhoFacC(kLev)
          ENDDO
         ENDDO
       ENDIF
      ENDIF
#endif /* ALLOW_ADDFLUID */

#ifdef ALLOW_FRICTION_HEATING
      IF ( addFrictionHeating ) THEN
        IF ( fluidIsAir ) THEN
C         conversion from in-situ Temp to Pot.Temp
          tmpFac = (atm_Po/rC(kLev))**atm_kappa
C         conversion from W/m^2/r_unit to K/s
          tmpFac = (tmpFac/atm_Cp) * mass2rUnit
        ELSE
C         conversion from W/m^2/r_unit to K/s
          tmpFac = recip_Cp * mass2rUnit
        ENDIF
        DO j=1,sNy
          DO i=1,sNx
            gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj)
     &         + frictionHeating(i,j,k,bi,bj)*tmpFac
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
          ENDDO
        ENDDO
      ENDIF
#endif /* ALLOW_FRICTION_HEATING */

      IF ( fluidIsAir .AND. atm_Rq.NE.zeroRL .AND. Nr.NE.1 ) THEN
C--   Compressible fluid: account for difference between moist and dry air
C     specific volume in Enthalpy equation (+ V.dP term), since only the
C     dry air part is accounted for in the (dry) Pot.Temp formulation.
C     Used centered averaging from interface to center (consistent with
C     conversion term in KE eq) and same discretisation ( [T*Q]_bar_k )
C     as for Theta_v in CALC_PHI_HYD

C     conversion from in-situ Temp to Pot.Temp
        tmpFac = (atm_Po/rC(kLev))**atm_kappa
C     conversion from W/kg to K/s
        tmpFac = tmpFac/atm_Cp
        km = kLev-1
        kc = kLev
        kp = kLev+1
        IF ( kLev.EQ.1 ) THEN
          DO j=1,sNy
           DO i=1,sNx
            tmpVar(i,j) = 0.
           ENDDO
          ENDDO
        ELSE
          delPI = atm_Cp*( (rC(km)/atm_Po)**atm_kappa
     &                   - (rC(kc)/atm_Po)**atm_kappa )
          DO j=1,sNy
           DO i=1,sNx
            tmpVar(i,j) = wVel(i,j,kc,bi,bj)*delPI*atm_Rq
     &                  *( theta(i,j,km,bi,bj)*salt(i,j,km,bi,bj)
     &                   + theta(i,j,kc,bi,bj)*salt(i,j,kc,bi,bj)
     &                   )*maskC(i,j,km,bi,bj)*0.25 _d 0
           ENDDO
          ENDDO
        ENDIF
        IF ( kLev.LT.Nr ) THEN
          delPI = atm_Cp*( (rC(kc)/atm_Po)**atm_kappa
     &                   - (rC(kp)/atm_Po)**atm_kappa )
          DO j=1,sNy
           DO i=1,sNx
            tmpVar(i,j) = tmpVar(i,j)
     &                  + wVel(i,j,kp,bi,bj)*delPI*atm_Rq
     &                  *( theta(i,j,kc,bi,bj)*salt(i,j,kc,bi,bj)
     &                   + theta(i,j,kp,bi,bj)*salt(i,j,kp,bi,bj)
     &                   )*maskC(i,j,kp,bi,bj)*0.25 _d 0
           ENDDO
          ENDDO
        ENDIF
        DO j=1,sNy
          DO i=1,sNx
            gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj)
     &         + tmpVar(i,j)*tmpFac
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
          ENDDO
        ENDDO
#ifdef ALLOW_DIAGNOSTICS
        IF ( useDiagnostics ) THEN
C     conversion to W/m^2
          tmpFac = rUnit2mass
          CALL DIAGNOSTICS_SCALE_FILL( tmpVar, tmpFac, 1,
     &                     'MoistCor', kc, 1, 3, bi,bj,myThid )
        ENDIF
#endif /* ALLOW_DIAGNOSTICS */
      ENDIF

C     Ocean: Add temperature surface forcing (e.g., heat-flux) in surface level
      IF ( kLev .EQ. kSurface ) THEN
       DO j=1,sNy
        DO i=1,sNx
          gT(i,j,kLev,bi,bj)=gT(i,j,kLev,bi,bj)
     &      +surfaceForcingT(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
        ENDDO
       ENDDO
      ELSEIF ( kSurface.EQ.-1 ) THEN
       DO j=1,sNy
        DO i=1,sNx
         IF ( kSurfC(i,j,bi,bj).EQ.kLev ) THEN
          gT(i,j,kLev,bi,bj)=gT(i,j,kLev,bi,bj)
     &      +surfaceForcingT(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF

      IF (linFSConserveTr) THEN
       DO j=1,sNy
        DO i=1,sNx
          IF (kLev .EQ. kSurfC(i,j,bi,bj)) THEN
            gT(i,j,kLev,bi,bj)=gT(i,j,kLev,bi,bj)
     &        +TsurfCor*recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
          ENDIF
        ENDDO
       ENDDO
      ENDIF

#ifdef SHORTWAVE_HEATING
C Penetrating SW radiation
c     IF ( usePenetratingSW ) THEN
       swfracb(1)=abs(rF(kLev))
       swfracb(2)=abs(rF(kLev+1))
       CALL SWFRAC(
     I             2, minusOne,
     U             swfracb,
     I             myTime, 1, myThid )
       kp1 = kLev+1
       IF (kLev.EQ.Nr) THEN
        kp1 = kLev
        swfracb(2)=0. _d 0
       ENDIF
       DO j=1,sNy
        DO i=1,sNx
         gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj)
     &   -Qsw(i,j,bi,bj)*(swfracb(1)*maskC(i,j,kLev,bi,bj)
     &                   -swfracb(2)*maskC(i,j,kp1, bi,bj))
     &    *recip_Cp*mass2rUnit
     &    *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
        ENDDO
       ENDDO
c     ENDIF
#endif

#ifdef ALLOW_FRAZIL
      IF ( useFRAZIL )
     &     CALL FRAZIL_TENDENCY_APPLY_T(
     U                 gT(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
#endif /* ALLOW_FRAZIL */

#ifdef ALLOW_SHELFICE
      IF ( useShelfIce )
     &     CALL SHELFICE_FORCING_T(
     U                   gT(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   myTime, 0, myThid )
#endif /* ALLOW_SHELFICE */

#ifdef ALLOW_ICEFRONT
      IF ( useICEFRONT )
     &     CALL ICEFRONT_TENDENCY_APPLY_T(
     U                   gT(1-OLx,1-OLy,kLev,bi,bj),
     I                   kLev, bi, bj, myTime, 0, myThid )
#endif /* ALLOW_ICEFRONT */

#ifdef ALLOW_SALT_PLUME
      IF ( useSALT_PLUME )
     &     CALL SALT_PLUME_TENDENCY_APPLY_T(
     U                     gT(1-OLx,1-OLy,kLev,bi,bj),
     I                     iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                     myTime, 0, myThid )
#endif /* ALLOW_SALT_PLUME */

#ifdef ALLOW_RBCS
      IF (useRBCS) THEN
        CALL RBCS_ADD_TENDENCY(
     U                 gT(1-OLx,1-OLy,kLev,bi,bj),
     I                 kLev, bi, bj, 1,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_RBCS */

#ifdef ALLOW_OBCS
      IF (useOBCS) THEN
        CALL OBCS_SPONGE_T(
     U                   gT(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_OBCS */

#ifdef ALLOW_BBL
      IF ( useBBL ) CALL BBL_TENDENCY_APPLY_T(
     U                       gT(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_BBL */

#ifdef ALLOW_MYPACKAGE
      IF ( useMYPACKAGE ) THEN
        CALL MYPACKAGE_TENDENCY_APPLY_T(
     U                 gT(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_MYPACKAGE */

C=== Yidongfang Si =====
      Nsquare_const = 1.0e-06
      gravity_const = 9.81e+00
      temp_Alpha = 2.0e-04

C      topog_slope = 4.912113091736391
      topog_slope = 4.d0
      pi_const = 3.141592653589793e+00
      cos_slope = cos(pi_const*topog_slope/180.)
      sin_slope = sin(pi_const*topog_slope/180.)



C     Advection by true vertical velocity (cross-slope component)
      DO j=0,sNy+1
       DO i=1,sNx
        gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj) 
     &  - 0.5*(uVel(i,j,kLev,bi,bj)+uVel(i+1,j,kLev,bi,bj))
     &  *Nsquare_const*sin_slope/gravity_const/temp_Alpha
       ENDDO
C     For zonal grid at the re-entry boundaries
        gT(sNx+1,j,kLev,bi,bj) = gT(sNx+1,j,kLev,bi,bj)
     &  - 0.5*(uVel(sNx+1,j,kLev,bi,bj)+uVel(1,j,kLev,bi,bj))
     &  *Nsquare_const*sin_slope/gravity_const/temp_Alpha
      ENDDO

C     Advection by true vertical velocity (slope-normal component)
      IF ( kLev.GT.1 ) THEN
       DO j=0,sNy+1
        DO i=1,sNx+1
         gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj) 
     &  - 0.5*(wVel(i,j,kLev,bi,bj)+wVel(i,j,kLev-1,bi,bj))
     &  *Nsquare_const*cos_slope/gravity_const/temp_Alpha
        ENDDO
       ENDDO
      ENDIF
C     Upper boundary grid cell
      DO j=0,sNy+1
       DO i=1,sNx+1
         gT(i,j,1,bi,bj) = gT(i,j,1,bi,bj)
     &  - 0.5*(0.+wVel(i,j,1,bi,bj))
     &  *Nsquare_const*cos_slope/gravity_const/temp_Alpha
       ENDDO
      ENDDO

C==========================

#endif /* USE_OLD_EXTERNAL_FORCING */
      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: EXTERNAL_FORCING_S
C     !INTERFACE:
      SUBROUTINE EXTERNAL_FORCING_S(
     I           iMin,iMax, jMin,jMax, bi,bj, kLev,
     I           myTime, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R EXTERNAL_FORCING_S
C     | o Contains problem specific forcing for merid velocity.
C     *==========================================================*
C     | Adds terms to gS for forcing by external sources
C     | e.g. fresh-water flux, climatalogical relaxation, etc ...
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "SURFACE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     bi,bj     :: Current tile indices
C     kLev      :: Current vertical level index
C     myTime    :: Current time in simulation
C     myThid    :: Thread Id number
      INTEGER iMin, iMax, jMin, jMax, kLev, bi, bj
      _RL myTime
      INTEGER myThid

#ifdef USE_OLD_EXTERNAL_FORCING
C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j       :: Loop counters
C     kSurface  :: index of surface level
      INTEGER i, j
      INTEGER kSurface
CEOP

      IF ( fluidIsAir ) THEN
       kSurface = 0
      ELSEIF ( usingZCoords .AND. useShelfIce ) THEN
       kSurface = -1
      ELSEIF ( usingPCoords ) THEN
       kSurface = Nr
      ELSE
       kSurface = 1
      ENDIF

C--   Forcing term
#ifdef ALLOW_AIM
      IF ( useAIM ) CALL AIM_TENDENCY_APPLY_S(
     U                       gS(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_AIM */

#ifdef ALLOW_ATM_PHYS
      IF ( useAtm_Phys ) CALL ATM_PHYS_TENDENCY_APPLY_S(
     U                       gS(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_ATM_PHYS */

#ifdef ALLOW_FIZHI
      IF ( useFIZHI ) CALL FIZHI_TENDENCY_APPLY_S(
     U                       gS(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_FIZHI */

#ifdef ALLOW_ADDFLUID
      IF ( selectAddFluid.NE.0 .AND. salt_addMass.NE.UNSET_RL ) THEN
       IF ( ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &      .OR. convertFW2Salt.EQ.-1. _d 0 ) THEN
         DO j=1,sNy
          DO i=1,sNx
            gS(i,j,kLev,bi,bj) = gS(i,j,kLev,bi,bj)
     &        + addMass(i,j,kLev,bi,bj)*mass2rUnit
     &          *( salt_addMass - salt(i,j,kLev,bi,bj) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
C    &          *recip_deepFac2C(kLev)*recip_rhoFacC(kLev)
          ENDDO
         ENDDO
       ELSE
         DO j=1,sNy
          DO i=1,sNx
            gS(i,j,kLev,bi,bj) = gS(i,j,kLev,bi,bj)
     &        + addMass(i,j,kLev,bi,bj)*mass2rUnit
     &          *( salt_addMass - sRef(kLev) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
C    &          *recip_deepFac2C(kLev)*recip_rhoFacC(kLev)
          ENDDO
         ENDDO
       ENDIF
      ENDIF
#endif /* ALLOW_ADDFLUID */

C     Ocean: Add salinity surface forcing (e.g., fresh-water) in surface level
      IF ( kLev .EQ. kSurface ) THEN
       DO j=1,sNy
        DO i=1,sNx
          gS(i,j,kLev,bi,bj)=gS(i,j,kLev,bi,bj)
     &      +surfaceForcingS(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
        ENDDO
       ENDDO
      ELSEIF ( kSurface.EQ.-1 ) THEN
       DO j=1,sNy
        DO i=1,sNx
         IF ( kSurfC(i,j,bi,bj).EQ.kLev ) THEN
          gS(i,j,kLev,bi,bj)=gS(i,j,kLev,bi,bj)
     &      +surfaceForcingS(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF

      IF (linFSConserveTr) THEN
       DO j=1,sNy
        DO i=1,sNx
          IF (kLev .EQ. kSurfC(i,j,bi,bj)) THEN
            gS(i,j,kLev,bi,bj)=gS(i,j,kLev,bi,bj)
     &        +SsurfCor*recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
          ENDIF
        ENDDO
       ENDDO
      ENDIF

#ifdef ALLOW_SHELFICE
      IF ( useShelfIce )
     &     CALL SHELFICE_FORCING_S(
     U                   gS(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   myTime, 0, myThid )
#endif /* ALLOW_SHELFICE */

#ifdef ALLOW_ICEFRONT
      IF ( useICEFRONT )
     &     CALL ICEFRONT_TENDENCY_APPLY_S(
     U                   gS(1-OLx,1-OLy,kLev,bi,bj),
     I                   kLev, bi, bj, myTime, 0, myThid )
#endif /* ALLOW_ICEFRONT */

#ifdef ALLOW_SALT_PLUME
      IF ( useSALT_PLUME )
     &     CALL SALT_PLUME_TENDENCY_APPLY_S(
     U                     gS(1-OLx,1-OLy,kLev,bi,bj),
     I                     iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                     myTime, 0, myThid )
#endif /* ALLOW_SALT_PLUME */

#ifdef ALLOW_RBCS
      IF (useRBCS) THEN
        CALL RBCS_ADD_TENDENCY(
     U                 gS(1-OLx,1-OLy,kLev,bi,bj),
     I                 kLev, bi, bj, 2,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_RBCS */

#ifdef ALLOW_OBCS
      IF (useOBCS) THEN
        CALL OBCS_SPONGE_S(
     U                   gS(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_OBCS */

#ifdef ALLOW_BBL
      IF ( useBBL ) CALL BBL_TENDENCY_APPLY_S(
     U                       gS(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_BBL */

#ifdef ALLOW_MYPACKAGE
      IF ( useMYPACKAGE ) THEN
        CALL MYPACKAGE_TENDENCY_APPLY_S(
     U                 gS(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_MYPACKAGE */

#endif /* USE_OLD_EXTERNAL_FORCING */
      RETURN
      END
