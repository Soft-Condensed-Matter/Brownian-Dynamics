!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%***************************************************************************%%
!%%***************************************************************************%%
!%%**  PROGRAM      BROWNIAN DYNAMICS                                       **%%
!%%**  AUTHOR       ALEXIS TORRES CARBAJAL - ALPIXELS                       **%%
!%%**  LICENSE      LGPL-V3                                                 **%%
!%%**                                                                       **%%
!%%**  INTERACTION  LENNARD-JONES                                           **%%
!%%**  ALGORITHM    ERMACK-MACKAMON                                         **%%
!%%**  DATE         OCTOBER 24, 2022                                        **%%
!%%**                                                                       **%%
!%%**  OBSERVATION  THIS CODE IS DEVELOPED FOR THE IV INTERNATIONAL SCHOOL  **%%
!%%**               OF ENGEENIERING OF MATTER OUT OF EQUILIBIRUM - 2022     **%%
!%%***************************************************************************%%
!%%***************************************************************************%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN MODULE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE BDVAR
 IMPLICIT NONE
 INTEGER, PARAMETER:: D       = KIND(1.0)   !NUMERICAL PRECISION
 INTEGER, PARAMETER:: MP      = 1400        !MAXIMUM NUMBER OF PARTICLES       
 INTEGER, PARAMETER:: OUTDISP = 10000       !SHOW DATA ON DISPLAY
 INTEGER, PARAMETER:: OUTSAMP = 1000        !TAKE A STATISTICAL SAMPLE
 INTEGER, PARAMETER:: NBINX   = 20000       !MAXIMUM NUMBER OF PARTICLES
 INTEGER, PARAMETER:: OUTNSAM = 10000       !TIME CORRELATION
 INTEGER, PARAMETER:: IT0     = 100         !TAKE A NEW t=0 FOR CORRELATION
 INTEGER, PARAMETER:: T0MAX   = 100000      !MAXIMUM NUMBER OF t=0
 INTEGER, PARAMETER:: TMAXD   = 100000      !MAXIMUM TIME OF OBSERVATION

 REAL(D), PARAMETER:: SIGMA   = 1.0         !PARTICLES DIAMETER
 REAL(D), PARAMETER:: MASS    = 1.0         !PARTICLES MASS
 REAL(D), PARAMETER:: H       = 0.00001     !TIME STEP
 REAL(D), PARAMETER:: PI      = ACOS(-1.0)  !PI NUMBER
 REAL(D), PARAMETER:: DBIN    = 0.01        !BIN WIDTH FOR STATISTICS
 REAL(D), PARAMETER:: TOBS    = 500.0        !OBSERVATION TIME

 INTEGER:: BDSTEP,NPART,IDUMM
 INTEGER:: SSTAT,IFRAME,NBIN,ISAM
 INTEGER:: HGR(NBINX),IBD,SWITCH
 INTEGER:: NTEL,T0,TT0,TMAX
 INTEGER:: NTIME(TMAXD),TIME0(T0MAX),DELT

 REAL(D):: RCUT,USHIFT,EPOT
 REAL(D):: CUT,PHI,RHOSTAR
 REAL(D):: RX(MP),RY(MP)
 REAL(D):: FX(MP),FY(MP)
 REAL(D):: BOXX,BOXY,ESTAR
 REAL(D):: IBOXX,IBOXY,TSTAR
 REAL(D):: MSDX(TMAXD),MSDY(TMAXD)
 REAL(D):: RX0(MP,TMAXD),RY0(MP,TMAXD)
END MODULE BDVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM BROWNIANDYNAMICS
 USE BDVAR
 IMPLICIT NONE

 CALL BEGIN                                 !BEGIN SIMULATION

 DO IBD=1,BDSTEP
    CALL INTEGRA                            !POSITIONS @ t + ðt
    CALL BDFORCE                            !FORCES @ t + ðt
    CALL DYNAMIC                            !SHOW DATA AND COMPUTE OBSERVABLES
 ENDDO

 CALL FINIS                                 !END SIMULATION

 STOP
END PROGRAM BROWNIANDYNAMICS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SUBROUTINES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE BEGIN
 IMPLICIT NONE

 CALL READINP                               !READ SIMULATION PARAMETERS
 CALL ICONFIG                               !BUILD INITIAL CONFIGURATION
 CALL SETUPDB                               !SET UP SIMULATION PARAMETERS
 CALL CORRELA                               !SET UP CORRELATION FUNCTIONS
 CALL BDFORCE                               !COMPUTE FORCE ON EACH PARTICLE
 CALL SHOWSIM                               !SHOW SIMULATION CONDITIONS

 OPEN(UNIT=12,FILE="BDStatus.dat")          !SAVE SIMULATION EVOLUTION
 RETURN
END SUBROUTINE BEGIN
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE READINP
 USE BDVAR
 IMPLICIT NONE

 OPEN(UNIT=10,FILE="BD.inp",STATUS="OLD")

 READ(10,*)NPART                            !NUMBER OF CELLS
 READ(10,*)PHI                              !PACKING FRACTION
 READ(10,*)TSTAR                            !SYSTEM TEMPERATURE
 READ(10,*)CUT                              !CUT-OFF RADIUS
 READ(10,*)BDSTEP                           !SIMULATION STEPS
 READ(10,*)SSTAT                            !BEGIN STATISTICS

 CLOSE(UNIT=10)

 RETURN 
END SUBROUTINE READINP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE ICONFIG
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: I
 REAL(D):: DEL

 RHOSTAR=PHI*4.0/PI                         !SYSTEM DENSITY NUMBER
 NPART=INT(SQRT(REAL(NPART)))*INT(SQRT(REAL(NPART)))
 BOXX=SQRT(REAL(NPART)*SIGMA**2/RHOSTAR)    !SIMULATION BOX LENGTH (SQUARE DOMAIN)
 DEL=SQRT(1.0/RHOSTAR)                      !MEAN SEPARATION DISTANCE
 
 RX(1)=(-BOXX + DEL)*0.5                    !FIRST PARTICLE X POSITION
 RY(1)=(-BOXX + DEL)*0.5                    !FIRST PARTICLE Y POSITION

 DO I=2,NPART                               !SIMPLE SQUARE LATTICE
    RX(I)=RX(I-1) + DEL
    RY(I)=RY(I-1)
    IF(RX(I) .GT. 0.5*BOXX)THEN
      RX(I)=RX(1)
      RY(I)=RY(I-1) + DEL
    ENDIF
 ENDDO

 OPEN(UNIT=11,FILE="BDSnap.xyz")            !INITIAL CONFIGURATION SNAPSHOT

 WRITE(11,*)NPART
 WRITE(11,*)'FRAME',1
 
 DO I=1,NPART
    WRITE(11,111)'C',RX(I),RY(I),0.0
 ENDDO

 CLOSE(UNIT=11)
 111 FORMAT(A,3(F15.8))

 RETURN
END SUBROUTINE ICONFIG
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE SETUPDB
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: I
 REAL(D):: RC6

 CALL SYSTEM_CLOCK(IDUMM)                   !(RANDOM!) SEED ACCORDING TIME MACHINE
 T0=0                                       !NUMBER OF to
 NBIN=INT(BOXX*0.5/DBIN)                    !NUMBER OF BIN IN SIMULATION
 RCUT=CUT*SIGMA                             !CUT-OFF RADIUS
 ISAM=0                                     !NUMBER OF STATISTICAL SAMPLES
 IDUMM=-IDUMM                               !SEED FOR RANDON NUMBER GENERATOR
 SWITCH=0                                   !SET UP CORRELATION FUNCTIONS
 RC6=(SIGMA/RCUT)**6                        !ATTRACTIVE POTENTIAL TERM AT RCUT
 USHIFT=4.0*ESTAR*RC6*(RC6 - 1.0)           !POTENTIAL SHIFT TERM 
 ESTAR=1.0/TSTAR                            !INTERACTION WELL DEPTH BETWEEEN PARTICLES
 IFRAME=0                                   !NUMBER OF FRAME FOR MOVIE
 BOXY=BOXX                                  !Y SIMULATION BOX LENGTH
 IBOXX=1.0/BOXX                             !X INVERSE SIMULATION LENGTH
 IBOXY=1.0/BOXY                             !Y INVERSE SIMULATION LENGTH

 DO I=1,NBIN
    HGR(I)=0                                !INITIALIZATION HISTOGRAM GR
 ENDDO

 RETURN
END SUBROUTINE SETUPDB
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE SHOWSIM
 USE BDVAR
 IMPLICIT NONE

 OPEN(UNIT=13,FILE="WazapBD.txt")

 EPOT=EPOT/REAL(NPART)

 WRITE(6,113)  RHOSTAR,NPART,PHI,BOXX,EPOT,BDSTEP
 WRITE(13,113) RHOSTAR,NPART,PHI,BOXX,EPOT,BDSTEP

 CLOSE(UNIT=13)

 113 FORMAT(1X,//,'***  SIMULATION PARAMETERS  *** ',// &
          'REDUCED DENSITY               ',F15.8,/ &
          'NUMBER OF PARTICLES           ',I10  ,/ &
          'PACKING                       ',F15.8,/ &
          'SIMULATION X BOX SIZE         ',F15.8,/ &
          'POTENTIAL ENERGY              ',F15.8,/ &
          'SIMULATION STEPS              ',I10  ,//)

 RETURN
END SUBROUTINE SHOWSIM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE BDFORCE
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: I,J
 REAL(D):: X,Y,DX,DY,RIJ
 REAL(D):: FORCE,FZA
 REAL(D):: BETAUR,BUR

 EPOT=0.0

 DO I=1,NPART                               !SET FORCES TO ZERO
    FX(I)=0.0
    FY(I)=0.0
 ENDDO

 DO I=1,NPART-1                             !FORCE COMPUTED WITH NEWTON THIRD LAW
    X=RX(I)                                 !X POSITION I-TH REFERENCE PARTICLE
    Y=RY(I)                                 !Y POSITION I-TH REFERENCE PARTICLE
    DO J=I+1,NPART
       DX=X - RX(J)                         !X SEPARATION DISTANCE WITH J-TH PARTICLE
       DY=Y - RY(J)                         !Y SEPARATION DISTANCE WITH J-TH PARTICLE

       DX=DX - ANINT(DX*IBOXX)*BOXX         !X MINIMUM IMAGE CONVENTION
       DY=DY - ANINT(DY*IBOXY)*BOXY         !Y MINIMUM IMAGE CONVENTION

       RIJ=SQRT(DX*DX + DY*DY)              !I-J SEPARATION DISTANCE
       IF(RIJ .LE. RCUT)THEN                !CUT-OFF RADIUS CONDITION
         BUR=BETAUR(SIGMA,ESTAR,USHIFT,RIJ) !POTENTIAL BY PAIRS
         EPOT=EPOT + BUR                    !POTENTIAL ENERGY

         FZA=FORCE(SIGMA,ESTAR,RIJ)         !FORCE BY PAIRS

         FX(I)=FX(I) + FZA*DX               !X FORCE OVER I-TH PARTICLE
         FY(I)=FY(I) + FZA*DY               !Y FORCE OVER I-TH PARTICLE

         FX(J)=FX(J) - FZA*DX               !X FORCE OVER J-TH PARTICLE
         FY(J)=FY(J) - FZA*DY               !Y FORCE OVER J-TH PARTICLE
       ENDIF
    ENDDO
 ENDDO

 RETURN 
END SUBROUTINE BDFORCE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE INTEGRA
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: I
 REAL(D):: SIGGMA,RDX,RDY
 REAL(D):: RXN,RYN
 REAL(D):: GASDEV

 SIGGMA=SQRT(2.0*H)

 DO I=1,NPART                               !NEW POSITIONS @ t + ðt
    RDX=SIGGMA*GASDEV(IDUMM)                !X RANDOM DISPLACEMENT 
    RDY=SIGGMA*GASDEV(IDUMM)                !Y RANDOM DISPLACEMENT

    RXN=RX(I) + FX(I)*H + RDX               !X NEW POSITION
    RYN=RY(I) + FY(I)*H + RDY               !Y NEW POSITION

    RX(I)=RXN! - ANINT(RXN*IBOXX)*BOXX       !X PERIODIC BOUNDARY CONDITION
    RY(I)=RYN! - ANINT(RYN*IBOXY)*BOXY       !Y PERIODIC BOUNDARY CONDITION
 ENDDO

 RETURN
END SUBROUTINE INTEGRA
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE DYNAMIC
 USE BDVAR
 IMPLICIT NONE

 EPOT=EPOT/REAL(NPART)                      !POTENTIAL ENERGY PER PARTICLE
 IF(IBD .EQ. SSTAT)SWITCH=1

 IF(MOD(IBD,OUTDISP) .EQ. 0)THEN
    WRITE(6,03)IBD,EPOT
    WRITE(12,03)IBD,EPOT
 ENDIF

 IF(IBD .GE. SSTAT)THEN
   IF(MOD(IBD,OUTSAMP) .EQ. 0)CALL SAMPLE
   IF(MOD(IBD,OUTNSAM) .EQ. 0)CALL CORRELA
 ENDIF

 03 FORMAT(I10,F10.5,I10)   
 RETURN 
END SUBROUTINE DYNAMIC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE SAMPLE
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: I,J,LL
 REAL(D):: X,Y,DX,DY,RIJ

 ISAM=ISAM+1

 DO I=1,NPART-1
    X=RX(I)
    Y=RY(I)
    DO J=I+1,NPART
       DX=X - RX(J)
       DY=Y - RY(J)

       DX=DX - ANINT(DX*IBOXX)*BOXX
       DY=DY - ANINT(DY*IBOXY)*BOXY

       RIJ=SQRT(DX*DX + DY*DY)
       LL=INT(RIJ/DBIN) + 1
       IF(LL .LE. NBIN)HGR(LL)=HGR(LL) + 2
    ENDDO
 ENDDO

 RETURN
END SUBROUTINE SAMPLE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE CORRELA
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: I,T
 REAL(D):: CTIME,DTIME
 REAL(D):: MSDT
 
 IF(SWITCH .EQ. 0)THEN
   NTEL=0
   DTIME=H*REAL(OUTNSAM)
   TMAX=INT(TOBS/DTIME)
   DO I=1,TMAX
      NTIME(I)=0

      MSDX(I)=0.0D0
      MSDY(I)=0.0D0
   ENDDO

 ELSEIF(SWITCH .EQ. 1)THEN
   NTEL=NTEL+1

   IF(MOD(NTEL,IT0) .EQ. 0)THEN
     T0=T0+1
     TT0=MOD(T0-1,T0MAX) + 1
     TIME0(TT0)=NTEL
     DO I=1,NPART
        RX0(I,TT0)=RX(I)
        RY0(I,TT0)=RY(I)
     ENDDO    
   ENDIF

   DO T=1,MIN(T0,T0MAX)
      DELT=NTEL - TIME0(T) + 1
      IF(DELT .LT. TMAX)THEN
        NTIME(DELT)=NTIME(DELT) + 1
        DO I=1,NPART
           MSDX(DELT)=MSDX(DELT) + (RX(I) - RX0(I,T))**2
           MSDY(DELT)=MSDY(DELT) + (RY(I) - RY0(I,T))**2
        ENDDO
      ENDIF
   ENDDO
 ELSEIF(SWITCH .EQ. 2)THEN
   OPEN(UNIT=16,FILE="MSD.dat")

   DO I=1,TMAX-1
      CTIME=H*REAL(OUTNSAM*(I-1))
      MSDX(I)=MSDX(I)/(REAL(NPART*NTIME(I)))
      MSDY(I)=MSDY(I)/(REAL(NPART*NTIME(I)))
      MSDT=MSDX(I) + MSDY(I)

      WRITE(16,116)CTIME,MSDX(I),MSDY(I),MSDT
   ENDDO

   CLOSE(UNIT=16)
 ENDIF

 116 FORMAT(4(F15.8))
 RETURN
END SUBROUTINE CORRELA
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE FINIS
 USE BDVAR
 IMPLICIT NONE

 SWITCH=2                                   !COMPUTE CORRELATION AVERAGES

 CALL CORRELA                               !COMPUTE CORRELATIONS AVERAGES
 CALL RDF                                   !COMPUTE RADIAL DISTRIBUTION FUNCTION
 CALL PIC                                   !FINAL CONFIGURATION SNAPSHOT

 CLOSE(UNIT=12)
 RETURN
END SUBROUTINE FINIS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE RDF
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: I
 REAL(D):: R,R1,DARE,GR

 OPEN(UNIT=14,FILE="BDGr.dat")

 DO I=2,NBIN-1
    R1=REAL(I-1)*DBIN
    DARE=2.0*PI*R1*DBIN*REAL(NPART)
    GR=REAL(HGR(I))/(RHOSTAR*DARE*REAL(ISAM))
    R=(REAL(I)-0.5)*DBIN
    WRITE(14,114)R,GR
 ENDDO

 CLOSE(UNIT=14)
 114 FORMAT(2(F15.8))
 RETURN
END SUBROUTINE RDF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE PIC
 USE BDVAR
 IMPLICIT NONE
 INTEGER:: I

 OPEN(UNIT=15,FILE='BDPic.xyz')
 WRITE(15,*)NPART
 WRITE(15,*)'FRAME',1

 DO I=1,NPART
    WRITE(15,115)'C',RX(I),RY(I),0.0
 ENDDO

 CLOSE(UNIT=15)
 115 FORMAT(A,3(F15.8))
 RETURN
END SUBROUTINE PIC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION BETAUR(SIGMA,ESTAR,USHIFT,RIJ)     !LENNARD-JONES ENERGY
 INTEGER, PARAMETER:: D = KIND(1.0)
 REAL(D):: BETAUR,SIGMA,ESTAR,RIJ
 REAL(D):: USHIFT,RIJ6

 RIJ6=(SIGMA/RIJ)**6
 BETAUR=4.0*ESTAR*RIJ6*(RIJ6 - 1.0) + USHIFT

 RETURN
END FUNCTION BETAUR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION FORCE(SIGMA,ESTAR,RIJ)             !LENNARD-JONES FORCE
 INTEGER, PARAMETER:: D = KIND(1.0)
 REAL(D):: FORCE,SIGMA,ESTAR,RIJ 
 REAL(D):: RIJ2,RIJ6

 RIJ6=(SIGMA/RIJ)**6
 RIJ2=(1.0/RIJ)**2
 FORCE=24.0*ESTAR*RIJ6*(2.0*RIJ6 - 1.0)*RIJ2

 RETURN
END FUNCTION FORCE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION RAN1(IDUMM)
 INTEGER, PARAMETER:: D = KIND(1.0)
 INTEGER:: IDUMM,IA,IM,IQ,IR,NTAB,NDIV
 REAL(D):: RAN1,AM,EPS,RNMX
 PARAMETER( IA=16807, IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
            NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2E-7,RNMX=1.0-EPS )
 INTEGER:: J,K,IV(NTAB),IY
 SAVE IV,IY
 DATA IV /NTAB*0/,IY /0/

 IF(IDUMM .LE. 0 .OR.  IY .EQ. 0)THEN
   IDUMM=MAX(-IDUMM,1)
   DO J=NTAB+8,1,-1
      K=IDUMM/IQ
      IDUMM=IA*(IDUMM-K*IQ)-IR*K
      IF(IDUMM .LT. 0)IDUMM=IDUMM + IM
      IF(J .LE. NTAB)IV(J)=IDUMM
   ENDDO
   IY=IV(1)
 ENDIF

 K=IDUMM/IQ
 IDUMM=IA*(IDUMM-K*IQ) - IR*K
 IF(IDUMM .LT. 0)IDUMM=IDUMM + IM
 J=1 + IY/NDIV
 IY=IV(J)
 IV(J)=IDUMM
 RAN1=MIN(AM*IY,RNMX)

 RETURN
END FUNCTION RAN1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION GASDEV(IDUMM)
 INTEGER, PARAMETER:: D = KIND(1.0)
 INTEGER:: IDUMM
 REAL(D):: GASDEV
 INTEGER:: ISET
 REAL(D):: FAC,GSET,RSQ,V1,V2,RAN1
 SAVE ISET,GSET
 DATA ISET/0/


 IF(IDUMM .LT. 0)ISET=0
 IF(ISET .EQ. 0)THEN
   13 CONTINUE
   V1=2.0*RAN1(IDUMM) - 1.0
   V2=2.0*RAN1(IDUMM) - 1.0
   RSQ=V1**2 + V2**2
   IF( RSQ .GE. 1.0 .OR. RSQ .EQ. 0.0)GOTO 13
   FAC=SQRT(-2.0*LOG(RSQ)/RSQ)
   GSET=V1*FAC
   GASDEV=V2*FAC
   ISET=1
 ELSE
   GASDEV=GSET
   ISET=0
 ENDIF

 RETURN
END FUNCTION GASDEV
