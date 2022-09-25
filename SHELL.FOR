C $STORAGE:2
      PROGRAM COMPOSITE_SHELL
C-----------------------------------------------------
      INCLUDE 'SHELL.BLK'
      CALL SET UP
      write(6,*)' Generating mesh data....'
      CALL MESH
      CALL VOLUME
      IF (ICASE.EQ.1)CALL STATIC_ANALYSIS
      IF (ICASE.EQ.2)CALL STABILITY_ANALYSIS
      IF (ICASE.EQ.3)CALL FAST_BUCKLING_ANALYSIS
C      IF (ICASE.EQ.3)CALL FREQUENCY_ANALYSIS
      CALL QUIT
      END

      SUBROUTINE QUIT
C--------------------
      CLOSE (8,STATUS='DELETE')
      CLOSE (30,STATUS='DELETE')
      CLOSE (40,STATUS='DELETE')
      STOP
      END

      BLOCKDATA
C---------------------------------------

      INCLUDE 'SHELL.BLK'
      DATA PI/3.14159265/
      DATA AB,WT/
     *   -.7745966692, 0., .7745966692, 0.,
     *    .5555555556, .8888888889, .5555555556, 0./

C    *   -.86113 63115 94053, -.33998 10435 84856,
C    *    .33998 10435 84856,  .86113 63115 94053,
C    *    .34785 48451 37454,  .65214 51548 62546,
C    *    .65214 51548 62546,  .34785 48451 37354/

      DATA MAXNDF,NDF,NMAT,NTERMS/8,8,0,0/
      DATA PRE,DP,DEL1/17*0./
      DATA AME/225*0./,AMASS/3200*0./
      END

      SUBROUTINE SETUP
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'
      CHARACTER*14 FILE
      FAIL=.FALSE.
      CALL RNULL(CORDS,200)
      CALL RNULL(SEG,150)
      CALL RNULL(AJUNC,11)
      CALL RNULL(BCOND,160)
      CALL  NULL(IVAR,400)
      CALL RNULL(ALDATA,990)
      CALL RNULL(ANDATA,500)
      CALL RNULL(DISPL,400)

      OPEN (30,FILE='E:WORK.F01',STATUS='UNKNOWN',
     *          FORM='UNFORMATTED')
      OPEN (40,FILE='E:WORK.F02',STATUS='UNKNOWN',
     *          FORM='UNFORMATTED',ACCESS='DIRECT',RECL=34)
C---------------------------READ DATA FROM 8 ----------
      WRITE(6,1)
1     FORMAT(' Jobname : ? ')
      READ(5,2)FILE
2     FORMAT(A8)
      OPEN (UNIT=8, FILE='A:'//FILE//'.DMP', STATUS='UNKNOWN',
     *     FORM='UNFORMATTED')
      READ(8)N0,NMIN,NMAX,MAXROOTS,ICASE,IPRINT,CRITERION
      READ(8)NSEGS,NHIST,JUNCS,WIRL
      READ(8)CORDS(1,1),CORDS(2,1)
      READ(8)SEG,RING,HIST,BCOND,IROT
      READ(8)(((ALAYER(I,J,K),I=1,17),J=1,20),K=1,3)
      READ(8)(((ALAYER(I,J,K),I=1,17),J=1,20),K=4,6)
      READ(8)(((ALAYER(I,J,K),I=1,17),J=1,20),K=7,9)
      READ(8)((ALAYER(I,J,10),I=1,17),J=1,20)
      OPEN (10,FILE='A:'//FILE//'.OUT',STATUS='UNKNOWN',
     *   FORM='FORMATTED')
C------------------------------------------------------
      RETURN
      END
      SUBROUTINE STATIC_ANALYSIS
C------------------------------------------------------------

      INCLUDE 'SHELL.BLK'
      VIBRATION=.FALSE.
      BUCKLE=.FALSE.
      STAB=.FALSE.
      N=0
      IDASH=1
      CALL ASSIGN
      NEQNS0=NEQNS
      MBAND0=6
      MAXS0=MAXSIZ
      CALL IMOVE(IVAR,IVAR0,MAXVAR)
      PRE=0.
      DO 100 ISTEP=1,NHIST
      MAXINC=HIST(ISTEP,2)
      IF(MAXINC.LE.0)MAXINC=1
      PMAX=HIST(ISTEP,1)
      DPRE=(PMAX-PRE)/MAXINC
      IF(ABS(DPRE).LE.0.1E-10)GOTO 100
      DO 90 I=1,MAXINC
      PRE=PRE+DPRE
      CALL SHL001
      WRITE(10,80)PRE,ITER,DMAX
      CALL RESULTS(0)
      CALL CHECK_FAILURE(1)
90    CONTINUE
      CALL RESULTS(1)
100   CONTINUE
80    FORMAT(//,' PRES. MULT= ',1PE12.3,'  ITERATIONS= ',
     *  I4,'  MAX. DISP= ',E12.3)
      RETURN
      END
      SUBROUTINE MESH
C------------------------------------------------------------

      INCLUDE 'SHELL.BLK'
      NODE=1
      AJUNC(1)=NODE
      NEL=1
      ANDATA(1,1)=0.
      DO 100 NSEG=1,NSEGS
      P1=SEG(12,NSEG)
      P2=SEG(13,NSEG)
      NELSEG=SEG(2,NSEG)
      IF(SEG(3,NSEG).EQ.1.)THEN
      FACTOR=1./SEG(2,NSEG)
      ELSE
      FACTOR=(1.-SEG(3,NSEG))/(1.-SEG(3,NSEG)**NELSEG)
      END IF
      NSTYPE=SEG(1,NSEG)
      GOTO(120,140,160),NSTYPE
      GOTO 100
120   CONTINUE
      RAD1=1.0E+20
      DX=FACTOR*(SEG(4,NSEG)-CORDS(1,NODE))
      DR=FACTOR*(SEG(5,NSEG)-CORDS(2,NODE))
      ALEN=SQRT(DX*DX+DR*DR)
      GRADP=(P2-P1)/ALEN
      ANG(1)=ANGLE(DX,DR)
      DO 130 NOD=1,NELSEG
      NODE=NODE+1
      CORDS(1,NODE)=CORDS(1,NODE-1)+DX
      CORDS(2,NODE)=CORDS(2,NODE-1)+DR
      SL=SQRT(DX*DX+DR*DR)
      ALDATA(1,NEL)=SL
      ALDATA(2,NEL)=ANG(1)
      ALDATA(3,NEL)=RAD1
      ANDATA(1,NODE)=ANDATA(1,NODE-1)+SL
      ALDATA(7,NEL)=P1
      P1=P1+GRADP*SL
      ALDATA(8,NEL)=P1
      MATERIAL(NEL)=NSEG
      NEL=NEL+1
      DX=DX*SEG(3,NSEG)
      DR=DR*SEG(3,NSEG)
130   CONTINUE
      IF(CORDS(2,NODE).LT.0.)CORDS(2,NODE)=0.
      GOTO 180
140   CONTINUE
      DX=CORDS(1,NODE)-SEG(6,NSEG)
      DR=CORDS(2,NODE)-SEG(7,NSEG)
      ANG(1)=ANGLE(DR,DX)
      RAD1=SQRT(DX*DX+DR*DR)
      DX=SEG(4,NSEG)-SEG(6,NSEG)
      DR=SEG(5,NSEG)-SEG(7,NSEG)
      ANG(2)=ANGLE(DR,DX)
      RAD2=SQRT(DX*DX+DR*DR)
      RAD1=.5*(RAD1+RAD2)
      DPHI=ANG(2)-ANG(1)
      IF(IROT(NSEG).EQ. 1.AND.DPHI.LT.0.)ANG(2)=ANG(2)+2.*PI
      IF(IROT(NSEG).EQ.-1.AND.DPHI.GT.0.)ANG(1)=ANG(1)+2.*PI
      DPHI=ANG(2)-ANG(1)
      ALEN=ABS(RAD1*DPHI)
      GRADP=(P2-P1)/ALEN
      DPHI=FACTOR*DPHI
      DO 150 NOD=1,NELSEG
      NODE=NODE+1
      ALDATA(2,NEL)=-ANG(1)
      AMEAN=ANG(1)+.5*DPHI
      RMEAN=RAD1
      IF(AMEAN.GT..5*PI.AND.AMEAN.LT.1.5*PI)RMEAN=-RAD1
      ANG(1)=ANG(1)+DPHI
      SL=ABS(RAD1*DPHI)
      DPHI=SEG(3,NSEG)*DPHI
      CORDS(1,NODE)=SEG(6,NSEG)+RAD1*SIN(ANG(1))
      CORDS(2,NODE)=SEG(7,NSEG)+RAD1*COS(ANG(1))
      ALDATA(1,NEL)=SL
      ALDATA(3,NEL)=RMEAN
      ANDATA(1,NODE)=ANDATA(1,NODE-1)+SL
      ALDATA(7,NEL)=P1
      P1=P1+GRADP*SL
      ALDATA(8,NEL)=P1
      MATERIAL(NEL)=NSEG
      NEL=NEL+1
150   CONTINUE
      IF(CORDS(2,NODE).LT.0.)CORDS(2,NODE)=0.
      GOTO 180
160   CONTINUE
C  GENERAL SHELL SEGMENT.    NOT TO BE USED NOW.
      RAD1=1.0E+20
      READ(5,333)ITITLE
333   FORMAT(20A4)
      ALEN=0.
      DO 170 NOD=1,NELSEG
      NODE=NODE+1
      READ(5,*)CORDS(1,NODE),CORDS(2,NODE)
      DX=CORDS(1,NODE)-CORDS(1,NODE-1)
      DR=CORDS(2,NODE)-CORDS(2,NODE-1)
      SL=SQRT(DX*DX+DR*DR)
      ANG(1)=ANGLE(DX,DR)
      ALDATA(1,NEL)=SL
      ALDATA(2,NEL)=ANG(1)
      ALDATA(3,NEL)=RAD1
      ANDATA(1,NODE)=ANDATA(1,NODE-1)+SL
      ALEN=ALEN+SL
      MATERIAL(NEL)=NSEG
      NEL=NEL+1
170   CONTINUE
      NODE=NODE-NELSEG
      GRADP=(P2-P1)/ALEN
      DO 175 NOD=1,NELSEG
      ALDATA(7,NODE)=P1
      SL=ALDATA(1,NODE)
      P1=P1+GRADP*SL
      ALDATA(8,NODE)=P1
      NODE=NODE+1
175   CONTINUE
180   CONTINUE
      AJUNC(NSEG+1)=NODE
100   CONTINUE
      MAXNOD=NODE
      MAXNEL=NEL-1


      NDIFF=MAXNDF-NDF
      IF(NDIFF.EQ.0)GOTO 421
      CALL RNULL(AW,17)
      DO  NEL=1,MAXNEL
      NN=MAXNDF+3
      DO  J=1,NDIFF
      NN=NN-1
      IREC=NDIFF*(NEL-1)+J
      WRITE(40,REC=IREC)(AW(I),I=1,NN)
      end do
      end do
421   CONTINUE


      WRITE(10,199)
199   FORMAT(//,' NODE     X            R           S')
      DO 200 NODE=1,MAXNOD
      WRITE(10,201)NODE,(CORDS(J,NODE),J=1,2),ANDATA(1,NODE)
200   CONTINUE
201   FORMAT(I5,5(1PE12.3))
      WRITE(10,202)
202   FORMAT(//,
     * ' ELEMENT  LENGTH    ANGLE(1)       RAD1         ',
     * 'P1            P2')
      DO 300 NEL=1,MAXNEL
      ANG(1)=ALDATA(2,NEL)*180./PI
      WRITE(10,301)NEL,ALDATA(1,NEL),ANG(1),ALDATA(3,NEL),
     * ALDATA(7,NEL),ALDATA(8,NEL)
300   CONTINUE
301   FORMAT(I5,7(1PE12.3),I8)
      WRITE(10,399)
399   FORMAT(//,20X,'RING LOAD MULTIPLIERS',/
     *  ' JNCTN NODE    AXIAL      CIRCUM.      RADIAL      MOMENT',//)
      DO 400 JUNC=1,JUNCS
      NODE=AJUNC(JUNC)
      WRITE(10,401)JUNC,NODE,(RING(J,JUNC),J=1,4)
400   CONTINUE
401   FORMAT(2I5,4(1PE12.3))
      RETURN
      END
      FUNCTION ANGLE(X,Y)
C------------------------------------------------------------

      PI=3.14159265
      IF(X.EQ.0.)THEN
      A=.5*PI*SIGN(1.,Y)
      ELSE
      A=ATAN(Y/X)
      IF(X.LT.0.)A=A+PI
      END IF
      if(a.lt.0.)a=a+2.*pi
      ANGLE=A
      RETURN
      END
      SUBROUTINE SOLVE(I1)
C------------------------------------------------------------

      INCLUDE 'SHELL.BLK'
      IF(I1.EQ.1)CALL ADDLOD
      DO 9 I=1,JUNCS
      DO 9 J=1,4
      K=BCOND(1,J,I)
      IF(K.EQ.1)GOTO 9
      NODE=AJUNC(I)
      VAR=BCOND(2,J,I)
      NEQN=IVAR(J,NODE)
      IF(NEQN.EQ.0)GOTO 9
      DO 8 K=2,MBAND
      L=NEQN-K+1
      IF(L.LE.0)GOTO 1
      M=L+(K-1)*NEQNS
      ALOAD(L)=ALOAD(L)-VAR*STIFF(M)
      STIFF(M)=0.
1     L=NEQN+K-1
      IF(L.GT.NEQNS)GOTO 8
      M=NEQN+(K-1)*NEQNS
      ALOAD(L)=ALOAD(L)-VAR*STIFF(M)
      STIFF(M)=0.
8     CONTINUE
      STIFF(NEQN)=1.
      ALOAD(NEQN)=VAR
9     CONTINUE

      IF(I1.EQ.2)THEN
      CALL DETERM(STIFF,MAXSIZ,DET,IDET,NEQNS,MBAND)

      ELSE IF (I1.EQ.3)THEN
C  PREDICTION OF BUCKLING PRESSURES BY EIGENVALUE ANALYSIS.
      CALL INVERSE_ITERATIONS(STIFF,AMASS,3200,NEQNS,MBAND,DISP,
     *  ALOAD,EIG,KOUNT,FAIL)
      IF(FAIL)THEN
      WRITE(10,*)' INVERSE_ITERATIONS FAILS TO FIND THE EIGENVALUE'
      CALL QUIT
      END IF

      PSTAR=PRE+EIG*DPRE
      EIGEN(1)=EIG

      ELSE IF (I1.EQ.4)THEN
C  VIBRATION FREQUENCIES
      RETURN

      ELSE
C  STATIC DISPLACEMENTS
      CALL MOVE(DISPL,DISPL0,MAXVAR)

      IRUN=1
      CALL LINEQN(STIFF,MAXSIZ,ALOAD,NEQNS,MBAND,IRUN,FAIL)
      IF(FAIL)THEN
      WRITE(10,101)
101   FORMAT(' POSSIBLE CAUSE OF FAILURE MAY BE DUE TO AXISYMMETRIC
     * COLLAPSE OF THE SHELL')
      CALL QUIT
      END IF
      ERROR=0.
      DMAX=0.
      DO 20 NODE=1,MAXNOD
      DO 20 I=1,4
      NEQN=IVAR(I,NODE)
      IF(NEQN.EQ.0)GOTO 20
      DISPL(I,NODE)=DISPL(I,NODE)+ALOAD(NEQN)
      IF(ABS(DMAX).GT.ABS(DISPL(I,NODE)))GOTO 20
      DMAX=DISPL(I,NODE)
      IF(ABS(DMAX).LT.1.0E-20)GOTO 20
      ERROR=ALOAD(NEQN)/DMAX
      ERROR=ABS(ERROR)
20    CONTINUE
      END IF

      RETURN
      END
      SUBROUTINE ASSIGN
C------------------------------------------------------------

      INCLUDE 'SHELL.BLK'
      DO 10 NODE=1,MAXNOD
      DO 10 J=1,4
      IVAR(J,NODE)=1
10    CONTINUE
      MBAND=8
      RAD=CORDS(2,1)
      IF(RAD.EQ.0..AND.N.EQ.0)IVAR(4,1)=0
      RAD=CORDS(2,MAXNOD)
      IF(RAD.EQ.0..AND.N.EQ.0)IVAR(4,MAXNOD)=0
      IF(N.EQ.0.AND.IDASH.EQ.1)THEN
      MBAND=6
      DO 20 NODE=1,MAXNOD
      IVAR(2,NODE)=0
20    CONTINUE
      END IF
      IF(N.EQ.0.AND.IDASH.EQ.2)THEN
      MBAND=2
      DO 30 NODE=1,MAXNOD
      IVAR(1,NODE)=0
      IVAR(3,NODE)=0
      IVAR(4,NODE)=0
30    CONTINUE
      END IF
      DO 35 I=1,JUNCS
      DO 35 J=1,4
      K=BCOND(1,J,I)
      IF(K.EQ.1)GOTO 35
      IF(BCOND(2,J,I).NE.0.)GOTO 35
      NODE=AJUNC(I)
      IVAR(J,NODE)=0
35    CONTINUE
      NEQNS=0
      DO 40 NODE=1,MAXNOD
      DO 40 I=1,4
      IF(IVAR(I,NODE).EQ.0)GOTO 40
      NEQNS=NEQNS+1
      IVAR(I,NODE)=NEQNS
40    CONTINUE
      MAXVAR=MAXNOD*4
      MAXSIZ=NEQNS*MBAND
      WRITE(10,50)N,NEQNS,MBAND,MAXSIZ
50    FORMAT(//,' WAVE,N=',I3,'  NEQNS=',I4,'  MBAND=',I2,
     *       '  MAXSIZ=',I6,/)
      RETURN
      END
      SUBROUTINE ADDLOD
C------------------------------------------------------------

      INCLUDE 'SHELL.BLK'
      DO 100 JUNC=1,JUNCS
      NODE=AJUNC(JUNC)
      DO 20 I=1,4
      NEQN=IVAR(I,NODE)
      IF(NEQN.EQ.0)GOTO 20
      ALOAD(NEQN)=ALOAD(NEQN)+PRE*RING(I,JUNC)
20    CONTINUE
100   CONTINUE
      RETURN
      END
      SUBROUTINE REDUCE
C------------------------------------------------------------

      INCLUDE 'SHELL.BLK'
      DO 10 I=1,NDIFF
      N1=MAXNDF-I
      NN=N1+1
      IF(SE(NN,NN).EQ.0.)GOTO 3
      DO 2 L=1,N1
      AW(2)=SE(L,NN)/SE(NN,NN)
      RE(L)=RE(L)-RE(NN)*AW(2)
      DO 2 M=L,N1
      AW(1)=SE(M,NN)/SE(NN,NN)
      SE(L,M)=SE(L,M)-SE(L,NN)*AW(1)
      AME(L,M)=AME(L,M)-AME(L,NN)*AW(1)-
     *  AME(M,NN)*AW(2)+AME(NN,NN)*AW(1)*AW(2)
2     CONTINUE
3     CONTINUE
      IREC=NDIFF*(NEL-1)+I
      WRITE(40,REC=IREC)(SE(J,NN),J=1,NN),RE(NN),DEL(NN)
10    CONTINUE
      RETURN
      END
      SUBROUTINE RECOVR
C------------------------------------------------------------

      INCLUDE 'SHELL.BLK'
      IF(NDIFF.EQ.0)RETURN
      DO 10 I=1,NDIFF
      N1=MAXNDF-I
      NN=N1+1
      IREC=NDIFF*(NEL-1)+I
      READ(40,REC=IREC)(SE(J,NN),J=1,NN),RE(NN),DEL(NN)
10    CONTINUE
      DO 20 I=1,NDIFF
      N1=NDF-1+I
      NN=N1+1
      DEL1(NN)=RE(NN)
      DO 15 J=1,N1
      DEL1(NN)=DEL1(NN)-SE(J,NN)*DEL1(J)
15    CONTINUE
      IF(SE(NN,NN).EQ.0.)SE(NN,NN)=1.
      DEL1(NN)=DEL1(NN)/SE(NN,NN)
      DEL(NN)=DEL(NN)+DEL1(NN)
20    CONTINUE
      RETURN
      END
      SUBROUTINE TRANSF
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'

      IF(NDIFF.GT.0)CALL REDUCE
      CALL TRANS(RE,ANG,AW,1,0)
      CALL TRANS(RE,ANG,AW,2,0)

      CALL TRANS(SE,ANG,AW,1,1)
      CALL TRANS(SE,ANG,AW,1,2)
      CALL TRANS(SE,ANG,AW,2,2)

      CALL TRANS(AME,ANG,AW,1,1)
      CALL TRANS(AME,ANG,AW,1,2)
      CALL TRANS(AME,ANG,AW,2,2)
      RETURN
      END
      SUBROUTINE TRANS(A,ANG,AW,I,J)
C------------------------------------------------------------
      REAL A(15,1),AW(100)
      REAL ANG(2)
      IF(J.EQ.0)THEN
      I1=1+(I-1)*4
      I2=I1+2
      AW(1)=A(I1,1)
      AW(2)=A(I2,1)
      CO=COS(ANG(I))
      SI=SIN(ANG(I))
      A(I1,1)=CO*AW(1)-SI*AW(2)
      A(I2,1)=SI*AW(1)+CO*AW(2)
      RETURN
      END IF
      N=0
      I1=4*(I-1)+1
      J1=4*(J-1)+1
      I2=I1+3
      J2=J1+3
      DO 10 K=I1,I2
      DO 10 L=J1,J2
      N=N+1
      IF(K.LE.L)THEN
      AW(N)=A(K,L)
      ELSE
      AW(N)=A(L,K)
      END IF
10    CONTINUE
      CI=COS(ANG(I))
      SI=SIN(ANG(I))
      CJ=COS(ANG(J))
      SJ=SIN(ANG(J))
      A(I1,J1)=CI*(AW(1)*CJ-AW(3)*SJ)-SI*(AW(9)*CJ-AW(11)*SJ)
      A(I1,J1+1)=CI*AW(2)-AW(10)*SI
      A(I1,J1+2)=CI*(AW(1)*SJ+AW(3)*CJ)-SI*(AW(9)*SJ+AW(11)*CJ)
      A(I1,J1+3)=AW(4)*CI-AW(12)*SI
      IF(I.NE.J)THEN
      A(I1+1,J1)=AW(5)*CJ-AW(7)*SJ
      A(I1+2,J1)=SI*(AW(1)*CJ-AW(3)*SJ)+CI*(AW(9)*CJ-AW(11)*SJ)
      A(I1+2,J1+1)=AW(2)*SI+AW(10)*CI
      A(I1+3,J1)=AW(13)*CJ-AW(15)*SJ
      A(I1+3,J1+2)=AW(13)*SJ+AW(15)*CJ
      END IF
      A(I1+1,J1+2)=AW(5)*SJ+AW(7)*CJ
      A(I1+2,J1+2)=SI*(AW(1)*SJ+AW(3)*CJ)+CI*(AW(9)*SJ+AW(11)*CJ)
      A(I1+2,J1+3)=AW(4)*SI+AW(12)*CI
      RETURN
      END
      SUBROUTINE MERGE
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'
      DO 10 I=1,4
      NW(I)=IVAR(I,NEL)
      NW(I+4)=IVAR(I,NEL+1)
10    CONTINUE
      DO 20 I=1,8
      IF(NW(I).EQ.0)GOTO 20
      ALOAD(NW(I))=ALOAD(NW(I))+RE(I)
      DO 30 J=I,8
      IF(NW(J).EQ.0)GOTO 30
      I1=MIN0(NW(I),NW(J))
      J1=MAX0(NW(I),NW(J))
      J1=J1-I1+1
      J1=I1+(J1-1)*NEQNS
      STIFF(J1)=STIFF(J1)+SE(I,J)
      AMASS(J1)=AMASS(J1)+AME(I,J)
30    CONTINUE
20    CONTINUE
      RETURN
      END
      SUBROUTINE SHAPE(S1,I1)
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'
      DO 1 I=1,7
      DO 1 J=1,MAXNDF
1     SHP(I,J)=0.
      R1=S1/SL
      R2=R1*R1
      R3=R1*R2

      IF(N.EQ.0.AND.IDASH.EQ.2)GOTO 10
C**********   U   ******************
      GOTO(11,14,12,14,13),NDIFF+1
11    SHP(1,1)=1.-R1
      SHP(1,5)=R1
      GOTO 14

12    SHP(1,1)=1.-3.*R1+2.*R2
      SHP(1,5)=-R1+2.*R2
      SHP(1,9)=4.*R1-4.*R2
      GOTO 14

13    SHP(1,1)=1.-5.5*R1+9.*R2-4.5*R3
      SHP(1,5)=R1-4.5*R2+4.5*R3
      SHP(1,9)=9.*R1-22.5*R2+13.5*R3
      SHP(1,10)=-4.5*R1+18.*R2-13.5*R3

14    CONTINUE
C**********   W   ******************
      SHP(3,3)=1.-3.*R2+2.*R3
      SHP(3,4)=-SL*(R1-2.*R2+R3)
      SHP(3,1)=-SHP(3,4)*RAD1I
      SHP(3,7)=3.*R2-2.*R3
      SHP(3,8)=-SL*(-R2+R3)
      SHP(3,5)=-SHP(3,8)*RAD1I

C**********  DU/DS  ******************
      GOTO(21,24,22,24,23),NDIFF+1
21    SHP(4,1)=-1./SL
      SHP(4,5)=1./SL
      GOTO 24

22    SHP(4,1)=(-3.+4.*R1)/SL
      SHP(4,5)=(-1.+4.*R1)/SL
      SHP(4,9)=(4.-8.*R1)/SL
      GOTO 24

23    SHP(4,1)=(-5.5+18.*R1-13.5*R2)/SL
      SHP(4,9)=(9.-45.*R1+40.5*R2)/SL
      SHP(4,10)=(-4.5+36.*R1-40.5*R2)/SL
      SHP(4,5)=(1.-9.*R1+13.5*R2)/SL

24    CONTINUE
C**********  DW/DS  *****************
      SHP(6,3)=6.*(R2-R1)/SL
      SHP(6,4)=-(1.-4.*R1+3.*R2)
      SHP(6,1)=-SHP(6,4)*RAD1I
      SHP(6,7)=-SHP(6,3)
      SHP(6,8)=-(-2.*R1+3.*R2)
      SHP(6,5)=-SHP(6,8)*RAD1I

C**********  D2W/DS2  ***************
      SHP(7,4)=-(-4.+6.*R1)/SL
      SHP(7,1)=-SHP(7,4)*RAD1I
      SHP(7,3)=(-6.+12.*R1)/(SL*SL)
      SHP(7,8)=-(-2.+6.*R1)/SL
      SHP(7,5)=-SHP(7,8)*RAD1I
      SHP(7,7)=-SHP(7,3)

10    IF(N.EQ.0.AND.IDASH.EQ.1)GOTO 20
C**********   V   ******************
      GOTO(31,34,32,34,33),NDIFF+1
31    SHP(2,2)=1.-R1
      SHP(2,6)=R1
      GOTO 34

32    SHP(2,2)=1.-3.*R1+2.*R2
      SHP(2,6)=-R1+2.*R2
      SHP(2,10)=4.*R1-4.*R2
      GOTO 34

33    SHP(2,2)=1.-5.5*R1+9.*R2-4.5*R3
      SHP(2,6)=R1-4.5*R2+4.5*R3
      SHP(2,11)=9.*R1-22.5*R2+13.5*R3
      SHP(2,12)=-4.5*R1+18.*R2-13.5*R3

34    CONTINUE
C**********  DV/DS  ******************
      GOTO(41,20,42,20,43),NDIFF+1
41    SHP(5,2)=-1./SL
      SHP(5,6)=1./SL
      GOTO 20

42    SHP(5,2)=(-3.+4.*R1)/SL
      SHP(5,6)=(-1.+4.*R1)/SL
      SHP(5,10)=(4.-8.*R1)/SL
      GOTO 20

43    SHP(5,2)=(-5.5+18.*R1-13.5*R2)/SL
      SHP(5,11)=(9.-45.*R1+40.5*R2)/SL
      SHP(5,12)=(-4.5+36.*R1-40.5*R2)/SL
      SHP(5,6)=(1.-9.*R1+13.5*R2)/SL

20    CONTINUE

      IF(I1.EQ.0)RETURN
      DO 60 I=1,7
      DE(I)=0.
      DO 60 J=1,MAXNDF
      DE(I)=DE(I)+SHP(I,J)*DEL(J)
60    CONTINUE
      U=DE(1)
      V=DE(2)
      W=DE(3)
      BETA=DE(6)-U*RAD1I

      RETURN
      END
      SUBROUTINE DMAT
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'
      NSEG=MATERIAL(NEL)
      IF(NSEG.EQ.NMAT)RETURN
      NMAT=NSEG
      LAYERS=SEG(8,NSEG)
      HMEAN=0.
      hdensity=0.
      DO 10 LAYER=1,LAYERS
      HMEAN=HMEAN+ALAYER(1,LAYER,NSEG)
      hdensity=hdensity+alayer(1,layer,nseg)*alayer(17,layer,nseg)
10    CONTINUE
      HMEAN=.5*HMEAN
      H2=-HMEAN
      DO 20 I=1,6
      DO 20 J=1,6
20    D(I,J)=0.
      DO 30 LAYER=1,LAYERS
      E11=ALAYER(2,LAYER,NSEG)
      E22=ALAYER(3,LAYER,NSEG)
      N21=ALAYER(4,LAYER,NSEG)
      G12=ALAYER(5,LAYER,NSEG)
      WANG=ALAYER(11,LAYER,NSEG)*PI/180.
      N12=N21*E11/E22
      DENO=(1.-N12*N21)
      Q11=E11/DENO
      Q22=E22/DENO
      Q12=N12*E22/DENO
      Q66=G12
      CS=COS(WANG)
      SN=SIN(WANG)
      C4=CS*CS*CS*CS
      S4=SN*SN*SN*SN
      CS2=CS*CS*SN*SN
      CS3=CS*SN*SN*SN
      SC3=SN*CS*CS*CS
      QB11=Q11*C4+2.*(Q12+2.*Q66)*CS2+Q22*S4
      QB12=(Q11+Q22-4.*Q66)*CS2+Q12*(C4+S4)
      QB16=-Q22*CS3+Q11*SC3-(Q12+2.*Q66)*(SC3-CS3)
      QB22=Q11*S4+2.*(Q12+2.*Q66)*CS2+Q22*C4
      QB26=-Q22*SC3+Q11*CS3+(Q12+2.*Q66)*(SC3-CS3)
      QB66=(Q11+Q22-2.*Q12)*CS2+Q66*(C4-2.*CS2+S4)
      H=ALAYER(1,LAYER,NSEG)
      H1=H2
      H2=H1+H
      H22=.5*(H2*H2-H1*H1)
      H23=(H2*H2*H2-H1*H1*H1)/3.
      D(1,1)=D(1,1)+QB11*H
      D(1,2)=D(1,2)+QB12*H
      D(2,2)=D(2,2)+QB22*H
      D(3,3)=D(3,3)+QB66*H
      D(1,4)=D(1,4)+QB11*H22
      D(1,5)=D(1,5)+QB12*H22
      D(2,5)=D(2,5)+QB22*H22
      D(3,6)=D(3,6)+QB66*H22
      D(4,4)=D(4,4)+QB11*H23
      D(4,5)=D(4,5)+QB12*H23
      D(5,5)=D(5,5)+QB22*H23
      D(6,6)=D(6,6)+QB66*H23
30    CONTINUE
      D(2,1)=D(1,2)
      D(2,4)=D(1,5)
      D(4,1)=D(1,4)
      D(4,2)=D(2,4)
      D(5,1)=D(1,5)
      D(5,2)=D(2,5)
      D(5,4)=D(4,5)
      D(6,3)=D(3,6)
      RETURN
      END
      SUBROUTINE LDISPL(I1)
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'

      DO 1010 I=1,4
      DEL(I)=DISPL(I,NEL)
      DEL(I+4)=DISPL(I,NEL+1)
      IF(I1.EQ.0)GOTO 1010
      DEL1(I)=DEL(I)-DISPL0(I,NEL)
      DEL1(I+4)=DEL(I+4)-DISPL0(I,NEL+1)
1010  CONTINUE

      CO=COS(ANG(1))
      SI=SIN(ANG(1))
      AW(1)=DEL(1)*CO+DEL(3)*SI
      AW(2)=-DEL(1)*SI+DEL(3)*CO
      DEL(1)=AW(1)
      DEL(3)=AW(2)

      IF(I1.EQ.0)GOTO 1020
      AW(1)=DEL1(1)*CO+DEL1(3)*SI
      AW(2)=-DEL1(1)*SI+DEL1(3)*CO
      DEL1(1)=AW(1)
      DEL1(3)=AW(2)

1020  CONTINUE
      CO=COS(ANG(2))
      SI=SIN(ANG(2))
      AW(1)=DEL(5)*CO+DEL(7)*SI
      AW(2)=-DEL(5)*SI+DEL(7)*CO
      DEL(5)=AW(1)
      DEL(7)=AW(2)

      IF(I1.EQ.0)GOTO 1030
      AW(1)=DEL1(5)*CO+DEL1(7)*SI
      AW(2)=-DEL1(5)*SI+DEL1(7)*CO
      DEL1(5)=AW(1)
      DEL1(7)=AW(2)

1030  CONTINUE
      IF(NDIFF.GT.0)CALL RECOVR
      RETURN
      END
      SUBROUTINE LDATA
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'
      SL=ALDATA(1,NEL)
      ANG(1)=ALDATA(2,NEL)
      RAD1=ALDATA(3,NEL)
      RAD1I=1.+1./RAD1
      RAD1I=RAD1I-1.
      ANG(2)=ANG(1)-SL*RAD1I
      CALL DMAT

      P1=ALDATA(7,NEL)*PRE
      P2=ALDATA(8,NEL+1)*PRE
      GRADP=(P2-P1)/SL
      RETURN
      END
      subroutine volume
      INCLUDE 'SHELL.BLK'
      vol=0.
      do 100 nel=1,maxnel
      call ldata
      call dmat
      do 200 it=1,4
      s1=.5*sl*(1.+ab(it))
      weight=wt(it)
      PHI=ANG(1)-S1*RAD1I
      IF(ABS(RAD1).LT.1.0E10)THEN
      RAD=CORDS(2,NEL)+RAD1*(COS(PHI)-COS(ANG(1)))
      ELSE
      RAD=CORDS(2,NEL)+S1*SIN(PHI)
      END IF
      WAIT=SL*RAD*WEIGHT*PI
      vol=vol+wait*hdensity
200   continue
100   continue
      write(6,300)vol
300   format(' Mass of the shell =',1pe12.3)
      write(10,300)vol
      return
      end
      SUBROUTINE SHL001
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'
      ITER=1
      N=0
      IDASH=1
      MBAND=6
      NEQNS=NEQNS0
      MAXSIZ=MAXS0
      CALL RNULL (DISPL,400)
      CALL IMOVE(IVAR0,IVAR,MAXVAR)
      BUCKLE=.FALSE.
      STAB=.FALSE.
      write(6,11)pre
11    format(' Pressure multiplier = ',1pe12.3)
1000  CONTINUE
      write(6,22)iter
22    format(' Iteration = ',I3)
      CALL RNULL(STIFF,MAXSIZ)
      CALL RNULL(ALOAD,NEQNS)

      DO 2000 NEL=1,MAXNEL
      CALL LDATA
      CALL LDISPL(0)
      CALL RNULL(SE,225)
      CALL RNULL(RE,15)

      DO 1200 IT=1,4
      S1=.5*SL*(1.+AB(IT))
      WEIGHT=WT(IT)
      PHI=ANG(1)-S1*RAD1I
      IF(ABS(RAD1).LT.1.0E10)THEN
      RAD=CORDS(2,NEL)+RAD1*(COS(PHI)-COS(ANG(1)))
      ELSE
      RAD=CORDS(2,NEL)+S1*SIN(PHI)
      END IF
      P=P1+S1*GRADP
      CALL SHAPE(S1,1)
      CALL STRAIN(1)
      CALL MULT1
1200  CONTINUE
      CALL TRANSF
      CALL MERGE
2000  CONTINUE
      CALL SOLVE(1)
      ITER=ITER+1

      IF(ITER.GT.30)THEN
      WRITE(10,111)
111   format(' NEWTON-RAPHSON ITERATIONS FAIL TO CONVERGE IN SHL001.',
     *  ' POSSIBLE REASON OF THE FAILURE MAY BE ',/,
     *  '      1.  AN AXI-SYMMETRIC COLLAPSE or',/,
     *  '      2.  YOUR LOAD STEP IS TOO LARGE.',//)
      CALL QUIT
      END IF

      IF(ERROR.GT.0.1E-2)GOTO 1000
      RETURN
      END
      SUBROUTINE RESULTS(I1)
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'
      IF(IPRINT.EQ.1)REWIND 30
      IF(I1.EQ.0)CALL RNULL(INDEX,49)
99    FORMAT('1 PRESSURE OR LOAD MULTIPLIER=',1PE12.3,//,
     * '   STATION        U           V           W          BETA    ',
     * '    N(FI)     N(THETA)   N(FI,THETA)    M(FI)  ',
     * '   M(THETA)  M(FI,THETA)')
      IF(I1.EQ.1)WRITE(10,99)PRE
100   FORMAT(1PE12.3,10E12.3)
      DO 2000 NEL=1,MAXNEL
      CALL LDATA
      CALL LDISPL(0)
c      DO 1000 S1=0.15*SL,1.01*SL,0.35*SL
      DO 1000 S1=0., 1.01*sl, 0.5*sl
      PHI=ANG(1)-S1*RAD1I
      IF(ABS(RAD1).LT.1.0E10)THEN
      RAD=CORDS(2,NEL)+RAD1*(COS(PHI)-COS(ANG(1)))
      ELSE
      RAD=CORDS(2,NEL)+S1*SIN(PHI)
      END IF
      P=P1+S1*GRADP
      CALL SHAPE(S1,1)
      STATION=ANDATA(1,NEL)+S1
      CALL STRAIN(1)
      IF(I1.EQ.1)WRITE(10,100)STATION,U,V,W,BETA,FRC
      CALL STRESS(I1)
1000  CONTINUE
2000  CONTINUE
      IF(IPRINT.NE.0)THEN
      REWIND 30
      IF (I1.EQ.1)WRITE(10,109)
109   FORMAT('1 STATION    LAYER     Z           EPS-1        EPS-2',
     *'    EPS-12       SIGMA-1     SIGMA-2      SIGMA-12',//)
110   CONTINUE
      READ(30,ERR=150)STATION,LAYER,Z,E1,E2,E12,S1,S2,S12
      WRITE(10,120)STATION,LAYER,Z,E1,E2,E12,S1,S2,S12
120   FORMAT(1PE12.3,I5,7E12.3)
      GOTO 110
150   CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE CHECK_FAILURE(J1)
C-----------------------------
      INCLUDE 'SHELL.BLK'
      REAL FINDEX(4)
      WRITE(10,101)
101   FORMAT('  MAXIMUM INDICES OF FAILURE',/
     *       '            MAX.STRS.    MAX.STRN   TSAI-HILL',
     *'    HOFFMAN  TSAI-WU.STRS  TSAI-WU.STRN    OWEN',/)
      WRITE(10,102)(INDEX(1,I),I=1,7)
102   FORMAT('  STATION',1PE12.3,6E12.3)
      WRITE(10,103)(INDEX(2,I),I=1,7)
103   FORMAT('  LAYER  ',7(6X,F6.1))
      WRITE(10,104)(INDEX(3,I),I=1,7)
104   FORMAT('    Z    ',1PE12.3,6E12.3)
      WRITE(10,105)(INDEX(4,I),I=1,7)
105   FORMAT(' SIGMA-1 ',1PE12.3,6E12.3)
      WRITE(10,106)(INDEX(5,I),I=1,7)
106   FORMAT(' SIGMA-2 ',1PE12.3,6E12.3)
      WRITE(10,107)(INDEX(6,I),I=1,7)
107   FORMAT(' SIGMA-12',1PE12.3,6E12.3)
      WRITE(10,108)(INDEX(7,I),I=1,7)
108   FORMAT(' INDEX   ',1PE12.3,6E12.3)
      I1=0
      IF(CRITERION.EQ.8)THEN
      J=0
      DO 200 I=1,7
      IF(INDEX(7,I).GT.1.)J=J+1
200   CONTINUE
      IF(J.EQ.7)THEN
      PLOWER=PRE-DPRE
      WRITE(10,201)PLOWER,PRE
      I1=1
      END IF
      ELSE
      IF(NTERMS.LT.4)THEN
      NTERMS=NTERMS+1
      FINDEX(NTERMS)=INDEX(7,CRITERION)-1.
      ELSE
      FINDEX(1)=FINDEX(2)
      FINDEX(2)=FINDEX(3)
      FINDEX(3)=FINDEX(4)
      FINDEX(4)=INDEX(7,CRITERION)-1.
      END IF
      IF(FINDEX(NTERMS).GT.0.)THEN
      CALL INTERPOLATE(PRE,FINDEX,NTERMS,DPRE,PLOWER,PUPPER)
      WRITE(10,201)PLOWER,PUPPER
      I1=1
      END IF
      END IF
      IF(I1.EQ.0)RETURN
201   FORMAT(' WARNING **** MATERIAL FAILURE ****
     * BETWEEN PRES. MULT.',1PE12.3,'  AND ', E12.3)
      IF(J1.EQ.1)THEN
      PRE=PLOWER
      CALL SHL001
      CALL RESULTS(1)
      write(6,301)pre
301   format(' Pressure multiplier at the failure point =',1pe12.3)
      CALL QUIT
      END IF
      RETURN
      END
      SUBROUTINE STRESS(I1)
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'
C  FROM LOCAL STRESSES TO FIBRE STRESSES
      NSEG=NMAT
      LAYERS=SEG(8,NSEG)
      HMEAN=0.
      DO 10 LAYER=1,LAYERS
      HMEAN=HMEAN+ALAYER(1,LAYER,NSEG)
10    CONTINUE
      HMEAN=.5*HMEAN
      H2=-HMEAN
      DO 30 LAYER=1,LAYERS
      H=ALAYER(1,LAYER,NSEG)
      H1=H2
      H2=H1+H
      WANG=ALAYER(11,LAYER,NSEG)*PI/180.
      CO=COS(WANG)
      SI=SIN(WANG)
      C2=CO*CO
      S2=SI*SI
      CS=CO*SI
      DO 30 IZ=1,2
      Z=H1+H*(IZ-1)
C   STRAINS IN MERIDIONAL AND GLOBAL DIRECTIONS.
      EX=EPS(1)-Z*EPS(4)
      EY=EPS(2)-Z*EPS(5)
      EXY=EPS(3)+Z*EPS(6)
C   STRAINS IN FIBRE DIRECTIONS
      E1=C2*EX+S2*EY+CS*EXY
      E2=S2*EX+C2*EY-CS*EXY
      E12=-CS*(EX-EY)+(C2-S2)*EXY

      E11=ALAYER(2,LAYER,NSEG)
      E22=ALAYER(3,LAYER,NSEG)
      N21=ALAYER(4,LAYER,NSEG)
      G12=ALAYER(5,LAYER,NSEG)
      N12=N21*E11/E22
      DENO=(1.-N12*N21)
      Q11=E11/DENO
      Q22=E22/DENO
      Q12=N12*E22/DENO
      Q66=G12

C  STRESSES IN FIBRE DIRECTIONS
      ST1=Q11*E1+Q12*E2
      ST2=Q12*E1+Q22*E2
      ST12=Q66*E12
      IF(I1.EQ.1.AND.IPRINT.EQ.1)
     *     WRITE(30)STATION,LAYER,Z,E1,E2,E12,ST1,ST2,ST12
      IF(I1.EQ.0)CALL CHECK(LAYER,Z,ST1,ST2,ST12,E1,E2,E12)
30    CONTINUE
      RETURN
      END
      SUBROUTINE CHECK(LAYER,Z,S1,S2,T12,E1,E2,E12)
C---------------------------------------------------
      INCLUDE 'SHELL.BLK'
      REAL FI(7)
      X1T=ALAYER(6,LAYER,NSEG)
      X1C=ALAYER(7,LAYER,NSEG)
      X2T=ALAYER(8,LAYER,NSEG)
      X2C=ALAYER(9,LAYER,NSEG)
      S12=ALAYER(10,LAYER,NSEG)
      E1T=ALAYER(12,LAYER,NSEG)
      E1C=ALAYER(13,LAYER,NSEG)
      E2T=ALAYER(14,LAYER,NSEG)
      E2C=ALAYER(15,LAYER,NSEG)
      ES12=ALAYER(16,LAYER,NSEG)

C  MAX. STRESS THEORY

      FI(1)=0.
      IF(S1.GT.0.)THEN
      FI(1)=AMAX1(FI(1),S1/X1T)
      ELSE
      FI(1)=AMAX1(FI(1),-S1/X1C)
      END IF

      IF(S2.GT.0.)THEN
      FI(1)=AMAX1(FI(1),S2/X2T)
      ELSE
      FI(1)=AMAX1(FI(1),-S2/X2C)
      END IF

      FI(1)=AMAX1(FI(1),ABS(T12)/S12)

C  MAX. STRAIN THEORY.

      FI(2)=0.
      IF(E1.GT.0.)THEN
      FI(2)=AMAX1(FI(2),E1/E1T)
      ELSE
      FI(2)=AMAX1(FI(2),-E1/E1C)
      END IF

      IF(E2.GT.0.)THEN
      FI(2)=AMAX1(FI(2),E2/E2T)
      ELSE
      FI(2)=AMAX1(FI(2),-E2/E2C)
      END IF

      FI(2)=AMAX1(FI(2),ABS(E12)/ES12)

C  TSAI-HILL STRESS THEORY

      FI(3)=0.
      IF(S1.GT.0.)THEN
      A=S1/X1T
      FI(3)=A*A
      ELSE
      A=S1/X1C
      FI(3)=A*A
      END IF
      IF(S2.GT.0.)THEN
      B=S2/X2T
      C=S2/X1T
      FI(3)=FI(3)+B*B
      ELSE
      B=S2/X2C
      C=S2/X1C
      FI(3)=FI(3)+B*B
      END IF
      FI(3)=FI(3)+T12*T12/(S12*S12)-A*C

C  HOFFMAN THEORY

      FI(4)=S1*(X1C-X1T)/(X1C*X1T)+S2*(X2C-X2T)/(X2C*X2T)
     *     +S1*S1/(X1T*X1C)+S2*S2/(X2T*X2C)+T12*T12/(S12*S12)
     *     -S1*S2/(X1T*X1C)

C  TSAI-WU STRESS THEORY

      FI(5)=S1*(X1C-X1T)/(X1C*X1T)+S2*(X2C-X2T)/(X2C*X2T)
     *     +S1*S1/(X1C*X1T)+S2*S2/(X2T*X2C)+T12*T12/(S12*S12)
     *     -S1*S2/SQRT(X1T*X1C*X2T*X2C)

C  TSAI-WU STRAIN THEORY

      F11=1./(X1C*X1T)
      F22=1./(X2C*X2T)
      F1=(X1C-X1T)*F11
      F2=(X2C-X2T)*F22
      F33=1./(S12*S12)
      F12=.5*SQRT(F11*F22)
      FI(6)=(F1*Q11+F2*Q12)*E1+(F1*Q12+F2*Q22)*E2
     *     +(F11*Q12*Q12+2.*F12*Q11*Q12+F22*Q12*Q12)*E1*E1
     *     +(F11*Q12*Q12+2.*F12*Q12*Q22+F22*Q22*Q22)*E2*E2
     *     +F33*Q66*Q66*E12*E12
     *+2.*(F11*Q11*Q12+F12*(Q11*Q22+Q12*Q12)+F22*Q12*Q22)*E1*E2

C  OWEN THEORY

      IF(S1.GT.0..AND.S2.GT.0.)FI(7)=
     *  S1*S1/(X1T*X1T)+S2*S2/(X2T*X2T)+T12*T12/(S12*S12)
      IF(S1.LT.0..AND.S2.LT.0.)FI(7)=
     *  S1*S1/(X1C*X1C)+S2*S2/(X2C*X2C)+T12*T12/(S12*S12)
      IF(S1.GT.0..AND.S2.LT.0.)FI(7)=
     *  S1*S1/(X1T*X1T)-S1*S2/(X1T*X2C)+S2*S2/(X2C*X2C)+
     *  T12*T12/(S12*S12)
      IF(S1.LT.0..AND.S2.GT.0.)FI(7)=
     *  S1*S1/(X1C*X1C)-S1*S2/(X1C*X2T)+S2*S2/(X2T*X2T)+
     *  T12*T12/(S12*S12)

      DO 100 I=1,7
      IF(FI(I).LT.INDEX(7,I))GOTO 100
      INDEX(1,I)=STATION
      INDEX(2,I)=FLOAT(LAYER)
      INDEX(3,I)=Z
      INDEX(4,I)=S1
      INDEX(5,I)=S2
      INDEX(6,I)=T12
      INDEX(7,I)=FI(I)
100   CONTINUE

      RETURN
      END
      SUBROUTINE LINEQN(S,MAXS,R,NSIZE,MBAND,IRUN,FAIL)
C----------------------------------------
      REAL S(MAXS),R(NSIZE)
      LOGICAL*1 FAIL
      IF(IRUN.EQ.2)GOTO 1000
      DO 790 N=1,NSIZE
      DO 780 L=2,MBAND
      NL=N+(L-1)*NSIZE
      IF (S(NL).EQ.0.)GOTO 780
      I=N+L-1
      IF(S(N).EQ.0.)GOTO 2000
      C=S(NL)/S(N)
      J=0
      DO 750 K=L,MBAND
      J=J+1
      IJ=I+(J-1)*NSIZE
      NK=N+(K-1)*NSIZE
750   S(IJ)=S(IJ)-C*S(NK)
      S(NL)=C
780   CONTINUE
790   CONTINUE
      IRUN=2
1000  CONTINUE
      DO 830 N=1,NSIZE
      DO 820 L=2,MBAND
      NL=N+(L-1)*NSIZE
      IF(S(NL).EQ.0.)GOTO 820
      I=N+L-1
      R(I)=R(I)-S(NL)*R(N)
820   CONTINUE
      IF(S(N).EQ.0.)GOTO 2000
830   R(N)=R(N)/S(N)
      DO 860 M=2,NSIZE
      N=NSIZE+1-M
      DO 850 L=2,MBAND
      NL=N+(L-1)*NSIZE
      IF(S(NL).EQ.0.)GOTO 850
      K=N+L-1
      R(N)=R(N)-S(NL)*R(K)
850   CONTINUE
860   CONTINUE
      RETURN
2000  CONTINUE
      FAIL=.TRUE.
      WRITE(6,*)' LINEQN FAILS'
      RETURN
      END
      SUBROUTINE DETERM(A,MAXA,D,ID,NEQ,IBAND)
      REAL A(MAXA),C
      IFAIL=1
      D=1.0
      ID=0
      TOL=1.0E-30
C
C     MATRIX REDUCTION BY GAUSS ELIMINATION.
C
      DO 100 I=1,NEQ
      DO 200 J=2,IBAND
      LOCA=(J-1)*NEQ+I
      IF(ABS(A(LOCA)).LT.TOL)GOTO 200
      K=I+J-1
      IF(ABS(A(I)).LT.TOL)GOTO 500
      C=A(LOCA)/A(I)
      L=0
      DO 300 M=J,IBAND
      L=L+1
      KL=(L-1)*NEQ+K
      IM=(M-1)*NEQ+I
300   A(KL)=A(KL)-C*A(IM)
      A(LOCA)=C
200   CONTINUE
100   CONTINUE

      DO 400 I=1,NEQ
      D=D*A(I)
410   CONTINUE
      IF(ABS(D).GE.1.)GOTO 430
420   CONTINUE
      IF(ABS(D).GE.0.1)GOTO 400
      D=D*10.0
      ID=ID-1
      GOTO 420
430   CONTINUE
      D=0.1*D
      ID=ID+1
      GOTO 410
400   CONTINUE
      RETURN
500   WRITE(6,1)I,A(I)
1     FORMAT(/,5X,I5,'TH ELEMENT ALONG THE LEADING DIAGONAL OF THE MATRI
     1X IN DETERM IS',E12.5)
      D=0.0
      STOP
      END
      SUBROUTINE NULL(I,N)
      INTEGER I(N)
      DO 1 J=1,N
1     I(J)=0
      RETURN
      END
      SUBROUTINE RNULL(A,N)
      REAL A(N)
      DO 1 I=1,N
1     A(I)=0.
      RETURN
      END
      SUBROUTINE MOVE(FROM,TO,N)
      REAL FROM(N),TO(N)
      DO 1 I=1,N
1     TO(I)=FROM(I)
      RETURN
      END
      SUBROUTINE MOVE8(FROM,TO,N)
      REAL FROM(N),TO(N)
      DO 1 I=1,N
1     TO(I)=FROM(I)
      RETURN
      END
      SUBROUTINE IMOVE(FROM,TO,N)
      INTEGER FROM(N),TO(N)
      DO 1 I=1,N
1     TO(I)=FROM(I)
      RETURN
      END
      SUBROUTINE SET_VARIABLES
C-----------------------------
      INCLUDE 'SHELL.BLK'

      IF(N0.EQ.0)GOTO 10
      N=N0
      IDASH=1
      CALL ASSIGN
      NEQNSN=NEQNS
      MBANDN=MBAND
      MAXSN=MAXSIZ
      IDASHN=IDASH
      CALL IMOVE(IVAR,IVARN,MAXVAR)

10    CONTINUE
      N=0
      IDASH=1
      CALL ASSIGN
      NEQNS0=NEQNS
      MBAND0=6
      MAXS0=MAXSIZ
      CALL IMOVE(IVAR,IVAR0,MAXVAR)
      RETURN
      END
      SUBROUTINE FAST_BUCKLING_ANALYSIS
C--------------------------------------
      INCLUDE 'SHELL.BLK'
      CALL SET_VARIABLES
      PRE=HIST(1,1)
      IF(PRE.EQ.0.)PRE=1.
      DPRE=0.1*PRE
      CALL BUCKLING_ANALYSIS
      RETURN
      END
      SUBROUTINE STABILITY_ANALYSIS
C----------------------------------
      INCLUDE 'SHELL.BLK'

      VIBRATION=.FALSE.
      CALL SET_VARIABLES
      PRE=0.
      DO 100 ISTEP=1,NHIST
      MAXINC=HIST(ISTEP,2)
      PMAX=HIST(ISTEP,1)
      IF(MAXINC.LE.0)MAXINC=1
      DPRE=(PMAX-PRE)/MAXINC
      IF(ABS(DPRE).LT.1.0E-10)GOTO 100
      NTERM=0
      DO 90 INC=1,MAXINC
      PRE=PRE+DPRE
      CALL SHL001
      BUCKLE=.TRUE.
      STAB=.TRUE.
      CALL SHL002
      WRITE(10,80)PRE,DET,IDET,N0,ITER,ERROR
      IF(NTERM.EQ.0)NDET=IDET
      IF(NTERM.LT.4)THEN
      NTERM=NTERM+1
      ELSE
      FUN(1)=FUN(2)
      FUN(2)=FUN(3)
      FUN(3)=FUN(4)
      END IF
      FUN(NTERM)=DET*10.0**(IDET-NDET)
      IF(ISTEP*INC.EQ.1)GOTO 89
      IF(DET*DET0.LE.0.)GOTO 110
89    CONTINUE
      DET0=DET
90    CONTINUE
100   CONTINUE
80    FORMAT(' PRES/MULT =',1PE12.3,' DET =',F7.4,
     *  ' 10**',I5,' WAVE NO. =',I2,' ITERATIONS=',I2,
     *' ERROR=',E12.3)
      WRITE(10,85)
85    FORMAT(' NO CHANGE IN SIGN OF THE DETERMINANT OF '
     *  'THE STABILITY MATRIX.')
      RETURN
110   CONTINUE

      CALL INTERPOLATE(PRE,FUN,NTERM,DPRE,PLOWER,PUPPER)
      WRITE(10,160)PLOWER,PUPPER,N0
160   FORMAT(//,' BUCKLING PRESSURE LIES IN BETWEEN ',1PE12.3,
     *  '  AND ', 1E10.3, 'WAVE NO.',I5)
      WRITE(10,195)
195   FORMAT(1X,60('--'),/)
      PRE=PLOWER
      DPRE=0.01*PRE
      CALL BUCKLING_ANALYSIS
      RETURN
      END
      SUBROUTINE INTERPOLATE(X0,FUN,NTERM,STEP,XL,XU)
      REAL FUN(NTERM)
C  X0=CURRENT VALUE OF X, I.E. AT FUN(NTERM)
C  THIS PROCEDURE USES BACKWARD-DIFFERENCE FORMULA TO INTERPOLATE AND
C  CALCULATE THE X VALUE AT THE ZERO OF THE FUNCTION BETWEEN NTERM-1 AND NTERM
C  STEP. XL AND XU ARE LOWER AND UPPER BOUND OF THE ZERO.
      Y0=FUN(NTERM-1)
      DO 150 A=-1.0,0.0001,.01
      Y=FUN(NTERM)
      IF(NTERM.GT.1)THEN
      D1=FUN(NTERM)-FUN(NTERM-1)
      Y=Y+A*D1
      END IF
      IF(NTERM.GT.2)THEN
      D1=FUN(NTERM)-2.*FUN(NTERM-1)+FUN(NTERM-2)
      Y=Y+A*(A+1.)*D1/2.
      END IF
      IF(NTERM.GT.3)THEN
      D1=FUN(NTERM)-3.*FUN(NTERM-1)+3.*FUN(NTERM-2)-FUN(NTERM-3)
      Y=Y+A*(A+1.)*(A+2.)*D1/6.
      END IF
      IF(Y*Y0.LE.0.)GOTO 170
150   CONTINUE
170   CONTINUE
      XU=X0+A*STEP
      XL=XU-.01*STEP
      RETURN
      END
      SUBROUTINE SHL002
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'

      NEQNS=NEQNSN
      MBAND=MBANDN
      MAXSIZ=MAXSN
      CALL IMOVE(IVARN,IVAR,MAXVAR)
      CALL RNULL(STIFF,MAXSIZ)
      CALL RNULL(AMASS,MAXSIZ)
      CALL RNULL(ALOAD,NEQNS)
      DO 2000 NEL=1,MAXNEL
      CALL RNULL(SE,225)
      CALL RNULL(AME,225)
      CALL LDATA
      IF(BUCKLE.AND.STAB)THEN
      CALL LDISPL(0)
      ELSE
      CALL LDISPL(1)
      END IF
      DO 1200 IT=1,4
      S1=.5*SL*(1.+AB(IT))
      WEIGHT=WT(IT)
      PHI=ANG(1)-S1*RAD1I
      IF(ABS(RAD1).LT.1.0E10)THEN
      RAD=CORDS(2,NEL)+RAD1*(COS(PHI)-COS(ANG(1)))
      ELSE
      RAD=CORDS(2,NEL)+S1*SIN(PHI)
      END IF
      P=P1+S1*GRADP
      DP=P*DPRE/PRE

      IF(BUCKLE)THEN
      N=0
      IDASH=1
      CALL SHAPE(S1,1)
      CALL STRAIN(1)
      END IF

      N=N0
      IDASH=IDASHN
      CALL SHAPE(S1,0)
      CALL STRAIN(0)

      CALL MULT1
1200  CONTINUE
      CALL TRANSF
      CALL MERGE
2000  CONTINUE
      IF(STAB)CALL SOLVE(2)
      IF(BUCKLE.AND..NOT.STAB)CALL SOLVE(3)
      RETURN
      END
      SUBROUTINE SHL003
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'

      N=N0
      IDASH=IDASHN
      NEQNS=NEQNSN
      MBAND=MBANDN
      MAXSIZ=MAXSN
      CALL IMOVE(IVARN,IVAR,MAXVAR)
      CALL RNULL(STIFF,MAXSIZ)
      CALL RNULL(AMASS,MAXSIZ)

      DO 2000 NEL=1,MAXNEL
      CALL RNULL(SE,225)
      CALL RNULL(AME,225)
      CALL LDATA
      CALL LDISPL(0)
      DO 1200 IT=1,4
      S1=.5*SL*(1.+AB(IT))
      WEIGHT=WT(IT)
      PHI=ANG(1)-S1*RAD1I
      IF(ABS(RAD1).LT.1.0E10)THEN
      RAD=CORDS(2,NEL)+RAD1*(COS(PHI)-COS(ANG(1)))
      ELSE
      RAD=CORDS(2,NEL)+S1*SIN(PHI)
      END IF
      CALL SHAPE(S1,0)
      CALL STRAIN(0)
      CALL MULT1
1200  CONTINUE
      CALL TRANSF
      CALL MERGE
2000  CONTINUE
      CALL SOLVE(4)
      RETURN
      END
      SUBROUTINE EIGENVECTOR
C--------------------------
      INCLUDE 'SHELL.BLK'
      WRITE(6,11)PRE,N0
11    FORMAT(' FINAL BUCKLING LOAD =',1PE12.3,'  WAVE No.',I3)

      WRITE(10,100)N0
100   FORMAT(//,' BUCKLING OR VIBRATION MODE --  WAVE NO. ',I3,//,
     * '   STATION        U           V           W          BETA    ',
     *  //)
200   FORMAT(1PE12.3,4E12.3)

      N=N0
      IDASHN=1
      NEQNS=NEQNSN
      MBAND=MBANDN
      MAXSIZ=MAXSN
      CALL IMOVE(IVARN,IVAR,MAXVAR)

      DO 2000 NODE=1,MAXNOD
      DO 1000 I=1,4
      NEQN=IVAR(I,NODE)
      AW(I)=0.
      IF(NEQN.NE.0)AW(I)=DISP(NEQN)
1000  CONTINUE
      WRITE(10,200)ANDATA(1,NODE),(AW(I),I=1,4)
2000  CONTINUE
      RETURN
      END
      SUBROUTINE INVERSE_ITERATIONS(A,B,MAXA,N,MBAND,X,Y,EIGEN,KOUNT,
     *  FAIL)
C--------------------------------------------------------------------
      REAL A(MAXA),B(MAXA),X(N),Y(N)
      LOGICAL*1 FAIL
      M=MBAND-1
      EIGEN=0.
      IRUN=1
      KOUNT=0
      DO 100 I=1,N
100   X(I)=1.
200   CONTINUE
      KOUNT=KOUNT+1
      DO 300 I=1,N
      J1=I-M
      J2=I+M
      IF(J1.LE.0)J1=1
      IF(J2.GT.N)J2=N
      Y(I)=0.
      DO 300 J=J1,J2
      IF(J.GE.I)THEN
      K=J-I+1
      K=I+(K-1)*N
      ELSE
      K=I-J+1
      K=J+(K-1)*N
      END IF
      Y(I)=Y(I)+B(K)*X(J)
300   CONTINUE
      CALL LINEQN(A,MAXA,Y,N,MBAND,IRUN,FAIL)
      IF(FAIL)RETURN
      YMAX=0.
      DO 400 I=1,N
      IF(ABS(YMAX).LT.ABS(Y(I)))YMAX=Y(I)
400   CONTINUE
      IF(ABS(YMAX).LT.1.0E-10)THEN
      WRITE(6,*)' NULL EIGENVECTOR FOUND DURING INVERSE_ITERATIONS'
      CALL QUIT
      END IF
      DO 500 I=1,N
500   X(I)=Y(I)/YMAX
      ERROR=ABS(1.-EIGEN*YMAX)
      EIGEN=1./YMAX
      IF(ERROR.LT.1.0E-6.OR.KOUNT.GT.100)RETURN
      GOTO 200
      END
      SUBROUTINE BUCKLING_ANALYSIS
C---------------------------------
      INCLUDE 'SHELL.BLK'

      CALL SHL001
      WRITE(10,95)PRE,ITER,ERROR
95    FORMAT(' PRES/LOAD MULT.= ',1PE12.3,' ITERATIONS=',I2,
     *  ' ERROR=',E12.3)
      PRE=PRE+DPRE
      CALL SHL001
      WRITE(10,95)PRE,ITER,ERROR
      PMIN=1.0E30
      PRE=PRE-DPRE
      IF(NMIN.LT.1)NMIN=1
      DO 140 N0=NMIN,NMAX
      BUCKLE=.TRUE.
      STAB=.FALSE.
      CALL SHL002
      EIGEN(1)=-EIGEN(1)
      PSTAR=PRE+DPRE*EIGEN(1)
      WRITE(10,141)N0,PSTAR
      WRITE(6,141)N0,PSTAR
      IF(PMIN.GT.PSTAR)THEN
      PMIN=PSTAR
      N0MIN=N0
      END IF
140   CONTINUE
141   FORMAT(10X,'WAVE NO. ',I2,' CRITICAL PRES/LOAD MULT.=',1PE12.3)
      N0=N0MIN
      PRE=PMIN
      DPRE=0.001*PRE
      CALL SHL001
      WRITE(10,95)PRE,ITER,ERROR
      CALL RESULTS(0)
      CALL CHECK_FAILURE(0)
      CALL RESULTS(1)
      PRE=PRE+DPRE
      CALL SHL001
      WRITE(10,95)PRE,ITER,ERROR
      PRE=PRE-DPRE
      BUCKLE=.TRUE.
      STAB=.FALSE.
      CALL SHL002
      CALL EIGENVECTOR
      RETURN
      END
      SUBROUTINE STRAIN(I1)
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'
      if(rad.eq.0.)then
      rn=0.
      si=0.
      co=0.
      else
      RN=FLOAT(N)/RAD
      SI=SIN(PHI)/RAD
      CO=COS(PHI)/RAD
      endif
      DO 100 I=1,MAXNDF
      DO 5 J=1,9
      BMAT(J,I)=0.
5     CONTINUE

      IF(N.EQ.0.AND.IDASH.EQ.2)GOTO 10
      BMAT(7,I)=SHP(6,I)-SHP(1,I)*RAD1I
10    IF(N.EQ.0.AND.IDASH.EQ.1)GOTO 20
      BMAT(8,I)=.5*(-RN*SHP(1,I)-SHP(5,I)-SI*SHP(2,I))
      BMAT(9,I)=-RN*SHP(3,I)-CO*SHP(2,I)

20    IF(N.EQ.0.AND.IDASH.EQ.2)GOTO 30
      BMAT(1,I)=SHP(4,I)+SHP(3,I)*RAD1I
      BMAT(2,I)=RN*SHP(2,I)+SI*SHP(1,I)+CO*SHP(3,I)
      BMAT(4,I)=SHP(7,I)-SHP(4,I)*RAD1I
      BMAT(5,I)=RN*BMAT(9,I)+SI*BMAT(7,I)
30    IF(N.EQ.0.AND.IDASH.EQ.1)GOTO 100
      BMAT(3,I)=-RN*SHP(1,I)+SHP(5,I)-SI*SHP(2,I)
      BMAT(6,I)=RN*BMAT(7,I)+SI*BMAT(9,I)+CO*SHP(5,I)
      BMAT(6,I)=2.*BMAT(6,I)
100   CONTINUE

      IF(VIBRATION)RETURN
      IF(I1.EQ.0)GOTO 510

      DO 400 I=1,9
      EPS(I)=0.
      EPS1(I)=0.
      DO 400 J=1,MAXNDF
      EPS(I)=EPS(I)+BMAT(I,J)*DEL(J)
      EPS1(I)=EPS1(I)+BMAT(I,J)*DEL1(J)
400   CONTINUE
      EPS(1)=EPS(1)+.5*EPS(7)**2
      EPS1(1)=EPS1(1)+.5*EPS1(7)*EPS1(7)

C      IF(BUCKLE.AND..NOT.STAB)THEN
C      DO 450 I=1,9
C      A=EPS(I)
C      EPS(I)=EPS1(I)
C      EPS1(I)=A
C450   CONTINUE
C      END IF


      DO 500 I=1,6
      FRC(I)=0.
      FRC1(I)=0.
      DO 490 J=1,6
      FRC(I)=FRC(I)+D(I,J)*EPS(J)
490   FRC1(I)=FRC1(I)+D(I,J)*EPS1(J)
      FRC1(I)=FRC1(I)-FRC(I)
500   CONTINUE

510   CONTINUE
      DO 410 I=1,MAXNDF
      BMAT(1,I)=BMAT(1,I)+EPS(7)*BMAT(7,I)
      BMAT(3,I)=BMAT(3,I)+EPS(7)*BMAT(9,I)
410   CONTINUE

      RETURN
      END
      SUBROUTINE MULT1
C------------------------------------------------------------
      INCLUDE 'SHELL.BLK'
      if(rad.eq.0.)then
      co=0.
      else
      CO=COS(PHI)/RAD
      end if
      WAIT=.5*SL*RAD*WEIGHT*PI
      IF(N.EQ.0)WAIT=2.*WAIT
      DO 10 I=1,MAXNDF
      IF(VIBRATION)GOTO 5
      AW(1)=0.
      DO 1 J=1,6
      AW(1)=AW(1)-FRC(J)*BMAT(J,I)
1     CONTINUE
      AW(1)=AW(1)+P*(SHP(3,I)+U*SHP(1,I)*RAD1I+V*CO*SHP(2,I)
     *     -W*(RAD1I+CO)*SHP(3,I))
      RE(I)=RE(I)+WAIT*AW(1)
5     CONTINUE

      DO 10 J=I,MAXNDF
      AW(1)=0.
      DO 111 K=1,6
      DO 111 L=1,6
      AW(1)=AW(1)+BMAT(K,I)*D(K,L)*BMAT(L,J)
111   CONTINUE

      IF (VIBRATION)GOTO 20

      AW(2)=BMAT(7,I)*BMAT(7,J)+BMAT(8,I)*BMAT(8,J)
      AW(3)=BMAT(8,I)*BMAT(8,J)+BMAT(9,I)*BMAT(9,J)
      AW(4)=BMAT(7,I)*BMAT(9,J)+BMAT(7,J)*BMAT(9,I)
      AW(5)=SHP(1,I)*SHP(1,J)*RAD1I+CO*SHP(2,I)*SHP(2,J)
     *              -(RAD1I+CO)*SHP(3,I)*SHP(3,J)
      AW(1)=AW(1)+FRC(1)*AW(2)+FRC(2)*AW(3)
     *         +FRC(3)*AW(4)-P*AW(5)
20    CONTINUE
      SE(I,J)=SE(I,J)+WAIT*AW(1)
      IF(VIBRATION)THEN
C  WARNING H IS THE THICKNESS OF THE SHELL WALL AND IT IS UNDEFINED AT
C  THIS STAGE.
      AME(I,J)=AME(I,J)+RO*H*WAIT*(SHP(1,I)*SHP(1,J)+
     *        SHP(2,I)*SHP(2,J)+SHP(3,I)*SHP(3,J))
      ELSE
      AME(I,J)=AME(I,J)-WAIT*(FRC1(1)*AW(2)+FRC1(2)*AW(3)
     *        +FRC1(3)*AW(4)-DP*AW(5))
      END IF

10    CONTINUE
      RETURN
      END
      SUBROUTINE FREQUENCY_ANALYSIS
C----------------------------------
      INCLUDE 'SHELL.BLK'
      RETURN
      END
