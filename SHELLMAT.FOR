C  $STORAGE:2
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

      SUBROUTINE LAGRANGE_INTERPOLATION(X,Y,X0,Y0,XVALUE,YVALUE,N)
      REAL X(10),Y(10)

		     ! USE LAGRANGIAN INTERPOLATING FORMULAE.
		     ! CALCULATE YVALUE FOR GIVEN XVALUE
		     ! X0,Y0 ARE NEW POINTS TO BE ADDED TO THE DATA POINTS
      IF(N.GT.10)STOP 'Error in Lagrange interpolation routine. N>10'
      do 200 i=1,n-1
      x(i)=x(i+1)
200   y(i)=y(i+1)
      X(N)=X0
      Y(N)=Y0
      YVALUE=0.
      DO 100 I=1,N
      IF(Y(I).EQ.0.)GOTO 100
      TERM=1.
      DO 50 J=1,N
      IF(J.EQ.I)GOTO 50
      XJI=X(J)-X(I)
      IF(XJI.NE.0.)TERM=TERM*(X(J)-XVALUE)/XJI
50    CONTINUE
      YVALUE=YVALUE+TERM*Y(I)
100   CONTINUE
      RETURN
      END
