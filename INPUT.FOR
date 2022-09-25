C$STORAGE:2
      PROGRAM COMPOSITE_SHELL_DATA_GENERATOR
C-------------------------------------------
      INCLUDE 'SHELL.BLK'
      CHARACTER JOBSES*14,JOBDMP*14,ANSWER*1,JOBNAME*7,JOBSES1*14
      WRITE(6,333)
333   FORMAT(' Jobname (1..7 CHARACTERS ) ? ')
      READ(5,444)JOBNAME
444   FORMAT(A7)
      JOBDMP='A:'//JOBNAME//'.DMP'
      OPEN (UNIT=8,FILE=JOBDMP,STATUS='unknown',FORM='UNFORMATTED')
      WRITE(6,1)
1     FORMAT(/,' Do you want to read data from the file created',/,
     * ' from the previous session ? ')
      READ(5,555)ANSWER
555   FORMAT(A)
      JOBSES='A:'//JOBNAME//'.SES'
      IF(ANSWER.EQ.'Y'.OR.ANSWER.EQ.'y')THEN
      JOBSES1='SHELL.SES'
      CALL READ_DATA(JOBSES,JOBSES1)
      WRITE(6,101)JOBSES1
101   FORMAT(/,' New session file ',a14,' has been created')
      ELSE
      CALL INPUT_DATA(JOBSES)
      END IF
      WRITE(8)N0,NMIN,NMAX,MAXROOTS,ICASE,IPRINT,CRITERION
      WRITE(8)NSEGS,NHIST,JUNCS,WIRL
      WRITE(8)CORDS(1,1),CORDS(2,1)
      WRITE(8)SEG,RING,HIST,BCOND,IROT
      WRITE(8)(((ALAYER(I,J,K),I=1,17),J=1,20),K=1,3)
      WRITE(8)(((ALAYER(I,J,K),I=1,17),J=1,20),K=4,6)
      WRITE(8)(((ALAYER(I,J,K),I=1,17),J=1,20),K=7,9)
      WRITE(8)((ALAYER(I,J,10),I=1,17),J=1,20)
      STOP
      END
      SUBROUTINE INPUT_DATA(JOBSES)
C------------------------------INTERACTIVE INPUT-------------

      INCLUDE 'SHELL.BLK'
      CHARACTER*8 CDISP(4),ANSWER*1,JOBSES*14
      DATA CDISP/'       u','       v','       w','Rotation'/
      OPEN (UNIT=7,FILE=JOBSES,STATUS='unknown',FORM='FORMATTED')

      WRITE(6,50)
50    FORMAT(T20,'Analysis (1:static, 2:stability, 3:',
     *'buckling, 4:fequency) ? ',T2,I5)
      READ(5,*)ICASE
      WRITE(7,50)ICASE

      IF(ICASE.NE.4)THEN
      WRITE(6,70)
70    FORMAT(T20,'Do you want to print strains & stresses (Y/N) ? '
     *,T2,5X,A1)
      READ(5,111)ANSWER
      WRITE(7,70)ANSWER
      IPRINT=0
      IF(ANSWER.EQ.'Y'.OR.ANSWER.EQ.'y')IPRINT=1
      WRITE(6,71)
71    FORMAT(T20,'Failure criterion ? (1:Max. stress, 2:Max. strain,'
     * ,/,T30,'3:Tsai-Hill, 4:Hoffman, 5:Tsai-Wu stress,',/,
     * T30,'6:Tsai-Wu strain, 7:Owen, 8:All) ? ',T2,I5)
      READ(5,*)CRITERION
      WRITE(7,71)CRITERION
      END IF


      IF(ICASE.NE.1)THEN
      WRITE(6,51)
51    FORMAT(T20,'Starting wave number (n) ? ',T2,I5)
      READ(5,*)N0
      WRITE(7,51)N0
      WRITE(6,52)
52    FORMAT(T20,'Minimum wave number ? ',T2,I5)
      READ(5,*)NMIN
      WRITE(7,52)NMIN
      WRITE(6,53)
53    FORMAT(T20,'Maximum wave number ? ',T2,I5)
      READ(5,*)NMAX
      WRITE(7,53)NMAX
      WRITE(6,54)
54    FORMAT(T20,'Number of eigenvalues per n ? ',T2,I5)
      READ(5,*)MAXROOTS
      WRITE(7,54)MAXROOTS
      END IF

      WRITE(6,1)
1     FORMAT(T20,'Number of segments ? ',T2,I5)
      READ(5,*)NSEGS
      WRITE(7,1)NSEGS

      WRITE(6,3)
3     FORMAT(T20,'Starting x co-ord. of the structure ? ',T2,
     * 1PE14.6)
      READ(5,*)CORDS(1,1)
      WRITE(6,3001)
3001  FORMAT(T20,'Starting r co-ord. of the structure ? ',T2,
     *    1PE14.6)
      READ(5,*)CORDS(2,1)
      WRITE(7,3)CORDS(1,1)
      WRITE(7,3001)CORDS(2,1)

      DO 100 I=1,NSEGS
      WRITE(6,4)I
      WRITE(7,4)I
4     FORMAT(T20,'************Segment No.',I2,'************')
      WRITE(6,5)
5     FORMAT(T20,'Type (1: cone/cylinder, 2: toroid/sphere) ? '
     * ,T2,F7.1)
      READ(5,*)SEG(1,I)
      WRITE(7,5)SEG(1,I)
      WRITE(6,6)
6     FORMAT(T20,'Number of elements in the segment ? ',T2,F7.1)
      READ(5,*)SEG(2,I)
      WRITE(7,6)SEG(2,I)
      WRITE(6,7)
7     FORMAT(T20,'Step ratio (>1 for incr. and <1 for decr) ? '
     *,T2,1PE14.6)
      READ(5,*)SEG(3,I)
      WRITE(7,7)SEG(3,I)
      WRITE(6,8)
8     FORMAT(T20,'x co-ord. of the end of the segment ? '
     *,T2,1PE14.6)
      READ(5,*)SEG(4,I)
      WRITE(6,8001)
8001  FORMAT(T20,'r co-ord. of the end of the segment ? '
     *,T2,1PE14.6)
      READ(5,*)SEG(5,I)
      WRITE(7,8)SEG(4,I)
      WRITE(7,8001)SEG(5,I)
      IF(SEG(1,I).EQ.1)GOTO 120
      WRITE(6,9)
9     FORMAT(T20,'x co-ord. of the centre of curvature ? '
     *,T2,1PE14.6)
      READ(5,*)SEG(6,I)
      WRITE(7,9)SEG(6,I)
      WRITE(6,9001)
9001  FORMAT(T20,'r co-ord. of the centre of curvature ? '
     *,T2,1PE14.6)
      READ(5,*)SEG(7,I)
      WRITE(7,9001)SEG(7,I)
      WRITE(6,10)
10    FORMAT(T20,'Arc direction (+1 if counter-clockwise, -1 if
     * clockwise) ?  ',T2,I5)
      READ(5,*)IROT(I)
      WRITE(7,10)IROT(I)
120   CONTINUE
      WRITE(6,11)
11    FORMAT(t20,'Is the composition of this segment new',
     *' (Y/N) ? ',t2,5X,A1)
      READ(5,111)ANSWER
      WRITE(7,11)ANSWER
111   FORMAT(A1)
      IF(ANSWER.EQ.'N'.OR.ANSWER.EQ.'n')THEN
      WRITE(6,222)
222   FORMAT(t20,'Same as which segment ? ',t2,I5)
      READ(5,*)NWORK
      WRITE(7,222)NWORK
      SEG(8,I)=SEG(8,NWORK)
      DO  J=1,SEG(8,I)
      DO  K=1,17
         ALAYER(K,J,I)=ALAYER(K,J,NWORK)
      end do
      end do
      GOTO 100
      END IF
      WRITE(6,12)
12    FORMAT(t20,'Number of layers ? ',t2,I5)
      READ(5,*)LAYERS
      WRITE(7,12)LAYERS
      SEG(8,I)=LAYERS
      DO 200 J=1,LAYERS
      WRITE(6,131)J
131   FORMAT(t20,'Layer : ',i2,' New properties (Y/N) ? ',t2,5x,A1)
      READ(5,111)ANSWER
      WRITE(7,131)J,ANSWER
      IF(ANSWER.EQ.'N'.OR.ANSWER.EQ.'n')THEN
      WRITE(6,132)
132   FORMAT(t20,'Same as which layer ? ',t2,I5)
      READ(5,*)LAYER
      WRITE(7,132)LAYER
      DO  K=1,17
         ALAYER(K,J,I)=ALAYER(K,LAYER,I)
      end do
      GOTO 199
      END IF
      WRITE(6,121)
121   FORMAT(T20,'Density ? ',T2,1PE14.6)
      READ(5,*)ALAYER(17,J,I)
      WRITE(7,121)ALAYER(17,J,I)
      WRITE(6,13)
13    FORMAT(t20,'Thickness ? ',t2,1PE14.6)
      READ(5,*)ALAYER(1,J,I)
      WRITE(7,13)ALAYER(1,J,I)
      WRITE(6,14)
14    FORMAT(t20,'Youngs Modulus, E11 ? ',t2,1PE14.6)
      READ(5,*)ALAYER(2,J,I)
      WRITE(7,14)ALAYER(2,J,I)
      WRITE(6,15)
15    FORMAT(t20,'Youngs Modulus, E22 ? ',t2,1PE14.6)
      READ(5,*)ALAYER(3,J,I)
      WRITE(7,15)ALAYER(3,J,I)
      WRITE(6,16)
16    FORMAT(T20,'Poissons ratio (nu12, smaller of the two) ? '
     * ,T2,1PE14.6)
      READ(5,*)ALAYER(4,J,I)
      WRITE(7,16)ALAYER(4,J,I)
      WRITE(6,17)
17    FORMAT(T20,'Shear Modulus, G12 ? ',T2,1PE14.6)
      READ(5,*)ALAYER(5,J,I)
      WRITE(7,17)ALAYER(5,J,I)
      WRITE(6,18)
18    FORMAT(T20,'Ultimate strength in tension, X1T ? ',T2,1PE14.6)
      READ(5,*)ALAYER(6,J,I)
      WRITE(7,18)ALAYER(6,J,I)
      WRITE(6,19)
19    FORMAT(T20,'Ultimate strength in compression, X1C ? ',
     * T2,1PE14.6)
      READ(5,*)ALAYER(7,J,I)
      WRITE(7,19)ALAYER(7,J,I)
      WRITE(6,20)
20    FORMAT(T20,'Ultimate strength in tension, X2T ? ',T2,1PE14.6)
      READ(5,*)ALAYER(8,J,I)
      WRITE(7,20)ALAYER(8,J,I)
      WRITE(6,21)
21    FORMAT(T20,'Ultimate strength in compression, X2C ? ',
     * T2,1PE14.6)
      READ(5,*)ALAYER(9,J,I)
      WRITE(7,21)ALAYER(9,J,I)
      WRITE(6,22)
22    FORMAT(T20,'Shear strength, S12 ? ',T2,1PE14.6)
      READ(5,*)ALAYER(10,J,I)
      WRITE(7,22)ALAYER(10,J,I)
      WRITE(6,60)
60    FORMAT(T20,'Ultimate strain in tension, E1T ? ',T2,1PE14.6)
      READ(5,*)ALAYER(12,J,I)
      WRITE(7,60)ALAYER(12,J,I)
      WRITE(6,61)
61    FORMAT(T20,'Ultimate strain in compression E1C ? ',T2,
     *  1PE14.6)
      READ(5,*)ALAYER(13,J,I)
      WRITE(7,61)ALAYER(13,J,I)
      WRITE(6,62)
62    FORMAT(T20,'Ultimate strain in tension, E2T ? ',T2,1PE14.6)
      READ(5,*)ALAYER(14,J,I)
      WRITE(7,62)ALAYER(14,J,I)
      WRITE(6,63)
63    FORMAT(T20,'Ultimate strain in compression, E2C ? ',T2,
     *  1PE14.6)
      READ(5,*)ALAYER(15,J,I)
      WRITE(7,63)ALAYER(15,J,I)
      WRITE(6,64)
64    FORMAT(T20,'Ultimate in-plane shear strain, E12 ? ',T2,
     *  1PE14.6)
      READ(5,*)ALAYER(16,J,I)
      WRITE(7,64)ALAYER(16,J,I)
199   CONTINUE
      WRITE(6,23)
23    FORMAT(T20,'Winding angle ? ',T2,1PE14.6)
      READ(5,*)ALAYER(11,J,I)
      WRITE(7,23)ALAYER(11,J,I)

200   CONTINUE
100   CONTINUE

      JUNCS=NSEGS+1
      IF(ICASE.EQ.4)GOTO 1000

      DO 300 I=1,NSEGS
      WRITE(6,24)I
24    FORMAT(T20,'*****Segment ',i2,'*****')
      WRITE(7,24)I
      WRITE(6,241)
241   FORMAT(T20,'Pressure multiplier (at start) ? ',T2,1PE14.6)
      READ(5,*)SEG(12,I)
      WRITE(7,241)SEG(12,I)
      WRITE(6,242)
242   FORMAT(T20,'Pressure multiplier (at finish) ? ',T2,1PE14.6)
      READ(5,*)SEG(13,I)
      WRITE(7,242)SEG(13,I)
300   CONTINUE

      DO 350 I=1,JUNCS
      WRITE(6,351)I
      WRITE(7,351)I
351   FORMAT(t20,'Ring load multipliers at junction ',i2)
      WRITE(6,352)
352   FORMAT(T20,'Axial load multiplier ? ',T2,1PE14.6)
      READ(5,*)RING(1,I)
      WRITE(7,352)RING(1,I)
      WRITE(6,353)
353   FORMAT(T20,'Circum. load multiplier ? ',T2,1PE14.6)
      READ(5,*)RING(2,I)
      WRITE(7,353)RING(2,I)
      WRITE(6,354)
354   FORMAT(T20,'Radial load multiplier ? ',T2,1PE14.6)
      READ(5,*)RING(3,I)
      WRITE(7,354)RING(3,I)
      WRITE(6,355)
355   FORMAT(T20,'Moment load multiplier ? ',T2,1PE14.6)
      READ(5,*)RING(4,I)
      WRITE(7,355)RING(4,I)
350   CONTINUE
      WRITE(6,25)
25    FORMAT(T20,'Number of Load/pressure history points ? ',T2,I5)
      READ(5,*)NHIST
      WRITE(7,25)NHIST
      IF(NHIST.NE.0)THEN
      DO 400 I=1,NHIST
      WRITE(6,26)I
26    FORMAT(T20,'Point ',I2,' load ',T2,1PE14.6)
      READ(5,*)HIST(I,1)
      WRITE(7,26)I,HIST(I,1)
      WRITE(6,261)
261   FORMAT(T27,'Increments ? ',T2,F7.1)
      READ(5,*)HIST(I,2)
      WRITE(7,261)HIST(I,2)
400   CONTINUE
      END IF
      WRITE(6,27)
27    FORMAT(T20,'Rotational speed (rad/s) of the shell about the '
     *,'axis ? ',T2,1PE14.6)
      READ(5,*)WIRL
      WRITE(7,27)WIRL

1000  CONTINUE

      DO 500 I=1,JUNCS
      DO  J=1,4
         BCOND(1,J,I)=1.
      end do
      WRITE(6,28)I
28    FORMAT(T20,'Are there any constraints at junction ',
     *  i2,' ? ',T2,5X,A1)
      READ(5,111)ANSWER
      WRITE(7,28)I,ANSWER
      IF(ANSWER.EQ.'Y'.OR.ANSWER.EQ.'y')THEN
      DO 502 J=1,4
      WRITE(6,29)CDISP(J)
29    FORMAT(T20,A8,' :( 1 free, 0 constrained) ? ',T2,I5)
      READ(5,*)K
      WRITE(7,29)CDISP(J),K
      IF(K.EQ.0)THEN
      BCOND(1,J,I)=0
      WRITE(6,30)cdisp(j)
30    FORMAT(T20,'Prescribed value of ',a8,' ? ',T2,1PE14.6)
      READ(5,*)BCOND(2,J,I)
      WRITE(7,30)cdisp(j),BCOND(2,J,I)
      END IF
502   CONTINUE
      END IF
500   CONTINUE

      END
      SUBROUTINE READ_DATA(JOBSES,JOBSES1)
C---------------------------READ DATA FROM 7 ----------

      INCLUDE 'SHELL.BLK'
      CHARACTER*8 CDISP(4),ANSWER*1,JOBSES*14,BLANK*1,JOBSES1*14
      DATA CDISP/'       u','       v','       w','Rotation'/

      OPEN (UNIT=7,FILE=JOBSES,STATUS='unknown',FORM='FORMATTED')
      OPEN (UNIT=9,FILE=JOBSES1,STATUS='unknown',FORM='FORMATTED')
1111  FORMAT(I6)
2222  FORMAT(E15.6)
3333  FORMAT(6X,A1)
4444  FORMAT(F8.1)
50    FORMAT(T20,'Analysis (1:static, 2:stability, 3:',
     *'buckling, 4:fequency) ? ',T2,I5)
      READ(7,1111)ICASE
      WRITE(6,50)ICASE
      WRITE(9,50)ICASE

70    FORMAT(T20,'Do you want to print strains & stresses (Y/N) ? '
     * ,T2,5X,A1)
      IF(ICASE.NE.4)THEN
      READ(7,3333)ANSWER
      WRITE(6,70)ANSWER
      WRITE(9,70)ANSWER
      IPRINT=0
      IF(ANSWER.EQ.'Y'.OR.ANSWER.EQ.'y')IPRINT=1
71    FORMAT(T20,'Failure criterion ? (1:Max. stress, 2:Max. strain,'
     * ,/,T30,'3:Tsai-Hill, 4:Hoffman, 5:Tsai-Wu stress,',/,
     * T30,'6:Tsai-Wu strain, 7:Owen, 8:All) ? ',T2,I5)
      READ(7,111)BLANK
      READ(7,111)BLANK
      READ(7,1111)CRITERION
      WRITE(6,71)CRITERION
      WRITE(9,71)CRITERION
      END IF

      IF(ICASE.NE.1)THEN
51    FORMAT(T20,'Starting wave number (n) ? ',T2,I5)
      READ(7,1111)N0
      WRITE(9,51)N0
      WRITE(6,51)N0
52    FORMAT(T20,'Minimum wave number ? ',T2,I5)
      READ(7,1111)NMIN
      WRITE(9,52)NMIN
      WRITE(6,52)NMIN
53    FORMAT(T20,'Maximum wave number ? ',T2,I5)
      READ(7,1111)NMAX
      WRITE(9,53)NMAX
      WRITE(6,53)NMAX
54    FORMAT(T20,'Number of eigenvalues per n ? ',T2,I5)
      READ(7,1111)MAXROOTS
      WRITE(9,54)MAXROOTS
      WRITE(6,54)MAXROOTS
      END IF

1     FORMAT(T20,'Number of segments ? ',T2,I5)
      READ(7,1111)NSEGS
      WRITE(9,1)NSEGS
      WRITE(6,1)NSEGS

3     FORMAT(T20,'Starting x co-ord. of the structure ? ',T2,
     *  1PE14.6)
      READ(7,2222)CORDS(1,1)
3001  FORMAT(T20,'Starting r co-ord. of the structure ? ',T2,
     *  1PE14.6)
      READ(7,2222)CORDS(2,1)
      WRITE(6,3)CORDS(1,1)
      WRITE(9,3)CORDS(1,1)
      WRITE(6,3001)CORDS(2,1)
      WRITE(9,3001)CORDS(2,1)

      DO 100 I=1,NSEGS
      READ(7,111)BLANK
      WRITE(6,4)I
      WRITE(9,4)I
4     FORMAT(T20,'************Segment No.',I2,'************')
5     FORMAT(T20,'Type (1: cone/cylinder, 2: toroid/sphere) ? '
     * ,T2,F7.1)
      READ(7,4444)SEG(1,I)
      WRITE(6,5)SEG(1,I)
      WRITE(9,5)SEG(1,I)
6     FORMAT(T20,'Number of elements in the segment ? ',T2,F7.1)
      READ(7,4444)SEG(2,I)
      WRITE(6,6)SEG(2,I)
      WRITE(9,6)SEG(2,I)
7     FORMAT(T20,'Step ratio (>1 for incr. and <1 for decr) ? '
     *,T2,1PE14.6)
      READ(7,2222)SEG(3,I)
      WRITE(6,7)SEG(3,I)
      WRITE(9,7)SEG(3,I)
8     FORMAT(T20,'x co-ord. of the end of the segment ? '
     *,T2,1PE14.6)
      READ(7,2222)SEG(4,I)
8001  FORMAT(T20,'r co-ord. of the end of the segment ? '
     *,T2,1PE14.6)
      READ(7,2222)SEG(5,I)
      WRITE(6,8)SEG(4,I)
      WRITE(9,8)SEG(4,I)
      WRITE(6,8001)SEG(5,I)
      WRITE(9,8001)SEG(5,I)
      IF(SEG(1,I).EQ.1)GOTO 120
9     FORMAT(T20,'x co-ord. of the centre of curvature ? '
     *,T2,1PE14.6)
      READ(7,2222)SEG(6,I)
      WRITE(6,9)SEG(6,I)
      WRITE(9,9)SEG(6,I)
9001  FORMAT(T20,'r co-ord. of the centre of curvature ? '
     *,T2,1PE14.6)
      READ(7,2222)SEG(7,I)
      WRITE(6,9001)SEG(7,I)
      WRITE(9,9001)SEG(7,I)
10    FORMAT(T20,'Arc direction (+1 if counter-clockwise, -1 if
     * clockwise) ?  ',T2,I5)
      READ(7,1111)IROT(I)
      WRITE(6,10)IROT(I)
      WRITE(9,10)IROT(I)
120   CONTINUE
11    FORMAT(t20,'Is the composition of this segment new ? '
     *  ,t2,5X,A1)
      READ(7,3333)ANSWER
      WRITE(6,11)ANSWER
      WRITE(9,11)ANSWER
111   FORMAT(A1)
      IF(ANSWER.EQ.'N'.OR.ANSWER.EQ.'n')THEN
222   FORMAT(t20,'Same as which segment ? ',t2,I5)
      READ(7,1111)NWORK
      WRITE(6,222)NWORK
      WRITE(9,222)NWORK
      SEG(8,I)=SEG(8,NWORK)
      DO  J=1,SEG(8,I)
      DO  K=1,17
         ALAYER(K,J,I)=ALAYER(K,J,NWORK)
      end do
      end do
      GOTO 100
      END IF
12    FORMAT(t20,'Number of layers ? ',t2,I5)
      READ(7,1111)LAYERS
      WRITE(6,12)LAYERS
      WRITE(9,12)LAYERS
      SEG(8,I)=LAYERS
      DO 200 J=1,LAYERS
131   FORMAT(t20,'Layer : ',i2,' New properties (Y/N) ? ',t2,5X,A1)
      READ(7,3333)ANSWER
      WRITE(6,131)J,ANSWER
      WRITE(9,131)J,ANSWER
      IF(ANSWER.EQ.'N'.OR.ANSWER.EQ.'n')THEN
132   FORMAT(t20,'Same as which layer ? ',t2,I5)
      READ(7,1111)LAYER
      WRITE(6,132)LAYER
      WRITE(9,132)LAYER
      DO  K=1,17
           ALAYER(K,J,I)=ALAYER(K,LAYER,I)
      end do
      GOTO 199
      END IF
121   FORMAT(T20,'Density ? ',T2,1PE14.6)
      READ(7,2222)ALAYER(17,J,I)
      WRITE(6,121)ALAYER(17,J,I)
      WRITE(9,121)ALAYER(17,J,I)
c     WRITE(7,121)ALAYER(17,J,I)
13    FORMAT(t20,'Thickness ? ',t2,1PE14.6)
      READ(7,2222)ALAYER(1,J,I)
      WRITE(6,13)ALAYER(1,J,I)
      WRITE(9,13)ALAYER(1,J,I)
14    FORMAT(t20,'Youngs Modulus, E11 ? ',t2,1PE14.6)
      READ(7,2222)ALAYER(2,J,I)
      WRITE(6,14)ALAYER(2,J,I)
      WRITE(9,14)ALAYER(2,J,I)
15    FORMAT(t20,'Youngs Modulus, E22 ? ',t2,1PE14.6)
      READ(7,2222)ALAYER(3,J,I)
      WRITE(6,15)ALAYER(3,J,I)
      WRITE(9,15)ALAYER(3,J,I)
16    FORMAT(T20,'Poissons ratio (nu12, smaller of the two) ? '
     * ,T2,1PE14.6)
      READ(7,2222)ALAYER(4,J,I)
      WRITE(6,16)ALAYER(4,J,I)
      WRITE(9,16)ALAYER(4,J,I)
17    FORMAT(T20,'Shear Modulus, G12 ? ',T2,1PE14.6)
      READ(7,2222)ALAYER(5,J,I)
      WRITE(6,17)ALAYER(5,J,I)
      WRITE(9,17)ALAYER(5,J,I)
18    FORMAT(T20,'Ultimate strength in tension, X1T ? ',T2,1PE14.6)
      READ(7,2222)ALAYER(6,J,I)
      WRITE(6,18)ALAYER(6,J,I)
      WRITE(9,18)ALAYER(6,J,I)
19    FORMAT(T20,'Ultimate strength in compression, X1C ? ',
     * T2,1PE14.6)
      READ(7,2222)ALAYER(7,J,I)
      WRITE(6,19)ALAYER(7,J,I)
      WRITE(9,19)ALAYER(7,J,I)
20    FORMAT(T20,'Ultimate strength in tension, X2T ? ',T2,1PE14.6)
      READ(7,2222)ALAYER(8,J,I)
      WRITE(6,20)ALAYER(8,J,I)
      WRITE(9,20)ALAYER(8,J,I)
21    FORMAT(T20,'Ultimate strength in compression, X2C ? ',
     * T2,1PE14.6)
      READ(7,2222)ALAYER(9,J,I)
      WRITE(6,21)ALAYER(9,J,I)
      WRITE(9,21)ALAYER(9,J,I)
22    FORMAT(T20,'Shear strength, S12 ? ',T2,1PE14.6)
      READ(7,2222)ALAYER(10,J,I)
      WRITE(6,22)ALAYER(10,J,I)
      WRITE(9,22)ALAYER(10,J,I)
60    FORMAT(T20,'Ultimate strain in tension, E1T ? ',T2,1PE14.6)
      READ(7,2222)ALAYER(12,J,I)
      WRITE(6,60)ALAYER(12,J,I)
      WRITE(9,60)ALAYER(12,J,I)
61    FORMAT(T20,'Ultimate strain in compression E1C ? ',T2,1PE14.6)
      READ(7,2222)ALAYER(13,J,I)
      WRITE(6,61)ALAYER(13,J,I)
      WRITE(9,61)ALAYER(13,J,I)
62    FORMAT(T20,'Ultimate strain in tension, E2T ? ',T2,1PE14.6)
      READ(7,2222)ALAYER(14,J,I)
      WRITE(6,62)ALAYER(14,J,I)
      WRITE(9,62)ALAYER(14,J,I)
63    FORMAT(T20,'Ultimate strain in compression, E2C ? ',T2,
     *  1PE14.6)
      READ(7,2222)ALAYER(15,J,I)
      WRITE(6,63)ALAYER(15,J,I)
      WRITE(9,63)ALAYER(15,J,I)
64    FORMAT(T20,'Ultimate in-plane shear strain, E12 ? ',T2,
     *  1PE14.6)
      READ(7,2222)ALAYER(16,J,I)
      WRITE(6,64)ALAYER(16,J,I)
      WRITE(9,64)ALAYER(16,J,I)
199   CONTINUE
23    FORMAT(T20,'Winding angle ? ',T2,1PE14.6)
      READ(7,2222)ALAYER(11,J,I)
      WRITE(6,23)ALAYER(11,J,I)
      WRITE(9,23)ALAYER(11,J,I)

200   CONTINUE
100   CONTINUE

      JUNCS=NSEGS+1
      IF(ICASE.EQ.4)GOTO 1000

      DO 300 I=1,NSEGS
      READ(7,111)BLANK
24    FORMAT(T20,'*****Segment ',i2,'*****')
      WRITE(6,24)I
      WRITE(9,24)I
241   FORMAT(T20,'Pressure multiplier (at start) ? ',T2,1PE14.6)
      READ(7,2222)SEG(12,I)
      WRITE(6,241)SEG(12,I)
      WRITE(9,241)SEG(12,I)
242   FORMAT(T20,'Pressure multiplier (at finish) ? ',T2,1PE14.6)
      READ(7,2222)SEG(13,I)
      WRITE(6,242)SEG(13,I)
      WRITE(9,242)SEG(13,I)
300   CONTINUE

      DO 350 I=1,JUNCS
      READ(7,111)BLANK
      WRITE(6,351)I
      WRITE(9,351)I
351   FORMAT(t20,'Ring load multipliers at junction ',i2)
352   FORMAT(T20,'Axial load multiplier ? ',T2,1PE14.6)
      READ(7,2222)RING(1,I)
      WRITE(6,352)RING(1,I)
      WRITE(9,352)RING(1,I)
353   FORMAT(T20,'Circum. load multiplier ? ',T2,1PE14.6)
      READ(7,2222)RING(2,I)
      WRITE(6,353)RING(2,I)
      WRITE(9,353)RING(2,I)
354   FORMAT(T20,'Radial load multiplier ? ',T2,1PE14.6)
      READ(7,2222)RING(3,I)
      WRITE(6,354)RING(3,I)
      WRITE(9,354)RING(3,I)
355   FORMAT(T20,'Moment load multiplier ? ',T2,1PE14.6)
      READ(7,2222)RING(4,I)
      WRITE(6,355)RING(4,I)
      WRITE(9,355)RING(4,I)
350   CONTINUE
25    FORMAT(T20,'Number of Load/pressure history points ? ',T2,I5)
      READ(7,1111)NHIST
      WRITE(6,25)NHIST
      WRITE(9,25)NHIST
      IF(NHIST.NE.0)THEN
      DO 400 I=1,NHIST
26    FORMAT(T20,'Point ',I2,' load ',T2,1PE14.6)
      READ(7,2222)HIST(I,1)
      WRITE(6,26)I,HIST(I,1)
      WRITE(9,26)I,HIST(I,1)
261   FORMAT(T27,'Increments ? ',T2,F7.1)
      READ(7,4444)HIST(I,2)
      WRITE(6,261)HIST(I,2)
      WRITE(9,261)HIST(I,2)
400   CONTINUE
      END IF
27    FORMAT(T20,'Rotational speed (rad/s) of the shell about the '
     *,'axis ? ',T2,1PE14.6)
      READ(7,2222)WIRL
      WRITE(6,27)WIRL
      WRITE(9,27)WIRL

1000  CONTINUE

      DO 500 I=1,JUNCS
      DO  J=1,4
         BCOND(1,J,I)=1.
      end do
28    FORMAT(T20,'Are there any constraints at junction ',
     *  i2,' ? ',T2,5X,A1)
      READ(7,3333)ANSWER
      WRITE(6,28)I,ANSWER
      WRITE(9,28)I,ANSWER
      IF(ANSWER.EQ.'Y'.OR.ANSWER.EQ.'y')THEN
      DO 502 J=1,4
29    FORMAT(T20,A8,' :( 1 free, 0 constrained) ? ',T2,I5)
      READ(7,1111)K
      WRITE(6,29)CDISP(J),K
      WRITE(9,29)CDISP(J),K
      IF(K.EQ.0)THEN
      BCOND(1,J,I)=0
30    FORMAT(T20,'Prescribed value of ',a8,' ? ',T2,1PE14.6)
      READ(7,2222)BCOND(2,J,I)
      WRITE(6,30)CDISP(J),BCOND(2,J,I)
      WRITE(9,30)CDISP(J),BCOND(2,J,I)
      END IF
502   CONTINUE
      END IF
500   CONTINUE

      END
