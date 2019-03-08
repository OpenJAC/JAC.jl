************************************************************************
*                                                                      *
      SUBROUTINE BUBBLE (JPOL,FAIL)
*                                                                      *
*   Reduces a circuit of order 2 , giving  delta  function and phase   *
*   factors.                                                           *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME,NAMSUB
      INTEGER ARR,TAB1
      LOGICAL FAIL,SUMVAR
*
      PARAMETER (
     :   MANGM = 60, M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/NAM/NAMSUB
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      DATA NAME/'BUBBLE'/
*
      NAMSUB = NAME
      K2 = 2
      K23 = 3
      I1 = 1
      I2 = 1
      IT1 = NPOINT(1)
      IT2 = NPOINT(2)
*
      IF (IT2 .EQ. ILAST) THEN
         IF (IT1 .NE. IFIRST) THEN
            IT2 = IT1
            IT1 = ILAST
         ENDIF
         I1 = -1
         K23 = 2
         I2 = 2
      ENDIF
*
      CALL PHASE (IT1,JDIAG,M4TRD)
      K = ABS ((3*ARR(IT2,1)+2*ARR(IT2,2)+ARR(IT2,3))/2)+1
      IF (K .NE. 4) CALL PHASE2 (JDIAG(IT2,K))
      IF (NBNODE .EQ. 2) RETURN
      IL1 = IL(IT2)+I1
      IT = IH(IL1)
      ARR(IT,K23) = ARR(IT1,K23)
      L = JDIAG(IT1,K23)
      L1 = JDIAG(IT,K23)
      JDIAG(IT,K23) = L
*
      IF (JPOL .NE. 1) THEN
         CALL DELTA (L,L1,FAIL)
         IF (FAIL) RETURN
      ELSE
         MP = MP-1
         KW(2,JWC) = L
         J6(J6C-1) = L
         J6(J6C) = L
         IF (K .EQ. 2) J8(J8C) = L
      ENDIF
*
      TAB1(L,I2) = IT
*
      IF (IT1 .NE. ILAST) THEN
         IF (IT2 .EQ. ILAST) THEN
            TAB1(L,1) = IH(2)
            IL1 = 2
            K2 = 1
         ENDIF
*
         DO 1 I = IL1,NBNODE
            IT = IH(I)
            IL(IT) = I-K2
            IH(I-K2) = IT
    1    CONTINUE
*
      ENDIF
*
      J9(J9C+1) = L
      J9C = J9C+2
      J9(J9C) = L
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE CHANGE (L,K)
*                                                                      *
*   Exchanges the free ends in either first or last triad of JDIAG.    *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INTEGER ARR,TAB1
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM,
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,    MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
*
      CALL PHASE (L,JDIAG,M4TRD)
      JP = JDIAG(L,K)
      JDIAG(L,K) = JDIAG(L,1)
      JDIAG(L,1) = JP
      JAR = ARR(L,K)
      ARR(L,K) = ARR(L,1)
      ARR(L,1) = JAR
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE CHKLP1 (FAIL)
*                                                                      *
*   This routine checks if  there are active triads with two identi-   *
*   cal  arguments.  This is a loop of  order 1 (a "lollypop").  The   *
*   other argument  must then be zero;  i.e. j1(j) = 1 in 2j+1 nota-   *
*   tion. Suppression of the loop introduces factors and phases. Two   *
*   Two triads  become inactive.  All this is  performed by invoking   *
*   ZERO with first argument 1.                                        *
*                                                                      *
*   Written by  Marcel Klapisch for correcting  an error detected by   *
*   Charlotte F Fischer.    This version includes Marcel's fix of 12   *
*   June 1992.                                                         *
*                                         Last revision: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INTEGER ARROW
      LOGICAL FAIL,FREE,TABS
      CHARACTER*6 NAME,NAMSUB
*
      PARAMETER (MANGM = 60,MTRIAD = 12,M2TRD = 2*MTRIAD)
*
      COMMON/COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /NAM/NAMSUB
     :      /TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
*
      DATA NAME/'CHKLP1'/
*
      NAMSUB = NAME
      NBTR1 = 2*(N-1)
      DO 10 L = 1,NBTR1
      IF (.NOT. TABS(L)) THEN
         JDIF = 0
         IF     (J23(L,1) .EQ. J23(L,2)) THEN
            JDIF = J23(L,3)
         ELSEIF (J23(L,1) .EQ. J23(L,3)) THEN
            JDIF = J23(L,2)
         ELSEIF (J23(L,2) .EQ. J23(L,3)) THEN
            JDIF = J23(L,1)
         ENDIF
         IF (JDIF .NE. 0) THEN
*
*   Putting the link to 0. ZERO changes NBTR
*
            FAIL = .FALSE.
            IF (J1(JDIF) .NE. 1.AND. .NOT. FREE(JDIF)) THEN
               FAIL = .TRUE.
               IF (IBUG3 .EQ. 1) WRITE (99,300) JDIF,J1(JDIF)
               RETURN
            ELSE
	       CALL ZERO (1,JDIF,FAIL)
               IF (FAIL) RETURN
            ENDIF
         ENDIF
      ENDIF
  10  CONTINUE
      IF (JDIF .NE. 0) CALL PRINTJ (NAME,4)
*
      RETURN
*
  300 FORMAT (1X,'JDIF = ',1I2,'; should be 0; J1(JDIF) = ',1I2,
     :           '; RECUP -> 0.')
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE CHVAR (JP,NBC,KBC,JT,JINV,NSUM)
*                                                                      *
*   Change  the  order  of  summation variable to be able to perform   *
*   separately the summations in GENSUM.                               *
*                                                                      *
*                                           Last update: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      LOGICAL JT
*
      DIMENSION JINV(NSUM),JP(NBC),JT(NSUM)
*
      KB = KBC+1
*
      IF (KB .LE. NBC) THEN
         DO 1 I = KB,NBC
            JK = JP(I)
            IF (JT(JK)) THEN
               KBC = KBC+1
               JP(I) = JP(KBC)
               JP(KBC) = JINV(JK)
            ENDIF
    1    CONTINUE
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE CUT1L (FAIL)
*                                                                      *
*   Cut  on  one  line, that  was  left as a free end in JDIAG. Puts   *
*   corresponding delta in J23.                                        *
*                                                                      *
*                                           Last update: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARR,TAB1
      LOGICAL FAIL,SUMVAR,FREE
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
*
      DATA NAME/'CUT1L '/
*
      IT = ITFREE(1)
      J0 = JDIAG(IT,1)
      CALL DELTA (J0,M,FAIL)
      IF (FAIL) GOTO 2
      CALL DELTA (JDIAG(IT,3),JDIAG(IT,2),FAIL)
      IF (FAIL) GOTO 2
      JDIAG(IT+1,3) = JDIAG(IT,3)
*
      IF (ARR(IT,2) .EQ. ARR(IT,3)) THEN
         ARR(IT+1,3) = 1
         ARR(IT-1,2) = -1
      ELSEIF (ARR(IT,2) .LT. ARR(IT,3)) THEN
         ARR(IT+1,3) = -1
         ARR(IT-1,2) = 1
      ENDIF
*
      J9C = J9C+1
      J9(J9C) = JDIAG(IT,3)
      J = 2
      CALL ZERO (J,J0,FAIL)
      IF (FAIL) GOTO 2
      IL1 = IL(IT+1)
*
      DO 1 I = IL1,NBNODE
         IT = IH(I)
         ILP = I-1
         IL(IT) = ILP
         IH(ILP) = IT
    1 CONTINUE
*
      NBNODE = NBNODE-1
*
    2 CALL PRINTJ (NAME,12)
      RETURN
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE CUT2L (FAIL)
*                                                                      *
*   Cut on two lines that were left as free ends in JDIAG. Puts cor-   *
*   responding delta in J23.                                           *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARR,TAB1,ARROW
      LOGICAL FAIL,TABS,SUMVAR
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,    MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      DATA NAME/'CUT2L '/
*
*
      IT1 = ITFREE(1)
      IT2 = ITFREE(2)
      JT1 = JDIAG(IT1,1)
      JT2 = JDIAG(IT2,1)
      CALL DELTA (JT1,JT2,FAIL)
      IF (FAIL) GOTO 1
      IF (ARR(IT1,1) .EQ. ARR(IT2,1)) CALL PHASE2 (JT1)
      ARR(IT2,1) = -ARR(IT1,1)
      JDIAG(IT2,1) = JT1
      TAB1(JT1,2) = IT2
      J9(J9C+1) = JT1
      J9C = J9C+2
      J9(J9C) = JT1
      CALL OTHERJ (0,JT1,L1,LC1,K1)
      CALL OTHERJ (0,JT2,L2,LC2,K2)
      J23(L2,LC2) = JT1
      LINE(JT1,K1) = L2
      LCOL(JT1,K1) = LC2
      ARROW(L2,LC2) = -ARROW(L1,LC1)
*
    1 CALL PRINTJ (NAME,12)
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE CUTNL (FAIL)
*                                                                      *
*   This subroutine  examines the case where there are more than two   *
*   free ends, but they are contiguous, so that the graph can be cut   *
*   without destroying the flat structure.                             *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARROW,ARR,TAB1
      LOGICAL TABS,SUMVAR,FAIL
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /KEEP/JKP(2,3),JARR(2,3),IT2,IT3,IT5
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      DATA NAME/'CUTNL '/
*
      NTF = ITFREE(NFREE)-ITFREE(1)
      IF (NTF .GT. NFREE) GOTO 8
      IT2 = ITFREE(1)
      IT3 = ITFREE(NFREE)
      IT1 = IT2-1
      IT4 = IT3+1
*
      IF (NTF .NE. NFREE) THEN
*
         JT = JDIAG(IT2,3)
         CALL DELTA (JT,JDIAG(IT3,2),FAIL)
*
         IF (FAIL) GOTO 9
*
         IF (ARR(IT2,3) .EQ. ARR(IT3,2)) THEN
            CALL PHASE2 (JT)
            ARR(IT2,3) = -ARR(IT2,3)
            ARR(IT1,2) = -ARR(IT1,2)
         ENDIF
*
         JDIAG(IT3,2) = JT
         JDIAG(IT4,3) = JT
         J9(J9C+1) = JT
         J9C = J9C+2
         J9(J9C) = JT
         NBTR = NBTR+NFREE
         IT5 = 0
*
      ELSE
*
         NFR = 0
*
         DO 3 IT5 = IT2,IT3
            NFR = NFR+1
            IF (ITFREE(NFR) .GT. IT5) GOTO 4
    3    CONTINUE
*
    4    JKP(1,1) = JDIAG(IT5,1)
         JARR(1,1) = -ARR(IT5,1)
         JKP(1,2) = JDIAG(IT2,3)
         JARR(1,2) = -ARR(IT2,3)
         JKP(1,3) = JDIAG(IT3,2)
         JARR(1,3) = -ARR(IT3,2)
*
         DO 5 J = 1,3
            JKP(2,J) = JDIAG(IT5,J)
            JARR(2,J) = ARR(IT5,J)
    5    CONTINUE
*
         JDIAG(IT5,2) = JDIAG(IT3,2)
         ARR(IT5,2) = ARR(IT3,2)
         JDIAG(IT5,3) = JDIAG(IT2,3)
         ARR(IT5,3) = ARR(IT2,3)
         ILP = IL(IT2)
         IL(IT5) = ILP
         IH(ILP) = IT5
         NBTR = NBTR+NFREE+2
         CALL PHASE (IT5,JDIAG,M4TRD)
         K = ABS ((3*ARR(IT5,1)+2*ARR(IT5,2)+ARR(IT5,3))/2+1)
         IF (K .NE. 4) CALL PHASE2 (JDIAG(IT5,K))
*
      ENDIF
*
      IL1 = IL(IT4)
*
      DO 7 I = IL1,NBNODE
         IT = IH(I)
         ILP = I-NFREE
         IL(IT) = ILP
         IH(ILP) = IT
    7 CONTINUE
*
      NBNODE = NBNODE-NFREE
      NFIN = 0
      GOTO 8
*
    9 FAIL = .TRUE.
    8 CALL PRINTJ (NAME,8)
*
      RETURN
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE DELTA (JA,JB,FAIL)
*                                                                      *
*   Test for delta(JA,JB). If they are summation variables, the sec-   *
*   ond  is  changed  into  the first everywhere. if they are fixed,   *
*   their value is checked, and fail put to .TRUE. if they differ.     *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      LOGICAL FAIL,CUT,SUMVAR,FREE
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /CUTDIG/CUT
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /DIM/J6CC,J7CC,J8CC,J9CC,JWCC,JDELC
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
*
      IF (IBUG3 .EQ. 1) WRITE (99,300) JA,SUMVAR(JA),JB,SUMVAR(JB)
      IF (SUMVAR(JA) .AND. SUMVAR(JB)) GOTO 2
      IF (FREE(JA) .OR. FREE(JB)) THEN
         JDEL = JDEL+1
         LDEL(JDEL,1) = JA
         LDEL(JDEL,2) = JB
         SUMVAR(JA) = .FALSE.
         SUMVAR(JB) = .FALSE.
         RETURN
      ENDIF
*
      IF (J1(JA) .NE. J1(JB)) FAIL = .TRUE.
      CUT = .TRUE.
      RETURN
*
    2 IF (J6C .NE. J6CC) THEN
         J61 = J6CC+1
*
         DO 3 I = J61,J6C
            IF (J6(I) .EQ. JB) J6(I) = JA
    3    CONTINUE
*
      ENDIF
*
      IF (J7C .NE. J7CC) THEN
         J71 = J7CC+1
*
         DO 5 I = J71,J7C
            IF (J7(I) .EQ. JB) J7(I) = JA
    5    CONTINUE
      ENDIF
*
      IF (J8C .NE. J8CC) THEN
         J81 = J8CC+1
*
         DO 7 I = J81,J8C
            IF (J8(I) .EQ. JB) J8(I) = JA
    7    CONTINUE
      ENDIF
*
      IF (J9C .NE. J9CC) THEN
         J91 = J9CC+1
*
         DO 9 I = J91,J9C
            IF (J9(I) .EQ. JB) J9(I) = JA
    9    CONTINUE
      ENDIF
*
      IF (JWC .NE. JWCC) THEN
         JW1 = JWCC+1
*
         DO 14 I = JW1,JWC
            DO 13 J = 1,6
               IF (KW(J,I) .EQ. JB) KW(J,I) = JA
   13       CONTINUE
   14    CONTINUE
      ENDIF
*
      IF (JDEL .NE. JDELC) THEN
         JDEL1 = JDELC+1
*
         DO 17 I = JDEL1,JDEL
            DO 16 J = 1,2
               IF (LDEL(I,J) .EQ. JB) LDEL(I,J) = JA
   16       CONTINUE
   17    CONTINUE
*
         SUMVAR(JB) = .FALSE.
      ENDIF
*
      RETURN
*
  300 FORMAT (/'From DELTA: JA = ',I2,L2,5X,'JB = ',I2,L2)
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE DIAGRM (JUMP)
*                                                                      *
*   This subroutine builds up a flat diagram from the triads J23 and   *
*   places them in JDIAG . Arrows  are in ARR (INTEGER). The diagram   *
*   is built so as to maximize the number of triads involved, within   *
*   a one-step-forward-check process. If the diagram does not inclu-   *
*   de all the NBTR triads, it will have 'free ends'. JDIAG has dim-   *
*   ension double that of  J23 , because the path may proceed either   *
*   way.                                                               *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARR,TAB1,ARROW
      LOGICAL TABS,FREE,SUMVAR
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /BUILD/IAL(M4TRD),IF1,IF2,NODE
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
*
      DATA NAME/'DIAGRM'/
*
      DATA NB/0/
*
*   Initialization
*
      IF (JUMP .GT. 2) GOTO 17
      IF (JUMP .LT. 2) NB = 0
    1 NB = NB+1
      IF (TABS(NB)) GOTO 1
      NODE = NBTR
      ILAST = NBTR
*
      DO 2 J = 1,3
         JDIAG(NODE,J) = J23(NB,J)
         ARR(NODE,J) = ARROW(NB,J)
    2 CONTINUE
*
      TABS(NB) = .TRUE.
*
      DO 15 I = 1,MP
         IAL(I) = 0
   15 CONTINUE
*
      IF1 = JDIAG(NODE,1)
      IF2 = JDIAG(NODE,3)
      IAL(IF1) = 1
      IAL(IF2) = 1
   17 NTIME = 0
      I1 = 1
      K1 = 1
      K2 = 2
      K3 = 3
    3 JB = JDIAG(NODE,K2)
      CALL OTHERJ (0,JB,L,LC,KP)
      CALL NEIBOR (LC,L1,L2)
*
*   Check consistency of triads
*
      IF (TABS(L)) THEN
         WRITE (*,300)
         STOP
      ENDIF
*
      CALL WAY (L,L1,L2,ICH,ND)
      NODE = NODE+I1
      TABS(L) = .TRUE.
      JDIAG(NODE,K3) = J23(L,LC)
      ARR(NODE,K3) = ARROW(L,LC)
      ICT = ICH*I1
*
      IF (ICH .LE. 0) THEN
         LP = L1
         L1 = L2
         L2 = LP
      ENDIF
*
      IF (ICT .LE. 0) CALL PHASE (L,J23,M2TRD)
      JDIAG(NODE,K1) = J23(L,L1)
      ARR(NODE,K1) = ARROW(L,L1)
      JDIAG(NODE,K2) = J23(L,L2)
      ARR(NODE,K2) = ARROW(L,L2)
      J = J23(L,L1)
      IAL(J) = IAL(J)+1
      J = J23(L,L2)
      IAL(J) = IAL(J)+1
      IF (ND .LT. 1) GOTO 3
      NTIME = NTIME+1
      ILAST = MAX (NODE,ILAST)
      IFIRST = MIN (NODE,NBTR)
      NBP = IAL(IF1)+IAL(IF2)
      IF ((NBP .GT. 3) .OR. (NTIME .GT. 1)) THEN
         NBNODE = ILAST-IFIRST+1
         NBTR = NBTR-NBNODE
*
*   Definition of free ends and other quantities.
*
         CALL INTAB
         CALL PRINTJ (NAME,12)
         GOTO 50
      ENDIF
*
      IF (NBP .GT. 2) THEN
         IF (IAL(IF1) .LE. IAL(IF2)) THEN
            JT = JDIAG(NBTR,1)
            JAR = ARR(NBTR,1)
            JDIAG(NBTR,1) = JDIAG(NBTR,3)
            ARR(NBTR,1) = ARR(NBTR,3)
            JDIAG(NBTR,3) = JT
            ARR(NBTR,3) = JAR
            CALL PHASE (NBTR,JDIAG,M4TRD)
         ENDIF
      ENDIF
*
      NODE = NBTR
      I1 = -1
      K2 = 3
      K3 = 2
      GOTO 3
*
   50 RETURN
*
  300 FORMAT ('DIAGRM: Flat graph impossible to build.')
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE DRACAH (I,J,K,L,M,N,RAC)
*                                                                      *
*   SUBROUTINE  to calculate Racah coefficients. The arguments I, J,   *
*   K, L, M, N should be twice their actual value. Works for integer   *
*   and  half-integer  values of  angular momenta. The routine makes   *
*   use of the GAM  array, thus  SUBROUTINE FACTT must be called be-   *
*   fore this routine is used.                                         *
*                                                                      *
*   Written by N S Scott                    Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      PARAMETER (MFACT = 500)
*
      COMMON/FACTS/GAM(MFACT)
*
      J1 = I+J+M
      J2 = K+L+M
      J3 = I+K+N
      J4 = J+L+N
      IF (((2*MAX (I,J,M)-J1) .GT. 0) .OR. (MOD (J1,2) .NE. 0)) GOTO 2
      IF (((2*MAX (K,L,M)-J2) .GT. 0) .OR. (MOD (J2,2) .NE. 0)) GOTO 2
      IF (((2*MAX (I,K,N)-J3) .GT. 0) .OR. (MOD (J3,2) .NE. 0)) GOTO 2
      IF (((2*MAX (J,L,N)-J4) .GT. 0) .OR. (MOD (J4,2) .NE. 0)) GOTO 2
      GOTO 1
   2  RAC = 0.0D 00
      RETURN
*
   1  CONTINUE
      J1 = J1/2
      J2 = J2/2
      J3 = J3/2
      J4 = J4/2
      J5 = (I+J+K+L)/2
      J6 = (I+L+M+N)/2
      J7 = (J+K+M+N)/2
      NUMIN = MAX (J1,J2,J3,J4)+1
      NUMAX = MIN (J5,J6,J7)+1
      RAC = 1.0D 00
      ICOUNT = 0
*
      IF (NUMIN .EQ. NUMAX) GOTO 4
      NUMIN = NUMIN+1
*
      DO 3 KK = NUMIN,NUMAX
         KI = NUMAX-ICOUNT
         RAC = 1.0D 00
     :       -(RAC*DBLE(KI*(J5-KI+2)*(J6-KI+2)*(J7-KI+2))/
     :      DBLE((KI-1-J1)*(KI-1-J2)*(KI-1-J3)*(KI-1-J4)))
         ICOUNT = ICOUNT+1
   3  CONTINUE
*
      NUMIN = NUMIN-1
   4  RAC = RAC*((-1.0D 00)**(J5+NUMIN+1))
     : *EXP( (GAM(NUMIN+1)-GAM(NUMIN-J1)
     : -GAM(NUMIN  -J2)-GAM(NUMIN  -J3)-GAM(NUMIN  -J4)-GAM(J5+2-NUMIN)
     : -GAM(J6+2-NUMIN)-GAM(J7+2-NUMIN))+((GAM(J1+1-I)+GAM(J1+1-J)
     : +GAM(J1+1-M)-GAM(J1+2)+GAM(J2+1-K)+GAM(J2+1-L)+GAM(J2+1-M)
     : -GAM(J2+2)+GAM(J3+1-I)+GAM(J3+1-K)+GAM(J3+1-N)-GAM(J3+2)
     : +GAM(J4+1-J)+GAM(J4+1-L)+GAM(J4+1-N)-GAM(J4+2))
     : *0.5D 00  ))
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE GENSUM (J6C,J7C,J8C,J9C,JWC,J6,J7,J8,J9,JW,JDEL,
     :                   LDEL,SUMVAR,MP,J6P,J7P,J8P,J9P,JWORD,NLSUM,
     :                   NBJ,NB6J,K6CP,K7CP,K8CP,K9CP,JSUM4,JSUM5,
     :                   JSUM6,INV6J,RECUP)
*                                                                      *
*   Carries  out the summation over coefficients defined by the arr-   *
*   ays J6, J7, J8, LDEL and  JW  to give RECUP. The entry is either   *
*   made from NJGRAF or directly  assuming that the arrays J6,...,JW   *
*   have already been determined  by  a previous entry to NJGRAf and   *
*   that the summation is required for another set of j values defi-   *
*   ned by the array J1. RECUP is the recoupling coefficient.          *
*                                                                      *
*   Call(s) to: [NJGRAF]: DRACAH, RDIAG.                               *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      LOGICAL LDIAG,NOEL,FREE,SUMVAR
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10,MFACT = 500)
*
      PARAMETER (EPSIL = 1.0D-10)
*
      DIMENSION MAT(MTRIAD,MTRIAD),NOEL(MTRIAD),
     :   MAXLP(MTRIAD),JSUM2(MTRIAD),JSUM3(MTRIAD),
     :   JSUM(2,M6J),JWTEST(M6J),WSTOR(M6J),IPAIR(2,2),LDIAG(MTRIAD)
      DIMENSION XJ1(MANGM),IST(6)
      DIMENSION J12(4,MTRIAD,MTRIAD)
      DIMENSION J6P(MANGMP),J7P(MANGMP),J8P(MANGMP),J9P(MANGMP),
     :   JWORD(6,M6J),NBJ(MSUM),NB6J(MSUM),K6CP(MSUM),
     :   K7CP(MSUM),K8CP(MSUM),K9CP(MSUM),JSUM6(MTRIAD),
     :   JSUM4(MTRIAD,M6J),JSUM5(MTRIAD,M6J),INV6J(M6J)
      DIMENSION J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :   J9(MANGMP),JW(6,M6J),LDEL(M6J,2),SUMVAR(MANGM)
*
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /FACTS/GAM(MFACT)
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
*
      DATA MXCSVR/4/
*
*   evaluates all terms in J6, J7, J8, J9, LDEL, JW which do not
*   involve a summation. The result is stored in RECUP and IASTOR
*
      IF (IBUG3 .EQ. 1) THEN
*
         DO 139 I = 1,M
            XJ1(I) = 0.5D 00*DBLE (J1(I)-1)
  139    CONTINUE
*
         WRITE (99,400) (XJ1(I),I = 1,M)
         WRITE (99,306) NLSUM
         WRITE (99,401)
      ENDIF
*
      MM = M+1
      J1(MM) = 1
*
*   Test delta functions
*
      J1(MM) = 1
      IF (JDEL .LE. 0) GOTO 180
*
      DO 141 I = 1,JDEL
         I1 = LDEL(I,1)
         I2 = LDEL(I,2)
         IF ((I1 .GT. MM) .OR. (I2 .GT. MM))THEN
            IF (I1.GT.MM) J1(I1) = J1(I2)
            IF (I2.GT.MM) J1(I2) = J1(I1)
         ELSE
            IF (J1(I1) .NE. J1(I2)) THEN
               RECUP = 0.0D 00
               RETURN
            ENDIF
         ENDIF
  141 CONTINUE
*
  180 RECUP = 1.0D 00
      IF (JWC .NE. 0) THEN
*
*   Multiply RECUP by all Racah coefficients which do not involve a
*   summation
*
         IF (IBUG3 .EQ. 1) WRITE (99,309)
*
         DO 7 I = 1,JWC
            IF (INV6J(I) .GT. 0) GOTO 7
            DO 3 J = 1,6
               I1 = JW(J,I)
               IST(J) = J1(I1) - 1
    3       CONTINUE
*
            CALL DRACAH (IST(1),IST(2),IST(3),IST(4),IST(5),IST(6),X1)
            IF (IBUG3 .EQ. 1) WRITE (99,305) (XJ1(JW(K,I)),K = 1,6),X1
            RECUP = RECUP*X1
*
    7    CONTINUE
*
      ENDIF
*
      SQR = 1.0D 00
*
      IF (J6C .NE. 0) THEN
         DO 12 I = 1,J6C
            I1 = J6(I)
            SQR = SQR*J1(I1)
   12     CONTINUE
      ENDIF
*
      SPR = 1.0D 00
*
      IF (J9C .NE. 0) THEN
         DO 144 I = 1,J9C
            I1 = J9(I)
            SPR = SPR*J1(I1)
  144    CONTINUE
      ENDIF
*
      RECUP = RECUP*SQRT (SQR/SPR)
      IF (ABS(RECUP) .LT. EPSIL) GOTO 145
      IASTOR = 0
*
      IF (J7C .NE. 0) THEN
         DO 17 I = 1,J7C
            I1 = J7(I)
            IASTOR = IASTOR + J1(I1) -1
   17    CONTINUE
      ENDIF
*
      IF (J8C .NE. 0) THEN
         DO 22 I = 1,J8C
            I1 = J8(I)
            IASTOR = IASTOR +2*(J1(I1)-1)
   22    CONTINUE
      ENDIF
*
      IF (NLSUM .LE. 0) THEN
         IASTOR = IASTOR/2
*
*   No summation involved. End of computation
*
         STOR1 = 1.0D 00
         STOR = 1.0D 00
         IF (MOD (IASTOR,2) .EQ. 1) RECUP = -RECUP
         IF (IBUG3 .EQ. 1) WRITE (99,303) RECUP
         RETURN
*
      ENDIF
*
*   Evaluation of the part involving summations.
*
      NFS = 0
      JWR = 0
      J6F = 0
      J7F = 0
      J8F = 0
      J9F = 0
      NPS = 0
   25 NPS = NPS+1
      IF (IBUG3 .EQ. 1) WRITE (99,302) NPS
*
*   Loop on the disconnected summations
*
      IAS = 0
      NSUM = NBJ(NPS)-NFS
      JWRD = NB6J(NPS)-JWR
      J6CP = K6CP(NPS)
      J7CP = K7CP(NPS)
      J8CP = K8CP(NPS)
      J9CP = K9CP(NPS)
*
*   The range of values of each summation variable is defined by
*   establishing a matrix of the links between variables.
*   MAT(I,J) contains:
*       I = J  Number of possible values of I due to triangular
*              relations with non-variables, i.e. constants.
*       I > J  Number of links between I and J through constants
*       I < J  Value of the constant, if the above is 1. If not,
*              these values are srored in J12(L,I,J) where there
*              is room for MXCSVR such values (L .LE. 4)
*
      DO 52 I = 1,NSUM
         DO 152 J = 1,NSUM
            MAT(I,J) = 0
  152    CONTINUE
   52 CONTINUE
*
      DO 66 I1 = 1,NSUM
         I1T = I1+NFS
         I2 = JSUM6(I1T)
         DO 65 I3 = 1,I2
            I = JSUM5(I1T,I3)
            J = JSUM4(I1T,I3)
            GOTO (54,55,56,57,58,59),J
*
*   The rows of the IPAIR arrays give limits of summation imposed
*
   54       IPAIR(1,1) = JWORD(2,I)
            IPAIR(1,2) = JWORD(5,I)
            IPAIR(2,1) = JWORD(3,I)
            IPAIR(2,2) = JWORD(6,I)
            GOTO 60
*
   55       IPAIR(1,1) = JWORD(1,I)
            IPAIR(1,2) = JWORD(5,I)
            IPAIR(2,1) = JWORD(4,I)
            IPAIR(2,2) = JWORD(6,I)
            GOTO 60
*
   56       IPAIR(1,1) = JWORD(1,I)
            IPAIR(1,2) = JWORD(6,I)
            IPAIR(2,1) = JWORD(4,I)
            IPAIR(2,2) = JWORD(5,I)
            GOTO 60
*
   57       IPAIR(1,1) = JWORD(2,I)
            IPAIR(1,2) = JWORD(6,I)
            IPAIR(2,1) = JWORD(3,I)
            IPAIR(2,2) = JWORD(5,I)
            GOTO 60
*
   58       IPAIR(1,1) = JWORD(1,I)
            IPAIR(1,2) = JWORD(2,I)
            IPAIR(2,1) = JWORD(3,I)
            IPAIR(2,2) = JWORD(4,I)
            GOTO 60
*
   59       IPAIR(1,1) = JWORD(1,I)
            IPAIR(1,2) = JWORD(3,I)
            IPAIR(2,1) = JWORD(2,I)
            IPAIR(2,2) = JWORD(4,I)
*
   60       DO 63 I4 = 1,2
               KM = 0
               DO 62 I5 = 1,2
                  IF (IPAIR(I4,I5) .GT. MP) KM = KM+1
   62          CONTINUE
*
               JJ1 = IPAIR(I4,1)
               JJ2 = IPAIR(I4,2)
               IF (KM .EQ. 1) GOTO 67
               IF (KM .GT. 1) GOTO 63
*
*   One variable linked to two constants. Fix the diagonal MAT(I,I)
*
              JT1 = J1(JJ1)-1
              JT2 = J1(JJ2)-1
              JMIN = ABS (JT1-JT2)
              JMAX = JT1+JT2
*
              IF (MAT(I1,I1) .GT. 1) THEN
*
*   If there are several couples of constants, take the more
*   stringent combination
*
                 JMIN = MAX (JMIN,JSUM(1,I1))
                 JMAX = MIN (JMAX,JSUM(2,I1))
                 IF (JMAX .GE. JMIN) THEN
                    JSUM(1,I1) = JMIN
                    JSUM(2,I1) = JMAX
                    MAT(I1,I1) = (JMAX-JMIN)/2+1
                    GOTO 63
                 ELSE
                    RECUP = 0.0D 00
                    GOTO 110
                 ENDIF
              ELSEIF (MAT(I1,I1) .LT. 1) THEN
*
*   First time
*
                  MAT(I1,I1) = (JMAX-JMIN)/2+1
                  JSUM(1,I1) = JMIN
                  JSUM(2,I1) = JMAX
               ENDIF
*
               GOTO 63
*
*   One variable linked to one constant and one variable  non diagonal
*   element
*
   67          JT1 = MIN (JJ1,JJ2)
               JT2 = MAX (JJ1,JJ2)-MP
               IF (JT2 .GT. I1) GOTO 63
               JT4 = J1(JT1)-1
               K = MAT(I1,JT2)
               IF (K .EQ. 0) GOTO 107
*
               DO 71 LL = 1,K
                  IF (JT4 .EQ. J12(LL,JT2,I1)) GOTO 63
   71          CONTINUE
*
  107          K = K+1
               IF (K .GT. MXCSVR) GOTO 63
               MAT(I1,JT2) = K
               J12(K,JT2,I1) = JT4
*
   63       CONTINUE
   65    CONTINUE
   66 CONTINUE
*
*   Reduce the diagonal elements by taking into account the non
*   diagonal elements, and keep the latter only if needed
*
  150 ICHAN = 0
*
      DO 74 I = 1,NSUM
         NOEL(I) = .TRUE.
         I1 = I-1
         IF (I1 .EQ. 0) GOTO 170
         DO 72  J = 1,I1
            IF ((MAT(I,J) .EQ. 0) .OR. (MAT(J,J) .EQ. 0)) GOTO 72
            IK1 = I
            IK2 = J
            CALL RDIAG (I,J,IK1,IK2,ICHAN,MAT,JSUM,J12)
            NOEL(I) = .FALSE.
   72    CONTINUE
  170    IF (I .EQ. NSUM) GOTO 74
         I2 = I+1
*
         DO 73 J = I2,NSUM
            IF ((MAT(J,I) .EQ. 0) .OR. (MAT(J,J) .EQ. 0)) GOTO 73
            IK1 = J
            IK2 = I
            CALL RDIAG (I,J,IK1,IK2,ICHAN,MAT,JSUM,J12)
   73    CONTINUE
   74 CONTINUE
*
      IF (ICHAN .NE. 0) GOTO 150
      GOTO 220
*
*   Carry out the summations.
*
  220 DO 230 I = 1,NSUM
         JSUM3(I) = 1
         LDIAG(I) = .FALSE.
         IF (MAT(I,I) .EQ. 1) LDIAG(I) = .TRUE.
  230 CONTINUE
*
      DO 231 I = 1,JWRD
         JWTEST(I) = 1
  231 CONTINUE
*
      STOR = 0.0D 00
      STOR1 = 1.0D 00
      NOLP = 0
      IP = 1
  240 NOLP = NOLP+1
*
*   Find the range of JSUM2(NOLP)
*   NOLP is the index  of the summation variable
*
      JMIN = JSUM(1,NOLP)
      JMAX = JSUM(2,NOLP)
      IF (NOEL(NOLP)) GOTO 241
      NO1 = NOLP-1
*
      DO 242 NJ = 1,NO1
         IF (MAT(NOLP,NJ) .EQ. 1) THEN
            JJ1 = MAT(NJ,NOLP)
            JJ2 = JSUM2(NJ)
            JMIN = MAX (JMIN,IABS(JJ2-JJ1))
            JMAX = MIN (JMAX,JJ1+JJ2)
         ELSEIF (MAT(NOLP,NJ) .GT. 1) THEN
            K = MAT(NOLP,NJ)
            JJ2 = JSUM2(NJ)
*
            DO 245 I = 1,K
            JJ1 = J12(I,NJ,NOLP)
            JMIN = MAX (JMIN,IABS(JJ2-JJ1))
            JMAX = MIN (JMAX,JJ1+JJ2)
  245       CONTINUE
*
         ENDIF
*
  242 CONTINUE
*
  241 JSUM2(NOLP) = JMIN
      MAXLP(NOLP) = JMAX
      IF (LDIAG(NOLP)) JSUM3(NOLP) = 0
      IF (NOLP .LT. NSUM) GOTO 240
*
      DO 260 JJ = JMIN,JMAX,2
         JSUM2(NSUM) = JJ
*
*   Determine which RACAH coefficients need re-evaluating and
*   set JWTEST appropriately
*
      DO 114 J = IP,NSUM
         IF (JSUM3(J) .LE. 0) GOTO 114
         I2 = JSUM6(J)
*
         DO 113 I1 = 1,I2
            I3 = JSUM5(J,I1)
            JWTEST(I3) = 1
  113    CONTINUE
  114 CONTINUE
*
      DO 98 J = 1,JWRD
         IF (JWTEST(J) .EQ. 0) GOTO 98
         JWJ = J+JWR
*
         DO 90 I = 1,6
            IF (JWORD(I,JWJ) .LE. MP) THEN
               I1 = JWORD(I,JWJ)
               IST(I) = J1(I1) - 1
            ELSE
               I1 = JWORD(I,JWJ)-MP-NFS
               IST(I) = JSUM2(I1)
            ENDIF
   90    CONTINUE
*
         CALL DRACAH (IST(1),IST(2),IST(3),IST(4),IST(5),IST(6),X1)
         WSTOR(J) = X1
         IF (IBUG3 .EQ. 1) THEN
            DO 99 I = 1,6
               XJ1(I) = 0.5D 00*DBLE (IST(I))
   99       CONTINUE
*
            WRITE (99,305) (XJ1(I), I = 1,6),X1
         ENDIF
   98 CONTINUE
*
*   Form product of Racah coefficients, (2J+1) factors and (-1)
*   factors in STOR1
*
      DO 126 I = 1,JWRD
         STOR1 = STOR1*WSTOR(I)
  126 CONTINUE
*
*   IASTOR contains the power of (-1) which is common to all terms
*
      IX2 = 0
      IJ6CP = 1
      IF (J6CP .NE. J6F) THEN
         JB = J6F+1
*
         DO 128 I = JB,J6CP
            I1 = J6P(I)-NFS
            IJ6CP = IJ6CP*(JSUM2(I1)+1)
  128    CONTINUE
      ENDIF
*
      IF (J9CP .NE. J9F) THEN
         JB = J9F+1
*
         DO 147 I = JB,J9CP
            I1 = J9P(I)-NFS
            IJ6CP = IJ6CP/(JSUM2(I1)+1)
  147    CONTINUE
      ENDIF
*
      STOR1 = STOR1*SQRT (DBLE (IJ6CP))
*
      IF (J7CP .NE. J7F) THEN
         JB = J7F+1
*
         DO 131 I = JB,J7CP
            I1 = J7P(I)-NFS
            IX2 = IX2 + JSUM2(I1)
  131    CONTINUE
      ENDIF
*
      IF (J8CP .NE. J8F) THEN
         JB = J8F+1
*
         DO 134 I = JB,J8CP
            I1 = J8P(I)-NFS
            IX2 = IX2 + 2*(JSUM2(I1))
  134    CONTINUE
      ENDIF
*
      IF (MOD(IX2,2) .EQ. 1) THEN
         IAS = -1
         IX2 = IX2+1
      ENDIF
*
      IX2 = IX2/2
*
*   Add term into STOR and reset STOR1 to 1 ready for next term
*
      IF (MOD(IX2,2) .EQ. 1) STOR1 = -STOR1
      STOR = STOR + STOR1
      STOR1 = 1.0D 00
      NSUM1 = NSUM-1
      IF (NSUM1 .EQ. 0) GOTO 260
*
      DO 261 IK = 1,NSUM1
         JSUM3(IK) = 0
  261 CONTINUE
*
      DO 262 IK = 1,JWRD
         JWTEST(IK) = 0
  262 CONTINUE
*
  260 CONTINUE
*
  250 NOLP = NOLP-1
*
      IF (NOLP .NE. 0) THEN
         IF (LDIAG(NOLP)) GOTO 250
         JSUM3(NOLP) = 1
         JSUM2(NOLP) = JSUM2(NOLP)+2
         IF (JSUM2(NOLP) .GT. MAXLP(NOLP)) GOTO 250
         IP = NOLP
*
*   Proceed to next variable
*
         GOTO 240
*
      ENDIF
*
      RECUP = RECUP*STOR
      IF (IBUG3 .EQ. 1) WRITE (99,307) NPS,STOR,RECUP
      IF (ABS(RECUP) .LT. EPSIL) GOTO 145
      JWR = JWRD+JWR
      NFS = NSUM+NFS
      J6F = J6CP
      J7F = J7CP
      J8F = J8CP
      J9F = J9CP
      IASTOR = IASTOR+IAS
*
*   Proceed to next sum
*
      IF (NPS .LT. NLSUM) GOTO 25
      IASTOR = IASTOR/2
      IF (MOD (IASTOR,2) .NE. 0) RECUP = -RECUP
      IF (IBUG3 .EQ. 1) WRITE (99,304) RECUP
  110 RETURN
*
*   No summations. Check that there are no inconsistencies. Then
*   multiply by (-1) factor and exit
*
  145 RECUP = 0.0D 00
      RETURN
*
  302 FORMAT (' Sum Nr.',I3)
  303 FORMAT (' No summation. Recoupling coefficient = ',G15.8)
  304 FORMAT (' Recoupling coefficient = ',G15.8)
  305 FORMAT (6F5.1,10X,G15.8)
  306 FORMAT (' Number of independent sums:',I3)
  307 FORMAT (' Sum Nr.',I2,' Sum value = ',G15.8,' RECUP = ',G15.8)
  309 FORMAT (' Not involving summation variable')
  400 FORMAT (//' Printout from SUBROUTINE GENSUM'
     :       //' Values of angular momenta in *REAL* FORMAT'
     :        /(14F5.1))
  401 FORMAT (/' Racah W functions (6J)'
     :       /' Arguments in *REAL* FORMAT',18X,'value')
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE INTAB
*                                                                      *
*   This SUBROUTINE called at the end of DIAGRM, fixes the arrays IH   *
*   and IL - so to speak hardware and logical addresses of triads in   *
*   JDIAG . Also  determines the number of free ends NFREE and their   *
*   location ITFREE.                                                   *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INTEGER ARR,TAB1
      LOGICAL FREE
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /BUILD/IAL(M4TRD),IF1,IF2,NODE
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
*
      DO 1 I = 1,M
         IAL(I) = 1
    1 CONTINUE
*
      DO 3 I = IFIRST,ILAST
         J = JDIAG(I,1)
         K = IAL(J)
         TAB1(J,K) = I
         IAL(J) = K+1
    3 CONTINUE
*
      IFR = IFIRST-1
*
      DO 4 I = IFIRST,ILAST
         IT = I-IFR
         IL(I) = IT
         IH(IT) = I
    4 CONTINUE
*
      J = JDIAG(IFIRST,3)
      K = IAL(J)
      IF (K .GT. 1) TAB1(J,2) = TAB1(J,1)
      TAB1(J,1) = IFIRST
      IAL(J) = 3
      J = JDIAG(ILAST,2)
      TAB1(J,2) = ILAST
      IAL(J) = 3
      NFREE = 0
*
      DO 7 I = IFIRST,ILAST
         J = JDIAG(I,1)
         IF (IAL(J) .NE. 3) THEN
            NFREE = NFREE+1
            ITT = ILAST+NFREE
            TAB1(J,2) = ITT
            IL(ITT) = NFREE*1000
            ITFREE(NFREE) = I
         ENDIF
    7 CONTINUE
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE LOLPOP (FAIL)
*                                                                      *
*   Reduces a loop with one line and one node in the flat graph.       *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME,NAMSUB
      INTEGER ARR,TAB1
      LOGICAL FAIL,SUMVAR,FREE
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      DIMENSION KP(3),KS(3)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /NAM/NAMSUB
*
      DATA NAME/'LOLPOP'/
      DATA KP/2,3,1/
      DATA KS/0,1,-1/
*
      NAMSUB = NAME
      I1 = NPOINT(1)
      K3 = 2
      IF (I1 .EQ. ILAST) K3 = 3
      L = JDIAG(I1,K3)
      CALL DELTA (L,MP,FAIL)
      IF (FAIL) RETURN
      K = KP(K3)
      IF (ARR(I1,K) .LT. 0) CALL PHASE2 (JDIAG(I1,K))
      K1 = KS(K3)
      IL1 = IL(I1)+K1
      I2 = IH(IL1)
      L1 = JDIAG(I2,1)
      CALL DELTA (L1,JDIAG(I2,K3),FAIL)
      IF (FAIL) RETURN
      IF (ARR(I2,K3) .EQ. K1) CALL PHASE2 (L1)
      IL2 = IL(I2)+K1
      I3 = IH(IL2)
      K2 = K3+K1
      JDIAG(I3,K2) = L1
      ARR(I3,K2) = ARR(I2,1)
      J9C = J9C+1
      J9(J9C) = L1
      J6C = J6C+1
      J6(J6C) = JDIAG(I1,1)
      IF (K3 .EQ. 3) RETURN
*
      DO 1 I = 3,NBNODE
         IT = IH(I)
         ILP = I-2
         IL(IT) = ILP
         IH(ILP) = IT
    1 CONTINUE
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE NEIBOR (LC,L1,L2)
*                                                                      *
*   Gives the positions of the other two arguments in the triad.       *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      IF (LC .LT. 2) THEN
         L1 = 2
         L2 = 3
      ELSEIF (LC .EQ. 2) THEN
         L1 = 3
         L2 = 1
      ELSE
         L1 = 1
         L2 = 2
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE NJGRAF (RECUP,IGEN,FAIL)
*                                                                      *
*   Program to calculate a general recoupling coefficient. This ver-   *
*   sion is slightly modified  (Colin T Johnson, Oxford University).   *
*   The changes are as follows:                                        *
*                                                                      *
*      1. PARAMETER IGEN has been included in the argument list:       *
*            IGEN =  0  normal call to NJGRAF                          *
*            IGEN = -1  GENSUM is not called                           *
*      2. The contents of COMMON blocks /ARGU/ and /SUMARG/ are used   *
*         by  GENSUM  to calculate the recoupling coefficient. These   *
*         COMMON  blocks  have  been removed from GENSUM. Their con-   *
*         tents are passed to  GENSUM  through the argument list in-   *
*         stead, so that NJGRAF can be called to set up formulae for   *
*         both the direct and exchange cases in COR and BREIT.         *
*      3. Extra  dimension  tests  have  been  included  in routines   *
*         NJGRAF, PRINTJ, SPRATE, VAR and ZERO. These are  discussed   *
*         below.                                                       *
*      4. An  extra routine  RDIAG  has been introduced to remove an   *
*         extended  DO loop from GENSUM, to conform with the FORTRAN   *
*         77 standard.                                                 *
*                                                                      *
*                                                                      *
*   Description of some COMMON blocks. A full discussion is given in   *
*   the NJGRAF program description (Bar-Shalom and Klapisch op cit).   *
*                                                                      *
*      COMMON block COUPLE                                             *
*                                                                      *
*         M                The total number of angular momentum val-   *
*                          ues in the initial and final states         *
*         N                The number of basic angular momentum val-   *
*                          ues that are coupled                        *
*         J1(I),           The angular momentum values stored as 2J+1  *
*            I = 1,M                                                   *
*         J2(I,J),         The position in the J1 array of the init-   *
*            I = 1,(N-1),  ial state triads                            *
*            J = 1,3                                                   *
*         J3(I,J),         The position in the J1 array of the final   *
*            I = 1,(N-1),  state triads                                *
*            J = 1,3                                                   *
*         FREE(I),         If FREE(I) = .TRUE., no reference is made   *
*            I = 1,M       to the value of J1(I) when establishing a   *
*                          formula in  NJGRAF .  GENSUM  may then be   *
*                          called  for  repeated  occurences of this   *
*                          formula  with  differing values of J1(I).   *
*                          If J1(I) does  not  vary between calls to   *
*                          GENSUM then FREE(I) should be set .FALSE.   *
*                          so that zero branches  can be removed be-   *
*                          fore the formula is established.            *
*                                                                      *
*      COMMON block DEBUG                                              *
*                                                                      *
*         IBUG1            Not used                                    *
*         IBUG2            Not used                                    *
*         IBUG3            Debug prints in NJGRAF and GENSUM if 1      *
*         IBUG4            Not used                                    *
*         IBUG5            Not used                                    *
*         IBUG6            Not used                                    *
*                                                                      *
*      COMMON block ARGU                                               *
*                                                                      *
*         J6C              The number of elements in the K6 array      *
*         J7C              The number of elements in the K7 array      *
*         J8C              The number of elements in the K8 array      *
*         J9C              The number of elements in the K9 array      *
*         JWC              The number of columns in the KW array       *
*         J6(I),           Each entry corresponds to a  factor  SQRT   *
*            I = 1,J6C     (2J+1) in RECUP. The value  of  J6  GIVES   *
*                          position in  J1  array where  J  value is   *
*                          found                                       *
*         J7(I),           Each entry corresponds to factor  (-1)**J   *
*            I = 1,J7C     in RECUP                                    *
*         J8(I),           Each entry corresponds to a factor (-1)**   *
*            I = 1,J8C     (2J) in RECUP                               *
*         J9(I),           Each entry corresponds to a factor (2J+1)   *
*            I = 1,J9C     **(-0.5) in RECUP                           *
*         KW(I,J),         Each column corresponds to a Racah coeff-   *
*            I = 1,6,      icient in RECUP                             *
*            J = 1,JWC                                                 *
*         JDEL             The number of delta functions               *
*         LDEL(I,J),       The arguments of the delta functions        *
*              J = 1,2                                                 *
*         SUMVAR(I)        .TRUE. for ang. mom. I (a summation vari-   *
*                          able                                        *
*         MP               The index of the last variable              *
*                                                                      *
*   The arrays  J6, J7, J8, J9 and  KW, Are evaluated by NJGRAF. The   *
*   summation over the variables in  J6, J7, J8, J9 and  KW, and the   *
*   evaluation of RECUP is carried out in GENSUM. GENSUM  can be re-   *
*   entered directly to evaluate different  recoupling  coefficients   *
*   with the same structure  by just  altering the numbers in the J1   *
*   array.                                                             *
*                                                                      *
*   This is the main program. It handles all the analysis of the re-   *
*   coupling  coefficient without referring explicitly to the values   *
*   of angular  momenta  which  are in J1(J),except for zero in case   *
*   FREE = .FALSE. . Like NJSYM it  prepares arrays of arguments for   *
*   phase factors, (2*J+1) factors and  6j-coefficients to be compu-   *
*   ted in GENSUM, which can also be called separately when only the   *
*   numerical values of angular momenta change. These variable angu-   *
*   lar momenta should be declared  FREE(J)  = .TRUE. , so  that the   *
*   formula prepared for GENSUM should be correct when J1 is not ze-   *
*   ro. FAIL will be TRUE when the  recoupling  coefficient  is zero   *
*   because of unsatisfied delta or other similar causes.              *
*                                                                      *
*   This version holds the array dimensions in parameter statements.   *
*   The dimensions are labelled:                                       *
*                                                                      *
*      MANGM  : Dimension of the J1 and FREE arrays in /COUPLE/, and   *
*               the  first  dimension of the LINE and LCOL arrays in   *
*               /TREE/. Also  the  dimension  of the SUMVAR array in   *
*               /ARGU/, AND OF THE INVER array in routine SPRATE. It   *
*               is tested for  M  on entry to  NJGRAF, and for MP in   *
*               routine SPRATE.                                        *
*      MTRIAD : Dimension of the  J2 and  J3 arrays in /COUPLE/. The   *
*               dimensions of these  arrays  are checked on entry to   *
*               NJGRAF in addition  MTRIAD sets the dimension of the   *
*               JSUM6 array and the first dimension of the JSUM4 and   *
*               JSUM5  arrays in /SUMARG/. Also gives the dimensions   *
*               of some  temporary working arrays in SPRATE and GEN-   *
*               SUM. In these  cases  mtriad sets the maximum number   *
*               of summation variables  in any particular sum, which   *
*               is tested in SPRATE.                                   *
*      M2TRD  : (=2*MTRIAD) Dimension of the J23 ,  ARROW  and  TABS   *
*               arrays in /TREE/. Also  the  dimension of the npoint   *
*               array in /GRAPH/.                                      *
*      M4TRD  : (=4*MTRIAD) Dimension of the  JDIAG,  ARR, IL and IH   *
*               arrays in /GRAPH/, and of the IAL array in /BUILD/.    *
*      M3MNGM : Dimension of the J6 array in /ARGU/, tested in SPRATE  *
*               Dimension of the J7 array in /ARGU/, tested in SPRATE  *
*               Dimension of the J8 array in /ARGU/, tested in SPRATE  *
*      MANGMP : Dimension of the J9 array in /ARGU/, tested in SPRATE  *
*               MANGMP also sets the dimension of the J6P,  J7P, J8P   *
*               and J9P arrays in /SUMARG/, And of the JNS  array in   *
*               routine VAR. The dimension of the JNS array is  tes-   *
*               ted in VAR.                                            *
*      M6J    : Dimension of the JW(or KW) and LDEL arrays in /ARGU/,  *
*               and of the JWORD and INV6J arrays in /SUMARG/.  Also   *
*               the second dimension of the  JSUM4 and  JSUM5 arrays   *
*               in /SUMARG/. In addition it  gives the dimensions of   *
*               a  number of  temporary  working  arrays in routines   *
*               SPRATE and GENSUM. M6J is tested in SPRATE.            *
*      MFACT  : The dimension of the factorial array GAM in /FACTS /.  *
*      MSUM   : Dimension of the NBJ, NB6J, K6CP, K7CP, K8CP and K9CP  *
*               arrays in /SUMARG/. MSUM is the  maximum  number  of   *
*               sums allowed, and is tested in routine SPRATE.         *
*      MTAB   : The dimension of the JTAB array in  routine  PRINTJ.   *
*               MTAB is tested in PRINTJ.                              *
*      MZERO  : Dimension of the JZERO array in /ZER/. MZERO is tes-   *
*               ted in routine ZERO.                                   *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME,NAMSUB
      LOGICAL FAIL,FIND,TABS,CUT,FREE,SUMVAR
      INTEGER ARROW,ARR,TAB1
*
      PARAMETER (MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :           MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :           MFACT = 500, M6J = 20, MSUM = 10)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /CONS/ZRO,HALF,TENTH,ONE,TWO,THREE,TEN
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
     :      /CUTDIG/CUT
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /FACTS /GAM(MFACT)
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /NAM/NAMSUB
     :      /SUMARG/J6P(MANGMP),J7P(MANGMP),J8P(MANGMP),J9P(MANGMP),
     :              JWORD(6,M6J),NLSUM,NBJ(MSUM),NB6J(MSUM),K6CP(MSUM),
     :              K7CP(MSUM),K8CP(MSUM),K9CP(MSUM),JSUM6(MTRIAD),
     :              JSUM4(MTRIAD,M6J),JSUM5(MTRIAD,M6J),INV6J(M6J)
     :      /TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
*
      DATA NAME/'NJGRAF'/
*
*   Debug printout
*
      IF (IBUG3 .EQ. 1) THEN
         WRITE (99,300)
         WRITE (99,301) M,N-1
         WRITE (99,302)
         WRITE (99,303) (J1(I),I = 1,M)
         WRITE (99,304) (FREE(I),I = 1,M)
         WRITE (99,305)
         WRITE (99,306) ((J2(I,J),J = 1,3),
     :                    (J3(I,J),J = 1,3),I = 1,N-1)
      ENDIF
*
*   Test the dimension of the J1 array
*
      IF (M+1 .GT. MANGM) THEN
         WRITE (*,307)
         WRITE (*,308) M+1,MANGM
         STOP
      ENDIF
*
*   Test the dimensions of the J2 and J3 arrays
*
      IF (N-1 .GT. MTRIAD) THEN
         WRITE (*,307)
         WRITE (*,309) N-1,MTRIAD
         STOP
      ENDIF
*
*   Initializations
*
      DO 502 I = N,MTRIAD
         DO 501 J = 1,3
            J2(I,J) = 0
            J3(I,J) = 0
  501    CONTINUE
  502 CONTINUE
*
      FAIL = .FALSE.
      J6C = 0
      J7C = 0
      J8C = 0
      J9C = 0
      JWC = 0
      JDEL = 0
      CALL SETDM
      NFIN = 0
      CUT = .FALSE.
*
*   Building up of the unstructured graph
*
      CALL SETTAB (FAIL)
*
*   Exit with RECUP set to zero if any angular momentum is
*   impossible
*
      M = M+1
      IF (FAIL) GOTO 7
*
      M = M-1
*
*   Locate and eliminate any zero angular momentum; simplify the
*   graph
*
      JF = 0
      JF1 = 0
      CALL ZERO (JF1,JF,FAIL)
      IF (FAIL) GOTO 7
*
      MP = M
      IF (NBTR .EQ. 0) GOTO 6
      JUMP = 1
*
    1 CALL CHKLP1 (FAIL)
      IF (FAIL) GOTO 7
*
*   Build a flat diagram out of the unstructured graph; several flat
*   diagrams may constitute the original graph, in which case there
*   are possible cuts; the flat diagrams will have free ends if cut
*
      CALL DIAGRM (JUMP)
      NFIN = MAX (0,NFREE-2)
*
      IF (NFIN .NE. 0) THEN
         JUMP = 3
*
*   Handling of free ends if a cut was found
*
         CALL CUTNL (FAIL)
         IF (FAIL) GOTO 7
      ELSE
         JUMP = 2
         IF (NFREE .EQ. 1) THEN
            CALL CUT1L (FAIL)
            IF (FAIL) GOTO 7
         ELSEIF (NFREE .GT. 1) THEN
            CALL CUT2L (FAIL)
            IF (FAIL) GOTO 7
         ENDIF
      ENDIF
*
      NBTR = NBTR+NFIN
      IF (NBTR .NE. 0) CUT = .TRUE.
*
*   Analysis of the flat diagram.
*   Closed circuits of increasing order NC are searched, analysed,
*   and taken out of the flat diagram, thus reducing the number of
*   nodes, NBNODE.
*
      NC = 0
   10 NC = NC+1
      CALL SEARCH (FIND)
      IF (.NOT. FIND) GOTO 10
      NCP = NC-2
      JPOL = 0
      IF ((M .EQ. MP) .AND. (NC.GT.3)) CALL SETDM
      IF (IPARTL .GT. 2) CALL POLYGN (JPOL)
      GOTO (11,12,13,14),NC
   11 CALL LOLPOP (FAIL)
      IF (FAIL) GOTO 7
      GOTO 15
   12 CALL BUBBLE (JPOL,FAIL)
      IF (FAIL) GOTO 7
      GOTO 15
   13 CALL TRIANG (FAIL)
      IF (FAIL) GOTO 7
      GOTO 15
   14 CALL SQUARE
   15 NBNODE = NBNODE-2
      IF (NBNODE .EQ. 0) GOTO 9
      IFIRST = IH(1)
      ILAST = IH(NBNODE)
*
*   PRINTJ is an all purpose printing SUBROUTINE called from many
*   places
*
      CALL PRINTJ (NAMSUB,8)
      IF (NBNODE .EQ. NFIN) GOTO 9
      NC = NCP
*
*   Proceed to other circuits of order NC-1
*
      GOTO 10
    9 IF (NBTR .EQ. 0) GOTO 6
      IF (JUMP .EQ. 3) CALL ORDTRI
*
*   At this stage, the flat diagram has been reduced to nodes
*   involving free ends. Proceed to build other flat diagrams
*   if necessary.
*
      GOTO 1
*
*   All parts of the original graph have been reduced.
*
    7 RECUP = ZRO
      M = M-1
      RETURN
    6 CALL PRINTJ (NAME,0)
*
*   Preparation of the results, and separation in several sums
*   if cuts have been detected, also in the flat diagram itself
*
      CALL SPRATE (M)
      M = M-1
*
*   GENSUM computes the numerical value of the recoupling
*   coefficient
*
      IF (IGEN .NE. -1) CALL GENSUM (J6C,J7C,J8C,J9C,JWC,J6,J7,J8,J9,KW,
     :                               JDEL,LDEL,SUMVAR,MP,J6P,J7P,J8P,
     :                               J9P,JWORD,NLSUM,NBJ,NB6J,K6CP,K7CP,
     :                               K8CP,K9CP,JSUM4,JSUM5,JSUM6,INV6J,
     :                               RECUP)
*
      RETURN
*
  300 FORMAT (//' ++++++++++ NJGRAF ++++++++++'/)
  301 FORMAT (' Total number of angular momenta (M) = ',1I3
     :      //' Number of triads in each of the ''left-hand'' and',
     :        ' ''right-hand'' states (N-1) = ',1I3)
  302 FORMAT (/' (2J+1)-value for each angular momentum:')
  303 FORMAT (1X,42I3)
  304 FORMAT (1X,42L3)
  305 FORMAT (/' ''Left-hand'' triads',10X,'''Right-hand'' triads')
  306 FORMAT (1X,3I3,19X,3I3)
  307 FORMAT (/' ***** Error in NJGRAF *****'/)
  308 FORMAT (' M+1 = ',1I3,', exceeds PARAMETER MANGM = ',1I3)
  309 FORMAT (' N-1 = ',1I3,', exceeds PARAMETER MTRIAD = ',1I3)
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE ORDTRI
*                                                                      *
*   This subroutine orders the triads which were left with free ends   *
*   as consequence of cutting,so that the new graph will start there.  *
*                                                                      *
*                                           Last update: 16 Aug 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARROW,ARR,TAB1
      LOGICAL SUMVAR,TABS
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /BUILD/IAL(M4TRD),IF1,IF2,NODE
     :      /KEEP/JKP(2,3),JARR(2,3),IT2,IT3,IT5
*
      DATA NAME/'ORDTRI'/
*
      DO 10 I = 1,MP
         IAL(I) = 0
   10 CONTINUE
*
      IF (NFIN .NE. 0) THEN
         NBT1 = NBTR-1
         NBT = NBT1+NFIN
         NBTT = NBT+1
         NB = 0
         GOTO 31
      ENDIF
*
      NF = NBTR-ITFREE(1)
*
      IF (IT5 .EQ. 0) THEN
         NBT1 = NBTR-1
         N0 = 0
         NFT = NFREE
         ISW = 2
         GOTO 100
      ENDIF
*
      NFT = IT5-IT2
      NM = NFT+NBTR+1
      NBT1 = NBTR
*
      DO 21 J = 1,3
         JDIAG(NBTR,J) = JKP(1,J)
         ARR(NBTR,J) = JARR(1,J)
   21 CONTINUE
*
      JT = JDIAG(NM,1)
      N0 = 0
      ISW = 1
      GOTO 100
*
   22 N0 = NFT
*
      DO 211 J = 1,3
         JDIAG(NM,J) = JKP(2,J)
         ARR(NM,J) = JARR(2,J)
  211 CONTINUE
*
      NBT1 = NBT1+1
      NFT = IT3-IT5
      ISW = 3
      GOTO 100
*
   24 NBT1 = K-NFT
*
   23 NODE = NBT1+NFT
      CALL CHANGE (NODE,2)
      GOTO 40
*
   31 DO 35 I = 1,NBNODE
         I1 = IH(I)
         IF (IL(I1) .GT. ILAST) GOTO 35
         I2 = NBT1+I
         IF (I1 .GT. NBTT) GOTO 33
         IF (I1 .EQ. I2) GOTO 32
         IF (IL(I2) .LE. NBNODE) GOTO 35
*
   33    DO 34 J = 1,3
            JDIAG(I2,J) = JDIAG(I1,J)
            ARR(I2,J) = ARR(I1,J)
   34    CONTINUE
*
         IL(I1) = ILAST+I
   32    NB = NB+1
         IL(I2) = 0
*
   35 CONTINUE
*
      IF (NB .NE. NFIN) GOTO 31
      NODE = NBT
   40 IF1 = JDIAG(NBTR,1)
      IF2 = JDIAG(NBTR,3)
*
      DO 51 I = NBTR,NODE
         DO 50 K = 1,3
            J = JDIAG(I,K)
            IAL(J) = IAL(J)+1
   50    CONTINUE
   51 CONTINUE
*
      ILAST = NODE
      CALL PRINTJ (NAME,8)
*
      RETURN
*
  100 IF (NF .LE. 0) THEN
         NFR = N0
         I1 = 1
      ELSE
         NFR = NFT+1
         I1 = -1
      ENDIF
*
      DO 4 I = 1,NFT
         IK = NFR+I1*I
         IT = ITFREE(IK)
         K = NBT1+IK
*
         DO 3 J = 1,3
            JDIAG(K,J) = JDIAG(IT,J)
            ARR(K,J) = ARR(IT,J)
    3    CONTINUE
*
    4 CONTINUE
*
      GOTO (22,23,24),ISW
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE OTHERJ (LIN,J,LO,LCO,K)
*                                                                      *
*   Gives the other triad where a given J occurs and its position.     *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INTEGER ARROW
      LOGICAL TABS
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
*
      LO = LINE(J,1)
      IF ((LO .EQ. LIN) .OR. (TABS(LO))) THEN
         K = 1
         LO = LINE(J,2)
         LCO = LCOL(J,2)
      ELSE
         K = 2
         LCO = LCOL(J,1)
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE PHASE (L,JM,NDIM)
*                                                                      *
*   Phase factor arising from non-cyclic permutation of arguments in   *
*   triad L. JM may be either J23 or JDIAG.                            *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      LOGICAL SUMVAR
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      DIMENSION JM(NDIM,3)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      J7(J7C+1) = JM(L,1)
      J7(J7C+2) = JM(L,2)
      J7C = J7C+3
      J7(J7C) = JM(L,3)
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE PHASE2 (J)
*                                                                      *
*   Adds a phase factor (-1)**2J .                                     *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      LOGICAL SUMVAR
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      J8C = J8C+1
      J8(J8C) = J
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE POLYGN (JPOL)
*                                                                      *
*   This routine reduces a circuit of arbitrary order NC. It exchan-   *
*   ges nodes on the flat diagram until the distance on the axis be-   *
*   tween nodes equals one. Each exchange introduces a summation va-   *
*   riable  and  a 6j-symbol. The circuit has a maximum of NPART = 2   *
*   disconnected parts on the axis.                                    *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARR,TAB1
      LOGICAL SUMVAR
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      DATA NAME/'POLYGN'/
*
      NC1 = NC+1
      NC2 = NC
      NBC = IPARTL-2
*
   10 DO 8 I = 1,NBC
         IT2 = NPOINT(NC1-I)
         IT1 = NPOINT(NC2-I)
         JB = JDIAG(IT1,1)
         JC = JDIAG(IT2,1)
         JDIAG(IT1,1) = JC
         JDIAG(IT2,1) = JB
         JAR = ARR(IT1,1)
         ARR(IT1,1) = ARR(IT2,1)
         ARR(IT2,1) = JAR
         JE = JDIAG(IT1,2)
         MP = MP+1
         SUMVAR(MP) = .TRUE.
         JDIAG(IT1,2) = MP
         JDIAG(IT2,3) = MP
*
         IF (TAB1(JB,1) .EQ. IT1) THEN
            TAB1(JB,1) = IT2
         ELSE
            TAB1(JB,2) = IT2
         ENDIF
*
         IF (TAB1(JC,1) .EQ. IT2) THEN
            TAB1(JC,1) = IT1
         ELSE
            TAB1(JC,2) = IT1
         ENDIF
*
         IF (ARR(IT1,2) .LE. 0) THEN
            CALL PHASE2 (JE)
            ARR(IT1,2) = 1
            ARR(IT2,3) = -1
         ENDIF
*
         JWC = JWC+1
         KW(1,JWC) = JB
         KW(2,JWC) = MP
         KW(3,JWC) = JE
         KW(4,JWC) = JC
         KW(5,JWC) = JDIAG(IT2,2)
         KW(6,JWC) = JDIAG(IT1,3)
         J6(J6C+1) = MP
         J6C = J6C+2
         J6(J6C) = MP
    8 CONTINUE
*
      NC = NC-NBC
*
      IF (NC .GT. 4) THEN
         NBC = IPARTS-2
         NC1 = IPARTS+1
         NC2 = IPARTS
         GOTO 10
      ENDIF
*
      IF (NPART .NE. 1) THEN
         NPOINT(3) = NPOINT(NC1)
         NPOINT(4) = NPOINT(NC1+1)
      ENDIF
*
      IF (NC .EQ. 2) JPOL = 1
      CALL PRINTJ (NAME,10)
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE PRINTJ (NAMES,JP)
*                                                                      *
*   This  SUBROUTINE  prints  intermediate  results in standard form   *
*   from wherever it is called.                                        *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*8 IBLANK,IFREE,IFR
      CHARACTER*6 NAMES,NSETTB
      CHARACTER*4 I6,I7,I8,I9,IJ1
      CHARACTER IM,IP,IS(3)
      INTEGER ARR,TAB1,ARROW
      LOGICAL TABS,SUMVAR,FREE
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10,MTAB = 30,MZERO = 20)
*
      DIMENSION IX(7),JTAB(MTAB,3)
*
      COMMON/TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /CONST/I6C,I7C,I8C,I9C,IDEL,IWC
     :      /ZER/NZERO,JZERO(MZERO)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
*
      EQUIVALENCE(I6C,IX(1))
*
      DATA IBLANK,IFREE,IP,IM/'        ','FREE END','+','-'/
      DATA NSETTB/'SETTAB'/
      DATA I6,I7,I8,I9,IJ1/'I6=','I7=','I8=','I9=','J1='/
*
      IF (IBUG3 .NE. 1) RETURN
      WRITE (99,1050) NAMES
*
*   Initialise variables
*
      I6C = 1
      I7C = 1
      I8C = 1
      I9C = 1
      IDEL = 1
      IWC = 1
*
      JUMP = JP
      IF (JUMP .EQ. 0) THEN
*
         DO 9 I = 1,7
            IX(I) = 1
    9    CONTINUE
*
         WRITE (99,1020) IJ1,(J1(I),I = 1,M)
      ENDIF
*
      IF (JUMP .LT. 8) GOTO 20
      WRITE (99,1000) NBNODE,NBTR,NFIN,IFIRST,ILAST,NFREE
      JUMP = JUMP-8
      WRITE (99,1001)
      K = 0
*
      DO 1 I = 1,NBNODE
         IT = IH(I)
         IFR = IBLANK
         JT = JDIAG(IT,1)
*
         IF ((TAB1(JT,2) .NE. IT) .OR. (JT .EQ. JDIAG(IFIRST,3))) THEN
            K = K+1
            IF (K .GT. MTAB) THEN
               WRITE (*,100) K,MTAB
               STOP
            ENDIF
            JTAB(K,1) = JT
            JTAB(K,2) = TAB1(JT,1)
            JTAB(K,3) = TAB1(JT,2)
         ENDIF
*
         IF (TAB1(JT,2) .GT. ILAST) IFR = IFREE
*
         DO 2 J = 1,3
            IS(J) = IP
            IF (ARR(IT,J) .LT. 1) IS(J) = IM
    2    CONTINUE
*
         WRITE (99,1002) (IS(J),J = 1,3)
         WRITE (99,1003) IL(IT),IT,IFR,(JDIAG(IT,J),J = 1,3)
*
    1 CONTINUE
*
      WRITE (99,1004)
      NTIME = 0
      JT = JDIAG(IFIRST,3)
      IF (JT .NE. JDIAG(ILAST,2)) THEN
         IF (TAB1(JT,2) .LT. 1000) GOTO 5
      ENDIF
    4 K = K+1
      IF (K .GT. MTAB) THEN
         WRITE (*,101) K,MTAB
         STOP
      ENDIF
      JTAB(K,1) = JT
      JTAB(K,2) = TAB1(JT,1)
      JTAB(K,3) = TAB1(JT,2)
    5 NTIME = NTIME+1
*
      IF (NTIME .NE. 2) THEN
         JT = JDIAG(ILAST,2)
         IF (TAB1(JT,2) .EQ. 1000) GOTO 4
      ENDIF
*
      WRITE (99,1005) ((JTAB(I,J),J = 1,3),I = 1,K)
      WRITE (99,1006) (I,SUMVAR(I),I = 1,MP)
   20 IF (JUMP .LT. 4) GOTO 30
      JUMP = JUMP-4
      NBTR1 = 2*N-2
      WRITE (99,1010) NBTR1
      K = 0
*
      DO 11 I = 1,NBTR1
         IF (TABS(I)) GOTO 11
         K = K+1
*
         DO 12 J = 1,3
            IS(J) = IP
            IF (ARROW(I,J) .LT. 1) IS(J) = IM
   12    CONTINUE
*
         WRITE (99,1012) (IS(J),J = 1,3)
         WRITE (99,1013) K,I,(J23(I,J),J = 1,3)
*
   11 CONTINUE
*
      WRITE (99,1014)
      MM = M
      IF (NAMES .NE. NSETTB) MM = M-1
      WRITE (99,1015) (I,(LINE(I,J),LCOL(I,J),J = 1,2),I = 1,MM)
*
   30 IF (JUMP .GE. 2) THEN
         JUMP = JUMP-2
         WRITE (99,1030) NC,NPART,IPARTL,IPARTS,ICROSS,
     :                    (NPOINT(I),I = 1,NC)
      ENDIF
*
      IF (JUMP .GE. 1) WRITE (99,1040) NZERO,(I,JZERO(I),I = 1,NZERO)
      IF (J6C .GE. I6C) WRITE (99,1020) I6,(J6(I),I = I6C,J6C)
      IF (J7C .GE. I7C) WRITE (99,1020) I7,(J7(I),I = I7C,J7C)
      IF (J8C .GE. I8C) WRITE (99,1020) I8,(J8(I),I = I8C,J8C)
      IF (J9C .GE. I9C) WRITE (99,1020) I9,(J9(I),I = I9C,J9C)
      IF (JDEL .GE. IDEL) WRITE (99,1021)
     :                    ((LDEL(I,J),J = 1,2),I = IDEL,JDEL)
      IF (JWC .GE. IWC) WRITE (99,1022)
     :                    ((KW(J,I),J = 1,6),I = IWC,JWC)
      I6C = J6C+1
      I7C = J7C+1
      I8C = J8C+1
      I9C = J9C+1
      IDEL = JDEL+1
      IWC = JWC+1
      RETURN
*
  100 FORMAT (' Dimension error in PRINTJ. K = ',I5,' MTAB = ',I5)
  101 FORMAT (' Dimension error IN PRINTJ. K = ',I5,' MTAB = ',I5)
 1000 FORMAT (/10X,'NBNODE = ',I3,10X,'NBTR = ',I3,10X,'NFIN = ',I3,
     :   /10X,'IFIRST = ',I3,10X,'ILAST = ',I3,9X,'NFREE = ',I3)
 1001 FORMAT (//7X,'IL',3X,'IH',14X,'JDIAG'//)
 1002 FORMAT (28X,3(A1,2X))
 1003 FORMAT (7X,I2,3X,I2,2X,A8,2X,3I3/)
 1004 FORMAT (/5X,'TAB1'/)
 1005 FORMAT (4(I3,1H),2X,I3,I5,5X))
 1006 FORMAT (/2X,'SUMVAR = ',15(I3,L1))
 1010 FORMAT (//10X,'J23',10X,'NBTR1 = ',I3//)
 1012 FORMAT (18X,3(A1,2X))
 1013 FORMAT (I9,I5,2X,3I3/)
 1014 FORMAT (/3X,'J  L1 K1  L2 K2')
 1015 FORMAT (4(I4,1H),I3,I3,I4,I3))
 1020 FORMAT (/3X,A4,3X,3(20I3/))
 1021 FORMAT (/3X,'DELTA = ',7(I5,I3))
 1022 FORMAT (/3X,'KW(ARG. OF 6J)',6I3)
 1030 FORMAT (//2X,'NC = ',I2,4X,'NPART = ',I2,4X,'IPARTL = ',I2,4X,
     :   'IPARTS = ',I2,4X,'ICROSS = ',I2,4X,/2X,'NPOINT = ',20I3)
 1040 FORMAT (//2X,'NZERO = ',I2,5X,12(I4,1H),I3))
 1050 FORMAT (///3X,'Printout after calling SUBROUTINE ',A7)
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE RDIAG (I,J,IK1,IK2,ICHAN,MAT,JSUM,J12)
*                                                                      *
*   Called by  GENSUM to establish the range of values of the summa-   *
*   tion variables.  This routine replaces an extended range do loop   *
*   in GENSUM, to conform with the FORTRAN 77 standard.                *
*                                                                      *
*                                           Last update: 02 Sep 1987   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      PARAMETER (MTRIAD = 12,M6J = 20)
*
      DIMENSION MAT(MTRIAD,MTRIAD),JMNP(5),JMXP(5)
      DIMENSION J12(4,MTRIAD,MTRIAD),JSUM(2,M6J)
*
      JMIN1 = 0
      JMAX1 = 1000
      K = MAT(IK1,IK2)
*
      DO 1 L1 = 1,K
*
         L3 = MAT(J,J)
         JJ1 = JSUM(1,J)
         JND = J12(L1,IK2,IK1)
         JMIN = 1000
         JMAX = 0
         JMNP(L1) = 0
         JMXP(L1) = 1000
*
         DO 2 L2 = 1,L3
*
            JMN = IABS(JND-JJ1)
            JMX = JND+JJ1
            JMIN = MIN (JMN,JMIN)
            JMAX = MAX (JMX,JMAX)
            JMNP(L1) = MAX (JMN,JMNP(L1))
            JMXP(L1) = MIN (JMX,JMXP(L1))
            JJ1 = JJ1+2
*
    2    CONTINUE
*
         JMIN1 = MAX (JMIN1,JMIN)
         JMAX1 = MIN (JMAX1,JMAX)
*
    1 CONTINUE
*
      IF (MAT(I,I) .EQ. 0) THEN
         JSUM(1,I) = JMIN1
         JSUM(2,I) = JMAX1
         MAT(I,I) = (JMAX1-JMIN1)/2+1
         ICHAN = ICHAN+1
         GOTO 3
      ENDIF
*
      IF (JSUM(1,I) .LT. JMIN1) THEN
         JSUM(1,I) = JMIN1
         ICHAN = ICHAN+1
      ENDIF
*
      IF (JSUM(2,I) .GT. JMAX1) THEN
         JSUM(2,I) = JMAX1
         ICHAN = ICHAN+1
      ENDIF
*
    3 K1 = 0
*
      DO 4 L1 = 1,K
         IF ((JMNP(L1) .LE. JSUM(1,I)) .AND.
     :       (JMXP(L1) .GE. JSUM(2,I))) GOTO 4
         K1 = K1+1
         J12(K1,IK2,IK1) = J12(L1,IK2,IK1)
    4 CONTINUE
*
      IF (K1 .NE. K) THEN
         MAT(IK1,IK2) = K1
         ICHAN = ICHAN+1
      ENDIF
*
      MAT(IK2,IK1) = J12(1,IK2,IK1)
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE SEARCH (FIND)
*                                                                      *
*   This  routine locates circuits or loops of order  NC. NPOINT(NC)   *
*   are the  indices  of the points (triads) pertaining to the first   *
*   such loop found.  NPART  is the number of separate parts (groups   *
*   of contiguous points) on  the  axis of the flat graph. IPARTS is   *
*   the number of points in the smallest  part. IPARTL is the number   *
*   of points in  the  largest  part.  The  SUBROUTINE finds all the   *
*   possible loops  of  order  3  and 4. For NC .GE. 5, it looks for   *
*   only those who are partitionned in NPART .LE. 2. which can even-   *
*   tually reduce to a loop of  order  4  without breaking the basic   *
*   structure of the flat graph. ICROSS = -1, if lines cross.          *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARR,TAB1
      LOGICAL FIND
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
*
      DATA NAME/'SEARCH'/
*
*   Initialization
*
      FIND = .FALSE.
      NCM1 = NC-1
      NCM = NC-2
      ICROSS = 0
*
*   First treat two cases that do not involve do loops:
*
*   1. One isolated point, either the first or the last
*
      NPART = 1
      IPARTL = NC-1
      IPARTS = 1
*
*   A. First
*
      I1 = IFIRST
      K3 = 3
      K2 = 2
  200 JA = JDIAG(I1,1)
      JC = JDIAG(I1,K3)
*
      IF (JA .EQ. JC) THEN
         IF (NC .GT. 1) THEN
            WRITE (*,300) I1,K3,JA,JC,NC
            STOP
         ENDIF
         NPOINT(1) = I1
         GOTO 900
      ENDIF
*
      I2 = TAB1(JA,K2)
      I3 = TAB1(JC,K2)
*
      IF (ABS(IL(I3)-IL(I2))-NCM .LT. 0) THEN
         WRITE (*,301) I2,I3,JA,JC,K2,NC
         STOP
      ENDIF
*
      IF (ABS(IL(I3)-IL(I2))-NCM .GT. 0) THEN
*
*   B. Last
*
         IF (I1 .NE. IFIRST) GOTO 250
         I1 = ILAST
         K3 = 2
         K2 = 1
         GOTO 200
      ENDIF
*
      IC = 1
      NPOINT(IC) = I1
      I20 = MIN (I2,I3)
      I21 = IL(I20)
      I31 = I21+NCM1
*
      DO 203 II = I21,I31
         IC = IC+1
         NPOINT(IC) = IH(II)
  203 CONTINUE
*
      IF (NC .LE. 2) THEN
         IF (JDIAG(IFIRST,1) .NE. JDIAG(ILAST,1))
     :       CALL PHASE (I1,JDIAG,M4TRD)
         GOTO 900
      ENDIF
*
      IF (I1 .NE. ILAST) THEN
         IT = I2
         JT = JDIAG(ILAST,2)
         K4 = 2
         I4 = ILAST
      ELSE
         IT = I3
         JT = JDIAG(IFIRST,3)
         K4 = 3
         I4 = IFIRST
      ENDIF
*
      IF (IT .EQ. I20) CALL PHASE (I1,JDIAG,M4TRD)
      IF ((JT .EQ.JA) .OR. (JT.EQ.JC)) CALL CHANGE (I4,K4)
      GOTO 900
*
*   2. Two isolated points,first and last
*
  250 IF (NC .EQ. 1) RETURN
      IF (NC .LE. 3) GOTO 100
      IPARTL = NC-2
      IPARTS = 1
      I1 = IFIRST
      I2 = ILAST
      JA = JDIAG(I1,1)
      JB = JDIAG(I1,3)
*
      IF (TAB1(JA,2) .NE. I2) THEN
         JA = JDIAG(I1,3)
         JB = JDIAG(I1,1)
         IF (TAB1(JA,2) .NE. I2) GOTO 100
      ENDIF
*
      IF (JA .EQ. JDIAG(I2,1)) THEN
         JC = JDIAG(I2,2)
      ELSE
         JC = JDIAG(ILAST,1)
      ENDIF
*
      I3 = TAB1(JB,2)
      I4 = TAB1(JC,1)
      IDIST = IL(I4)-IL(I3)
*
      IF (ABS(IDIST)-(NCM-1) .LT. 0) THEN
         WRITE (*,302) I3,I4,JB,JC,IDIST,NC
         STOP
      ENDIF
      IF (ABS(IDIST)-(NCM-1) .EQ. 0) THEN
         NPOINT(1) = ILAST
         NPOINT(2) = IFIRST
         ICROSS = SIGN (1,IDIST)
         IC = 2
         I20 = MIN (I3,I4)
         I21 = IL(I20)
         I31 = I21+NCM
*
         DO 261 II = I21,I31
            IC = IC+1
            NPOINT(IC) = IH(II)
  261    CONTINUE
*
         IF (JA .EQ. JDIAG(IFIRST,1)) CALL CHANGE (IFIRST,3)
         IF (JA .EQ. JDIAG(ILAST,1)) CALL CHANGE (ILAST,2)
         GOTO 900
      ENDIF
*
*   First general case: all points in one group
*
  100 NPART = 1
      IPARTS = 0
      IPARTL = NC
      K3 = 1
*
      DO 101 IN = 1,NBNODE
         I = IH(IN)
  108    JA = JDIAG(I,K3)
         IF (I .NE. TAB1(JA,2))THEN
            I2 = TAB1(JA,2)
*
            IF (IL(I2)-IN-NCM1 .LT. 0) THEN
               WRITE (*,303) IN,I,I2,IL(I2),JA,NC
               STOP
            ENDIF
            IF (IL(I2)-IN-NCM1 .EQ. 0) THEN
               I21 = IL(I2)
               IC = 0
*
               DO 103 II = IN,I21
                  IC = IC+1
                  NPOINT(IC) = IH(II)
  103          CONTINUE
*
               IF (JA .EQ. JDIAG(IFIRST,3)) CALL CHANGE (IFIRST,3)
               IF (JA .EQ. JDIAG(ILAST,2)) CALL CHANGE (ILAST,2)
               GOTO 900
            ENDIF
         ENDIF
*
         IF (IN .EQ. 1) THEN
            IF (K3 .NE. 3) THEN
               K3 = 3
               GOTO 108
            ELSE
               K3 = 1
            ENDIF
         ENDIF
*
  101 CONTINUE
*
*   Search did not find loop NC .LE. 3
*
      IF (NC .LE. 3) RETURN
*
*   General case of loop partitionned in 2 groups. DO loop
*   on IPARTS
*
      NPART = 2
      NC2 = NC/2
      K3 = 1
      K2 = 1
*
      DO 400 IPS = 2,NC2
         JPS = IPS-1
         NBN = NBNODE-JPS
*
         DO 401 I1 = 1,NBN
            I = IH(I1)
            I2 = IH(I1+JPS)
  402       JA = JDIAG(I,K3)
            JD = JDIAG(I2,K2)
*
            IF (I .EQ. TAB1(JA,1)) THEN
               II2 = TAB1(JD,2)
               II1 = TAB1(JA,2)
            ELSE
               II1 = TAB1(JA,1)
               II2 = TAB1(JD,1)
            ENDIF
*
            IDIST = IL(II1)-IL(II2)
*
            IF (ABS (IDIST)-(NCM-JPS) .LT. 0) THEN
               WRITE (*,304) JPS,I1,I,I2,JA,JD,II1,II2,IDIST,NC
               STOP
            ENDIF
            IF (ABS (IDIST)-(NCM-JPS) .GT. 0) GOTO 420
            ICROSS = SIGN (1,IDIST)
            IC = 0
            I21 = IL(I2)
*
            DO 410 II = I1,I21
               IC = IC+1
               NPOINT(IC) = IH(II)
  410       CONTINUE
*
            I20 = MIN (II1,II2)
            I30 = MAX (II1,II2)
            I21 = IL(I20)
            I31 = IL(I30)
*
            DO 411 II = I21,I31
               IC = IC+1
               NPOINT(IC) = IH(II)
  411       CONTINUE
*
            IPARTS = IPS
            IPARTL = NC-IPS
            IF ((JDIAG(IFIRST,3) .EQ. JA) .OR.
     :          (JDIAG(IFIRST,3) .EQ. JD)) CALL CHANGE (IFIRST,3)
            IF ((JDIAG(ILAST,2) .EQ. JA) .OR.
     :          (JDIAG(ILAST,2) .EQ. JD)) CALL CHANGE (ILAST,2)
            GOTO 900
*
  420       IF (I1 .EQ. 1) THEN
               IF (K3 .EQ. 3) THEN
                  K3 = 1
                  GOTO 401
               ELSE
                  K3 = 3
                  GOTO 402
               ENDIF
            ENDIF
*
            IF (I2 .EQ. ILAST) THEN
               IF (K2 .NE. 2) THEN
                  K2 = 2
                  GOTO 402
               ENDIF
            ENDIF
*
  401    CONTINUE
  400 CONTINUE
*
*   SEARCH did not find circuit of order NC
*
      RETURN
*
*   Loop found
*
  900 FIND = .TRUE.
      CALL PRINTJ (NAME,10)
*
      RETURN
*
*   Error printout
*
  300 FORMAT (' Error in SEARCH. I1,K3,JA,JC,NC = ',5I5)
  301 FORMAT (' Error in SEARCH. I2,I3,JA,JC,K2,NC = ',6I5)
  302 FORMAT (' Error in SEARCH. I3,I4,JB,JC,IDIST,NC = ',6I5)
  303 FORMAT (' Error in SEARCH. IN,I,I2,IL(I2),JA,NC = ',6I5)
  304 FORMAT (' Error in SEARCH. JPS,I1,I,I2,JA,JD,II1,II2,IDIST,NC = '
     :       ,10I5)
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE SETDM
*                                                                      *
*   Sets dimensions of arrays.                                         *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      LOGICAL SUMVAR
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /DIM/J6CC,J7CC,J8CC,J9CC,JWCC,JDELC
*
      JWCC = JWC
      JDELC = JDEL
      J6CC = J6C
      J7CC = J7C
      J8CC = J8C
      J9CC = J9C
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE SETTAB (FAIL)
*                                                                      *
*   Builds up the unstructured graph. Sets the array J23, containing   *
*   the  two lists of original triads J2 and J3, and the correspond-   *
*   ing arrows  on the  angular  momenta lines. Also establishes the   *
*   numerical and phase factors  connecting  recoupling  coefficient   *
*   and graphs, according to Yutsis, Levinson, and Vanagas. For this   *
*   purpose determines the total J.                                    *
*                                                                      *
*                                           Last update: 16 Ocy 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARROW
      LOGICAL FAIL,TABS,FREE,SUMVAR
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
     :      /BUILD/IAL(M4TRD),IF1,IF2,NODE
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      DATA NAME/'SETTAB'/
*
      DO 502 I = 1,M2TRD
         DO 501 J = 1,3
            J23(I,J) = 0
  501    CONTINUE
  502 CONTINUE
*
      IPR = N-1
      NBTR = IPR+IPR
*
      DO 4 I = 1,IPR
         DO 5 J = 1,2
            J23(I,J) = J2(I,J)
            ARROW(I,J) = 1
    5    CONTINUE
         TABS(I) = .FALSE.
         J23(I,3) = J2(I,3)
         ARROW(I,3) = -1
    4 CONTINUE
*
      IPR1 = IPR+1
*
      DO 7 I = IPR1,NBTR
         II = I-IPR
         DO 6 J = 1,2
            J23(I,J) = J3(II,J)
            ARROW(I,J) = -1
    6    CONTINUE
         TABS(I) = .FALSE.
         J23(I,3) = J3(II,3)
         ARROW(I,3) = 1
    7 CONTINUE
*
      DO 11 J = 1,NBTR
         J8(J) = J23(J,1)
   11 CONTINUE
*
      J8C = NBTR+IPR
      NB1 = NBTR+1
*
      DO 12 J = NB1,J8C
         I = J-IPR
         J8(J) = J23(I,3)
   12 CONTINUE
*
      J6C = NBTR
*
      DO 13 J = 1,J6C
         J6(J) = J23(J,3)
   13 CONTINUE
*
      DO 10 I = 1,M
         SUMVAR(I) = .FALSE.
         IAL(I) = 1
   10 CONTINUE
*
      DO 9 I = 1,NBTR
         DO 8 J = 1,3
            JI = J23(I,J)
            K = IAL(JI)
            LINE(JI,K) = I
            LCOL(JI,K) = J
            IAL(JI) = K+1
    8    CONTINUE
    9 CONTINUE
*
      IT = 0
*
      DO 18 I = 1,NBTR
*
         JT = J23(I,3)
*
         IF (IAL(JT) .EQ. 3) THEN
*
            CALL OTHERJ (I,JT,L,LC,K)
            IF (LC .EQ. 3) GOTO 19
*
         ELSE
*
            IF (IT .EQ. 1) THEN
               CALL DELTA (JT1,JT,FAIL)
               IF (FAIL) GOTO 20
               K = LINE(JT,1)
               KC = LCOL(JT,1)
               LINE(JT1,2) = K
               LCOL(JT1,2) = KC
               LINE(JT,2) = LINE(JT1,1)
               LCOL(JT,2) = LCOL(JT1,1)
               J23(K,KC) = JT1
               IAL(JT) = 1
               GOTO 19
            ENDIF
*
            JT1 = JT
            IT = 1
*
         ENDIF
*
   18 CONTINUE
*
   19 J9(J9C+1) = JT
      J9C = J9C+2
      J9(J9C) = JT
*
   20 CALL PRINTJ (NAME,4)
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE SPRATE (M)
*                                                                      *
*   This  subroutine  prepares  the  information to be transfered to   *
*   GENSUM for numerical evaluation.                                   *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME
      CHARACTER*5 NME
      LOGICAL SUM6J,T6J,JT,JS,SUMVAR,CUT
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      DIMENSION JTEM4(MTRIAD,M6J),JTEM5(MTRIAD,M6J),JTEM6(MTRIAD),
     :   NSUM6J(M6J),J6SUM(M6J)
      DIMENSION SUM6J(M6J),T6J(M6J),JT(MTRIAD),JS(MTRIAD),
     :   INVER(MANGM),JNSUM(MTRIAD),JINV(MTRIAD),N6JN(M6J),IN6J(M6J),
     :   JSUMT(M6J,6)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),JW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /CUTDIG/CUT
     :      /DIM/J6CC,J7CC,J8CC,J9CC,JWCC,JDELC
     :      /SUMARG/J6P(MANGMP),J7P(MANGMP),J8P(MANGMP),J9P(MANGMP),
     :              JWORD(6,M6J),NLSUM,NBJ(MSUM),NB6J(MSUM),K6CP(MSUM),
     :              K7CP(MSUM),K8CP(MSUM),K9CP(MSUM),JSUM6(MTRIAD),
     :              JSUM4(MTRIAD,M6J),JSUM5(MTRIAD,M6J),INV6J(M6J)
*
*   Test that array dimensions have not been exceeded.
*
      IF (MP .GT. MANGM) THEN
         NMX = MANGM
         NPX = MP
         NAME = 'MANGM '
         NME  = 'MP   '
      ELSEIF (JWC .GT. M6J) THEN
         NMX = M6J
         NPX = JWC
         NAME = 'M6J   '
         NME  = 'JWC  '
      ELSEIF (J6C .GT. M3MNGM) THEN
         NMX = M3MNGM
         NPX = J6C
         NAME = 'M3MNGM'
         NME  = 'J6C  '
      ELSEIF (J7C .GT. M3MNGM) THEN
         NMX = M3MNGM
         NPX = J7C
         NAME = 'M3MNGM'
         NME  = 'J7C  '
      ELSEIF (J8C .GT. M3MNGM) THEN
         NMX = M3MNGM
         NPX = J8C
         NAME = 'M3MNGM'
         NME  = 'J8C  '
      ELSE
         IF (J9C .LE. MANGMP) GOTO 54
         NMX = MANGMP
         NPX = J9C
         NAME = 'MANGMP'
         NME  = 'J9C  '
      ENDIF
*
   60 WRITE (*,300) NAME,NME,NPX,NMX
      STOP
*
*   Determination of effective summation variables and their
*   relationships with 6j coefficients.
*
   54 DO 2 I = 1,JWC
         INV6J(I) = 0
         SUM6J(I) = .FALSE.
    2 CONTINUE
*
      NSUM = 0
      NLSUM = 0
      IF (MP .EQ. M) RETURN
      M1 = M+1
*
      DO 1 I = M1,MP
         IF (SUMVAR(I)) THEN
            NSUM = NSUM+1
            JSUM6(NSUM) = 0
            INVER(I) = NSUM
         ENDIF
    1 CONTINUE
*
      IF (NSUM .EQ. 0) RETURN
*
      IF (NSUM .GT. MTRIAD) THEN
         NMX = MTRIAD
         NPX = NSUM
         NAME = 'MTRIAD'
         NME  = 'NSUM '
         GOTO 60
      ENDIF
*
      KT = 0
*
      DO 4 I = 1,JWC
         DO 5 J = 1,6
            IK = JW(J,I)
            IF (.NOT. SUMVAR(IK)) GOTO 5
*
            IF (.NOT. SUM6J(I)) THEN
               SUM6J(I) = .TRUE.
               KT = KT+1
               J6SUM(KT) = 0
               NSUM6J(KT) = I
               INV6J(I) = KT
            ENDIF
*
            ISK = INVER(IK)
            I2 = JSUM6(ISK)+1
            JSUM6(ISK) = I2
            JSUM4(ISK,I2) = J
            JSUM5(ISK,I2) = KT
            I3 = J6SUM(KT)+1
            J6SUM(KT) = I3
            JSUMT(KT,I3) = ISK
    5    CONTINUE
    4 CONTINUE
*
      CALL VAR (J6,J6P,J6C,J6CP,J6CC,SUMVAR,MP,M,INVER)
      CALL VAR (J7,J7P,J7C,J7CP,J7CC,SUMVAR,MP,M,INVER)
      CALL VAR (J8,J8P,J8C,J8CP,J8CC,SUMVAR,MP,M,INVER)
      CALL VAR (J9,J9P,J9C,J9CP,J9CC,SUMVAR,MP,M,INVER)
*
      IF (.NOT. CUT) THEN
         NLSUM = 1
         NBJ(1) = NSUM
         NB6J(1) = KT
         K6CP(1) = J6CP
         K7CP(1) = J7CP
         K8CP(1) = J8CP
         K9CP(1) = J9CP
*
         DO 21 I = 1,KT
            I1 = NSUM6J(I)
            DO 22 J = 1,6
               JWORD(J,I) = JW(J,I1)
   22       CONTINUE
   21    CONTINUE
*
         DO 80 I = 1,NSUM
            ISU = JSUM6(I)
            DO 81 J = 1,ISU
               I1 = JSUM5(I,J)
               J1 = JSUM4(I,J)
               JWORD(J1,I1) = MP+I
   81       CONTINUE
   80    CONTINUE
*
         RETURN
      ENDIF
*
*   Separation of variables and sums in case a cut was detected.
*
      K6C = 0
      K7C = 0
      K8C = 0
      K9C = 0
      NJ = 0
      N6J = 0
*
      DO 9 I = 1,KT
         T6J(I) = .FALSE.
    9 CONTINUE
*
      DO 7 I = 1,NSUM
         JT(I) = .FALSE.
         JS(I) = .FALSE.
    7 CONTINUE
*
      J = 1
*
   10 NJ = NJ+1
      JNSUM(NJ) = J
      JINV(J) = NJ
      JT(J) = .TRUE.
   18 JS(J) = .TRUE.
      JS6 = JSUM6(J)
*
      DO 11 I = 1,JS6
         I6J = JSUM5(J,I)
*
         IF (.NOT. T6J(I6J)) THEN
            T6J(I6J) = .TRUE.
            N6J = N6J+1
            N6JN(N6J) = NSUM6J(I6J)
            IN6J(I6J) = N6J
         ENDIF
*
         J6J = J6SUM(I6J)
*
         DO 12 K = 1,J6J
            JK = JSUMT(I6J,K)
            IF (.NOT. JT(JK)) THEN
               NJ = NJ+1
               JNSUM(NJ) = JK
               JINV(JK) = NJ
               JT(JK) = .TRUE.
            ENDIF
   12    CONTINUE
*
   11 CONTINUE
*
      DO 13 JJ = 1,NSUM
         J = JJ
         IF ((.NOT. JS(JJ)) .AND. JT(JJ)) GOTO 18
   13 CONTINUE
*
      NLSUM = NLSUM+1
*
      IF (NLSUM .GT. MSUM) THEN
         NMX = MSUM
         NPX = NLSUM
         NAME = 'MSUM  '
         NME = 'NLSUM'
         GOTO 60
      ENDIF
      NBJ(NLSUM) = NJ
      NB6J(NLSUM) = N6J
*
      IF (J6CP .NE. 0) CALL CHVAR (J6P,J6CP,K6C,JT,JINV,NSUM)
      K6CP(NLSUM) = K6C
      IF (J7CP .NE. 0) CALL CHVAR (J7P,J7CP,K7C,JT,JINV,NSUM)
      K7CP(NLSUM) = K7C
      IF (J8CP .NE. 0) CALL CHVAR (J8P,J8CP,K8C,JT,JINV,NSUM)
      K8CP(NLSUM) = K8C
      IF (J9CP .NE. 0) CALL CHVAR (J9P,J9CP,K9C,JT,JINV,NSUM)
      K9CP(NLSUM) = K9C
*
      IF (NJ .NE. NSUM) THEN
         DO 16 JJ = 1,NSUM
            J = JJ
            IF (.NOT. JT(JJ)) GOTO 10
   16    CONTINUE
      ENDIF
*
      DO 26 I = 1,KT
         I1 = N6JN(I)
         DO 27 J = 1,6
            JWORD(J,I) = JW(J,I1)
   27    CONTINUE
   26 CONTINUE
*
      DO 28 I = 1,NSUM
         IK = JNSUM(I)
         I2 = JSUM6(IK)
         JTEM6(I) = I2
         DO 29 J = 1,I2
            JTEM4(I,J) = JSUM4(IK,J)
            K = JSUM5(IK,J)
            JTEM5(I,J) = IN6J(K)
   29    CONTINUE
   28 CONTINUE
*
      DO 40 I = 1,NSUM
         I2 = JTEM6(I)
         JSUM6(I) = I2
         DO 41 J = 1,I2
            I1 = JTEM5(I,J)
            J1 = JTEM4(I,J)
            JSUM4(I,J) = J1
            JSUM5(I,J) = I1
            JWORD(J1,I1) = I+MP
   41    CONTINUE
   40 CONTINUE
*
      RETURN
*
  300 FORMAT (' Dimension error for ',A6
     :       /2X,A5,' = ',I5,' is out of allowed range',I4)
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE SQUARE
*                                                                      *
*   Reduces  a  circuit  of  order 4 in the two cases which are left   *
*   over by POLYGN, namely two disconnected groups of two points and   *
*   one group of two points plus  the  two  ends of the axis. In the   *
*   latter, the end of the axis is transferred  to the beginning. In   *
*   this  process,  one summation variable and two  6j  symbols  are   *
*   introduced.                                                        *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME,NAMSUB
      INTEGER ARR,TAB1
      LOGICAL SUMVAR
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /NAM/NAMSUB
*
      DATA NAME/'SQUARE'/
*
      NAMSUB = NAME
      MP = MP+1
      SUMVAR(MP) = .TRUE.
      IT1 = NPOINT(1)
      IT2 = NPOINT(2)
*
      IF (ICROSS .EQ. 1) THEN
         IT3 = NPOINT(3)
         IT4 = NPOINT(4)
         K23 = 3
         K32 = 2
      ELSE
         IT3 = NPOINT(4)
         IT4 = NPOINT(3)
         K23 = 2
         K32 = 3
      ENDIF
*
      L4 = JDIAG(IT2,1)
*
      IF (ARR(IT2,1) .LE. 0) THEN
         CALL PHASE2 (L4)
         ARR(IT2,1) = 1
         ARR(IT3,1) = -1
      ENDIF
*
      L2 = JDIAG(IT1,1)
      IF (ARR(IT1,1) .GT. 0) CALL PHASE2 (L2)
      JWC = JWC+1
      KW(1,JWC) = L4
      KW(2,JWC) = L2
      KW(3,JWC) = JDIAG(IT2,2)
      JJ1 = JDIAG(IT1,3)
      KW(4,JWC) = JJ1
      KW(5,JWC) = MP
      KW(6,JWC) = JDIAG(IT1,2)
      IF (ARR(IT1,2) .LT. 0) CALL PHASE2 (JDIAG(IT1,2))
      JWC = JWC+1
      KW(1,JWC) = L4
      KW(2,JWC) = L2
      JJ3 = JDIAG(IT3,K23)
      JJ2 = JDIAG(IT4,K32)
      KW(3,JWC) = JJ3
      KW(4,JWC) = JJ2
      KW(5,JWC) = MP
      KW(6,JWC) = JDIAG(IT3,K32)
      IF (ARR(IT3,K32) .LT. 0) CALL PHASE2 (JDIAG(IT3,K32))
      J6(J6C+1) = MP
      J6C = J6C+2
      J6(J6C) = MP
*
      IF (NPART .EQ. 1) THEN
         ITMIN = IT2
         ITMAX = IT3
      ELSE
         ITMIN = MIN (IT2,IT3)
         ITMAX = MAX (IT2,IT3)
      ENDIF
      ITMN = MIN (IT1,IT4)
      ITMX = MAX (IT1,IT4)
*
      TAB1(MP,1) = ITMIN
      TAB1(MP,2) = ITMAX
      JDIAG(IT2,1) = MP
      JDIAG(IT3,1) = MP
      JDIAG(IT2,3) = JJ1
      ARR(IT2,3) = ARR(IT1,3)
      JDIAG(IT3,K32) = JJ2
      ARR(IT3,K32) = ARR(IT4,K32)
*
      IF (ICROSS .EQ. 1) THEN
         J7(J7C+1) = L2
         J7(J7C+2) = L4
         CALL PHASE2 (L4)
         J7C = J7C+3
         J7(J7C) = MP
      ELSE
         CALL PHASE2 (JJ2)
      ENDIF
*
      ITLL = IL(ITMN)
      ITHL = IL(ITMX)
*
      DO 5 I = ITLL+1,ITHL-1
         IT = IH(I)
         ILP = I-1
         IL(IT) = ILP
         IH(ILP) = IT
    5 CONTINUE
      IF (ITHL .NE. NBNODE) THEN
         DO 6 I = ITHL+1,NBNODE
            IT = IH(I)
            ILP = I-2
            IL(IT) = ILP
            IH(ILP) = IT
    6    CONTINUE
      ENDIF
*
      IF (NPART .NE. 2) THEN
         TAB1(JJ1,1) = IH(1)
         TAB1(JJ1,2) = IH(NBNODE-2)
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE TRDEL (JJ1,JJ2,JJ3,NBN,FAIL)
*                                                                      *
*   Test for triangular delta. If not satisfied FAIL = .TRUE. .        *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      LOGICAL FAIL,SUMVAR,CUT,FREE
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /CUTDIG/CUT
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
*
      IF (SUMVAR(JJ1) .OR. SUMVAR(JJ2) .OR. SUMVAR(JJ3)) RETURN
      IF (NBN .GT. 4) CUT = .TRUE.
      IF ((.NOT. FREE(JJ1)) .AND. (.NOT. FREE(JJ2)) .AND.
     :    (.NOT.FREE(JJ3))) THEN
         I1 = J1(JJ1)
         I2 = J1(JJ2)
         I3 = J1(JJ3)
         IF ((I1 .LT. (ABS (I2-I3)+1)) .OR. (I1 .GT. (I2+I3-1)))
     :      FAIL = .TRUE.
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE TRIANG (FAIL)
*                                                                      *
*   Reduces  a triangle having one apex at either end of the axis of   *
*   the flat  diagram.  This introduces one 6j symbol and some phase   *
*   factors.                                                           *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME,NAMSUB
      INTEGER ARR,TAB1
      LOGICAL FAIL,SUMVAR
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /NAM/NAMSUB
*
      DATA NAME/'TRIANG'/
*
      NAMSUB = NAME
      IT1 = NPOINT(1)
      IT2 = NPOINT(2)
      IT3 = NPOINT(3)
      JWC = JWC+1
      KW(1,JWC) = JDIAG(IT3,2)
      KW(2,JWC) = JDIAG(IT2,3)
      KW(3,JWC) = JDIAG(IT3,1)
      IF (ARR(IT3,1) .GT. 0) CALL PHASE2 (KW(3,JWC))
      KW(4,JWC) = JDIAG(IT2,1)
      IF (ARR(IT2,1) .LT. 0) CALL PHASE2 (KW(4,JWC))
      K23 = 3
      IF (IT1 .EQ. IFIRST) K23 = 2
      KW(5,JWC) = JDIAG(IT1,K23)
      KW(6,JWC) = JDIAG(IT3,3)
      CALL TRDEL (KW(1,JWC),KW(2,JWC),KW(5,JWC),NBNODE,FAIL)
      IF (FAIL) GOTO 15
      IF (ARR(IT3,3) .GT. 0) CALL PHASE2 (KW(6,JWC))
      JT1 = KW(5,JWC)
      JDIAG(IT3,1) = JT1
      JDIAG(IT3,3) = KW(2,JWC)
      ARR(IT3,1) = ARR(IT1,K23)
      ARR(IT3,3) = ARR(IT2,3)
*
      IF (IT1 .NE. IFIRST) THEN
         TAB1(JT1,1) = IT3
         TAB1(JT1,2) = IH(NBNODE-1)
         K12 = 1
      ELSE
         TAB1(JT1,1) = IH(2)
         TAB1(JT1,2) = IT3
         K12 = 2
      ENDIF
*
      IL3 = IL(IT3)
*
      IF (IT1 .NE. ILAST) THEN
         IL2 = IL(IT2)-1
*
         DO 2 I = 2,IL2
            IT = IH(I)
            ILP = I-1
            IL(IT) = ILP
            IH(ILP) = IT
    2    CONTINUE
      ENDIF
*
      DO 1 I = IL3,NBNODE
         IT = IH(I)
         ILP = I-K12
         IL(IT) = ILP
         IH(ILP) = IT
    1 CONTINUE
*
   15 RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE VAR (JN,JNS,JNC,JNSC,JBC,SUMVAR,MP,M,INVER)
*                                                                      *
*   Test  for  variable  character and put in JNS if yes, and JN now   *
*   contains 0.                                                        *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      LOGICAL SUMVAR(MP)
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      DIMENSION JN(JNC),JNS(MANGMP),INVER(MP)
*
      JNSC = 0
      IF (JBC .NE. JNC) THEN
         JBBC = JBC+1
*
         DO 1 I = JBBC,JNC
            I1 = JN(I)
            IF (SUMVAR(I1)) THEN
               JNSC = JNSC+1
               IF (JNSC .GT. MANGMP) THEN
                  WRITE (*,300) JNSC,MANGMP
                  STOP
               ENDIF
               J = INVER(I1)
               JNS(JNSC) = J
               JN(I) = M
            ENDIF
    1    CONTINUE
      ENDIF
*
      RETURN
*
  300 FORMAT (' Dimension error in VAR. JNSC = ',I5,' MANGMP = ',I5)
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE WAY (L,KA,KB,ICH,NB)
*                                                                      *
*   Tests  one  step  forward  if  the way is free. First and second   *
*   arguments are interchanged or not according to ICH = -1, or +1.    *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INTEGER ARROW
      LOGICAL TABS
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /BUILD/IAL(M4TRD),IF1,IF2,NODE
*
      K1 = J23(L,KA)
      K2 = J23(L,KB)
      NB = IAL(K1)+IAL(K2)-1
      IF (NB) 3,2,8
    2 NB1 = IAL(K1)-IAL(K2)
      IF (NB1) 9,8,8
    3 CALL OTHERJ (L,K1,L1,LC1,LA)
      CALL OTHERJ (L,K2,L2,LC2,LB)
      CALL NEIBOR (LC1,I1,I2)
      CALL NEIBOR (LC2,I3,I4)
      JI1 = J23(L1,I1)
      JI2 = J23(L1,I2)
      JI3 = J23(L2,I3)
      JI4 = J23(L2,I4)
      IA = IAL(JI1)+IAL(JI2)
      IB = IAL(JI3)+IAL(JI4)
      NBP = IB+IA+1
      NBM = IB-IA
      GOTO (8,4,5,4,6),NBP
    4 IF (NBM) 9,8,8
    5 IF (NBM) 9,6,8
    6 IF ((JI3 .EQ. IF1) .OR. (JI3 .EQ. IF2) .OR.
     :    (JI4 .EQ. IF1) .OR. (JI4 .EQ. IF2)) GOTO 9
    8 ICH = 1
      GOTO 10
    9 ICH = -1
   10 RETURN
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE ZERO (J,JZ,FAIL)
*                                                                      *
*   Suppresses  one  line  and  two  nodes of the unstructured graph   *
*   introduces  zeros in the triads  J23. As a consequence the other   *
*   two arguments of the triad are put equal. If there was already a   *
*   zero in the triad which is changed, it is a special case.          *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARROW
      LOGICAL FAIL,TABS,FREE,SUMVAR,CUT,NOCUT
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10,MZERO = 20)
*
      COMMON/ZER/NZERO,JZERO(MZERO)
     :      /CUTDIG/CUT
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
     :      /KEEP/JKP(2,3),JARR(2,3),IT2,IT3,IT5
     :      /BUILD/IAL(M4TRD),IF1,IF2,NODE
     :      /TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      DATA NAME/'ZERO  '/
*
      NOCUT = .FALSE.
      NZERO = 0
*
      IF (J .GE. 1) THEN
         CALL OTHERJ (0,JZ,LIN,LC,K1)
         I = NZERO
         GOTO 8
      ENDIF
*
      DO 11 I = 1,M
         IF ((J1(I) .NE. 1) .OR. FREE(I) .OR. (IAL(I).LE.1)) GOTO 11
         NZERO = NZERO+1
         IF (NZERO .GT. MZERO) THEN
            WRITE (*,300) NZERO,MZERO
            STOP
         ENDIF
         JZERO(NZERO) = I
   11 CONTINUE
*
      NOCUT = .TRUE.
      M = M+1
      J1(M) = 1
      SUMVAR(M) = .FALSE.
      FREE(M) = .FALSE.
      IF (NZERO .EQ. 0) GOTO 7
      CALL PRINTJ (NAME,1)
      I = 0
    1 I = I+1
      JZ = JZERO(I)
      J = 0
   13 J = J+1
      LIN = LINE(JZ,J)
      IF (TABS(LIN)) GOTO 2
      LC = LCOL(JZ,J)
    8 CALL NEIBOR (LC,L1,L2)
      JJ1 = J23(LIN,L1)
      JJ2 = J23(LIN,L2)
*
      IF (JJ1 .EQ. JJ2) THEN
         J6C = J6C+1
         J6(J6C) = JJ1
         LO1=LIN
         LO2=LIN
         LCO1=L1
         LCO2=L2
         GOTO 10
      ENDIF
*
      CALL DELTA (JJ1,JJ2,FAIL)
      IF (FAIL) GOTO 7
*
      IF ((J1(JJ1) .NE. 1) .AND. (J1(JJ2) .NE. 1)) GOTO 15
      IF (J1(JJ1) .LT. J1(JJ2)) GOTO 15
      IF (J1(JJ1) .GT. J1(JJ2)) GOTO 19
*
      IF (NZERO .NE. 0) THEN
         DO 17 JJX = I,NZERO
            JJZ = JZERO(JJX)
            IF (JJ1 .EQ. JJZ) GOTO 15
            IF (JJ2 .EQ. JJZ) GOTO 19
   17    CONTINUE
      ENDIF
*
      GOTO 15
*
   19 JJZ = JJ2
      JJ2 = JJ1
      JJ1 = JJZ
*
   15 CALL OTHERJ (LIN,JJ1,LO1,LCO1,K1)
      CALL OTHERJ (LIN,JJ2,LO2,LCO2,K2)
      J9C = J9C+1
      J9(J9C) = JJ1
      J23(LO2,LCO2) = JJ1
      LINE(JJ1,K1) = LO2
      LCOL(JJ1,K1) = LCO2
*
   10 IF     (ARROW(LIN,L1) .LT. ARROW(LIN,L2)) THEN
         CALL PHASE2 (JJ1)
      ELSEIF (ARROW(LIN,L1) .EQ. ARROW(LIN,L2)) THEN
         ARROW(LO1,LCO1) = 1
         ARROW(LO2,LCO2) = -1
      ENDIF
*
      TABS(LIN) = .TRUE.
      NBTR = NBTR-1
      IF (NBTR .EQ. 0) GOTO 7
      IF (LO1 .EQ. LO2) THEN
         L = 6-LCO1-LCO2
         JT = J23(LO1,L)
         IF ((J1(JT) .EQ. 1) .AND. (.NOT.FREE(JT))) GOTO 2
         CALL DELTA (JT,M,FAIL)
         IF (FAIL) GOTO 7
         NZERO = NZERO+1
         JZERO(NZERO) = JT
      ENDIF
    2 IF (J .EQ. 1) GOTO 13
*
      IF (NBTR .NE. 0) THEN
         IF (I .LT. NZERO) GOTO 1
      ENDIF
*
    7 CALL PRINTJ (NAME,4)
      IF (NOCUT) CUT = .FALSE.
*
      RETURN
*
  300 FORMAT (' Dimension error in ZERO. NZERO = ',I5,' MZERO = ',I5)
*
      END
