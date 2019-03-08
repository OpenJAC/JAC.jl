************************************************************************
*
* This file contains several procedures from the GRASP-92 package 
* to enable the calculation of angular coefficients. Since parts of this 
* code date back to the late sixties, this code is (almost) impossible
* to modify; this has hampared the implementation of a more efficient
* computation of angular coefficients in the past. 
*
* No attempt has been made to 'improve' the readability of this code; 
* a few modifications need to be done as are clearly indicated beloww.
* This file just list all required procedures (as a list of external 
* Fortran-77 procedures without any explicit interface) in alphabetic 
* order. -- In a long term, we intent to replace this file and to 
* implement a more efficient scheme using the features of Fortran-90.
*
************************************************************************
*                                                                      *
      SUBROUTINE BREID (JA,JB,JA1,IPCA,JB1)
*                                                                      *
*   Computes closed shell contributions - aaaa and exchange only.      *
*                                                                      *
*   Call(s) to: [LIB92]: CLRX, CXK, ITRIG, TALK, SNRC.                 *
*                                                                      *
*                                           LAST UPDATE: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      DIMENSION CONE(7,20),JS(4),KAPS(4),KS(4),S(12)
*
      COMMON/BCORE/ICORE(149)
     :      /CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /M1/NQ1(149),NQ2(149)
     :      /M3/JLIST(149),KLIST(149),NPEEL,NCORE
     :      /ORB4/NP(149),NAK(149)
*
      PARAMETER (NUMAX = 20)
*
*   1.0  Initialization
*
      IF (IPCA .EQ. 2) THEN
         IA1 = KLIST(JA1)
      ELSE
         IA1 = JLIST(JA1)
      ENDIF
      IB1 = KLIST(JB1)
*
      ISG = 1
      IF (JA .EQ. JB) THEN
         IF ((ICORE(IA1) .NE. 0) .AND. (ICORE(IB1) .NE. 0)) THEN
            IF (JA .GT. 1) RETURN
            ISG = -1
         ENDIF
      ENDIF
*
      JS(1) = IA1
      JS(2) = IB1
      JS(3) = IA1
      JS(4) = IB1
      NQS1 = NQ1(IA1)
      NQS2 = NQ2(IB1)
      DO 1 I = 1,4
         KAPS(I) = 2*NAK(JS(I))
         KS(I) = IABS (KAPS(I))
    1 CONTINUE
      CONST = NQS1*NQS2
      IF (IBUG2 .NE. 0) WRITE (99,300) IA1,IB1
*
*   2.0  Set range of tensor indices
*
      CALL SNRC (JS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
      IF (IBUG2 .NE. 0) WRITE (99,301) ND1,ND2,NE1,NE2,IBRD,IBRE
      IF (IA1 .NE. IB1) GOTO 3
*
*   3.0 Calculate aaaa interaction
*
      DO 2 N = 1,ND2
         NU = ND1+2*(N-1)
         K = NU
         IF (MOD (K,2) .NE. 1) RETURN
         KAP1 = KAPS(1)/2
         GAM = CLRX (KAP1,NU,KAP1)
         DKSKS = KS(1)*KS(1)
         DNUNU1 = NU*(NU+1)
         COEF = CONST*TWO*DKSKS*GAM*GAM/DNUNU1
         IF (IBUG2 .NE. 0) WRITE (99,302) NU,GAM,COEF
         ITYPE = ISG*4
         CALL TALK (JA,JB,NU,IA1,IA1,IA1,IA1,ITYPE,COEF)
    2 CONTINUE
      RETURN
*
*   Calculate exchange interactions
*
    3 CONTINUE
      IF (IBRE .LT. 0) RETURN
      IF (NE2 .GT. NUMAX) THEN
         WRITE (*,304)
         STOP
      ENDIF
*
      DO 4 N = 1,NE2
         DO 4 MU = 1,7
            CONE(MU,N) = ZERO
    4 CONTINUE
*
      PROC = -CONST/DBLE (KS(1)*KS(2))
*
*   Negative sign arises from Pauli phase factor
*
      DO 10 N = 1,NE2
         NU = NE1+2*(N-1)
         K = NU
         IP = (KS(1)-KS(2))/2+K
         IPP = IP+1
         IF (NU .EQ. 0) GOTO 8
         KK = K+K+1
         IF (ITRIG(KS(1),KS(2),KK) .EQ. 0) GOTO 6
         PROD = PROC
         IF (MOD (IP,2) .NE. 0) PROD = -PROD
         CALL CXK (S,JS,KAPS,NU,K,IBRE,2)
         IF (IBUG2 .NE. 0) WRITE (99,303) PROD,(S(MU),MU = 1,3)
         DO 5 MU = 1,3
            CONE (MU,N) = CONE(MU,N)+PROD*S(MU)
    5    CONTINUE
*
    6    K = NU-1
         KK = K+K+1
         IF (ITRIG(KS(1),KS(2),KK) .EQ. 0) GOTO 8
         PROD = PROC
         IF (MOD (IPP,2) .NE. 0) PROD = -PROD
         CALL CXK (S,JS,KAPS,NU,K,IBRE,2)
         IF (IBUG2 .NE. 0) WRITE (99,303) PROD,(S(MU),MU = 1,3)
         DO 7 MU = 1,3
            CONE(MU,N) = CONE(MU,N)+PROD*S(MU)
    7    CONTINUE
*
    8    IF (N .EQ. NE2) GOTO 11
         K = NU+1
         KK = K+K+1
         PROD = PROC
         IF (MOD (IPP,2) .NE. 0) PROD = -PROD
         CALL CXK (S,JS,KAPS,NU,K,IBRE,2)
         IF (IBUG2 .NE. 0) WRITE (99,303) PROD,(S(MU),MU = 1,7)
         DO 9 MU = 1,7
            CONE (MU,N) = CONE(MU,N)+PROD*S(MU)
    9    CONTINUE
   10 CONTINUE
*
*   4.0  Output results
*
   11 CONTINUE
*
      DO 12 N = 1,NE2
         NU = NE1+2*(N-1)
         ITYPE = ISG*5
         CALL TALK (JA,JB,NU,IB1,IA1,IB1,IA1,ITYPE,CONE(1,N))
         CALL TALK (JA,JB,NU,IA1,IB1,IB1,IA1,ITYPE,CONE(2,N))
         CALL TALK (JA,JB,NU,IA1,IB1,IA1,IB1,ITYPE,CONE(3,N))
         IF (N .EQ. NE2) GOTO 12
         NUP1 = NU+1
         ITYPE = ISG*6
         CALL TALK (JA,JB,NUP1,IA1,IB1,IA1,IB1,ITYPE,CONE(4,N))
         CALL TALK (JA,JB,NUP1,IB1,IA1,IB1,IA1,ITYPE,CONE(5,N))
         CALL TALK (JA,JB,NUP1,IA1,IB1,IB1,IA1,ITYPE,CONE(6,N))
         CALL TALK (JA,JB,NUP1,IB1,IA1,IA1,IB1,ITYPE,CONE(7,N))
   12 CONTINUE
      RETURN
*
  300 FORMAT ('BREID: orbitals ',2I3)
  301 FORMAT (2X,'ND1 ND2 NE1 NE2 IBRD IBRE ',6I5)
  302 FORMAT (2X,'aaaa contribution: NU,GAM,COEF',I5,2(3X,1PD15.8))
  303 FORMAT (2X,'PROD = ',1PD15.8
     :       /' S',7D15.8)
  304 FORMAT ('BREID: Dimension error for NUMAX.')
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE BREIT (JA,JB,JA1,JB1,JA2,JB2)
*                                                                      *
*   Computes  the  coefficients  appearing in EQS. 5, 8, 9 AND 10 OF   *
*   I P Grant and B J McKenzie,  J Phys B 13 (1980) 2671--2681.  The   *
*   coefficients for each choice of  orbitals JA1, JB1, JA2, and JB2   *
*   depend on two further parameters  NU and K;  there are IMU inte-   *
*   grals for each such choice, where:                                 *
*                                                                      *
*                  IMU = 4          TYPE = 1                           *
*                        8                 2                           *
*                        1                 3                           *
*                        1                 4                           *
*                        3                 5                           *
*                        4                 6                           *
*                                                                      *
*   See the paper cited above for details.                             *
*                                                                      *
*   Call(s) to: [LIB92]: GENSUM, ITRIG, KNJ, LTAB, MODJ23, MUMDAD,     *
*                        NJGRAF, OCON, SETJ                            *
*               [RCI92]: CXK, SNRC, TALK.                              *
*                                                                      *
*                                           Last update: 14 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INTEGER PNTRIQ
*
      PARAMETER (
     :   MANGM=60,M3MNGM=3*MANGM,MANGMP=2*(MANGM/3),
     :   MTRIAD=12,M2TRD=2*MTRIAD,M4TRD=4*MTRIAD,
     :   M6J=20,MSUM=10)
*
      LOGICAL FREE,DSUMVR,ESUMVR,FAILD,FAILE
*
      DIMENSION COND(12,20),CONE(12,20),S(12)
      DIMENSION IS(4),KAPS(4),KS(4),NQS(4),ILS(4),LLS(4),IT1(4),IROWS(4)
      DIMENSION JS(4)
*
      DIMENSION JD6(M3MNGM),JD7(M3MNGM),JD8(M3MNGM),
     : JD9(MANGMP),KDW(6,M6J),LDDEL(M6J,2),DSUMVR(MANGM)
      DIMENSION JD6P(MANGMP),JD7P(MANGMP),JD8P(MANGMP),JD9P(MANGMP),
     : JDWORD(6,M6J),
     : NDBJ(MSUM),NDB6J(MSUM),KD6CP(MSUM),KD7CP(MSUM),KD8CP(MSUM),
     : KD9CP(MSUM),JDSUM6(MTRIAD),JDSUM4(MTRIAD,M6J),JDSUM5(MTRIAD,M6J),
     : INVD6J(M6J)
      DIMENSION JE6(M3MNGM),JE7(M3MNGM),JE8(M3MNGM),
     : JE9(MANGMP),KEW(6,M6J),LEDEL(M6J,2),ESUMVR(MANGM)
      DIMENSION JE6P(MANGMP),JE7P(MANGMP),JE8P(MANGMP),JE9P(MANGMP),
     : JEWORD(6,M6J),
     : NEBJ(MSUM),NEB6J(MSUM),KE6CP(MSUM),KE7CP(MSUM),KE8CP(MSUM),
     : KE9CP(MSUM),JESUM6(MTRIAD),JESUM4(MTRIAD,M6J),JESUM5(MTRIAD,M6J),
     : INVE6J(M6J)
*
      COMMON/COUPLE/MJA,NJA,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /L1/JBQ1(3,149),JBQ2(3,149),JTQ1(3),JTQ2(3)
     :      /L2/J2S(MTRIAD,3),J3S(MTRIAD,3)
     :      /M1/NQ1(149),NQ2(149)
     :      /M2/JJQ1(3,149),JJQ2(3,149)
     :      /M3/JLIST(149),KLIST(149),NPEEL,NCORE
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(149),NAK(149)
     :      /TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
*
      PARAMETER (EPS = 1.0D-10)
      PARAMETER (NUMAX = 20)
*
*   1.0  Initialize pointers and flags and set any
*        tables required.
*
*        In this segment, the array IS points to the
*        full list of orbitals, the array JS to the
*        array JLIST of peel orbital pointerS.
*
*   1.1  Initialization
*
      JS(1) = JA1
      JS(2) = JB1
      JS(3) = JA2
      JS(4) = JB2
      DO 1 I = 1,4
        IS(I) = JLIST(JS(I))
        KAPS(I) = 2*NAK(IS(I))
        KS(I) = ABS (KAPS(I))
    1 CONTINUE
      IA1 = IS(1)
      IB1 = IS(2)
      IA2 = IS(3)
      IB2 = IS(4)
      NQS(1) = NQ1(IA1)
      NQS(2) = NQ1(IB1)
      NQS(3) = NQ2(IA2)
      NQS(4) = NQ2(IB2)
*
      KJ23 = 0
      ISNJ = 0
*
      FAILD = .FALSE.
      FAILE = .FALSE.
      NBRJ = 3*NPEEL + 7
      DO 2 I = 1,(NBRJ-1)
        FREE(I) = .FALSE.
    2 CONTINUE
      FREE(NBRJ) = .TRUE.
*
*   2.0  Set quantum numbers of spectator shells.
*
      DO 4 J = 1,NW
        DO 3 K = 1,3
          JBQ1(K,J) = 0
          JBQ2(K,J) = 0
    3   CONTINUE
    4 CONTINUE
*
      DO 8 JJ = 1,NPEEL
        J = JLIST(JJ)
        IF ((J .NE. IA1) .AND. (J .NE. IB1)) THEN
          DO 5 K = 1,3
            JBQ1(K,J) = JJQ1(K,J)
    5     CONTINUE
        ENDIF
        IF ((J .NE. IA2) .AND. (J .NE. IB2)) THEN
          DO 6 K = 1,3
            JBQ2(K,J) = JJQ2(K,J)
    6     CONTINUE
        ENDIF
*
*   2.1  Examine spectator shells for orthogonality
*
        IF ((J .NE. IA1) .AND. (J .NE. IB1) .AND.
     :      (J .NE. IA2) .AND. (J .NE. IB2)) THEN
          DO 7 K = 1,3
            IF (JBQ1(K,J) .NE. JBQ2(K,J) ) GOTO  98
    7     CONTINUE
        ENDIF
    8 CONTINUE
*
*   3.0  Start main calculation
*        Begin with common factors
*
      CONST = OCON (IA1,IB1,IA2,IB2)
      IF (IBUG2 .NE. 0) WRITE (99,307) CONST
*
*   3.1  Set range of tensor index NU
*
      IF (IBUG2 .NE. 0) WRITE (99,302) IA1,IB1,IA2,IB2
      CALL SNRC (IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
      IF (IBUG2 .NE. 0) WRITE (99,303) ND1,ND2,NE1,NE2,IBRD,IBRE
      IF ((IBRD .LT. 0) .AND. (IBRE .LT. 0)) RETURN
      IF ((ND2 .GT. NUMAX) .OR. (NE2 .GT. NUMAX)) THEN
         KK = MAX (ND2,NE2)
         WRITE (*,301) KK
         STOP
      ENDIF
      IF (IBRD .GE. 0) THEN
        DO 10 N = 1,ND2
          DO 9 MU = 1,12
            COND(MU,N) = 0.0D 00
    9     CONTINUE
   10   CONTINUE
      ENDIF
      IF (IBRE .GE. 0) THEN
        DO 12 N = 1,NE2
          DO 11 MU = 1,12
            CONE(MU,N) = 0.0D 00
   11     CONTINUE
   12   CONTINUE
      ENDIF
*
*   3.2  Set parameters of summation over parent
*        (barred) terms in Eq. 2 (loc cit). The array
*        IROWS is formed to point to the list of
*        allowed parents of active shells in the
*        array NTAB
*
      CALL LTAB (IS,NQS,KS,IROWS)
*
      DO 13 I = 1,4
        II = IROWS(I)
        LLS(I) = ITAB(II)
        ILS(I) = JTAB(II)
   13 CONTINUE
*
*   4.0  Sum over all parent terms permitted by
*        angular momentum and seniority selection rules
*
      LLS1 = LLS(1)
      IF (LLS1 .NE. 1) FREE(JA1) = .TRUE.
      LLS2 = LLS(2)
      LLS3 = LLS(3)
      LLS4 = LLS(4)
*
      LS2 = ILS(2)
      DO 29 LB1 = 1,LLS2
        LS2 = LS2+3
        IT1(2) = NTAB(LS2)
        IT12 = IT1(2)
        IT2 = KS(2)
        IT3 = JJQ1(3,IB1)
        IF (ITRIG (IT12,IT2,IT3) .EQ. 0) GOTO 29
        IF (ABS (NTAB(LS2-2)-JJQ1(1,IB1) ) .NE. 1) GOTO 29
*
        LS1 = ILS(1)
        DO 28 LA1 = 1,LLS1
          LS1 = LS1+3
          IT1(1) = NTAB(LS1)
          IT11 = IT1(1)
          IT2 = KS(1)
          IF (IA1 .EQ. IB1) THEN
*
*   Treat IA1 .EQ. IB1 as a special case
*
            IT3 = IT1(2)
            IF (ITRIG (IT11,IT2,IT3) .EQ. 0) GOTO 28
            IF (ABS (NTAB(LS1-2)-NTAB(LS2-2)) .NE. 1) GOTO 28
            IF (LLS2 .NE. 1) FREE(NBRJ-8) = .TRUE.
          ELSE
            IT3 = JJQ1(3,IA1)
            IF (ITRIG (IT11,IT2,IT3) .EQ. 0) GOTO 28
            IF (ABS (NTAB(LS1-2)-JJQ1(1,IA1)) .NE. 1) GOTO 28
            IF (LLS2 .NE. 1) FREE(JB1) = .TRUE.
          ENDIF
*
          LS4 = ILS(4)
          DO 27 LB2 = 1,LLS4
            LS4 = LS4+3
            IT1(4) = NTAB(LS4)
            IT14 = IT1(4)
            IT2 = KS(4)
            IT3 = JJQ2(3,IB2)
            IF (ITRIG(IT14,IT2,IT3) .EQ. 0) GOTO 27
            IF (ABS (NTAB(LS4-2)-JJQ2(1,IB2)) .NE. 1) GOTO 27
*
            LS3 = ILS(3)
            DO 26 LA2 = 1,LLS3
              LS3 = LS3+3
              IT1(3) = NTAB(LS3)
              IT13 = IT1(3)
              IT2 = KS(3)
              IF (IA2 .EQ. IB2) THEN
*
*   TREAT IA2 .EQ. IB2 as a special case
*
                IT3 = IT1(4)
                IF (LLS4 .NE. 1) FREE(NBRJ-6) = .TRUE.
                IF (ITRIG (IT13,IT2,IT3) .EQ. 0) GOTO 26
                IF (ABS (NTAB(LS3-2)-NTAB(LS4-2)) .NE. 1) GOTO 26
              ELSE
                IT3 = JJQ2(3,IA2)
                IF (ITRIG (IT13,IT2,IT3) .EQ. 0) GOTO 26
                IF (ABS (NTAB(LS3-2)-JJQ2(1,IA2)) .NE. 1) GOTO 26
              ENDIF
*
*   At this point the current parent has been completely defined,
*   and its quantum numbers can now be set.  The JTQ arrays must
*   be set if IA1 .EQ. IB1 or IA2 .EQ. IB2. The matrix element should be
*   diagonal in barred quantum numbers.
*
              DO 15 K = 1,3
                JBQ1(K,IA1) = NTAB(LS1+K-3)
                JBQ2(K,IA2) = NTAB(LS3+K-3)
                JTQ1(K) = 0
                IF (IB1 .EQ. IA1) THEN
                  JTQ1(K) = NTAB(LS2+K-3)
                ELSE
                  JBQ1(K,IB1) = NTAB(LS2+K-3)
                ENDIF
                JTQ2(K) = 0
                IF (IB2 .EQ. IA2) THEN
                  JTQ2(K) = NTAB(LS4+K-3)
                ELSE
                  JBQ2(K,IB2) = NTAB(LS4+K-3)
                ENDIF
                DO 14 KK = 1,4
                  IF (JBQ1(K,IS(KK)) .NE. JBQ2(K,IS(KK))) GOTO 26
   14           CONTINUE
   15         CONTINUE
*
*   4.1 Evaluate product of 4 CFPs
*
              CALL MUMDAD (IS,KAPS,PROD)
              IF (ABS (PROD) .LT. EPS) GOTO 26
*
*    4.2  Set arrays for defining the recoupling
*         coefficient
*
              CALL SETJ (IS,JS,KS,NPEEL,KJ23)
*
              IF (ISNJ .EQ. 0) THEN
*
********************* N J G R A F   V E R S I O N **********************
*
*     Set up the arrays and variables for the direct case.
*
                IF (IBRD.GE.0) THEN
                  CALL NJGRAF (RECUP,-1,FAILD)
                  ISNJ = 1
                  IF (.NOT. FAILD) THEN
                  CALL KNJ(JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,JD9,
     :                     KDW,JDDEL,LDDEL,DSUMVR,MDP,
     :                     JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,
     :                     NDB6J,KD6CP,KD7CP,KD8CP,KD9CP,
     :                     JDSUM4,JDSUM5,JDSUM6,INVD6J)
                  ENDIF
                ENDIF
*
*   Set up the arrays and variables for the exchange case.
*
                IF (IBRE .GE. 0) THEN
                  CALL MODJ23
                  CALL NJGRAF (RECUP,-1,FAILE)
                  ISNJ = 2
                  IF (.NOT. FAILE) THEN
                  CALL KNJ(JE6C,JE7C,JE8C,JE9C,JEWC,JE6,JE7,JE8,JE9,
     :                     KEW,JEDEL,LEDEL,ESUMVR,MEP,
     :                     JE6P,JE7P,JE8P,JE9P,JEWORD,NELSUM,NEBJ,
     :                     NEB6J,KE6CP,KE7CP,KE8CP,KE9CP,
     :                     JESUM4,JESUM5,JESUM6,INVE6J)
                  ENDIF
                ENDIF
*
              IF (FAILD .AND. FAILE) GOTO 30
*
              ENDIF
*
*   4.3.1 Summation for direct terms
*
              IF ((IBRD .GE. 0). AND. (.NOT. FAILD)) THEN
*
                IMUD = 4
                IF (IBRD .GT. 1) IMUD = 1
                NCODE = 0
                DO 20 N = 1,ND2
                  NU = ND1+2*(N-1)
                  NUD = NU+NU+1
*
                  IF (NU .NE. 0) THEN
*
                    IF ((ITRIG (KS(1),KS(3),NUD) .NE. 0) .AND.
     :                  (ITRIG (KS(2),KS(4),NUD) .NE. 0)) THEN
*
                    K = NU
                    J1(MJA) = NUD
                    CALL GENSUM (JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,
     :                JD9,KDW,JDDEL,LDDEL,DSUMVR,MDP,
     :                JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,NDB6J,
     :                KD6CP,KD7CP,KD8CP,KD9CP,JDSUM4,JDSUM5,JDSUM6,
     :                INVD6J,X)
                    IF (IBUG2 .NE. 0) WRITE (99,304) NU,K,X
*
                    IF (ABS (X) .GE. EPS) THEN
                      X = X*PROD
                      CALL CXK(S,IS,KAPS,NU,K,IBRD,1)
                      IF (IBUG2 .NE. 0) WRITE (99,305)
     :                                  (S(III),III = 1,IMUD)
                      DO 16 MU = 1,IMUD
                        COND(MU,N) = COND(MU,N)+X*S(MU)
   16                 CONTINUE
                    ENDIF
*
                  ENDIF
*
*   K = NU-1
*
                  IF (IBRD .GT. 1) GOTO 20
*
                  K = NU-1
*
                  IF (NCODE .EQ. N) THEN
                    X=XCODE
                  ELSE
                    ITKMO = NUD-2
                    IF (ITRIG (KS(1),KS(3),ITKMO) .EQ. 0) GOTO 18
                    IF (ITRIG (KS(2),KS(4),ITKMO) .EQ. 0) GOTO 18
                    J1(MJA) = ITKMO
                    CALL GENSUM(JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,
     :                JD9,KDW,JDDEL,LDDEL,DSUMVR,MDP,
     :                JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,NDB6J,
     :                KD6CP,KD7CP,KD8CP,KD9CP,JDSUM4,JDSUM5,JDSUM6,
     :                INVD6J,X)
                  ENDIF
*
                  IF (IBUG2 .NE. 0) WRITE (99,304) NU,K,X
*
                  IF (ABS (X) .GE. EPS) THEN
                    X = X*PROD
                    CALL CXK (S,IS,KAPS,NU,K,IBRD,1)
                    IF (IBUG2 .NE. 0) WRITE (99,305)
     :                                (S(III),III = 1,4)
                    DO 17 MU = 1,4
                      COND (MU,N) = COND(MU,N)+X*S(MU)
   17               CONTINUE
                  ENDIF
*
                ENDIF
*
*   K = NU+1
*
   18           IF ((IBRD .GT. 1) .OR. (N .EQ. ND2)) GOTO 20
*
                NCODE = N+1
                XCODE = 0.0D 00
                ITKMO = NUD+2
*
                IF ((ITRIG(KS(1),KS(3),ITKMO) .NE. 0) .AND.
     :              (ITRIG(KS(2),KS(4),ITKMO) .NE. 0)) THEN
                  K = NU+1
                  J1(MJA) = ITKMO
                  CALL GENSUM (JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,
     :              JD9,KDW,JDDEL,LDDEL,DSUMVR,MDP,
     :              JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,NDB6J,
     :              KD6CP,KD7CP,KD8CP,KD9CP,JDSUM4,JDSUM5,JDSUM6,
     :              INVD6J,X)
                    XCODE = X
                    IF (IBUG2 .NE. 0) WRITE (99,304) NU,K,X
*
                    IF (ABS (X) .GE. EPS) THEN
                      X = X*PROD
                      CALL CXK (S,IS,KAPS,NU,K,IBRD,1)
                      IF (IBUG2 .NE. 0) WRITE (99,305)
     :                                  (S(III),III = 1,12)
                      DO 19 MU = 1,12
                        COND(MU,N) = COND(MU,N)+X*S(MU)
   19                 CONTINUE
                    ENDIF
*
                  ENDIF
*
   20           CONTINUE
*
              ENDIF
*
*   4.3.2 Summation for exchange terms
*
              IF ((IBRE .GE. 0) .AND. (.NOT. FAILE)) THEN
*
                NCODE = 0
*
                DO 25 N = 1,NE2
                  IMUE = 4
                  IF (IBRE .EQ. 2) IMUE = 1
                  IF (IBRE .EQ. 4) IMUE = 3
                  NU = NE1+2*(N-1)
                  NUD = NU+NU+1
*
                  IF (NU .NE. 0) THEN
*
                    IF ((ITRIG(KS(1),KS(4),NUD) .NE. 0) .AND.
     :                  (ITRIG(KS(2),KS(3),NUD) .NE. 0)) THEN
                      K = NU
                      J1(MJA) = NUD
                      CALL GENSUM(JE6C,JE7C,JE8C,JE9C,JEWC,JE6,JE7,JE8,
     :                  JE9,KEW,JEDEL,LEDEL,ESUMVR,MEP,
     :                  JE6P,JE7P,JE8P,JE9P,JEWORD,NELSUM,NEBJ,NEB6J,
     :                  KE6CP,KE7CP,KE8CP,KE9CP,JESUM4,JESUM5,JESUM6,
     :                  INVE6J,X)
                      IF (IBUG2 .NE. 0) WRITE (99,306) NU,K,X
*
                      IF (ABS (X) .GE. EPS) THEN
                        X = X*PROD
                        CALL CXK (S,IS,KAPS,NU,K,IBRE,2)
                        IF (IBUG2 .NE. 0) WRITE (99,305)
     :                                    (S(III),III = 1,IMUE)
                        DO 21 MU = 1,IMUE
                          CONE(MU,N) = CONE(MU,N)+X*S(MU)
   21                   CONTINUE
                      ENDIF
*
                    ENDIF
*
*   K = NU-1
*
                    IF (IBRE .EQ. 2) GOTO 25
*
                    IMUE = 4
                    IF (IBRE .EQ. 4) IMUE = 3
                    K = NU-1
*
                    IF (NCODE .EQ. N) THEN
                      X=XCODE
                    ELSE
                      ITKMO = NUD-2
                      IF (ITRIG (KS(1),KS(4),ITKMO) .EQ. 0) GOTO 23
                      IF (ITRIG (KS(2),KS(3),ITKMO) .EQ. 0) GOTO 23
                      J1(MJA) = ITKMO
                      CALL GENSUM (JE6C,JE7C,JE8C,JE9C,JEWC,JE6,JE7,JE8,
     :                  JE9,KEW,JEDEL,LEDEL,ESUMVR,MEP,
     :                  JE6P,JE7P,JE8P,JE9P,JEWORD,NELSUM,NEBJ,NEB6J,
     :                  KE6CP,KE7CP,KE8CP,KE9CP,JESUM4,JESUM5,JESUM6,
     :                  INVE6J,X)
                    ENDIF
*
                    IF (IBUG2 .NE. 0) WRITE (99,306) NU,K,X
*
                    IF (ABS (X) .GE. EPS) THEN
                      X = X*PROD
                      CALL CXK (S,IS,KAPS,NU,K,IBRE,2)
                      IF (IBUG2 .NE. 0) WRITE (99,305)
     :                                  (S(III),III = 1,IMUE)
                      DO 22 MU = 1,IMUE
                        CONE(MU,N) = CONE(MU,N)+X*S(MU)
   22                 CONTINUE
                    ENDIF
*
                  ENDIF
*
*   K = NU+1
*
   23             IF ((IBRE .EQ. 2) .OR. (N .EQ. NE2)) GOTO 25
*
                  NCODE = N+1
                  XCODE = 0.0D 00
                  IMUE = 12
                  IF (IBRE .EQ. 4) IMUE = 7
                  ITKMO = NUD+2
*
                  IF ((ITRIG (KS(1),KS(4),ITKMO) .NE. 0) .AND.
     :                (ITRIG (KS(2),KS(3),ITKMO) .NE. 0)) THEN
                    K = NU+1
                    J1(MJA) = ITKMO
                      CALL GENSUM (JE6C,JE7C,JE8C,JE9C,JEWC,JE6,JE7,JE8,
     :                  JE9,KEW,JEDEL,LEDEL,ESUMVR,MEP,
     :                  JE6P,JE7P,JE8P,JE9P,JEWORD,NELSUM,NEBJ,NEB6J,
     :                  KE6CP,KE7CP,KE8CP,KE9CP,JESUM4,JESUM5,JESUM6,
     :                  INVE6J,X)
                    XCODE = X
                    IF (IBUG2 .NE. 0) WRITE (99,306) NU,K,X
*
                    IF (ABS (X) .GE. EPS) THEN
                      X = X*PROD
                      CALL CXK (S,IS,KAPS,NU,K,IBRE,2)
                      IF (IBUG2 .NE. 0) WRITE (99,305)
     :                                  (S(III),III = 1,IMUE)
                      DO 24 MU = 1,IMUE
                        CONE (MU,N) = CONE(MU,N)+X*S(MU)
   24                 CONTINUE
                    ENDIF
*
                  ENDIF
*
   25           CONTINUE
*
              ENDIF
*
   26       CONTINUE
   27     CONTINUE
   28   CONTINUE
   29 CONTINUE
*
*   4.4 Insert outside factors
*
   30 IF (IBRD .GE. 0) THEN
        PRODD = KS(1)*KS(4)
        PRODD = CONST/SQRT (PRODD)
        IF ((IA1 .EQ. IB1) .AND. (IA2 .EQ. IB2))
     :     PRODD = 0.5D 00*PRODD
        DO 32 N = 1,ND2
          DO 31 MU = 1,12
            COND(MU,N) = COND(MU,N)*PRODD
   31     CONTINUE
   32   CONTINUE
      ENDIF
*
      IF (IBRE .GE. 0) THEN
        PRODE = KS(1)*KS(3)
        PRODE = -CONST/SQRT(PRODE)
        DO 34 N = 1,NE2
          DO 33 MU = 1,12
            CONE(MU,N) = CONE(MU,N)*PRODE
   33     CONTINUE
   34   CONTINUE
      ENDIF
*
*   5.0 Output results
*
      IF (IBRD .GE. 0) THEN
        DO 35 N = 1,ND2
          NU = ND1+2*(N-1)
          ITYPE = 1
          IF (IBRD .EQ. 2) ITYPE = 3
          IF (IBRD .EQ. 3) ITYPE = 4
          CALL TALK (JA,JB,NU,IA1,IA2,IB1,IB2,ITYPE,COND(1,N))
          IF (IBRD .GT. 1) GOTO 35
          CALL TALK (JA,JB,NU,IA2,IA1,IB2,IB1,ITYPE,COND(2,N))
          CALL TALK (JA,JB,NU,IA1,IA2,IB2,IB1,ITYPE,COND(3,N))
          CALL TALK (JA,JB,NU,IA2,IA1,IB1,IB2,ITYPE,COND(4,N))
          IF (N .EQ. ND2) GOTO 35
          NUP1 = NU+1
          ITYPE = 2
          CALL TALK (JA,JB,NUP1,IA1,IA2,IB1,IB2,ITYPE,COND(5,N))
          CALL TALK (JA,JB,NUP1,IB1,IB2,IA1,IA2,ITYPE,COND(6,N))
          CALL TALK (JA,JB,NUP1,IA2,IA1,IB2,IB1,ITYPE,COND(7,N))
          CALL TALK (JA,JB,NUP1,IB2,IB1,IA2,IA1,ITYPE,COND(8,N))
          CALL TALK (JA,JB,NUP1,IA1,IA2,IB2,IB1,ITYPE,COND(9,N))
          CALL TALK (JA,JB,NUP1,IB2,IB1,IA1,IA2,ITYPE,COND(10,N))
          CALL TALK (JA,JB,NUP1,IA2,IA1,IB1,IB2,ITYPE,COND(11,N))
          CALL TALK (JA,JB,NUP1,IB1,IB2,IA2,IA1,ITYPE,COND(12,N))
   35   CONTINUE
      ENDIF
*
      IF (IBRE .LT. 0) RETURN
*
      DO 36 N = 1,NE2
        NU = NE1+2*(N-1)
        IF (IBRE .NE. 4) THEN
          ITYPE = 1
          IF (IBRE .EQ. 2) ITYPE = 3
          CALL TALK (JA,JB,NU,IA1,IB2,IB1,IA2,ITYPE,CONE(1,N))
          IF (IBRE .EQ. 2) GOTO 36
          CALL TALK (JA,JB,NU,IB2,IA1,IA2,IB1,ITYPE,CONE(2,N))
          CALL TALK (JA,JB,NU,IA1,IB2,IA2,IB1,ITYPE,CONE(3,N))
          CALL TALK (JA,JB,NU,IB2,IA1,IB1,IA2,ITYPE,CONE(4,N))
          IF (N .EQ. NE2) GOTO 36
          NUP1 = NU+1
          ITYPE = 2
          CALL TALK (JA,JB,NUP1,IA1,IB2,IB1,IA2,ITYPE,CONE(5,N))
          CALL TALK (JA,JB,NUP1,IB1,IA2,IA1,IB2,ITYPE,CONE(6,N))
          CALL TALK (JA,JB,NUP1,IB2,IA1,IA2,IB1,ITYPE,CONE(7,N))
          CALL TALK (JA,JB,NUP1,IA2,IB1,IB2,IA1,ITYPE,CONE(8,N))
          CALL TALK (JA,JB,NUP1,IA1,IB2,IA2,IB1,ITYPE,CONE(9,N))
          CALL TALK (JA,JB,NUP1,IA2,IB1,IA1,IB2,ITYPE,CONE(10,N))
          CALL TALK (JA,JB,NUP1,IB2,IA1,IB1,IA2,ITYPE,CONE(11,N))
          CALL TALK (JA,JB,NUP1,IB1,IA2,IB2,IA1,ITYPE,CONE(12,N))
        ELSE
          ITYPE = 5
          CALL TALK (JA,JB,NU,IB1,IA1,IB1,IA1,ITYPE,CONE(1,N))
          CALL TALK (JA,JB,NU,IA1,IB1,IB1,IA1,ITYPE,CONE(2,N))
          CALL TALK (JA,JB,NU,IA1,IB1,IA1,IB1,ITYPE,CONE(3,N))
          IF (N .EQ. NE2) GOTO 36
          NUP1 = NU+1
          ITYPE = 6
          CALL TALK (JA,JB,NUP1,IA1,IB1,IA1,IB1,ITYPE,CONE(4,N))
          CALL TALK (JA,JB,NUP1,IB1,IA1,IB1,IA1,ITYPE,CONE(5,N))
          CALL TALK (JA,JB,NUP1,IA1,IB1,IB1,IA1,ITYPE,CONE(6,N))
          CALL TALK (JA,JB,NUP1,IB1,IA1,IA1,IB1,ITYPE,CONE(7,N))
        ENDIF
   36 CONTINUE
*
      RETURN
*
*   6.0 Fault diagnostic prints
*
   98 IF (IBUG2 .NE. 0) WRITE (99,300)
      RETURN
*
  300 FORMAT ('BREIT: Spectator quantum numbers not diagonal for',
     :        ' non-interacting shells')
  301 FORMAT ('BREIT: Increase second dimension of arrays',
     :        ' COND(MU,N) and CONE(MU,N) to the new value of NUMAX,'
     :       /' (at least ',1I3,').')
  302 FORMAT ('BREIT: Subshells ',4I5)
  303 FORMAT ('  ND1 ND2 NE1 NE2 IBRD IBRE',6I5)
  304 FORMAT ('  Direct NU K recoupling coef ',2I5,1P,D20.9)
  305 FORMAT (' S',1P,8D15.7)
  306 FORMAT ('  Exchange NU K recoupling coef ',2I5,1P,D20.9)
  307 FORMAT ('  Statistical factor ',1P,D20.9)
*
      END
************************************************************************
*                                                                      *
      FUNCTION CLRX (KAPPAA,K,KAPPAB)
*                                                                      *
*   The value of CLRX is the 3-j symbol:                               *
*                                                                      *
*                    ( JA        K        JB  )                        *
*                    ( 1/2       0       -1/2 )                        *
*                                                                      *
*   The  K'S are kappa angular quantum numbers. The formula is taken   *
*   from D M Brink and G R Satchler, <Angular Momentum>, second edi-   *
*   tion (Oxford: Clarendon press, 1968), p 138.   The logarithms of   *
*   the first  MFACT  factorials must be available in  COMMON/FACTS/   *
*   for this program to function correctly. Note that  N!  is stored   *
*   in FACT(N+1)                                                       *
*                                                                      *
*   No subroutines called.                                             *
*                                                                      *
*   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      PARAMETER (MFACT = 500)
*
      COMMON/FACTS/GAM(MFACT)
*
*   Determine the absolute values of the kappas
*
      KA = ABS (KAPPAA)
      KB = ABS (KAPPAB)
*
*   Perform the triangularity check
*
      IF ((ABS (KA-KB) .LE. K) .AND. (KA+KB-1 .GE. K)) THEN
*
*   Triangularity satisfied; compute the 3j coefficient
*
*   Begin with the logarithm of the square of the leading term
*
         EXPTRM = -LOG (DBLE (KA*KB))
*
*   Compute the logarithm of the square root of the leading term
*   and the factorial part that doesn't depend on the parity of
*   KA+KB+K (the delta factor)
*
         KAPKB = KA+KB
         KABKP = KAPKB+K
         KAMKB = KA-KB
         KBMKA = KB-KA
         EXPTRM = 0.5D 00
     :           *(EXPTRM+GAM(KAPKB-K  )+GAM(KAMKB+K+1)
     :                   +GAM(KBMKA+K+1)-GAM(KABKP  +1) )
*
*   The remainder depends on the parity of KA+KB+K
*
         IF (MOD (KABKP,2) .EQ. 0) THEN
*
*   Computation for even parity case
*
*   Include the phase factor: a minus sign if necessary
*
            IF (MOD (3*KABKP/2,2) .EQ. 0) THEN
               CLRX =  1.0D 00
            ELSE
               CLRX = -1.0D 00
            ENDIF
*
*   Include the contribution from the factorials
*
            EXPTRM = EXPTRM+GAM((KABKP  +2)/2)-GAM((KAPKB-K  )/2)
     :                     -GAM((KAMKB+K+2)/2)-GAM((KBMKA+K+2)/2)
*
         ELSE
*
*   Computation for odd parity case
*
*   Include the phase factor: a minus sign if necessary
*
            IF (MOD ((3*KABKP-1)/2,2) .EQ. 0) THEN
               CLRX =  1.0D 00
            ELSE
               CLRX = -1.0D 00
            ENDIF
*
*   Include the contribution from the factorials
*
            EXPTRM = EXPTRM+GAM((KABKP  +1)/2)-GAM((KAPKB-K+1)/2)
     :                     -GAM((KAMKB+K+1)/2)-GAM((KBMKA+K+1)/2)
*
         ENDIF
*
*   Final assembly
*
         CLRX = CLRX*EXP (EXPTRM)
*
      ELSE
*
*   Triangularity violated; set the coefficient to zero
*
         CLRX = 0.0D 00
*
      ENDIF
*
      RETURN
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE COR (JA,JB,JA1,JB1,JA2,JB2)
*                                                                      *
*   Computes  the  MCP  coefficients.  Equation numbers are those of   *
*   Computer Phys Commun 5 (1973) 263                                  *
*                                                                      *
*   Call(s) to: [LIB92]: CRE, ITRIG, LTAB, MODJ23, MUMDAD, OCON,       *
*                        SETJ, SKRC, SPEAK, KNJ.                       *
*               [NJGRAF]: NJGRAF, GENSUM.                              *
*                                                                      *
*                                           Last update: 02 Nov 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INTEGER PNTRIQ
*
      PARAMETER (
     :   MANGM=60,M3MNGM=3*MANGM,MANGMP=2*(MANGM/3),
     :   MTRIAD=12,
     :   M6J=20,MSUM=10)
*
      PARAMETER (EPS=1.0D-10)
*
      LOGICAL FREE,DSUMVR,ESUMVR,FAILD,FAILE
*
      PARAMETER (IDIM = 11)
      DIMENSION COND(IDIM),CONE(IDIM)
*
      DIMENSION KAPS(4),KS(4),NQS(4),ILS(4),LLS(4),IROWS(4)
      DIMENSION IS(4),JS(4)
*
      DIMENSION JD6(M3MNGM),JD7(M3MNGM),JD8(M3MNGM),
     : JD9(MANGMP),KDW(6,M6J),LDDEL(M6J,2),DSUMVR(MANGM)
      DIMENSION JD6P(MANGMP),JD7P(MANGMP),JD8P(MANGMP),JD9P(MANGMP),
     : JDWORD(6,M6J),
     : NDBJ(MSUM),NDB6J(MSUM),KD6CP(MSUM),KD7CP(MSUM),KD8CP(MSUM),
     : KD9CP(MSUM),JDSUM6(MTRIAD),JDSUM4(MTRIAD,M6J),JDSUM5(MTRIAD,M6J),
     : INVD6J(M6J)
      DIMENSION JE6(M3MNGM),JE7(M3MNGM),JE8(M3MNGM),
     : JE9(MANGMP),KEW(6,M6J),LEDEL(M6J,2),ESUMVR(MANGM)
      DIMENSION JE6P(MANGMP),JE7P(MANGMP),JE8P(MANGMP),JE9P(MANGMP),
     : JEWORD(6,M6J),
     : NEBJ(MSUM),NEB6J(MSUM),KE6CP(MSUM),KE7CP(MSUM),KE8CP(MSUM),
     : KE9CP(MSUM),JESUM6(MTRIAD),JESUM4(MTRIAD,M6J),JESUM5(MTRIAD,M6J),
     : INVE6J(M6J)
*
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
     :      /COUPLE/MJA,NJA,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /L1/JBQ1(3,149),JBQ2(3,149),JTQ1(3),JTQ2(3)
     :      /M1/NQ1(149),NQ2(149)
     :      /M2/JJQ1(3,149),JJQ2(3,149)
     :      /M3/JLIST(149),KLIST(149),NPEEL,NCORE
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(149),NAK(149)
     :      /TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
*
*   1.0  Initialize pointers and flags and set any
*        tables required.
*
*        In this segment, the array IS points to the
*        full list of orbitals, the array JS to the
*        array JLIST of peel orbital pointers.
*
*   1.1  Initialization
*
      JS(1) = JA1
      JS(2) = JB1
      JS(3) = JA2
      JS(4) = JB2
      DO 1 I = 1,4
        IS(I) = JLIST(JS(I))
        KAPS(I) = 2*NAK(IS(I))
        KS(I) = ABS (KAPS(I))
    1 CONTINUE
      IA1 = IS(1)
      IB1 = IS(2)
      IA2 = IS(3)
      IB2 = IS(4)
*
      KJ23 = 0
      ISNJ = 0
      FAILD = .FALSE.
      FAILE = .FALSE.
*
*   Initialize arrays
*
      DO 3 J = 1,NW
         DO 2 K = 1,3
            JBQ1(K,J) = 0
            JBQ2(K,J) = 0
    2    CONTINUE
    3 CONTINUE
*
      NBRJ = 3*NPEEL+7
      DO 4 I = 1,NBRJ
         FREE(I) = .FALSE.
    4 CONTINUE
*
*   2.0 Set tables of quantum numbers of spectator
*       shells.
*
      DO 8 JJ = 1,NPEEL
        J = JLIST(JJ)
        IF ((J .NE. IA1) .AND. (J .NE. IB1)) THEN
          DO 5 K = 1,3
            JBQ1(K,J) = JJQ1(K,J)
    5     CONTINUE
        ENDIF
        IF ((J .NE. IA2) .AND. (J .NE. IB2)) THEN
          DO 6 K = 1,3
            JBQ2(K,J) = JJQ2(K,J)
    6     CONTINUE
        ENDIF
*
*   2.1 Examine quantum numbers of spectator
*       shells for orthogonality and
*       exit if found.
*
        IF ((J .NE. IA1) .AND. (J .NE. IB1) .AND.
     :      (J .NE. IA2) .AND. (J .NE. IB2)) THEN
          DO 7 K = 1,3
            IF (JBQ1(K,J) .NE. JBQ2(K,J)) THEN
               IF (IBUG2 .NE. 0) WRITE (99,300)
               GOTO 41
            ENDIF
    7     CONTINUE
        ENDIF
    8 CONTINUE
*
*   3.1  Set range of the parameter k for Coulomb
*        integrals.
*        Terminate run if buffer store dimension
*        IDIM is too small.
*
      CALL SKRC (IS,KAPS,KS,KD1,KD2,KE1,KE2)
      IF ((KD2 .EQ. 0) .AND. (KE2 .EQ. 0)) GOTO 41
      IF ((KD2. GT. IDIM) .OR. (KE2 .GT. IDIM)) THEN
         KK = MAX (KE2,KD2)
         WRITE (*,301) KK
         STOP
      ENDIF
*
      IF (KD2 .NE. 0) THEN
         DO 9 K = 1,KD2
            COND(K) = ZERO
    9    CONTINUE
      ENDIF
      IF (KE2 .NE. 0) THEN
         DO 10 K = 1,KE2
            CONE(K) = ZERO
   10    CONTINUE
      ENDIF
*
      NQS(1) = NQ1(IA1)
      NQS(2) = NQ1(IB1)
      NQS(3) = NQ2(IA2)
      NQS(4) = NQ2(IB2)
*
*   3.3  Set parameters of summation over parent
*        (barred) terms in Eq. (5). The array IROWS
*        is formed to point at the list of allowed
*        parents of active shells in the array
*        NTAB.
*
      CALL LTAB (IS,NQS,KS,IROWS)
*
      DO 17 I = 1,4
         II = IROWS(I)
         LLS(I) = ITAB(II)
         ILS(I) = JTAB(II)
   17 CONTINUE
*
*   4.0  Sum contributions over all parent terms
*        permitted by angular momentum and seniority
*        selection rules.
*
      LLS1 = LLS(1)
      IF (LLS1 .NE. 1) FREE (JA1) = .TRUE.
      LLS2 = LLS(2)
      LLS3 = LLS(3)
      LLS4 = LLS(4)
*
      LS2 = ILS(2)
      DO 38 LB1 = 1,LLS2
        LS2 = LS2+3
        IT12 = NTAB(LS2)
        IT2 = KS(2)
        IT3 = JJQ1(3,IB1)
        IF (ITRIG (IT12,IT2,IT3) .EQ. 0) GOTO 38
        IF (ABS (NTAB(LS2-2)-JJQ1(1,IB1)) .NE. 1) GOTO 38
*
        LS1 = ILS(1)
        DO 37 LA1 = 1,LLS1
          LS1 = LS1+3
          IT11 = NTAB(LS1)
          IT2 = KS(1)
*
          IF (IA1 .EQ. IB1) THEN
*
*   Treat IA1 .EQ. IB1 as special case.
*
            IT3 = IT12
            IF (ITRIG (IT11,IT2,IT3) .EQ. 0) GOTO 37
            IF (ABS (NTAB(LS1-2)-NTAB(LS2-2)) .NE. 1) GOTO 37
            IF (LLS2 .NE. 1) FREE(NBRJ-8) = .TRUE.
          ELSE
            IT3 = JJQ1(3,IA1)
            IF (ITRIG (IT11,IT2,IT3) .EQ. 0) GOTO 37
            IF (ABS (NTAB(LS1-2)-JJQ1(1,IA1)) .NE. 1) GOTO 37
            IF (LLS2 .NE. 1) FREE(JB1) = .TRUE.
          ENDIF
*
          LS4 = ILS(4)
          DO 36 LB2 = 1,LLS4
            LS4 = LS4+3
            IT14 = NTAB(LS4)
            IT2 = KS(4)
            IT3 = JJQ2(3,IB2)
            IF (ITRIG (IT14,IT2,IT3) .EQ. 0) GOTO 36
            IF (ABS (NTAB(LS4-2)-JJQ2(1,IB2)) .NE. 1) GOTO 36
*
            LS3 = ILS(3)
            DO 35 LA2 = 1,LLS3
              LS3 = LS3+3
              IT13 = NTAB(LS3)
              IT2 = KS(3)
*
              IF (IA2 .EQ. IB2) THEN
*
*   Treat IA2 .EQ. IB2 as special case.
*
                IT3 = IT14
                IF (LLS4 .NE. 1) FREE(NBRJ-6) = .TRUE.
                IF (ITRIG (IT13,IT2,IT3) .EQ. 0) GOTO 35
                IF (ABS (NTAB(LS3-2)-NTAB(LS4-2)) .NE. 1) GOTO 35
              ELSE
                IT3 = JJQ2(3,IA2)
                IF (ITRIG (IT13,IT2,IT3) .EQ. 0) GOTO 35
                IF (ABS (NTAB(LS3-2)-JJQ2(1,IA2)) .NE. 1) GOTO 35
              ENDIF
*
*   At this point the current parent has been completely defined,
*   and its quantum numbers can now be set. The JTQ arrays must
*   be set if IA1 .EQ. IB1 or IA2 .EQ. IB2. The matrix element
*   should be diagonal in barred quantum numbers.
*
              DO 27 K = 1,3
                JBQ1(K,IA1) = NTAB(LS1+K-3)
                JBQ2(K,IA2) = NTAB(LS3+K-3)
                JTQ1(K) = 0
                IF (IB1 .EQ. IA1) THEN
                  JTQ1(K) = NTAB(LS2+K-3)
                ELSE
                  JBQ1(K,IB1) = NTAB(LS2+K-3)
                ENDIF
                JTQ2(K) = 0
                IF (IB2 .EQ. IA2) THEN
                  JTQ2(K) = NTAB(LS4+K-3)
                ELSE
                  JBQ2(K,IB2) = NTAB(LS4+K-3)
                ENDIF
                DO 26 KK = 1,4
                  IF (JBQ1(K,IS(KK)) .NE. JBQ2(K,IS(KK))) GOTO 35
   26           CONTINUE
   27         CONTINUE
*
*   4.1  Evaluate product of 4 CFPs, Eq. (5).
*
              CALL MUMDAD (IS,KAPS,PROD)
              IF (ABS (PROD) .LT. EPS) GOTO 35
*
*   4.2  Set arrays for defining the recoupling
*        coefficient.
*
              CALL SETJ (IS,JS,KS,NPEEL,KJ23)
*
              IF (ISNJ .EQ. 0) THEN
*
********************* N J G R A F   V e r s i o n **********************
*
*   Set up the arrays and variables for the direct case.
*   J1(NBRJ) ( = J1(MJA) ) is set to (2*KD1+1) so that NJGRAF is
*   called correctly.
*
                IF (KD2 .NE. 0) THEN
                  IF (KD2 .GT. 1) FREE(NBRJ) = .TRUE.
                  J1(NBRJ) = KD1+KD1+1
                  CALL NJGRAF (RECUP,-1,FAILD)
                  ISNJ = 1
                  IF (.NOT. FAILD) THEN
                  CALL KNJ (JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,JD9,
     :                     KDW,JDDEL,LDDEL,DSUMVR,MDP,
     :                     JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,
     :                     NDB6J,KD6CP,KD7CP,KD8CP,KD9CP,
     :                     JDSUM4,JDSUM5,JDSUM6,INVD6J)
                  ENDIF
                ENDIF
*
*
*   Set up the arrays and variables for the exchange case.
*   J1(NBRJ) ( = J1(MJA) ) is set to (2*KE1+1) so that NJGRAF is
*   called correctly.
*
                IF (KE2 .NE. 0) THEN
                  CALL MODJ23
                  FREE(NBRJ) = .FALSE.
                  IF (KE2 .GT. 1) FREE(NBRJ) = .TRUE.
                  J1(NBRJ) = KE1+KE1+1
                  CALL NJGRAF (RECUP,-1,FAILE)
                  ISNJ = 2
                  IF (.NOT. FAILE) THEN
                  CALL KNJ (JE6C,JE7C,JE8C,JE9C,JEWC,JE6,JE7,JE8,JE9,
     :                     KEW,JEDEL,LEDEL,ESUMVR,MEP,
     :                     JE6P,JE7P,JE8P,JE9P,JEWORD,NELSUM,NEBJ,
     :                     NEB6J,KE6CP,KE7CP,KE8CP,KE9CP,
     :                     JESUM4,JESUM5,JESUM6,INVE6J)
                  ENDIF
                ENDIF
              ENDIF
*
*   4.3  Calculate AD, Eq. (6),
*        without the phase factor.
*
              IF ((KD2 .NE. 0) .AND. (.NOT. FAILD)) THEN
                KK = KD1-2
                DO 30 K = 1,KD2
                  KK = KK+2
                  J1(MJA) = KK+KK+1
                  CALL GENSUM (JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,JD9,
     +            KDW,JDDEL,LDDEL,DSUMVR,MDP,
     +            JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,NDB6J,
     +            KD6CP,KD7CP,KD8CP,KD9CP,JDSUM4,JDSUM5,JDSUM6,INVD6J,
     +            X)
                  IF (ABS (X) .GE. EPS) COND(K)=COND(K)+X*PROD
   30           CONTINUE
              ENDIF
*
*   4.4  Calculate AE, Eq. (6),
*        without the phase factor.
*
              IF ((KE2 .NE. 0) .AND. (.NOT. FAILE)) THEN
                KK = KE1-2
                DO 34 K = 1,KE2
                  KK = KK+2
                  J1(MJA) = KK+KK+1
                  CALL GENSUM (JE6C,JE7C,JE8C,JE9C,JEWC,JE6,JE7,JE8,JE9,
     +            KEW,JEDEL,LEDEL,ESUMVR,MEP,
     +            JE6P,JE7P,JE8P,JE9P,JEWORD,NELSUM,NEBJ,NEB6J,
     +            KE6CP,KE7CP,KE8CP,KE9CP,JESUM4,JESUM5,JESUM6,INVE6J,
     +            Y)
                  IF (ABS(Y) .GE. EPS) CONE(K) = CONE(K)+Y*PROD
   34           CONTINUE
              ENDIF
*
              IF (FAILD .AND. FAILE) GOTO 500
*
   35       CONTINUE
   36     CONTINUE
   37   CONTINUE
   38 CONTINUE
*
*   4.5  Insert factors independent of barred
*        quantum numbers.
*        Output results
*
*        Begin with common statistical factors, Eq. (5).
*
  500 CONST = OCON (IA1,IB1,IA2,IB2)
*
      KAP1 = NAK(IA1)
      KAP2 = NAK(IB1)
      KAP3 = NAK(IA2)
      KAP4 = NAK(IB2)
*
*   4.6  Compute products of reduced matrix
*        elements, Eq. (7).
*        CRED for direct terms
*        CREE for exchange terms
*
      IF (KD2 .NE. 0) THEN
         PRODD = CONST/SQRT (DBLE (KS(1)*KS(4)))
         IF (MOD (KD1,2) .NE. 0) PRODD = -PRODD
         IF ((IA1 .EQ. IB1) .AND. (IA2 .EQ. IB2)) PRODD = PRODD*HALF
         KK = KD1-2
         DO 39 K = 1,KD2
            KK = KK+2
            CRED = CRE (KAP1,KK,KAP3)*CRE (KAP2,KK,KAP4)
            X = PRODD*COND(K)*CRED
            IF (ABS (X) .GE. EPS)
     :         CALL SPEAK (JA,JB,IA1,IB1,IA2,IB2,KK,X)
   39    CONTINUE
      ENDIF
*
      IF (KE2 .NE. 0) THEN
         PRODE = CONST/SQRT(DBLE (KS(1)*KS(3)))
         IF (MOD (KE1,2) .NE. 0) PRODE = -PRODE
         PRODE = -PRODE
         KK = KE1-2
         DO 40 K = 1,KE2
            KK = KK+2
            CREE = CRE (KAP1,KK,KAP4)*CRE (KAP2,KK,KAP3)
            X = PRODE*CONE(K)*CREE
            IF (ABS (X) .GE. EPS)
     :         CALL SPEAK (JA,JB,IA1,IB1,IB2,IA2,KK,X)
   40    CONTINUE
      ENDIF
*
   41 RETURN
*
  300 FORMAT('COR: Spectator quantum numbers not diagonal for',
     :       ' non-interacting shells')
  301 FORMAT('COR: Dimension error: reset PARAMETER IDIM to at least ',
     :       1I2,' and recompile.')
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE CORD (JA,JB,JA1,IPCA,JB1)
*                                                                      *
*   Computes the MCP coefficients for contributions involving closed   *
*   shells.  The  standard formulae are given in I P Grant, Advances   *
*   in Physics  19 (1970) 747, Eq. (8.33).  In this segment JA1, JB1   *
*   point to the JLIST array, IA1, IB1 to the full list of orbitals.   *
*                                                                      *
*   Call(s) to: CLRX, SPEAK.                                           *
*                                                                      *
*                                           Last update: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      PARAMETER (EPS = 1.0D-10)
*
      COMMON/ORB4/NP(149),NAK(149)
     :      /M1/NQ1(149),NQ2(149)
     :      /M3/JLIST(149),KLIST(149),NPEEL,NCORE
*
*   Set quantum numbers required.
*
      IF (IPCA .EQ. 2) THEN
         IA1 = KLIST(JA1)
      ELSE
         IA1 = JLIST(JA1)
      ENDIF
      IB1 = KLIST(JB1)
*
*   Force IA1 to be greater than IB1
*
      IF (IA1 .GT. IB1) THEN
         NS = IA1
         IA1 = IB1
         IB1 = NS
      ENDIF
*
      KAP1 = NAK(IA1)
      J1 = IABS (KAP1)
      NQS1 = NQ1(IA1)
*
      IF (IA1 .EQ. IB1) THEN
*
*   Case when IA1 .EQ. IB1
*
         X = DBLE (NQS1*(NQS1-1)/2)
         CALL SPEAK (JA,JB,IA1,IB1,IA1,IB1,0,X)
         NUMAX = J1+J1-2
         IF (NUMAX  .LE.  0) RETURN
         CONST = DBLE (NQS1*NQS1/2)
         DO 1 NU = 2,NUMAX,2
            GAM = CLRX (KAP1,NU,KAP1)
            X = -CONST*GAM*GAM
            IF (ABS (X) .GE. EPS)
     :         CALL SPEAK (JA,JB,IA1,IB1,IA1,IB1,NU,X)
    1    CONTINUE
*
*   Case when IA1 .NE. IB1
*
      ELSE
*
         KAP2 = NAK(IB1)
         J2 = ABS (KAP2)
         NQS2 = NQ1(IB1)
         CONST = DBLE (NQS1*NQS2)
         CALL SPEAK (JA,JB,IA1,IB1,IA1,IB1,0,CONST)
         NUMIN = ABS (J1-J2)
         NUMAX = J1+J2-1
         IF (KAP1*KAP2 .LT. 0) NUMIN = NUMIN+1
         DO 2 NU = NUMIN,NUMAX,2
            GAM = CLRX (KAP1,NU,KAP2)
            X = -CONST*GAM*GAM
            IF (ABS (X) .GE. EPS)
     :         CALL SPEAK (JA,JB,IA1,IB1,IB1,IA1,NU,X)
    2    CONTINUE
*
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      FUNCTION CRE (KAP1,K,KAP2)
*                                                                      *
*   Computes the relativistic reduced matrix element                   *
*                                                                      *
*                         (j1 || C(K) || j2),                          *
*                                                                      *
*   Eq. (5.15) of I P Grant, Advances in Physics 19 (1970) 762. KAP1,  *
*   KAP2 are the kappa values corresponding to j1, j2.  The triangle   *
*   conditions are tested by the routine CLRX.                         *
*                                                                      *
*   Call(s) to: [LIB92] CLRX.                                          *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      K1 = ABS (KAP1)
      DK1K2 = DBLE (4*K1*IABS (KAP2))
      CRE = SQRT (DK1K2)*CLRX (KAP1,K,KAP2)
      IF (MOD (K1,2) .EQ. 1) CRE  = -CRE
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE CXK (S,IS,KAPS,NU,K,IBR,IEX)
*                                                                      *
*   Computes  the  coefficients of radial integrals in the expansion   *
*   of the effective interaction strength: X(K,IA1,IB1,IA2,IB2).       *
*                                                                      *
*   Input variables:                                                   *
*                                                                      *
*      IS  : Orbital labels                                            *
*      KAPS: Values of 2*kappa                                         *
*      NU  : Order of radial integral                                  *
*      K   : Index of tensor operator                                  *
*      IEX : 1 for direct, 2 for exchange terms                        *
*      IBR : Classifies type of radial integral.There are 4 distinct   *
*            cases:                                                    *
*            IBR = 1 A. All states distinct                            *
*                    B. ((IA .EQ. IB) .AND. (IC .NE. ID)), or          *
*                       ((IA .NE. IB) .AND. (IC .EQ. ID))              *
*                    These give 12 distinct  radial  integrals, with   *
*                    values of K and NU limited only by angular mom-   *
*                    entum and parity                                  *
*            IBR = 2 ((IA .EQ. IC) .AND. (IB .NE. ID)) or              *
*                    ((IA .NE. IC) .AND. (IB .EQ. ID))                 *
*                    This case gives one non-zero integral when K =    *
*                    NU is ODD                                         *
*            IBR = 3 ((IA .EQ. IC) .AND. (IB .EQ. ID)) AND             *
*                    (IA .NE. IB)                                      *
*                    Integrals of magnetic F-type when K = NU is odd   *
*            IBR = 4 ((IA .EQ. ID) .AND. (IB .EQ. IC)) gives 3  mag-   *
*                    netic G-type integrals and  four  H-TYPE  inte-   *
*                    grals                                             *
*                                                                      *
*   Output:                                                            *
*                                                                      *
*      S   : Coefficients S(MU) MU = 1,12                              *
*                                                                      *
*                                                                      *
*   Call(s) to: [LIB92] CRE.                                           *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      DIMENSION IS(4),KAPS(4),S(12)
*
*   1.0  Initialization
*
      DO 1 MU = 1,12
         S(MU) = 0.0D 00
    1 CONTINUE
*
      IA = IS(1)
      IB = IS(2)
      IC = IS(3)
      ID = IS(4)
      KA = KAPS(1)/2
      KB = KAPS(2)/2
      KC = KAPS(3)/2
      KD = KAPS(4)/2
      IF (IEX .NE. 2) GOTO 2
      KK = KD
      IK = ID
      KD = KC
      ID = IC
      KC = KK
      IC = IK
    2 GOTO (3,8,11,12),IBR
      GOTO 17
*
*   2.0  IBR = 1 --- The general case
*
    3 CONTINUE
      IF (NU-K) 7,4,6
*
*   2.1  NU = K .GT. 0
*
    4 CONTINUE
      S(1) = -(KA+KC)*(KD+KB)
      IF (K .EQ. 0) GOTO 16
      D = K*(K+1)
      H = CRE (KA,K,KC)*CRE(KB,K,KD)
      IF (MOD (K,2) .NE. 0) H = -H
      S(1) = S(1)*H/D
      DO 5 MU = 2,4
         S(MU) = S(1)
    5 CONTINUE
      RETURN
*
*   2.2  NU = K+1
*
    6 CONTINUE
      DK1 = KC-KA
      DK2 = KD-KB
      FK = K
      GK = K+1
      G1 = DK1-GK
      G2 = DK1+GK
      G3 = DK2-GK
      G4 = DK2+GK
      KK = K+K+1
      H = CRE (KA,K,KC)*CRE(KB,K,KD)
      IF (MOD (K,2) .NE. 0) H = -H
      A = H*FK/GK/DBLE (KK*(KK+2))
      S(1) = A*G1*G3
      S(2) = A*G2*G4
      S(3) = A*G1*G4
      S(4) = A*G2*G3
      RETURN
*
*   2.2  NU = K-1
*
    7 CONTINUE
      DK1 = KC-KA
      DK2 = KD-KB
      FK = K
      GK = K+1
      F1 = DK1-FK
      F2 = DK1+FK
      F3 = DK2-FK
      F4 = DK2+FK
      G1 = DK1-GK
      G2 = DK1+GK
      G3 = DK2-GK
      G4 = DK2+GK
      KK = K+K+1
      H = CRE (KA,K,KC)*CRE(KB,K,KD)
      IF (MOD (K,2) .NE. 0) H = -H
      A = H*GK/FK/DBLE (KK*(KK-2))
      S(1) = A*F2*F4
      S(2) = A*F1*F3
      S(3) = A*F2*F3
      S(4) = A*F1*F4
      B = H/DBLE (KK*KK)
      S(5) = B*F2*G3
      S(6) = B*F4*G1
      S(7) = B*F1*G4
      S(8) = B*F3*G2
      S(9) = B*F2*G4
      S(10) = B*F3*G1
      S(11) = B*F1*G3
      S(12) = B*F4*G2
      RETURN
*
*   3.0  IBR = 2  Degenerate case: only one non-zero R-integral
*
    8 CONTINUE
      IF ((IA .EQ. IC) .AND. (IB .NE. ID)) GOTO 10
      IF ((IA .NE. IC) .AND. (IB .EQ. ID)) GOTO 9
      GOTO 17
*
    9 IK = IB
      IB = IA
      IA = IK
      IK = ID
      ID = IC
      IC = IK
*
      KK = KB
      KB = KA
      KA = KK
      KK = KD
      KD = KC
      KC = KK
*
   10 IF (MOD (K,2) .NE. 1) RETURN
      DK = K*(K+1)
      H = CRE (KA,K,KC)*CRE(KB,K,KD)/DK
      S(1) = H*DBLE (4*KA*(KB+KD))
      RETURN
*
*   4.0  IBR = 3. Direct magnetic F-integrals
*
   11 CONTINUE
      IF ((IA .NE. IC) .OR. (IB .NE. ID)) GOTO 17
      IF (MOD (K,2) .NE. 1) RETURN
      DK = K*(K+1)
      H = CRE(KA,K,KA)*CRE(KB,K,KB)/DK
      S(1) = H*DBLE (16*KA*KB)
      RETURN
*
*   5.0   IBR = 4. Exchange magnetic G- and H-integrals
*
   12 CONTINUE
      IF ((IA .NE. ID) .OR. (IB .NE. IC) )GOTO 17
      IF (NU-K) 15,13,14
*
*   5.1  NU = K
*
   13 CONTINUE
      S(1) = DBLE (KA+KB)*CRE(KA,K,KB)
      IP = ABS (KA)-ABS (KB)+K+1
      S(1) = S(1)*S(1)/DBLE(K*(K+1))
      IF (MOD (IP,2) .NE. 0) S(1) = -S(1)
      S(3) = S(1)
      S(2) = S(1)+S(1)
      RETURN
*
*   5.2  NU = K+1
*
   14 CONTINUE
      DK = KB-KA
      GK = K+1
      FK = K
      G1 = DK+GK
      G2 = DK-GK
      KK = K+K+1
      H = CRE (KA,K,KB)**2
      IF (KA*KB .LT. 0) H = -H
      A = H*FK/GK/DBLE(KK*(KK+2))
      S(1) = -A*G1*G1
      S(2) = -2.0D 00*A*G1*G2
      S(3) = -A*G2*G2
      RETURN
*
*   5.3  NU = K-1
*
   15 CONTINUE
      DK = KB-KA
      FK = K
      GK = K+1
      F1 = DK+FK
      F2 = DK-FK
      G1 = DK+GK
      G2 = DK-GK
      KK = K+K+1
      H = CRE (KA,K,KB)**2
      IF (KA*KB .LT. 0) H = -H
      A = H*GK/FK/DBLE (KK*(KK-2))
      S(1) = -A*F2*F2
      S(2) = -2.0D 00*A*F1*F2
      S(3) = -A*F1*F1
      B = H/DBLE (KK*KK)
      B = B+B
      S(4) = -B*F1*G2
*     S(5) = S(4)
      S(5) = -B*F2*G1
      S(6) = -B*F1*G1
      S(7) = -B*F2*G2
      RETURN
*
*   6.0  Special cases and errors
*
*   Illegal zero value of K in Type 1
*
   16 WRITE (*,300) IS(1),IS(2),IS(3),IS(4),NU,IBR,IEX
      STOP
*
*   Illegal combination of states in Type 3 or 4
*
   17 WRITE (*,301) IBR,IS(1),IS(2),IS(3),IS(4),NU,K,IEX
      STOP
*
  300 FORMAT ('CXK: Illegal value K = 0 -'
     :       /1X,4I3,2X,I3,2X,2I2)
  301 FORMAT ('CXK: Type ',I2,'-'
     :       /1X,I2,3X,4I3,2X,2I3,2X,I2)
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE FIXJ (JA1,JA2,KA,IS,KS,NS,KJ23)
*                                                                      *
*   Sets up the arrays J1, J2, J3 required by the recoupling package   *
*   NJSYM.                                                             *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      PARAMETER (MANGM=60,MTRIAD=12)
*
      LOGICAL FREE
*
      DIMENSION IS(2),KS(2)
*
      COMMON/M0/JJC1(149),JJC2(149)
     :      /M2/JJQ1(3,149),JJQ2(3,149)
     :      /M3/JLIST(149),KLIST(149),NPEEL,NCORE
     :      /COUPLE/MJA,NJA,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
     :      /L1/JBQ1(3,149),JBQ2(3,149),JTQ1(3),JTQ2(3)
*
*   Set up the J2 and J3 arrays
*
      NM1 = NS-1
      IF (KJ23 .EQ. 1) GOTO 5
      NS1 = NS+1
      N2 = NS+NS
      N3 = N2+NS
*
      J2(1,1) = N3+2
      J2(1,2) = N3+3
      J2(1,3) = N3+1
      J2(2,1) = JA1
      J2(2,2) = N3+1
      J2(2,3) = N3-1
*
      J3(1,1) = JA2
      J3(1,2) = N3+2
      J3(1,3) = N3
*
      IF (NS .EQ. 1) GOTO 3
*
      DO 1 JW = 1,NM1
         JJ = JW+2
         J2(JJ,1) = NS+JW-1
         J2(JJ,2) = JW+1
         J2(JJ,3) = NS+JW
*
         JK = JW+1
         J3(JK,1) = N2+JW-2
         J3(JK,2) = JW+1
         J3(JK,3) = N2+JW-1
    1 CONTINUE
*
      J2(3,1) = 1
      IF (JA1 .EQ. 1) J2(3,1) = N3-1
*
      J3(2,1) = 1
      IF (JA2 .EQ. 1) J3(2,1) = N3
*
      J2(NS1,3) = N2-1
*
      J3(NS1,1) = N3-2
      J3(NS1,2) = N3+3
      J3(NS1,3) = N2-1
*
      IF (JA1 .EQ. 1) GOTO 2
      JAF1 = JA1+1
      J2(JAF1,2) = N3-1
*
    2 IF (JA2 .EQ. 1) GOTO 4
      J3(JA2,2) = N3
*
      IF (NS .GT. 1) GOTO 4
    3 J3(2,1) = N3
      J3(2,2) = N3+3
      J3(2,3) = N3-1
*
    4 CONTINUE
*
*   Set the J1 array
*
    5 CONTINUE
      II = 0
*
      DO 6 JW = 1,NS
         IJ = JLIST(JW)
         II = II+1
         J1(II) = JBQ2(3,IJ)
    6 CONTINUE
*
      IF (NS .EQ. 1) GOTO 9
*
      DO 7 JW = 1,NM1
         II = II+1
         J1(II) = JJC1(JW)
    7 CONTINUE
*
      DO 8 JW = 1,NM1
         II = II+1
         J1(II) = JJC2(JW)
    8 CONTINUE
*
    9 CONTINUE
      II = II+1
      IJ = IS(1)
      J1(II) = JJQ1(3,IJ)
      J1(II+2) = KS(1)
      II = II+1
      IJ = IS(2)
      J1(II) = JJQ2(3,IJ)
      J1(II+2) = KS(2)
*
      II = II+3
      J1(II) = KA+KA+1
      MJA = II
      NJA = NS+2
*
      RETURN
*
      END
************************************************************************
*                                                                      *
      FUNCTION IROW1 (NELC,KSI)
*                                                                      *
*   Locate the row position of configuration j(**n) in table NTAB.     *
*                                                                      *
*                                           Last update: 14 Oct 1992   *
*                                                                      *
************************************************************************
*
      IF ((NELC .LE. 0) .OR. (NELC .GT. KSI)) THEN
         WRITE (*,300) NELC,KSI
         STOP
      ENDIF
*
      KQ1 = NELC-1
      KQ2 = KSI-KQ1
      KQL = MIN (KQ1,KQ2)+1
      IF (KQL .EQ. 1) THEN
         IROW1 = 1
      ELSE
         IROW1 = (KSI*(KSI-2))/8+KQL
      ENDIF
*
      RETURN
*
  300 FORMAT ('IROW1: ',I3,' electrons in shell with 2j+1 = ',I3)
*
      END
************************************************************************
*                                                                      *
      FUNCTION ITRIG (I1,I2,I3)
*                                                                      *
*   The  triangular delta. Input: Values of 2*J+1; Output: 1, IF J'S   *
*   form a triangle; 0, otherwise.                                     *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      I4 = I2-I3
      IF ((I1 .GE. (ABS (I4)+1)) .AND. ((I1 .LE. (I2+I3-1)))) THEN
         ITRIG = 1
      ELSE
         ITRIG = 0
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE KNJ (JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,JD9,KDW,
     :            JDDEL,LDDEL,DSUMVR,MDP,
     :            JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,NDB6J,
     :            KD6CP,KD7CP,KD8CP,KD9CP,JDSUM4,JDSUM5,JDSUM6,INVD6J)
*                                                                      *
*   This routine stores data for future calls to GENSUM.               *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      PARAMETER (
     :   MANGM=60,M3MNGM=3*MANGM,MANGMP=2*(MANGM/3),
     :   MTRIAD=12,
     :   M6J=20,MSUM=10)
*
      LOGICAL SUMVAR,DSUMVR
*
      DIMENSION JD6(M3MNGM),JD7(M3MNGM),JD8(M3MNGM),
     :          JD9(MANGMP),KDW(6,M6J),LDDEL(M6J,2),DSUMVR(MANGM)
      DIMENSION JD6P(MANGMP),JD7P(MANGMP),JD8P(MANGMP),JD9P(MANGMP),
     :          JDWORD(6,M6J),NDBJ(MSUM),
     :          NDB6J(MSUM),KD6CP(MSUM),KD7CP(MSUM),KD8CP(MSUM),
     :          KD9CP(MSUM),JDSUM6(MTRIAD),JDSUM4(MTRIAD,M6J),
     :          JDSUM5(MTRIAD,M6J),INVD6J(M6J)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :       J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),MP
     :      /SUMARG/J6P(MANGMP),J7P(MANGMP),J8P(MANGMP),J9P(MANGMP),
     :       JWORD(6,M6J),NLSUM,NBJ(MSUM),NB6J(MSUM),K6CP(MSUM),
     :       K7CP(MSUM),K8CP(MSUM),K9CP(MSUM),JSUM6(MTRIAD),
     :       JSUM4(MTRIAD,M6J),JSUM5(MTRIAD,M6J),INV6J(M6J)
*
      JD6C = J6C
      JD7C = J7C
      JD8C = J8C
      JD9C = J9C
      JDWC = JWC
      JDDEL = JDEL
      MDP = MP
      NDLSUM = NLSUM
      IF (J6C .NE. 0) THEN
         DO 1 I = 1,J6C
            JD6(I) = J6(I)
  1      CONTINUE
      ENDIF
      IF (J7C .NE. 0) THEN
         DO 2 I = 1,J7C
            JD7(I) = J7(I)
  2      CONTINUE
      ENDIF
      IF (J8C .NE. 0) THEN
         DO 3 I = 1,J8C
            JD8(I) = J8(I)
  3      CONTINUE
      ENDIF
      IF (J9C .NE. 0) THEN
         DO 4 I = 1,J9C
            JD9(I) = J9(I)
  4      CONTINUE
      ENDIF
      IF (JWC .NE. 0) THEN
         DO 5 I = 1,6
            DO 5 J = 1,JWC
               KDW(I,J) = KW(I,J)
  5      CONTINUE
         DO 6 I = 1,JWC
            INVD6J(I) = INV6J(I)
  6      CONTINUE
      ENDIF
      IF (JDEL .NE. 0) THEN
         DO 7 I = 1,2
            DO 7 J = 1,JDEL
               LDDEL(J,I) = LDEL(J,I)
  7      CONTINUE
      ENDIF
      IF (MP .NE. 0) THEN
         DO 8 I = 1,MP
            DSUMVR(I) = SUMVAR(I)
  8      CONTINUE
      ENDIF
      IF (NLSUM .NE. 0) THEN
         DO 9 I = 1,NLSUM
            NDBJ(I) = NBJ(I)
            NDB6J(I) = NB6J(I)
            KD6CP(I) = K6CP(I)
            KD7CP(I) = K7CP(I)
            KD8CP(I) = K8CP(I)
            KD9CP(I) = K9CP(I)
  9      CONTINUE
      ENDIF
      DO 10 I = 1,MANGMP
         JD6P(I) = J6P(I)
         JD7P(I) = J7P(I)
         JD8P(I) = J8P(I)
         JD9P(I) = J9P(I)
  10  CONTINUE
      DO 11 I = 1,MTRIAD
         JDSUM6(I) = JSUM6(I)
         DO 11 J = 1,M6J
            JDSUM4(I,J) = JSUM4(I,J)
            JDSUM5(I,J) = JSUM5(I,J)
  11  CONTINUE
      DO 12 I = 1,6
         DO 12 J = 1,M6J
            JDWORD(I,J) = JWORD(I,J)
  12  CONTINUE
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE LTAB (IS,NQS,KS,IROWS)
*                                                                      *
*   locates rows of possible parents of active shell states for acc-   *
*   essing  NTAB. It is assumed that empty shells have been elimina-   *
*   ted from consideration by SUBROUTINE RKCO.                         *
*                                                                      *
*                                           Last update: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      DIMENSION IS(4),NQS(4),KS(4),IROWS(4),KQ(4)
*
*old  COMMON/TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
*new
      nrows = 31
*new
*
      IF (IS(1) .EQ. IS(2)) NQS(1) = NQS(2)-1
      IF (IS(3) .EQ. IS(4)) NQS(3) = NQS(4)-1
*
      DO 1 I = 1,4
*
*   Check that input data are consistent
*
         IF ((NQS(I) .LE. 0) .OR. (NQS(I) .GT. KS(I))) THEN
            WRITE (*,300) NQS(I),IS(I),KS(I)
            STOP
         ENDIF
*
         KQ1 = NQS(I)-1
         KQ2 = KS(I)-KQ1
         KQ(I) = MIN (KQ1,KQ2)+1
         IF (KQ(I) .NE. 1) THEN
            IROWS(I) = (KS(I)*(KS(I)-2))/8+KQ(I)
         ELSE
            IROWS(I) = 1
         ENDIF
*
         IF (IROWS(I) .GT. NROWS) THEN
            WRITE (*,301)
            STOP
         ENDIF
*
    1 CONTINUE
*
      RETURN
*
  300 FORMAT ('LTAB: ',1I3,' Electrons in shell ',1I3,
     :        ' with 2j+1 = ',1I3)
  301 FORMAT ('LTAB: Extend COMMON block TERMS')
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE MODJ23
*                                                                      *
*   Restores  COMMON  block  /COUPLE/ from saved values for exchange   *
*   case.                                                              *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      PARAMETER (MANGM = 60,MTRIAD = 12)
*
      LOGICAL FREE
*
      COMMON/L2/J2S(MTRIAD,3),J3S(MTRIAD,3)
     :      /COUPLE/MJA,NJA,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
*
      NS2 = NJA-1
      DO 2 J = 1,3
         DO 1 I = 1,NS2
            J2(I,J) = J2S(I,J)
            J3(I,J) = J3S(I,J)
    1    CONTINUE
    2 CONTINUE
*
      I = J3(1,3)
      J3(1,3) = J2(1,1)
      J2(1,1) = I
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE MUMDAD (IS,KAPS,X)
*                                                                      *
*   Evaluate the product of 4 CFPs.                                    *
*                                                                      *
*   Call(s) to: [LIB92]: CFP.                                          *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      PARAMETER (EPS = 1.0D-10)
*
      DIMENSION IS(2,2),KAPS(2,2)
*
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /M1/NQ1(149),NQ2(149)
     :      /M2/JJQ1(3,149),JJQ2(3,149)
     :      /L1/JBQ1(3,149),JBQ2(3,149),JTQ1(3),JTQ2(3)
*
      X = 1.0D 00
*
*   First index
*
      LOCK = KAPS(1,1)
      IF (ABS (LOCK) .EQ. 2) GOTO 4
      II = IS(1,1)
      NEL = NQ1(II)
      IVP = JBQ1(1,II)
      IWP = JBQ1(2,II)
      IJP = JBQ1(3,II)-1
*
*   IA1 .NE. IB1 and IA2 .NE. IB2; use JJQ array.
*
      IF (IS(1,1) .EQ. IS(2,1)) GOTO 1
      IVD = JJQ1(1,II)
      IWD = JJQ1(2,II)
      IJD = JJQ1(3,II)-1
      GOTO 2
*
*   IA1 .EQ. IB1 or IA2 .EQ. IB2; JTQ array needed.
*
    1 NEL = NEL-1
      IVD = JTQ1(1)
      IWD = JTQ1(2)
      IJD = JTQ1(3)-1
    2 CALL CFP (LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C)
      IF (IBUG2 .NE. 0)
     :   WRITE (99,300) LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C
      IF (ABS (C) .LT. EPS) GOTO 17
      X = X*C
*
    4 LOCK = KAPS(2,1)
      IF (IABS (LOCK) .EQ. 2) GOTO 8
      II = IS(2,1)
      NEL = NQ1(II)
      IVD = JJQ1(1,II)
      IWD = JJQ1(2,II)
      IJD = JJQ1(3,II)-1
      IF (IS(1,1) .EQ. IS(2,1)) GOTO 5
      IVP = JBQ1(1,II)
      IWP = JBQ1(2,II)
      IJP = JBQ1(3,II)-1
      GOTO 6
    5 IVP = JTQ1(1)
      IWP = JTQ1(2)
      IJP = JTQ1(3)-1
    6 CALL CFP (LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C)
      IF (IBUG2 .NE. 0)
     :   WRITE (99,300) LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C
      IF (ABS (C) .LT. EPS) GOTO 17
      X = X*C
    8 CONTINUE
*
*   Second index
*
      LOCK = KAPS(1,2)
      IF (ABS (LOCK) .EQ. 2) GOTO 12
      II = IS(1,2)
      NEL = NQ2(II)
      IVP = JBQ2(1,II)
      IWP = JBQ2(2,II)
      IJP = JBQ2(3,II)-1
*
*   IA1 .NE. IB1 and IA2 .NE. IB2; use JJQ array.
*
      IF (IS(1,2) .EQ. IS(2,2)) GOTO 9
      IVD = JJQ2(1,II)
      IWD = JJQ2(2,II)
      IJD = JJQ2(3,II)-1
      GOTO 10
*
*   IA1 .EQ. IB1 or IA2 .EQ. IB2; JTQ array needed.
*
    9 NEL = NEL-1
      IVD = JTQ2(1)
      IWD = JTQ2(2)
      IJD = JTQ2(3)-1
   10 CALL CFP (LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C)
      IF (IBUG2 .NE. 0)
     :   WRITE (99,300) LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C
      IF (ABS (C) .LT. EPS) GOTO 17
      X = X*C
*
   12 LOCK = KAPS(2,2)
      IF (ABS (LOCK) .EQ. 2) GOTO 16
      II = IS(2,2)
      NEL = NQ2(II)
      IVD = JJQ2(1,II)
      IWD = JJQ2(2,II)
      IJD = JJQ2(3,II)-1
      IF (IS(1,2) .EQ. IS(2,2)) GOTO 13
      IVP = JBQ2(1,II)
      IWP = JBQ2(2,II)
      IJP = JBQ2(3,II)-1
      GOTO 14
   13 IVP = JTQ2(1)
      IWP = JTQ2(2)
      IJP = JTQ2(3)-1
   14 CALL CFP (LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C)
      IF (IBUG2 .NE. 0)
     :   WRITE (99,300) LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C
      IF (ABS (C) .LT. EPS) GOTO 17
      X = X*C
   16 CONTINUE
      RETURN
*
   17 X = 0.0D 00
      RETURN
*
  300 FORMAT ('MUMDAD: CFP ',I3,I4,I7,2I4,I7,2I4,1P,1D19.12)
*
      END
************************************************************************
*                                                                      *
      FUNCTION OCON (IA1,IB1,IA2,IB2)
*                                                                      *
*   Evaluates the  multiplicative statistical  factor. It is assumed   *
*   that states are ordered so that IA1 .LE. IB1, IA2 .LE. IB2.        *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      COMMON/M1/NQ1(149),NQ2(149)
*
      WA = DBLE (NQ1(IA1)*NQ1(IB1))
      IF (IA1 .EQ. IB1) WA = WA-DBLE (NQ1(IA1))
      WB = DBLE (NQ2(IA2)*NQ2(IB2))
      IF (IA2 .EQ. IB2) WB = WB-DBLE (NQ2(IB2))
      WC = WA*WB
      OCON = SQRT (WC)
*
*   Set phase factor (-1)**(DELTA P)
*
      LRD1 = MIN (IA2,IB2)+1
      LRD2 = MAX (IA2,IB2)
      IF (LRD1 .GT. LRD2) THEN
         IDR = 0
      ELSE
         IDR = 1
         DO 1 K = LRD1,LRD2
            IDR = IDR+NQ2(K)
    1    CONTINUE
      ENDIF
*
      LLD1 = MIN (IA1,IB1)+1
      LLD2 = MAX (IA1,IB1)
      IF (LLD1 .GT. LLD2) THEN
         IDL = 0
      ELSE
         IDL = 1
         DO 2 K = LLD1,LLD2
            IDL = IDL+NQ1(K)
    2    CONTINUE
      ENDIF
*
      IPHAS = IDR-IDL
      IF (MOD (IPHAS,2) .NE. 0) OCON = -OCON
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE RKCO (JA,JB,COR,CORD,INCOR)
*                                                                      *
*   Configurations JA, JB. Analyse the tables of quantum numbers set   *
*   in the COMMON  blocks M0 , M1, M2, M3  to determine all possible   *
*   sets of interacting  orbitals which give a non-vanishing Coulomb   *
*   matrix element,  and  initiates the calculation of coefficients.   *
*   The following conventions are in force: (1) labels 1, 2 refer to   *
*   left, right sides of matrix element respectively;   (2) pointers   *
*   JA1, JB1, JA2, JB2 point to the JLIST array of active  orbitals;   *
*   IA1, IB1, IA2, IB2 point to the complete list of orbitals.         *
*                                                                      *
*   Call(s) to: [LIB92]: COR, CORD, ISPAR, ITJPO, SETQNA, VIJOUT.      *
*                                                                      *
*                                           Last update: 02 Nov 1992   *
*                                                                      *
************************************************************************
*
      INTEGER PNTRIQ
*
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /DUMX/JLIS(149),JC1S(149),JC2S(149)
     :      /M0/JJC1(149),JJC2(149)
     :      /M1/NQ1(149),NQ2(149)
     :      /M2/JJQ1(3,149),JJQ2(3,149)
     :      /M3/JLIST(149),KLIST(149),NPEEL,NCORE
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(149),NAK(149)
*
*   The Hamiltonian is an even scalar operator
*
      IF (ITJPO (JA) .NE. ITJPO (JB)) RETURN
      IF (ISPAR (JA) .NE. ISPAR (JB)) RETURN
*
      CALL SETQNA (JA,JB)
      IF (IBUG2 .EQ. 1) CALL VIJOUT (JA,JB)
*
*   1.0 Analyse peel shell interactions
*
*   1.1 Analyse electron distribution in peel. (The full procedure is
*       needed only if the number of peel orbitals NPEEL .GE. 2)
*
      IF (NW .LT. 1) THEN
         PRINT *, 'RKCO: No subshells.'
         STOP
      ENDIF
      IF (NPEEL .EQ. 0) GOTO 48
      IF (NPEEL .EQ. 1) GOTO 43
*
*   Find differences in occupations, NDQ, for each peel orbital in
*   turn and use to set up labels of active orbitals maintaining the
*   convention JA1 .LE. JB1, JA2 .LE. JB2.
*
      IDQ = 0
      JA1 = 0
      JB1 = 0
      JA2 = 0
      JB2 = 0
      DO 10 JW = 1,NPEEL
         J = JLIST(JW)
         NDQ = NQ1(J) - NQ2(J)
         IF (IABS (NDQ) .GT. 2) RETURN
         IF (NDQ .LT. 0) GOTO 5
         IF (NDQ-1) 10,1,4
    1    IF (JA1 .GT. 0) GOTO 2
         JA1 = JW
         GOTO 3
    2    JB1 = JW
    3    IDQ = IDQ+1
         GOTO 10
    4    JA1 = JW
         IDQ = IDQ+2
         GOTO 10
    5    IF (NDQ+1) 9,6,10
    6    IF (JA2 .GT. 0) GOTO 7
         JA2 = JW
         GOTO 8
    7    JB2 = JW
    8    IDQ = IDQ+1
         GOTO 10
    9    JA2 = JW
         IDQ = IDQ+2
   10 CONTINUE
*
*   1.2 Calculate coefficients for all possible sets of active shells.
*
*   There are 4 cases, depending on the value of IDQ, the sum of the
*   absolute differences NDQ:
*
*   1.2.1 IDQ .GT. 4: matrix element null
*
      IF (IDQ .GT. 4) RETURN
      IF (IDQ .EQ. 4) GOTO 12
      IF (IDQ .NE. 2) GOTO 11
      KLAST = 1
      GOTO 16
   11 IF (IDQ .NE. 0) GOTO 54
      IF (JA .EQ. JB) GOTO 43
      KLAST = NPEEL
      GOTO 16
*
*   1.2.2 IDQ .EQ. 4: matrix element uniquely defined
*
   12 IF (JB1 .NE. 0) GOTO 13
      JB1 = JA1
   13 IF (JB2 .NE. 0) GOTO 14
      JB2 = JA2
   14 IF (IBUG2 .NE. 0) WRITE (99,301) JA1,JB1,JA2,JB2
      CALL COR (JA,JB,JA1,JB1,JA2,JB2)
      RETURN
*
*   1.2.3 IDQ .EQ. 2: One orbital fixed each side include all
*                     possible spectators.
*
*   Also IDQ .EQ. 0 for a matrix element off-diagonal in coupling
*   only. Must sum over all pairs of orbitals excluding core-core
*   terms
*
   16 DO 42 KWA = 1,KLAST
         IF (IDQ .EQ. 2) GOTO 17
         JA1 = KWA
         JA2 = KWA
   17    JT1 = JA1
         JT2 = JA2
         IT1 = JLIST(JA1)
         IT2 = JLIST(JA2)
         DO 25 KW = KWA,NPEEL
            K1 = JLIST(KW)
            IF (NQ1(K1)*NQ2(K1) .EQ. 0) GOTO 25
            JB1 = KW
            JB2 = KW
            JA1 = JT1
            JA2 = JT2
*
*   Interchange JA1 and JB1 and/or JA2 and JB2 if necessary
*
            IF (JA1-JB1) 20,19,18
   18       JT3 = JB1
            JB1 = JA1
            JA1 = JT3
            GOTO 20
   19       IB1 = JLIST(JB1)
            IF (NQ1(IB1) .LE. 1) GOTO 25
   20       IF (JA2-JB2) 23,22,21
   21       JT3 = JB2
            JB2 = JA2
            JA2 = JT3
            GOTO 23
   22       IB2 = JLIST(JB2)
            IF (NQ2(IB2) .LE. 1) GOTO 25
   23       IF (IBUG2 .NE. 0) WRITE (99,301) JA1,JB1,JA2,JB2
            CALL COR (JA,JB,JA1,JB1,JA2,JB2)
   25    CONTINUE
         IF ((IDQ .EQ. 0) .AND. (NCORE .EQ. 0)) GOTO 42
         IF ((NCORE .EQ. 0) .OR. (NAK(IT1) .NE. NAK(IT2))) RETURN
*
*   This section calculates the terms arising from active electrons
*   which are in closed shells
*
         NPEELM = NPEEL-1
         DO 26 I = 1,NPEEL
            JLIS(I) = JLIST(I)
   26    CONTINUE
         DO 27 I = 1,NPEELM
            JC1S(I) = JJC1(I)
            JC2S(I) = JJC2(I)
   27    CONTINUE
         DO 41 KW = 1,NCORE
            IJW = KLIST(KW)
            DO 28 I = 1,NPEEL
               IJ = JLIST(I)
               IF (IJW .LT. IJ) GOTO 29
   28       CONTINUE
            I = NPEEL+1
            GOTO 31
   29       IM = NPEEL-I+1
            DO 30 II = 1,IM
               JLIST(NPEEL+2-II) = JLIST(NPEEL+1-II)
               IF (NPEEL .EQ. II) GOTO 31
               JJC1(NPEEL+1-II) = JJC1(NPEEL-II)
               JJC2(NPEEL+1-II) = JJC2(NPEEL-II)
   30       CONTINUE
   31       CONTINUE
            IF (I .LT. 3) GOTO 32
            JJC1(I-1) = JJC1(I-2)
            JJC2(I-1) = JJC2(I-2)
            GOTO 33
   32       I1 = JLIST(1)
            JJC1(1) = JJQ1(3,I1)
            JJC2(1) = JJQ2(3,I1)
   33       JLIST(I) = IJW
            JA1 = JT1
            IF (JT1 .GE. I) JA1 = JA1+1
            JB1 = I
            JA2 = JT2
            IF (JT2 .GE. I) JA2 = JA2+1
            JB2 = I
            IF (JA1-JB1) 35,35,34
   34       JT3 = JB1
            JB1 = JA1
            JA1 = JT3
   35       CONTINUE
            IF (JA2-JB2) 37,37,36
   36       JT3 = JB2
            JB2 = JA2
            JA2 = JT3
   37       CONTINUE
            NPEEL = NPEEL+1
            IF (IBUG2 .NE. 0) THEN
               NPEELM = NPEEL-1
               WRITE (99,302) JA1,JB1,JA2,JB2,KW,KLIST(KW)
               WRITE (99,303) (JLIST(I),I = 1,NPEEL)
               WRITE (99,304) (JJC1(I),I = 1,NPEELM)
               WRITE (99,305) (JJC2(I),I = 1,NPEELM)
            ENDIF
            CALL COR (JA,JB,JA1,JB1,JA2,JB2)
            NPEEL = NPEEL-1
            NPEELM = NPEEL-1
            DO 39 I = 1,NPEEL
               JLIST(I) = JLIS(I)
   39       CONTINUE
            DO 40 I = 1,NPEELM
               JJC1(I)  = JC1S(I)
               JJC2(I)  = JC2S(I)
   40       CONTINUE
   41    CONTINUE
   42 CONTINUE
      RETURN
*
*   1.2.4 IDQ .EQ. 0 - diagonal case. Include all pairs with
*         JA1 = JA2, JB1 = JB2.
*
   43 DO 47 KW1 = 1,NPEEL
         K1 = JLIST(KW1)
         JB1 = KW1
         JB2 = KW1
         DO 46 KW2 = 1,KW1
            JA1 = KW2
            IF (JA1 .NE. JB1) GOTO 44
            IF (NQ1(K1) .LE. 1) GOTO 46
   44       JA2 = JA1
            IF (IBUG2 .NE. 0) WRITE (99,301) JA1,JB1,JA2,JB2
            CALL COR (JA,JB,JA1,JB1,JA2,JB2)
   46    CONTINUE
   47 CONTINUE
   48 IF (INCOR .LT. 1) RETURN
      IF (NCORE .EQ. 0) RETURN
*
*   2.0 The diagonal case. deal with contributions from core orbitals
*       if INCOR .EQ. 1.
*
      DO 53 KW1 = 1,NCORE
         JB1 = KW1
         JB2 = KW1
*
*   2.1 Calculate contribution from core/core terms
*
         IPCA = 2
         DO 50 KW2 = 1,KW1
            JA1 = KW2
            JA2 = KW2
            IF (IBUG2 .NE. 0) WRITE (99,301) JA1,JB1,JA1,JB1
            CALL CORD (JA,JB,JA1,IPCA,JB1)
   50    CONTINUE
*
*   2.2 Calculate contribution from peel/core terms
*
         IF (NPEEL .EQ. 0) GOTO 53
         IPCA = 1
         DO 52 KW2 = 1,NPEEL
            JA1 = KW2
            JA2 = KW2
            IF (IBUG2 .NE. 0) WRITE (99,301) JA1,JB1,JA1,JB1
            CALL CORD (JA,JB,JA1,IPCA,JB1)
   52    CONTINUE
   53 CONTINUE
      RETURN
*
*   3.0 Diagnostic print - NW .LT. 1
*
   54 WRITE (*,300)
      STOP
*
  300 FORMAT ('RKCO: Error.')
  301 FORMAT ('From RKCO:'
     :         /10X,' JA1 = ',I3,4X,' JB1 = ',I3,4X,' JA2 = ',I3,
     :              ' JB2 = ',I3)
  302 FORMAT ('From RKCO:'
     :   /10X,' JA1 = ',I3,4X,' JB1 = ',I3,4X,' JA2 = ',I3,
     :        ' JB2 = ',I3,' K2  = ',I3,   ' KW  = ',I3)
  303 FORMAT (1X,'JLIST : ',25I4)
  304 FORMAT (1X,'JJC1  : ',25I4)
  305 FORMAT (1X,'JJC2  : ',25I4)
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE SETJ (IS,JS,KS,NS,KJ23)
*                                                                      *
*   Sets the tables required by  the recoupling  coefficient package   *
*   NJGRAF. This routine loads  the COMMON block /COUPLE/ with para-   *
*   meters for the first call  of NJGRAF involving direct integrals.   *
*   Subsequent exchange calls  of NJGRAF must be preceeded by a call   *
*   of MODJ23 to restore these arrays to their correct initial state.  *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      PARAMETER (MANGM = 60,MTRIAD = 12)
*
      LOGICAL FREE
*
      DIMENSION IS(2,2),JS(2,2),KS(2,2)
*
      COMMON/L1/JBQ1(3,149),JBQ2(3,149),JTQ1(3),JTQ2(3)
     :      /L2/J2S(MTRIAD,3),J3S(MTRIAD,3)
     :      /M0/JJC1(149),JJC2(149)
     :      /M2/JJQ1(3,149),JJQ2(3,149)
     :      /M3/JLIST(149),KLIST(149),NPEEL,NCORE
     :      /COUPLE/MJA,NJA,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
*
*   1.0  Set J1 array
*
      II = 0
      DO 1 IJ = 1,NS
         I = JLIST(IJ)
         II = II+1
         J1(II) = JBQ1(3,I)
    1 CONTINUE
      IF (NS .EQ. 1) GOTO 4
      NS1 = NS-1
      DO 2 I = 1,NS1
         II = II+1
         J1(II) = JJC1(I)
    2 CONTINUE
      DO 3 I = 1,NS1
         II = II+1
         J1(II) = JJC2(I)
    3 CONTINUE
    4 CONTINUE
      DO 5 I = 1,2
         II = II+1
         IJ = IS(I,1)
         J1(II) = JJQ1(3,IJ)
         IF ((I .EQ. 1) .AND. (IS(1,1) .EQ. IS(2,1))) J1(II) = JTQ1(3)
         J1(II+4) = KS(I,1)
    5 CONTINUE
      DO 6 I = 1,2
         II = II+1
         IJ = IS(I,2)
         J1(II) = JJQ2(3,IJ)
         IF ((I .EQ. 1) .AND. (IS(1,2) .EQ. IS(2,2))) J1(II) = JTQ2(3)
         J1(II+4) = KS(I,2)
    6 CONTINUE
*
*   2.0  Set J2, J3 arrays if not already available
*
      NS2 = MAX (4,NS+2)
      IF (KJ23 .GT. 0) GOTO 14
*
      DO 7 I = 4,NS2
         J2(I,1) = NS+I-4
         J2(I,2) = I-2
         J2(I,3) = NS+I-3
         J3(I,1) = J2(I,1)+NS-1
         J3(I,2) = I-2
         J3(I,3) = J2(I,3)+NS-1
    7 CONTINUE
      J2(4,1) = 1
      J3(4,1) = 1
*
*   At this stage, the entries in rows corresponding to active
*   shells are set incorrectly.
*
*   3.0  Set rows 1 through 3
*
      NS3 = 3*NS
      J2(1,1) = NS3+5
      J2(1,2) = NS3+7
      J2(1,3) = NS3+3
      J2(2,1) = JS(1,1)
      J2(2,2) = NS3+3
      J2(2,3) = NS3-1
      J2(3,1) = JS(2,1)
      J2(3,2) = NS3+4
      J2(3,3) = NS3
*
      J3(1,1) = NS3+7
      J3(1,2) = NS3+4
      J3(1,3) = NS3+6
      J3(2,1) = JS(1,2)
      J3(2,2) = NS3+5
      J3(2,3) = NS3+1
      J3(3,1) = JS(2,2)
      J3(3,2) = NS3+6
      J3(3,3) = NS3+2
*
*   4.0  Set remaining resultants
*
      IJ1 = JS(1,1)
      IJ2 = JS(2,1)
      IF (IJ2 .GT. 1) J2(IJ2+2,2) = J2(3,3)
      IF (IJ2 .EQ. 1) J2(4,1) = J2(3,3)
      IF (IJ1 .NE. IJ2) GOTO 8
      J2(3,1) = J2(2,3)
      GOTO 9
*
    8 IF (IJ1 .GT. 1) J2(IJ1+2,2) = J2(2,3)
      IF (IJ1 .EQ. 1) J2(4,1) = J2(2,3)
*
    9 IJ1 = JS(1,2)
      IJ2 = JS(2,2)
      IF (IJ2 .GT. 1) J3(IJ2+2,2) = J3(3,3)
      IF (IJ2 .EQ. 1) J3(4,1) = J3(3,3)
      IF (IJ1 .NE. IJ2) GOTO 10
      J3(3,1) = J3(2,3)
      GOTO 11
*
   10 IF (IJ1 .GT. 1) J3(IJ1+2,2) = J3(2,3)
      IF (IJ1 .EQ. 1) J3(4,1) = J3(2,3)
*
*   All arrays now set. Put up flag KJ23.
*
   11 KJ23 = 1
      MJA = NS3+7
      NJA = NS+3
*
*   5.0  Save J2, J3 and return
*
      DO 13 J = 1,3
         DO 12 I = 1,NS2
            J2S(I,J) = J2(I,J)
            J3S(I,J) = J3(I,J)
   12    CONTINUE
   13 CONTINUE
      RETURN
*
*   6.0  Reset J2, J3 from buffers if KJ23 has been set
*
   14 DO 16 J = 1,3
         DO 15 I = 1,NS2
            J2(I,J) = J2S(I,J)
            J3(I,J) = J3S(I,J)
   15    CONTINUE
   16 CONTINUE
      RETURN
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE SETQNA (JA,JB)
*                                                                      *
*   This generates the  arrays  defining  the quantum numbers of the   *
*   states involved in the  matrix  element  linking  configurations   *
*   labelled by JA, JB.                                                *
*                                                                      *
*   Call(s) to: [LIB92]: ICHOP, IQ, JCUP, JQS.                         *
*                                                                      *
*                                           Last update: 30 Oct 1987   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INTEGER PNTRIQ
*
      COMMON/M0/JJC1(149),JJC2(149)
     :      /M1/NQ1(149),NQ2(149)
     :      /M2/JJQ1(3,149),JJQ2(3,149)
     :      /M3/JLIST(149),KLIST(149),NPEEL,NCORE
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(149),NAK(149)
*
*   List parameters defining all shells in both configurations, whether
*   participating or not
*
      DO 2 J = 1,NW
         NQ1(J) = IQ (J,JA)
         NQ2(J) = IQ (J,JB)
         DO 1 K = 1,3
            JJQ1(K,J) = JQS (K,J,JA)
            JJQ2(K,J) = JQS (K,J,JB)
    1    CONTINUE
    2 CONTINUE
*
*   Define coupling schemes: set JLIST array to define those shells
*   which are open in either configuration, and KLIST array to locate
*   the rest. Exclude shells which are empty in both configurations
*
      NPEEL = 0
      NCORE = 0
      DO 3 J = 1,NW
         IF ((ICHOP (J,JA) .NE. -1) .OR. (ICHOP (J,JB) .NE. -1)) THEN
            IF ((ICHOP (J,JA) .EQ. 1) .AND. (ICHOP (J,JB) .EQ. 1))
     :         THEN
               NCORE = NCORE+1
               KLIST(NCORE) = J
            ELSE
               NPEEL = NPEEL+1
               JLIST(NPEEL) = J
            ENDIF
         ENDIF
    3 CONTINUE
*
*   Return if not more than one shell is open
*
      IF (NPEEL .LE. 1) RETURN
*
*   Set arrays of coupling angular momenta interpolating closed
*   shells where necessary. Left hand side first ...
*
      JCNT = 1
      JCNTOP = 0
      JW1 = JLIST(1)
      JW2 = JLIST(2)
      IF (ICHOP (JW1,JA) .NE. 0) THEN
         JJC1(1) = JQS (3,JW2,JA)
         IF (ICHOP (JW2,JA) .EQ. 0) JCNTOP = 1
      ELSE
         JCNTOP = 1
         IF (ICHOP (JW2,JA) .EQ. 0) THEN
            JJC1(1) = JCUP (JCNT,JA)
            JCNT = JCNT+1
         ELSE
            JJC1(1) = JQS (3,JW1,JA)
         ENDIF
      ENDIF
*
      DO 4 J = 3,NPEEL
         JW = JLIST(J)
         IF (ICHOP (JW,JA) .NE. 0) THEN
            JJC1(J-1) = JJC1(J-2)
         ELSE
            IF (JCNTOP .NE. 0) THEN
               JJC1(J-1) = JCUP (JCNT,JA)
               JCNT = JCNT+1
            ELSE
               JJC1(J-1) = JQS (3,JW,JA)
            ENDIF
            JCNTOP = JCNTOP+1
         ENDIF
    4 CONTINUE
*
*   ... and repeat for right hand side
*
      JCNT = 1
      JCNTOP = 0
      JW1 = JLIST(1)
      JW2 = JLIST(2)
      IF (ICHOP (JW1,JB) .NE. 0) THEN
         JJC2(1) = JQS (3,JW2,JB)
         IF (ICHOP (JW2,JB) .EQ. 0) JCNTOP = 1
      ELSE
         JCNTOP = 1
         IF (ICHOP (JW2,JB) .EQ. 0) THEN
            JJC2(1) = JCUP (JCNT,JB)
            JCNT = JCNT+1
         ELSE
            JJC2(1) = JQS (3,JW1,JB)
         ENDIF
      ENDIF
*
      DO 5 J = 3,NPEEL
         JW = JLIST(J)
         IF (ICHOP (JW,JB) .NE. 0) THEN
            JJC2(J-1) = JJC2(J-2)
         ELSE
            IF (JCNTOP .NE. 0) THEN
               JJC2(J-1) = JCUP (JCNT,JB)
               JCNT = JCNT+1
            ELSE
               JJC2(J-1) = JQS (3,JW,JB)
            ENDIF
            JCNTOP = JCNTOP+1
         ENDIF
    5 CONTINUE
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE SKRC (IS,KAPS,KS,KD1,KD2,KE1,KE2)
*                                                                      *
*   Determines the range of the tensor rank k for Coulomb integral.    *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      DIMENSION IS(4),KAPS(4),KS(4)
*
      KD2 = 0
      KE2 = 0
*
*   Direct terms --- KD1 = minimum k , KD2 = number of terms
*
      ISD1 = 1
      IF (KAPS(1)*KAPS(3) .LT. 0) ISD1 = -1
      ISD2 = 1
      IF (KAPS(2)*KAPS(4) .LT. 0) ISD2 = -1
      KD1A = ABS (KS(1)-KS(3))
      IF (ISD1 .LT. 0) KD1A = KD1A+2
      KD1B = ABS (KS(2)-KS(4))
      IF (ISD2 .LT. 0) KD1B = KD1B+2
      IF (MOD ((KD1A-KD1B)/2,2) .NE. 0) GOTO 1
      KD2A = KS(1)+KS(3)-2
      IF (ISD1 .GT. 0) KD2A = KD2A-2
      KD2B = KS(2)+KS(4)-2
      IF (ISD2 .GT. 0) KD2B = KD2B-2
      KD1 = MAX (KD1A,KD1B)/2
      KD2 = MIN (KD2A,KD2B)/2
      KD2 = (KD2-KD1)/2+1
*
*   Exchange terms --- KE1 = minimum k , KE2 = number of terms
*
    1 CONTINUE
      IF ((IS(1) .EQ. IS(2)) .OR. (IS(3) .EQ. IS(4))) RETURN
      ISE1 = 1
      IF (KAPS(1)*KAPS(4) .LT. 0) ISE1 = -1
      ISE2 = 1
      IF (KAPS(2)*KAPS(3) .LT. 0) ISE2 = -1
      KE1A = ABS (KS(1)-KS(4))
      IF (ISE1 .LT. 0) KE1A = KE1A+2
      KE1B = ABS (KS(2)-KS(3))
      IF (ISE2 .LT. 0) KE1B = KE1B+2
      IF (MOD ((KE1A-KE1B)/2,2) .NE. 0) RETURN
      KE2A = KS(1)+KS(4)-2
      IF (ISE1 .GT. 0) KE2A = KE2A-2
      KE2B = KS(2)+KS(3)-2
      IF (ISE2 .GT. 0) KE2B = KE2B-2
      KE1 = MAX (KE1A,KE1B)/2
      KE2 = MIN (KE2A,KE2B)/2
      KE2 = (KE2-KE1)/2+1
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE SNRC (IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
*                                                                      *
*   Determines the range of tensor rank NU for direct/exchange terms,  *
*   and classifies the types of radial integral.                       *
*                                                                      *
*   Input variables:                                                   *
*                                                                      *
*      IS      : Orbital labels                                        *
*      KAPS    : Values of 2*kappa                                     *
*      KS      : Values of 2*J+1                                       *
*                                                                      *
*   Outputs:                                                           *
*                                                                      *
*      ND1/NE1 : Lowest NU value for direct/exchange types             *
*      ND2/NE2 : Corresponding number of contributing NU values: NU    *
*                = ND1, ND1+2 ,..., ND1+2*(ND2-1) etc                  *
*      IBRD    : Classify types of  radial  integrals  contributing;   *
*      IBRE      negative value implies null contribution              *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      DIMENSION IS(4),KAPS(4),KS(4)
*
      ND2 = 0
      NE2 = 0
*
*   2.0  Form limits for direct terms
*
      IAC = 1
      IF ((KAPS(1)*KAPS(3)) .LT. 0) IAC = -1
      IAD = 1
      IF ((KAPS(2)*KAPS(4)) .LT. 0) IAD = -1
      ND1 = ABS (KS(1)-KS(3))/2-1
      IF (IAC .EQ. -1) ND1 = ND1+1
      IF (ND1 .EQ. -1) ND1 = 1
      ND1A = ABS (KS(2)-KS(4))/2-1
      IF (IAD .EQ. -1) ND1A = ND1A+1
      IF (ND1A .EQ. -1) ND1A = 1
      IF (MOD (ND1-ND1A,2) .EQ. 0) GOTO 1
      IBRD = -1
      GOTO 2
    1 ND2 = ABS (KS(1)+KS(3))/2
      IF (IAC .EQ. -1) ND2 = ND2+1
      ND2A = ABS (KS(2)+KS(4))/2
      IF (IAD .EQ. -1) ND2A = ND2A+1
      ND1 = MAX (ND1,ND1A)
      ND2 = MIN (ND2,ND2A)
      ND2 = (ND2-ND1)/2+1
*
*   2.1  Identify type of radial integrals
*
      IBRD = 1
      IF (IS(1) .EQ. IS(3) .AND. IS(2) .NE. IS(4)) IBRD = 2
      IF (IS(1) .NE. IS(3) .AND. IS(2) .EQ. IS(4)) IBRD = 2
      IF (IS(1) .EQ. IS(3) .AND. IS(2) .EQ. IS(4)) IBRD = 3
*
*   3.0  Form limits for exchange terms
*
    2 IF ((IS(1) .NE. IS(2)) .AND. (IS(3) .NE. IS(4))) GOTO 3
      IBRE = -1
      RETURN
    3 CONTINUE
      IAC = 1
      IF ((KAPS(1)*KAPS(4)) .LT. 0) IAC = -1
      IAD = 1
      IF ((KAPS(2)*KAPS(3)) .LT. 0) IAD = -1
      NE1 = IABS (KS(1)-KS(4))/2-1
      IF (IAC .EQ. -1) NE1 = NE1+1
      IF (NE1 .EQ. -1) NE1 = 1
      NE1A = ABS (KS(2)-KS(3))/2-1
      IF (IAD .EQ. -1) NE1A = NE1A+1
      IF (NE1A .EQ. -1) NE1A = 1
      IF (MOD (NE1-NE1A,2) .EQ. 0) GOTO 4
      IBRE = -1
      RETURN
*
    4 NE2 = ABS (KS(1)+KS(4))/2
      IF (IAC .EQ. -1) NE2 = NE2+1
      NE2A = ABS (KS(2)+KS(3))/2
      IF (IAD .EQ. -1) NE2A = NE2A+1
      NE1 = MAX (NE1,NE1A)
      NE2 = MIN (NE2,NE2A)
      NE2 = (NE2-NE1)/2+1
*
*   3.1  Identify type of radial integrals
*
      IBRE = 1
      IF ((IS(1) .EQ. IS(4)) .AND. (IS(2) .NE. IS(3))) IBRE = 2
      IF ((IS(1) .NE. IS(4)) .AND. (IS(2) .EQ. IS(3))) IBRE = 2
      IF ((IS(1) .EQ. IS(3)) .AND. (IS(2) .EQ. IS(4))) IBRE = 4
      RETURN
*
      END
************************************************************************
*                                                                      *
      block data term
*                                                                      *
*   Taken from GRASP2.                                                 *
*                                                                      *
*   Table of subshell quantum numbers.    Symmetry  of the table for   *
*   particle/hole configurations is used to compress it.               *
*                                                                      *
*                                          Last updated: 30 Sep 1992   *
*                                                                      *
************************************************************************
*
      COMMON/TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
*
      DATA NROWS/31/
*
*   A row is defined by a subshell angular momentum and an occupation
*   number
*
*   Each entry ITAB gives the number of terms in a row
*
*   Each entry JTAB gives the starting location -1 of the first triad
*   in a row
*
*   Each triad in NTAB is (v,w,2J+1); here v is the seniority,
*   w resolves any degeneracy in the seniority scheme, and J is the
*   subshell total angular momentum
*
*   Empty subshell or full subshell
*
      DATA (ITAB(I),I =   1,  1)/
     :  1/
      DATA (JTAB(I),I =   1,  1)/
     :  0/
      DATA (NTAB(I),I =   1,  3)/
     :  0,0, 1/
*
*   s, p-   (j = 1/2)
*
      DATA (ITAB(I),I =   2,  2)/
     :  1/
      DATA (JTAB(I),I =   2,  2)/
     :  3/
      DATA (NTAB(I),I =   4,  6)/
     :  1,0, 2/
*
*   p, d-   (j = 3/2)
*
      DATA (ITAB(I),I =   3,  4)/
     :  1,2/
      DATA (JTAB(I),I =   3,  4)/
     :    6,  9/
      DATA (NTAB(I),I =   7, 15)/
     :  1,0, 4,
     :  0,0, 1,
     :  2,0, 5/
*
*  d, f-   (j = 5/2)
*
      DATA (ITAB(I),I =   5,  7)/
     :  1,3,3/
      DATA (JTAB(I),I =   5,  7)/
     :   15, 18, 27/
      DATA (NTAB(I),I =  16, 36)/
     :  1,0, 6,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9,
     :  1,0, 6,
     :  3,0, 4, 3,0,10/
*
*   f, g-   (j = 7/2)
*
      DATA (ITAB(I),I =   8, 11)/
     :  1,4,6,8/
      DATA (JTAB(I),I =   8, 11)/
     :   36, 39, 51, 69/
      DATA (NTAB(I),I =  37, 93)/
     :  1,0, 8,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9, 2,0,13,
     :  1,0, 8,
     :  3,0, 4, 3,0, 6, 3,0,10, 3,0,12, 3,0,16,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9, 2,0,13,
     :  4,0, 5, 4,0, 9, 4,0,11, 4,0,17/
*
*   g, h-   (j = 9/2)
*
      DATA (ITAB(I),I =  12, 16)/
     :  1,5,10,18,20/
      DATA (JTAB(I),I =  12, 16)/
     :   93, 96, 111,141,195/
      DATA (NTAB(I),I =  94,255)/
     :  1,0,10,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9, 2,0,13, 2,0,17,
     :  1,0,10,
     :  3,0, 4, 3,0, 6, 3,0, 8, 3,0,10, 3,0,12, 3,0,14, 3,0,16, 3,0,18,
     :  3,0,22,
     :  0,0, 1,
     :  2,0,5, 2,0,9, 2,0,13, 2,0,17,
     :  4,0, 1, 4,0, 5, 4,0, 7, 4,0, 9, 4,1, 9, 4,0,11, 4,0,13, 4,1,13,
     :  4,0,15, 4,0,17, 4,0,19, 4,0,21, 4,0,25,
     :  1,0,10,
     :  3,0, 4, 3,0, 6, 3,0, 8, 3,0,10, 3,0,12, 3,0,14, 3,0,16, 3,0,18,
     :  3,0,22,
     :  5,0, 2, 5,0, 6, 5,0, 8, 5,0,10, 5,0,12, 5,0,14, 5,0,16, 5,0,18,
     :  5,0,20, 5,0,26/
*
*   h, i-   (j = 11/2)
*
*   First two rows only
*
      DATA (ITAB(I),I =  17, 18)/
     :  1,6/
      DATA (JTAB(I),I =  17, 19)/
     :  255,258,277/
      DATA (NTAB(I),I = 256,276)/
     :  1,0,12,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9, 2,0,13, 2,0,17, 2,0,21/
*
*   i, k-   (j = 13/2)
*
*   First two rows only
*
      DATA (ITAB(I),I =  23, 24)/
     :  1,7/
      DATA (JTAB(I),I =  23, 25)/
     :  276,279,301/
      DATA (NTAB(I),I = 277,300)/
     :  1,0,14,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9, 2,0,13, 2,0,17, 2,0,21, 2,0,25/
*
*   k, l-   (j = 15/2)
*
*   First two rows only
*
      DATA (ITAB(I),I =  30, 31)/
     :  1,8/
      DATA (JTAB(I),I =  30, 32)/
     :  300,303,328/
      DATA (NTAB(I),I = 301,327)/
     :  1,0,16,
     :  0,0, 1,
     :  2,0, 5, 2,0, 9, 2,0,13, 2,0,17, 2,0,21, 2,0,25, 2,0,29/
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE TNSRJJ (KA,IOPAR,JA,JB,IA1,IA2,VSHELL)
*                                                                      *
*   The  main  program for evaluating the reduced matrix elements of   *
*   a one particle operator for configurations in jj-coupling.         *
*                                                                      *
*   Call(s) to: [LIB92]: CFP, FIXJ, GENSUM, ICHOP, IROW1, ISPAR,       *
*                        ITJPO, ITRIG, SETQNA, VIJOUT.                 *
*               [NJGRAF]: NJGRAF.                                      *
*                                                                      *
*                                           Last update: 02 Nov 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      INTEGER PNTRIQ
*
      PARAMETER (
     :   MANGM=60,M3MNGM=3*MANGM,MANGMP=2*(MANGM/3),
     :   MTRIAD=12,
     :   M6J=20,MSUM=10)
*
      PARAMETER (EPS=1.0D-10)
*
      LOGICAL FREE,SUMVAR,FAIL
*
      DIMENSION VSHELL(149)
      DIMENSION IS(2),KS(2)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :       J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),MP
     :      /CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
     :      /COUPLE/MJA,NJA,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
     :      /DUMX/JLIS(149),JC1S(149),JC2S(149)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /L1/JBQ1(3,149),JBQ2(3,149),JTQ1(3),JTQ2(3)
     :      /M0/JJC1(149),JJC2(149)
     :      /M1/NQ1(149),NQ2(149)
     :      /M2/JJQ1(3,149),JJQ2(3,149)
     :      /M3/JLIST(149),KLIST(149),NPEEL,NCORE
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(149),NAK(149)
      COMMON/SUMARG/J6P(MANGMP),J7P(MANGMP),J8P(MANGMP),J9P(MANGMP),
     :              JWORD(6,M6J),NLSUM,NBJ(MSUM),NB6J(MSUM),K6CP(MSUM),
     :              K7CP(MSUM),K8CP(MSUM),K9CP(MSUM),JSUM6(MTRIAD),
     :              JSUM4(MTRIAD,M6J),JSUM5(MTRIAD,M6J),INV6J(M6J)
     :      /TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
*
      IA1 = 0
      KK = KA+KA+1
      IF (ITRIG (ITJPO (JA),ITJPO (JB),KK) .EQ. 0) RETURN
      IF ((IOPAR .NE. 0) .AND. (ISPAR (JA)*ISPAR (JB)*IOPAR .NE. 1))
     :   RETURN
*
      CALL SETQNA (JA,JB)
      IF (IBUG4 .NE. 0) CALL VIJOUT (JA,JB)
*
      DO 1 IJ = 1,NW
         VSHELL(IJ) = ZERO
    1 CONTINUE
*
*   Analyse peel shell interactions
*
      IDQ = 0
      JA1 = 0
      JA2 = 0
*
      IF (NPEEL .NE. 0) THEN
*
        DO 2 JW = 1,NPEEL
          IJ = JLIST(JW)
          NDQ = NQ1(IJ)-NQ2(IJ)
          IF (ABS (NDQ) .GT. 1) GOTO 39
          IF (NDQ .GT. 0) THEN
            JA1 = JW
            IDQ = IDQ+1
          ELSEIF (NDQ .LT. 0) THEN
            JA2 = JW
            IDQ = IDQ+1
          ENDIF
    2   CONTINUE
*
        IF (IDQ .GT. 2) GOTO 39
*
*   Evaluate the array VSHELL
*
*   Then there are two possibilities IDQ = 0 or IDQ = 2
*   if IDQ = 0, then loop over all shells by index ISH
*   if IDQ = 2, then one orbital fixed on each side
*
        NS = NPEEL
      ENDIF
*
      IF (IDQ .EQ. 2) GOTO 19
*
*   Loop over shells when IDQ = 0
*
      ISH = 0
      IF (NPEEL .EQ. 0) GOTO 9
      DO 7 I = 1,NPEEL
    7    JLIS(I) = JLIST(I)
      IF (NPEEL .EQ. 1) GOTO 9
      NPEELM = NPEEL-1
      DO 8 I = 1,NPEELM
         JC1S(I) = JJC1(I)
    8    JC2S(I) = JJC2(I)
*
*   If ISH .GT. NW, then loop is over and return
*
    9 ISH = ISH+1
      IF (ISH .GT. NW) RETURN
      IF (ICHOP (ISH,JA) .EQ. -1) GOTO 9
      IF (IBUG6 .NE. 0) WRITE (99,308) ISH
      IF (ICHOP (ISH,JA) .EQ. 0) GOTO 16
*
*   Case one --- the ISH-th shell is in the core or in the peel and
*   closed for both sides
*
      I = 1
      IF (NPEEL.EQ.0) GOTO 15
      DO 10 I = 1,NPEEL
        IJ = JLIST(I)
        IF (ISH .LT. IJ) GOTO 11
   10 CONTINUE
      I = NPEEL+1
      GOTO 13
   11 IM = NPEEL-I+1
      DO 12 II = 1,IM
         JLIST(NPEEL+2-II) = JLIST(NPEEL+1-II)
         IF (NPEEL.EQ.II) GOTO 13
         JJC1(NPEEL+1-II) = JJC1(NPEEL-II)
         JJC2(NPEEL+1-II) = JJC2(NPEEL-II)
   12 CONTINUE
   13 CONTINUE
      IF (I .LT. 3) GOTO 14
      JJC1(I-1) = JJC1(I-2)
      JJC2(I-1) = JJC2(I-2)
      GOTO 15
   14 I1 = JLIST(1)
      JJC1(1) = JJQ1(3,I1)
      JJC2(1) = JJQ2(3,I1)
   15 JLIST(I) = ISH
      JA1 = I
      JA2 = I
      NS = NPEEL+1
      GOTO 19
*
*   Case two --- the ISH-th shell is in the peel and open for either
*   side
*
   16 NS = NPEEL
      DO 17 JW = 1,NPEEL
        NX = ISH-JLIST(JW)
        IF (NX.EQ.0) GOTO 18
   17 CONTINUE
   18 JA1 = JW
      JA2 = JW
*
*   Main computation
*
*     JA1, JA2 are the indices of interacting shells in JLIST
*     IA1, IA2 are the indices of interacting shells in NW
*
   19 IA1 = JLIST(JA1)
      IA2 = JLIST(JA2)
      KS1 = 2*ABS (NAK(IA1))
      KS2 = 2*ABS (NAK(IA2))
*
*   Check triangular condition for the active shells
*
      IF (ITRIG (KS1,KS2,KK).EQ.1) GOTO 99
      IF (IDQ .EQ. 2) RETURN
      GOTO 100
*
*   Set tables of quantum numbers of non-interacting spectator shells
*
   99 CONTINUE
*
      DO 26 JW = 1,NS
        IJ = JLIST(JW)
        IF (IJ .EQ. IA1) GOTO 23
        DO 22 K = 1,3
          JBQ1(K,IJ) = JJQ1(K,IJ)
   22   CONTINUE
*
   23   IF (IJ .EQ. IA2) GOTO 26
        DO 24 K = 1,3
          JBQ2(K,IJ) = JJQ2(K,IJ)
   24   CONTINUE
        IF ((IJ .EQ. IA1) .OR. (IJ .EQ. IA2)) GOTO 26
        DO 25 K = 1,3
          IF (JBQ1(K,IJ) .NE. JBQ2(K,IJ)) GOTO 40
   25   CONTINUE
   26 CONTINUE
*
*   Loop over parent states
*
      IS(1) = IA1
      IS(2) = IA2
      KS(1) = KS1
      KS(2) = KS2
      VAL = ZERO
      KJ23 = 0
      IX = 0
      FAIL = .FALSE.
*
      NELCTS = NQ2(IA2)
      L2 = IROW1(NELCTS,KS2)
      LLS2 = ITAB(L2)
      LS2 = JTAB(L2)
*
      DO 34 LB = 1,LLS2
        LS2 = LS2+3
        IT1 = NTAB(LS2)
        IT2 = KS2
        IT3 = JJQ2(3,IA2)
        IF (ITRIG (IT1,IT2,IT3) .EQ. 0) GOTO 34
        IF (ABS (NTAB(LS2-2)-JJQ2(1,IA2)) .NE. 1) GOTO 34
        DO 27 K = 1,3
          JBQ2(K,IA2) = NTAB(LS2+K-3)
   27   CONTINUE
*
        NELCTS = NQ1(IA1)
        L1 = IROW1(NELCTS,KS1)
        LLS1 = ITAB(L1)
        LS1 = JTAB(L1)
*
        DO 33 LA = 1,LLS1
          LS1 = LS1+3
          IT1 = NTAB(LS1)
          IT2 = KS1
          IT3 = JJQ1(3,IA1)
          IF (ITRIG (IT1,IT2,IT3) .EQ. 0) GOTO 33
          IF (ABS (NTAB(LS1-2)-JJQ1(1,IA1)) .NE. 1) GOTO 33
          DO 28 K = 1,3
           JBQ1(K,IA1) = NTAB(LS1+K-3)
   28     CONTINUE
*
          DO 20 K = 1,3
            IF (JBQ1(K,IA1) .NE. JBQ2(K,IA1)) GOTO 33
            IF (JBQ1(K,IA2) .NE. JBQ2(K,IA2)) GOTO 33
   20     CONTINUE
*
*   Parent shells now defined
*
          CALL FIXJ (JA1,JA2,KA,IS,KS,NS,KJ23)
          KJ23 = 1
*
          IF (IBUG6 .NE. 0) THEN
            MN1 = MJA
            NS1 = NJA-1
            WRITE (99,302)
            WRITE (99,303) (J1(J),J = 1,MN1)
            WRITE (99,304)
            DO 30 JW = 1,NS1
               WRITE (99,305) (J2(JW,K),K = 1,3),(J3(JW,K),K = 1,3)
   30       CONTINUE
          ENDIF
*
*   Evaluate recoupling coefficient
*
          IF (IX .EQ. 0) THEN
            DO 500 I = 1,MJA
              FREE(I) = .FALSE.
  500       CONTINUE
            IF (LLS2 .NE. 1) FREE(JA2) = .TRUE.
            CALL NJGRAF (RECUPS,-1,FAIL)
            IX = 1
            IF (FAIL) GOTO 501
          ENDIF
          CALL GENSUM (J6C,J7C,J8C,J9C,JWC,J6,J7,J8,J9,KW,JDEL,
     :                LDEL,SUMVAR,MP,
     :                J6P,J7P,J8P,J9P,JWORD,NLSUM,NBJ,NB6J,
     :                K6CP,K7CP,K8CP,K9CP,JSUM4,JSUM5,JSUM6,INV6J,
     :                RECUPS)
          IF (IBUG6 .NE. 0) WRITE (99,307) RECUPS
          IF (ABS(RECUPS) .LT. EPS) GOTO 33
*
*   Evaluates 2 CFPs
*
          IF (KS1 .EQ. 2) GOTO 31
          II = IA1
          NEL = NQ1(II)
          IVP = JBQ1(1,II)
          IWP = JBQ1(2,II)
          IJP = JBQ1(3,II)-1
          IVD = JJQ1(1,II)
          IWD = JJQ1(2,II)
          IJD = JJQ1(3,II)-1
          CALL CFP (KS1,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C)
          IF (IBUG6 .NE. 0) WRITE (99,306) KS1,NEL,IJD,IVD,IWD,IJP,IVP,
     :              IWP,C
          IF (ABS(C) .LT. EPS) GOTO 33
          RECUPS = RECUPS*C
*
31        IF (KS2 .EQ. 2) GOTO 32
          II = IA2
          NEL = NQ2(II)
          IVD = JJQ2(1,II)
          IWD = JJQ2(2,II)
          IJD = JJQ2(3,II)-1
          IVP = JBQ2(1,II)
          IWP = JBQ2(2,II)
          IJP = JBQ2(3,II)-1
          CALL CFP (KS2,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C)
          IF (IBUG6 .NE. 0) WRITE (99,306) KS2,NEL,IJD,IVD,IWD,
     :                                              IJP,IVP,IWP,C
          IF (ABS(C) .LT. EPS) GOTO 33
          RECUPS = RECUPS*C
*
   32     CONTINUE
          VAL = VAL+RECUPS
   33   CONTINUE
   34 CONTINUE
*
*   End of loop over parent states
*
  501 IF (IDQ .EQ. 2) GOTO 37
*
*   IDQ = 0 CASE
*
      VSHELL(ISH) = VAL*DBLE (NQ1(IA1))
*
*   Loop over all shells when IDQ = 0
*
100   CONTINUE
      IF (NPEEL .EQ. 0) GOTO 9
      DO 35 I = 1,NPEEL
   35    JLIST(I) = JLIS(I)
      IF (NPEEL .EQ. 1) GOTO 9
      NPEELM = NPEEL-1
      DO 36 I = 1,NPEELM
         JJC1(I)  = JC1S(I)
   36    JJC2(I)  = JC2S(I)
      GOTO 9
*
*   IDQ = 2 Case
*
*       Permutation factor for IDQ = 2
*
   37 CONTINUE
      VAL = VAL*SQRT (DBLE (NQ1(IA1)*NQ2(IA2)))
      LLD1 = MIN (IA1,IA2)+1
      LLD2 = MAX (IA1,IA2)
      IDL = 1
      IF (IA1 .LT. IA2) IDL = 0
      DO 38 K = LLD1,LLD2
        IDL = IDL+NQ1(K)
   38 CONTINUE
      IF (MOD(IDL,2) .NE. 0) VAL = -VAL
      VSHELL(1) = VAL
      RETURN
*
   39 IF (IBUG6 .NE. 0) WRITE (99,300)
      RETURN
   40 IF (IBUG6 .NE. 0) WRITE (99,301)
      RETURN
*
  300 FORMAT (' One side has more than one interacting electron')
  301 FORMAT (' Spectator quantum numbers not diagonal for non-interact'
     :   ,'ing shells')
  302 FORMAT (/' J1')
  303 FORMAT (24I5)
  304 FORMAT (' J2                   J3')
  305 FORMAT (3I5,I10,2I5)
  306 FORMAT(' CFP  ',I3,I4,I7,2I4,I7,2I4,1P,D20.9)
  307 FORMAT(/' Recoupling coefficient = ',1P,D19.12)
  308 FORMAT(//' ISH = ',I3)
*
      END
************************************************************************
*                                                                      *
      SUBROUTINE VIJOUT (JA,JB)
*                                                                      *
*   Prints  out tables of configurational quantum numbers defined by   *
*   SETQNA for current matrix element.                                 *
*                                                                      *
*                                           Last update: 14 Oct 1992   *
*                                                                      *
************************************************************************
*
      CHARACTER*2 IC,IJ
*
      DIMENSION IJ(2),JC(2),IC(2)
*
      COMMON/M0/JJC1(149),JJC2(149)
     :      /M1/NQ1(149),NQ2(149)
     :      /M2/JJQ1(3,149),JJQ2(3,149)
     :      /M3/JLIST(149),KLIST(149),NPEEL,NCORE
*
      DATA IJ/'/2','  '/
*
      IF (NPEEL .GT. 0) THEN
*
*   Identify CSFs
*
         WRITE (99,300) JA,JB
*
*   Print active shell quantum numbers from JLIST table
*
         WRITE (99,301)
*
         DO 1 J = 1,NPEEL
*
            JW = JLIST(J)
*
            JC(1) = JJQ1(3,JW)-1
            IF (MOD (JC(1),2) .EQ. 1) THEN
               IC(1) = IJ(1)
            ELSE
               JC(1) = JC(1)/2
               IC(1) = IJ(2)
            ENDIF
*
            JC(2) = JJQ2(3,JW)-1
*
            IF (MOD (JC(2),2) .EQ. 1) THEN
               IC(2) = IJ(1)
            ELSE
               JC(2) = JC(2)/2
               IC(2) = IJ(2)
            ENDIF
*
            WRITE (99,302) JW,NQ1(JW),JJQ1(1,JW),JJQ1(2,JW),JC(1),IC(1)
     :                       ,NQ2(JW),JJQ2(1,JW),JJQ2(2,JW),JC(2),IC(2)
    1    CONTINUE
*
*   Print coupling angular momenta if NPEEL .GE. 2
*
         IF (NPEEL .GT. 2) THEN
*
            WRITE (99,303)
            DO 2 J = 2,NPEEL
*
               JC(1) = JJC1(J-1)-1
*
               IF (MOD (JC(1),2) .EQ. 1) THEN
                  IC(1) = IJ(1)
               ELSE
                  JC(1) = JC(1)/2
                  IC(1) = IJ(2)
               ENDIF
*
               JC(2) = JJC2(J-1)-1
               IF (MOD (JC(2),2) .EQ. 1) THEN
                  IC(2) = IJ(1)
               ELSE
                  JC(2) = JC(2)/2
                  IC(2) = IJ(2)
               ENDIF
*
               WRITE (99,304) (JC(I),IC(I),I = 1,2)
    2       CONTINUE
         ENDIF
*
      ENDIF
*
      WRITE (99,305) NCORE
*
      RETURN
*
  300 FORMAT (/'From VIJOUT: CSF ',1I2,35X,'CSF ',1I2)
  301 FORMAT (3X,'subshell',4X,'q',4X,'v',2X,'w',2X,'J',
     :                     19X,'q',4X,'v',2X,'w',2X,'J')
  302 FORMAT (7X,I3,I6,I5,2I3,A2,15X,I3,I5,2I3,A2)
  303 FORMAT (' coupling schemes:')
  304 FORMAT (14X,I2,A2,27X,I2,A2)
  305 FORMAT (' there are ',I3,' inactive closed shells.'/)
*
      END
