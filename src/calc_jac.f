      PROGRAM calc_enthalpy
C
C*****DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END DOUBLE PRECISION
C*****SINGLE PRECISION
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END SINGLE PRECISION
      PARAMETER ( LENIWK = 40000, LENRWK = 400000, LENCWK = 2000,
     1            LENSYM = 16)
      DIMENSION IWORK (LENIWK), RWORK (LENRWK)
      REAL(8), ALLOCATABLE, DIMENSION(:) :: X
      REAL(8), ALLOCATABLE, DIMENSION(:) :: Y
      REAL(8), ALLOCATABLE, DIMENSION(:) :: Q
      CHARACTER CWORK(LENCWK)*(LENSYM)
      DATA LIN/5/, LOUT/6/, LINKCK/25/, LSAVE/7/, LIGN/9/, LREST/10/
C
C     LIN    = Unit number for Keyword input
C     LOUT   = Unit number for text output to terminal
C     LIGN   = Unit number for text output file
C     LSAVE  = Unit number for binary output file
C     LINKCK = Unit number for CHENKIN linking file
C     LREST  = Unit number for binary restart file
C     LENIWK = Length of integer work array
C     LENRWK = Length of real work array
C     LENCWK = Length of character work array
C     LENSYM = Length of a character string in character work array
C     IWORK  = Integer work array
C     RWORK  = Real work array
C     CWORK  = Character work array
C
C*****vms
C      SET I/O UNITS AND OPEN FILES. OPERATING SYSTEM IS vms VMS.
C      OPEN (LINKCK, STATUS='OLD', FORM='UNFORMATTED')
C      OPEN (LSAVE, STATUS='NEW', FORM='UNFORMATTED')
C      OPEN (LOUT, STATUS='NEW', FORM='FORMATTED')
C      OPEN (LIGN, STATUS='NEW', FORM='FORMATTED')
C      OPEN (LIN, STATUS='OLD', FORM='FORMATTED')
C      INQUIRE (FILE='restart', EXIST=LEXIST)
C      IF (LEXIST) OPEN (LREST,STATUS='OLD',FORM='UNFORMATTED')
C*****END vms
C
C*****unix
      OPEN (LINKCK, FORM='UNFORMATTED', FILE='input/cklink')
      OPEN (LIN, FORM='FORMATTED', FILE='input/datasheet')
      OPEN (LIGN, FORM='FORMATTED', FILE = 'output/rop')
      OPEN (LOUT, FORM='FORMATTED', FILE='output/terminalout')
C*****END unix
C
C     INITIALIZE WORK ARRAY BY inp FILE
C
      CALL CKINIT(LENIWK, LENRWK, LENCWK, LINKCK, LOUT,
     1            IWORK, RWORK, CWORK)
      CALL CKINDX(IWORK, RWORK, MM, KK, II, NFIT)
      ALLOCATE (X(KK), Y(KK), Q(II))

      WRITE (LIGN, *) "Time(s)  HBMS(kJ/g)  HBML(kJ/mole)"
C
C     READ EACH ROW IN DATASHEET
C
      DO 
          READ (LIN, *, END=999) TIME, P, T, (X(I), I = 1, KK)
C
C     CONVERT X TO Y
C
          CALL CKXTY  (X, IWORK, RWORK, Y)
C
C     CALCULATE ENTHALPY
C
          CALL CKHBMS (T, Y, IWORK, RWORK, HBMS)   ! mass units
          CALL CKHBML (T, X, IWORK, RWORK, HBML) ! molar units
C
C     CALCULATE ROP
C
          CALL CKQYP  (P, T, Y, IWORK, RWORK, Q)
C
C     PRINT OUT ENTHALPY
C
          WRITE (LIGN, *) TIME, (Q(I), I = 1, II)
      END DO
C     
999   CONTINUE
      STOP
      END