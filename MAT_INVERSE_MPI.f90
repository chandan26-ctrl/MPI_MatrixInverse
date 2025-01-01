!****************************************************************************
! *   FILE         = MAT_INVERSE_MPI.F90                                    *
! *   DESCRIPTION  = LU DECOMPSITION METHOD                                 *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *
!****************************************************************************
      PROGRAM MAT_INVERSE_MPI
      IMPLICIT NONE
      INCLUDE 'mpif.h'

      INTEGER, PARAMETER :: N =10, MASTER = 0, FROM_MASTER = 1, FROM_SLAVE = 2
      INTEGER :: NNODES,MY_ID,NSLAVES,SOURCE,NP,MSGTYPE
      INTEGER :: COLS,NCOL,EXTRA, OFFSET,I,J,K,T,LED,IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      REAL :: A(N,N),L(N,N),U(N,N),D(N,N), COEFF ,Y(N),SUMM

      CALL MPI_INIT( IERR )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, MY_ID, IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NNODES, IERR )
      NSLAVES = NNODES-1
      
! *************************** INITIALIZATION IN MASTER PROCESSOR *************************************
      IF (MY_ID .EQ. MASTER) THEN

      L=0.0
      U=0.0
      D=0.0

!     READ MAT A AND B 
      OPEN(1,FILE='amatrix.dat',STATUS='OLD')


       DO I=1,N
        DO J=1,N
          READ(1,*)A(I,J)
        END DO
        END DO

      CLOSE(1)


!     LU DECOMPOSITION OF MAT A
      DO K=1,N-1
        DO I=K+1,N
         COEFF=A(I,K)/A(K,K)
         L(I,K) = COEFF
          DO J=K+1,N
            A(I,J) = A(I,J)-COEFF*A(K,J)
          END DO
        END DO
      END DO

      DO I=1,N
        L(I,I) = 1.0
        D(I,I) = 1.0
      END DO

      DO J=1,N
        DO I=1,J
          U(I,J) = A(I,J)
        END DO
      END DO


!     SEND MATRIX DATA TO SLAVE PROCESSORS
        NCOL = N/NSLAVES
        EXTRA = MOD(N, NSLAVES)
        OFFSET = 1
        MSGTYPE = FROM_MASTER
      DO NP=1, NSLAVES
        IF (NP .LE. EXTRA) THEN
             COLS=NCOL+1
        ELSE
             COLS=NCOL
        END IF
             WRITE(*,*)'   SENDING',COLS,' COLS TO PROCESSOR',NP
        CALL MPI_SEND (OFFSET,1,  MPI_INTEGER, NP, MSGTYPE, MPI_COMM_WORLD, IERR)
        CALL MPI_SEND (COLS,  1,  MPI_INTEGER, NP, MSGTYPE, MPI_COMM_WORLD, IERR)
        CALL MPI_SEND( L,   N*N,  MPI_REAL,    NP, MSGTYPE, MPI_COMM_WORLD, IERR)
        CALL MPI_SEND( U,   N*N,  MPI_REAL,    NP, MSGTYPE, MPI_COMM_WORLD, IERR)
        CALL MPI_SEND (D(1,OFFSET),COLS*N, MPI_REAL, NP, MSGTYPE, MPI_COMM_WORLD, IERR)

          OFFSET=OFFSET+COLS

      END DO

!     RECEIVE RESULTS FROM SLAVE PROCESSORS
        MSGTYPE = FROM_SLAVE
      DO I=1, NSLAVES
          SOURCE = I
        CALL MPI_RECV( OFFSET, 1, MPI_INTEGER, SOURCE, MSGTYPE, MPI_COMM_WORLD, STATUS, IERR)
        CALL MPI_RECV( COLS, 1, MPI_INTEGER, SOURCE, MSGTYPE, MPI_COMM_WORLD, STATUS, IERR)
        CALL MPI_RECV( D(1,OFFSET), COLS*N, MPI_REAL,SOURCE, MSGTYPE, MPI_COMM_WORLD, STATUS, IERR)

      END DO

!     PRINT RESULT

      DO I=1,N
        PRINT*, (D(I,J), J = 1, N)
      END DO

      END IF

    

! *************************** COMPUTATION IN SLAVE PROCESSORS*************************************
      IF (MY_ID > MASTER) THEN
!     RECEIVE MATRIX DATA FROM THE MASTER PROCESSOR
        MSGTYPE = FROM_MASTER
        CALL MPI_RECV( OFFSET, 1, MPI_INTEGER, MASTER, MSGTYPE, MPI_COMM_WORLD, STATUS, IERR )
        CALL MPI_RECV( COLS, 1, MPI_INTEGER, MASTER, MSGTYPE, MPI_COMM_WORLD, STATUS, IERR ) 
        CALL MPI_RECV( L, N*N, MPI_REAL, MASTER, MSGTYPE, MPI_COMM_WORLD, STATUS, IERR )
        CALL MPI_RECV( U, N*N, MPI_REAL, MASTER, MSGTYPE, MPI_COMM_WORLD, STATUS, IERR )
        CALL MPI_RECV( D, COLS*N, MPI_REAL, MASTER, MSGTYPE, MPI_COMM_WORLD, STATUS, IERR )

      
      DO T=1, COLS
        DO I=1,N
           Y(I)=0.0
        END DO

!  SOLVE LY=D USING THE FORWARD SUBSTITUTION
        DO I=1,N
           SUMM =0.0
          DO J=1,N
             SUMM=SUMM-L(I,J)*Y(J)
          END DO
             Y(I)=(D(I,T)+SUMM)/L(I,I)
        END DO

! SOLVE UD=Y USING THE BACK SUBSTITUTION
        DO I=1,N
           D(I,T)=0.0
        END DO
        DO I=N,1,-1
           SUMM=0.0
          DO J=N,1,-1
             SUMM= SUMM-U(I,J)*D(J,T)
          END DO
             D(I,T)=(Y(I)+SUMM)/U(I,I)
        END DO

      END DO
      
!    SEND RESULTS BACK TO THE MASTER PROCESSOR
        MSGTYPE = FROM_SLAVE
        CALL MPI_SEND( OFFSET, 1, MPI_INTEGER, MASTER, MSGTYPE, MPI_COMM_WORLD, IERR)
        CALL MPI_SEND( COLS, 1, MPI_INTEGER, MASTER, MSGTYPE,  MPI_COMM_WORLD, IERR )
        CALL MPI_SEND( D, COLS*N, MPI_REAL, MASTER, MSGTYPE, MPI_COMM_WORLD, IERR )
      END IF
      CALL MPI_FINALIZE(IERR)
      END

