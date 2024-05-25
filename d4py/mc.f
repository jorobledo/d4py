C   BUGS CORREGIDOS EN LA RUTINA ANGUEMP (ERROR AL CAMBIAR EL PASO DE INTEGRACION
c   Y EN LA LLAMADA) 27/7/07
C   Relaci-n de los espectros que se registran:
C   ESPS(I,J,K):     Espectro acumulante saliente para la energ-a I--sima
C                    el orden de scattering J-simo, y el Q K-simo. El
C                    ESPS(I,10,K) contiene los -rdenes de scattering 10 o
C                    mayor, y el ESPS(I,11,K), la suma de K=2...10 (el
C                    m-ltiple total).
C   EPSPS2(I,J,K):   Acumula los cuadrados del anterior.
C   RESPS(I,J,K):    Valor medio del ESPS.
C   REPSPS2(I,J,K):  Error de ESPS.
C   ESP1C(I,J):      Scattering simple en el can para la energ-a I--sima
C                    y el Q K-simo.
C   EPSP1C2(I,J,K):  Acumula los cuadrados del anterior.
C   RESP1C(I,J):     Valor medio de ESP1C
C   REPSP1C2(I,J,K): Error de ESP1C
C   RESPT(I,J):      Scattering total para la energ-a I--sima y el Q K-simo.
C   REPSPT2(I,J,K):  Acumula los cuadrados del anterior.
C   ESP1NAS(I,J):    Scattering simple sin atenuaci-n.
C   RESP1NAS(I,J):   Valor medio del anterior.


      PROGRAM MCD48

      implicit none
      character*30 name
      integer ::  j, N, i
      integer, allocatable :: dseed(:)
      real(8), dimension(:), allocatable :: a1, a2, a3
      real(8), dimension(:,:), allocatable :: a4
      real(8) :: Uout, RAND1, RAND2
      real(8), dimension(1000):: RANDOM
      name = 'prueba'

      call random_init(.true., .true.)
      do N=1,100
      call RANDOM_NUMBER(RAND1)
      call RANDOM_NUMBER(RAND2)
      call INELASTICO_GENERALIZADO(10.0D0, 0.1D0, RAND1, RAND2, Uout)
      end do
      END

C---------------------------------------------------------------------
      ! MODULE MC

      ! private
      ! public :: MC2022, BOX_MULLER

C --------------------------------------------------------------------------

      SUBROUTINE MC2022(OUTPUT,
     #NCIC,NT,ANGDET0,DA,R1,R2,H,DH,AH,E0,WCO,
     #UTAB,SETOT,SEM,PSCAT,SETOTC,SEMC,PSCATC, IMOD,
     #ANGDET, RESPS, RESP1NAS,RESP1C, RESPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*30, intent(in) :: OUTPUT
      integer, intent(in) :: NCIC,NT
      INTEGER :: IRND2
cf2py integer optional, intent(in) :: ncic=5,nt=1000
      real(8), intent(in)::ANGDET0,DA,R1,R2,H,DH,AH,E0,WCO
      real(8), intent(in)::UTAB,SETOT,SEM,PSCAT,SETOTC,SEMC,PSCATC
cf2py real(8) optional,intent(in)::ANGDET0=0,DA=0.88,R1=0.3,R2=0.29
cf2py real(8) optional,intent(in)::H=6.0,DH=0.1,AH=5.0,E0=0.3243,WCO=1e-06
cf2py real(8) optional,intent(in)::UTAB=-0.482,SETOT=4.03,SEM=0.30,PSCAT=0.96,SETOTC=4.89,SEMC=0.37,PSCATC=0.78
      parameter (NQ=200, NDQE=400)
      INTEGER, intent(in) :: IMOD
C----- IMOD selecciona el modelo. 1 es elástico, 4 es inelástico generalizado.
cf2py INTEGER, optional, intent(in) :: IMOD=1
      REAL(8), dimension(NQ), intent(out) :: ANGDET,RESP1C,
     #RESP1NAS,RESPT
      REAL(8), dimension(NQ,11), intent(out) :: RESPS
C----- NQ es el número de puntos de Q deseados
C----- NDQE es el número de puntos de Mcd8dq.dat
C----- NINP es el número de puntos de Mcd8inp.dat
      CHARACTER*30 HORA,HORA1*8,HORA2*8
      DIMENSION POS(2,3),D(2,3),A(2,2),RANDOM(1000),
     #ESPS(NQ,11),ESPS2(NQ,11),RESPS2(NQ,11),
     #ESP1C(NQ),ESP1C2(NQ),RESP1C2(NQ),
     #RESPTSC(NQ),RESPT2(NQ),ESP1NAS(NQ),
     #ESPSINC(NQ),RESPSINC(NQ),NEVENT(10,2),QPRO(NQ),RQPRO(NQ),
     #QEMP(NDQE),DQCAN(NDQE),DQMUE(NDQE),DANGEMP(181),Q(NQ)
      COMMON /DATOSN/EIN,U,THETA,ALFA,DIS(8)
      COMMON /DATOSSAM/SE,SE0,SM,PS
      COMMON /DATOSCAN/SEC,SE0C,SMC,PSC
      COMMON /NEUTRON/VDIR(3),VPOS(3)
      COMMON /GEOMUEST/RADIO(4)
      COMMON /DATOSC/AC,BKTAUC,SIGBC
      COMMON /TEMPERATURA/XKBT
      DATA PI/3.1415926535897932D0/,HQS2M/2.07219093D-3/

      RADIO(1)=R1
      RADIO(2)=R2
      RADIO(3)=0.D0
      RADIO(4)=0.D0
C---------------------------------------------------------------------
C     Para g77 Linux
11    j=time()
      call srand(j)
      dseed=10000*rand()

C     Generador de semilla Lahey Fortran
c11    DSEED=RRAND()*10000.
C---------------------------------------------------------------------
C     Para LINUX ALPHA
c11    CALL RANDOM_SEED
c      CALL RANDOM_NUMBER(DSEED)
c      DSEED=DSEED*10000
c 11     DSEED=1000

C---------------------------------------------------------------------
C     Para Microsoft Fortran
C11    CALL GETTIM(I,J,K,L)
C      DSEED=L+K*10000./60.
C---------------------------------------------------------------------
C Se definen las variables IMED e ICAN.
C IMED: medio en el que se encuentra el neutron.
C       IMED=0 : neutron en muestra
C       IMED=1 : neutron en can
C ICAN: tipo de neutron
C       ICAN=0 : neutron que solo interactuo dentro de la muestra
C       ICAN=1 : neutron que interactuo por lo menos una vez en el can
C---------------------------------------------------------------------

      IF (DSEED.LE.0.D0) GOTO 11
C---------------------------------------------------------------------
C Inicializa espectros acumulantes y de salida
      DO 10 J=1,NQ
           ESP1C(J)=0.D0
           ESP1C2(J)=0.D0
           RESP1C(J)=0.D0
           RESP1C2(J)=0.D0
           RESPT(J)=0.D0
           RESPT2(J)=0.D0
           ESP1NAS(J)=0.D0
           RESP1NAS(J)=0.D0
           ESPSINC(J)=0.D0
           RESPSINC(J)=0.D0
           QPRO(J)=0.D0
           DO 10 K=1,11
             ESPS(J,K)=0.D0
             ESPS2(J,K)=0.D0
             RESPS(J,K)=0.D0
10           RESPS2(J,K)=0.D0

C Inicializa vector NEVENT: numero de eventos en funcion del numero
C de scatterings J y del tipo de neutron (de can o de muestra)

      DO 12 J=1,10
        DO 12 K=1,2
12        NEVENT(J,K)=0


      NTOT=0
      NOUT=0
      NHIS=0
      TRANSTOT=0.D0

C Lectura de secciones eficaces totales de muestra y can, microscopicas
C y macroscopicas, y probabilidades de scattering (frente a absorcion).
C Lectura de la distribucion angular de los neutrones dispersados por el can

c     Aca habria que incluir la nueva XS de muestra y can. Se podría
c     remplazar por un valor constante y ver agregarlo como parámetro.
!       OPEN(1,FILE='Mcd8inp.dat',STATUS='OLD')
!       READ(1,*)
!       DO 2 I=1,NINP
! 2     READ(1,*)UTAB(I),SETOT(I),SEM(I),PSCAT(I),SETOTC(I),SEMC(I),
!      #PSCATC(I)
!       CLOSE(1)
      OPEN(1,FILE='generated_files/Mcd8dq.dat',STATUS='OLD')
      READ(1,*)
      READ(1,997)DQMA0,DQMA1,DQMA2
      READ(1,997)DQCA0,DQCA1,DQCA2
      DO 3 I=1,NDQE
3     READ(1,*)QEMP(I),DQMUE(I),DQCAN(I)
      CLOSE(1)
997   format(1x,3(1pe14.7,1x))


C Lee parametros generales

      EIN=E0
      U=UTAB
      XKBT=0.0253D0
      AC=50.94D0
      BKTAUC=0.0253
      SIGBC=4.99D0

C Inicializa vector de Angulos
      DO 13 J=1,NQ
      ANGDET(J)=ANGDET0+(J-1)*DA
13    Q(J)=2.D0*DSQRT(E0/2.07219093D-3)*DSIN(ANGDET(J)/2*PI/180.D0)

C Inicializa indice del vector de numeros random, y llama a la rutina random

      KRAND=0
      HORA1='       '
      HORA2='       '

      CALL GGUBFS(DSEED,RANDOM)
      hora=ctime(time())
      HORA1=HORA(12:19)
      HORA2=HORA1

C Comienza el gran loop de NCIC ciclos de NT neutrones
      call random_init(.true., .true.)

      DO 1000 NUMCIC=1,NCIC
      DO 200 N=1,NT
         NTOT=NTOT+1
         KRAND=KRAND+1
         IF (KRAND.GT.1000) THEN          ! SI SE ACABAN LOS RANDOM GENERA
          CALL GGUBFS(DSEED,RANDOM)       ! NUEVOS
          KRAND=1
         END IF
      IRAND1=MOD(KRAND,1000)+1
      IRAND=MOD(IRAND1,1000)+1
      CALL RCERO(RADIO(1),AH,POS,A,RANDOM(IRAND1),RANDOM(IRAND))
      D(1,1)=DSIN(A(2,1)/180.D0*PI)*DCOS(A(2,2)/180.D0*PI)
      D(1,2)=DSIN(A(2,1)/180.D0*PI)*DSIN(A(2,2)/180.D0*PI)
      D(1,3)=DCOS(A(2,1)/180.D0*PI)
      DO 20,I=1,3
20    D(2,I)=D(1,I)

      NS=1                             ! NUMERO DE SCATTERINGS DEL NEUTRON
      WGT=1.D0
      ICAN=0

c     Defino los parámetros iniciales

      X=UTAB
      SE=SETOT
      SM=SEM
      PS=PSCAT
      SEX=SETOTC
      SMC=SEMC
      PSC=PSCATC
      SE0=SE
      SE0C=SEC

C Comienza loop de scattering multiple

100      CONTINUE
            IRAND=MOD(IRAND,1000)+1
C Inicializa vectores direccion y posicion del neutron (auxiliares)
            VDIR(1)=D(2,1)
            VDIR(2)=D(2,2)
            VDIR(3)=D(2,3)
            VPOS(1)=POS(1,1)
            VPOS(2)=POS(1,2)
            VPOS(3)=POS(1,3)
C Calcula distancias de vuelo en distintos medios segun posicion y direccion
C actual
            CALL DISTAN(1,H,DH)
C Sortea la distancia de vuelo y actualiza el peso
            CALL LAMBDA(XL,RANDOM(IRAND),WGT,ICAN,IMED,TRANS)
            IF (NS.EQ.1) THEN
                NHIS=NHIS+1
                TRANSTOT=TRANSTOT+TRANS
            END IF

            IRAND=MOD(IRAND,1000)+1
C Llama a la ruleta rusa. Si el neutron muere llama uno nuevo
            CALL RULRUS(RANDOM(IRAND),WGT,WCO,IFLAG)
            IF (IFLAG.EQ.0) GOTO 200
C Actualiza posicion del neutron
            CALL POSNUE(POS,XL,D)

C  Calcula la contribucion 1/S(E0) d2S/dOdE (E,E) exp[-s(E) D]

         VPOS(1)=POS(2,1)
         VPOS(2)=POS(2,2)
         VPOS(3)=POS(2,3)
         VDIR(1)=0.D0                     ! Detectores en el plano yz
         LL=NS
         IF (NS.GT.10) LL=10

C Caso del neutron dentro de la muestra  ---------------------------+
      IF (IMED.EQ.0) THEN
C Genera distribucion ang. de scattering de la muestra para la ener.
         IF (IMOD.EQ.1) THEN
          CALL GENDANG(U,QEMP,DQMUE,DQMA0,DQMA1,DQMA2,DANGEMP)
          CALL NORMDQ(U,QEMP,DQMUE,DQMA0,DQMA1,DQMA2,DQNOR)
         END IF
         IF (ICAN.EQ.0) THEN
           IF (NS.LE.10) NEVENT(LL,1)=NEVENT(LL,1)+1
         ELSE
           IF (NS.LE.10) NEVENT(LL,2)=NEVENT(LL,2)+1
         END IF
      DO I=1,NQ
         ALFA=ANGDET(I)/180.D0*PI           ! Angulo aparente
         VDIR(2)=DSIN(ALFA)                 ! Direcciones al detector
         VDIR(3)=DCOS(ALFA)                 ! del angulo dado
         THETA=DACOS(VDIR(2)*D(2,2)+VDIR(3)*D(2,3))  !Ang. real de int.
         CALL DISTAN(1,H,DH)                   ! Distancias a recorrer
         IF (IMOD.EQ.1)
     +   CALL CONTRANG(1,QEMP,DQMUE,DQNOR,DQMA0,DQMA1,DQMA2,VAL)

         ESPS(I,LL)=ESPS(I,LL)+VAL*WGT
         ESPS2(I,LL)=ESPS2(I,LL)+(VAL*WGT)**2.D0
         IF (LL.EQ.1)  THEN

         IF (IMOD.EQ.1)
     +   CALL CONTRANG(0,QEMP,DQMUE,DQNOR,DQMA0,DQMA1,DQMA2,VAL)
         ESP1NAS(I)=ESP1NAS(I)+VAL*WGT
         END IF
      END DO
      END IF
C Fin del caso del neutron dentro de la muestra --------+
C Caso del neutron en el can             ---------------+
      IF (IMED.EQ.1) THEN
C Genera distribucion angular de scattering del can para la energia dada
         CALL GENDANG(U,QEMP,DQCAN,DQCA0,DQCA1,DQCA2,DANGEMP)
         CALL NORMDQ(U,QEMP,DQCAN,DQCA0,DQCA1,DQCA2,DQNOR)
         IF (NS.LE.10) NEVENT(LL,2)=NEVENT(LL,2)+1

      DO I=1,NQ
         ALFA=ANGDET(I)/180.D0*PI           ! Angulo aparente
         VDIR(2)=DSIN(ALFA)                 ! Direcciones al detector
         VDIR(3)=DCOS(ALFA)                 ! del angulo dado
         THETA=DACOS(VDIR(2)*D(2,2)+VDIR(3)*D(2,3))  !Ang. real de int.
         CALL DISTAN(1,H,DH)                   ! Distancias a recorrer
C         CALL SMPSNU(UINFC,USUPC,VAL,1,IMED)
         CALL CONTRANG(1,QEMP,DQCAN,DQNOR,DQCA0,DQCA1,DQCA2,VAL)
         IF (LL.EQ.1) THEN   !LL es indicador de multiple scattering. Si LL=1 es single scattering
           !write(12,*) VAL, WGT  , LL
           ESP1C(I)=ESP1C(I)+VAL*WGT
           ESP1C2(I)=ESP1C2(I)+(VAL*WGT)**2.D0
         ELSE        ! caso multiple scattering
           ESPS(I,LL)=ESPS(I,LL)+VAL*WGT
           !write(12,*) ESPS(I,LL)+VAL*WGT , LL
           ESPS2(I,LL)=ESPS2(I,LL)+(VAL*WGT)**2.D0
c           ESP1C(I)=ESP1C(I)+VAL*WGT            ! ME INTERESA TAMBIÉN EL MULTIPLE SCATTERING CAN INCLUIDO EN SCATTERING CAN
c           ESP1C2(I)=ESP1C2(I)+(VAL*WGT)**2.D0
         END IF
      END DO
      END IF
C Fin del caso del neutron en el can     ----------+
C Comienza a analizar un nuevo scattering
      NS=NS+1
C Sortea nueva energia y angulo para el caso de neutron en la muestra

      IRAND=MOD(IRAND,1000)+1
      U1=U
      IF(IMED.EQ.0) THEN
            UC=UCENT
            US=USUP
            UI=UINF
      ELSE
            UC=UCENTC
            US=USUPC
            UI=UINFC
      END IF
c            CALL ENERGIA(U,IFLAG,IMED,UC,US,UI,RANDOM(IRAND))


      IF (IFLAG.EQ.0) THEN
            NOUT=NOUT+1
            GOTO 200
      END IF

      IRAND1=MOD(IRAND,1000)+1
      IRAND=MOD(IRAND1,1000)+1
      if (IMOD.EQ.1) then
      CALL ANGUEMP(A,POS,D,DANGEMP,RANDOM(IRAND1),RANDOM(IRAND))
      X=UTAB
      end if
      if (IMOD.EQ.2) then
      call ANGULOS(A,POS,D,U1,U,IMED,RANDOM(IRAND1),RANDOM(IRAND))
      X=U
      end if
      if (IMOD.EQ.4) then
      call RANDOM_NUMBER(RAND1)
      call RANDOM_NUMBER(RAND2)
      call INELASTICO_GENERALIZADO(HQS2M,XKBT, RAND1, RAND2, U) ! updates energy to inelastic energy.
      X=U
      end if
c      CALL VALINTERP(U)
      SE=SETOT
      SM=SEM
      PS=PSCAT
      SEX=SETOTC
      SMC=SEMC
      PSC=PSCATC
      SE0=SE
      SE0C=SEC

         GOTO 100
200    CONTINUE

C Calcula valores medios y errores

      DO 40 K=1,NQ
         IF (NTOT.GT.0) RESP1C(K)=ESP1C(K)/DBLE(NTOT)
         IF (NTOT.GT.0) RESP1NAS(K)=ESP1NAS(K)/DBLE(NTOT)
C         IF (NTOT.GT.0) RESPSINC(K)=ESPSINC(K)/DBLE(NTOT)
         IF (NTOT.GT.0) RQPRO(K)=QPRO(K)/DBLE(NTOT)
         IF (NHIS.GT.0) TPRO=TRANSTOT/NHIS
         IF (NTOT.GT.1)
     #RESP1C2(K)=DSQRT(DABS(ESP1C2(K)/DBLE(NTOT)-(ESP1C(K)/
     #DBLE(NTOT))**2.D0)/DBLE(NTOT-1))
      DO 40 LL=1,10
         NN=NEVENT(LL,1)+NEVENT(LL,2)
         IF (NTOT.GT.0) RESPS(K,LL)=ESPS(K,LL)/DBLE(NTOT)
40       IF (NN.GT.1)
     #RESPS2(K,LL)=DSQRT(DABS(ESPS2(K,LL)/DBLE(NN)-(ESPS(K,LL)/
     #DBLE(NN))**2.D0)/DBLE(NN-1))

C Suma el total del scattering multiple y su error

      DO 42 K=1,NQ
      RESPS(K,11)=0.D0
      RESPS2(K,11)=0.D0
      DO 41 LL=2,10
      RESPS(K,11)=RESPS(K,11)+RESPS(K,LL)
41    RESPS2(K,11)=RESPS2(K,11)+RESPS2(K,LL)**2
42    RESPS2(K,11)=DSQRT(RESPS2(K,11))

C Suma el scattering total y su error

      DO 51 K=1,NQ
      RESPT(K)=RESPS(K,1)+RESP1C(K)+RESPS(K,11) ! RESPS(K,11) ES SCATTERING TOTAL
      RESPTSC(K)=RESPS(K,1)+RESPS(K,11)
51    RESPT2(K)=DSQRT(RESPS2(K,1)**2+RESP1C2(K)**2+
     #            RESPS2(K,11)**2)


C Escribe salidas

      OPEN (1,FILE = trim(OUTPUT) // ".log")
      WRITE (1,*)'NEUTRONES TOTALES=        ',NTOT
      WRITE (1,*)'NEUTRONES FUERA DE RANGO= ',NOUT
      WRITE (1,*)'TRANSMISION PROMEDIO ............ ',TPRO
      WRITE(1,'(10I8)')NEVENT
      CLOSE(1)
      OPEN (1,FILE = trim(OUTPUT) // ".out")

      WRITE(1,998)
998   FORMAT('#    ANGULO     SINGLE       S/ATEN
     #CAN        MULT     IPLE     TOTAL')
      DO 31 K=1,NQ
      FMUL=0.D0
      AT=0.D0
      IF(RESPT(K).NE.0.D0) FMUL=RESPS(K,1)/RESPT(K)
      IF(RESP1NAS(K).NE.0.D0) AT=RESPS(K,1)/RESP1NAS(K)
31    WRITE(1,999)ANGDET(K),RESPS(K,1),RESP1NAS(K),RESP1C(K),
     #            RESPS(K,11),RESPT(K)
      CLOSE(1)

999   FORMAT(1X,12(1PE12.5,1X))

      hora=ctime(time())
      HORA1=HORA(12:19)
      HORA2=HORA1

1000  CONTINUE
      END SUBROUTINE MC2022

C---------------------------------------------------------------------
C---------------------------------------------------------------------
C---------------------------------------------------------------------
c                 COMIENZO DE SUBRUTINAS AUXILIARES
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C---------------------------------------------------------------------

C--------------------------------------------------------------------
C     Calcula la posici-n sobre la cara anterior de la muestra, por
C     donde ingresa el neutr-n, y el -ngulo de entrada del mismo.
C--------------------------------------------------------------------
      SUBROUTINE RCERO(RADIO,EH,R,A,RAND1,RAND2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(2,3),A(2,2)
      A(1,1)=0.D0
      A(1,2)=0.D0
C---- Aqu- se calcula la posici-n de entrada
      R(1,1)=EH*(RAND2-.5D0)
      R(1,2)=RADIO*(2.D0*RAND1-1.D0)
      R(1,3)=-DSQRT(RADIO**2.D0-R(1,2)**2.D0)
C---- Aqu- se calcula el -ngulo de entrada
      A(2,1)=0.D0
      A(2,2)=0.D0
      RETURN
      END

C--------------------------------------------------------------------
C     Calcula la distancia que debe volar el neutr-n desde la posici-n
C  POS, en la direcci-n definida por el versor D(I)
C
C  Caso IEJE=1 (eje del cilindro = x)
C                                  x
C                                 |
C                              ___|___
C                             /   |   \RADIO
C                            |\___T___/|
C               Neutrones    |    |    |
C              ----------->  |    |____|__________
C                            |    /    |          z
C                            |   /     |
C                            | _/_____ |
C                            |//      \|
C                             \_______/
C                            /
C                           y
C  Caso IEJE=2 (eje del cilindro = z)
C                                  x
C                                 |
C                                 |
C                                -|---
C                               / |  /\
C               Neutrones      |  | /  \
C              ----------->    |  || __|__________
C                              |  /|   |          z
C                              | / |   |
C                              |/  |   |
C                              /   |   |
C                             / \___\ /
C                            /
C                           y
C--------------------------------------------------------------------
      SUBROUTINE DISTAN(IEJE,H,DH)
      IMPLICIT REAL *8 (A-H,O-Z)
      IMPLICIT INTEGER *4 (I-N)
      REAL(8) :: H,DH
      DIMENSION ZZ(8)
      COMMON /GEOMUEST/RADIO(4)
      COMMON /NEUTRON/VDIR(3),VPOS(3)
      COMMON /DATOSN/EIN,U,THETA,ALFA,DIS(8)
      DO 1 I=1,8
1     DIS(I)=0.D0

      IF (IEJE.EQ.1) THEN
        VPOSEJE=VPOS(1)
        VPOSR1=VPOS(2)
        VPOSR2=VPOS(3)
        VDIREJE=VDIR(1)
        VDIRR1=VDIR(2)
        VDIRR2=VDIR(3)
      ELSE IF (IEJE.EQ.2) THEN
        VPOSEJE=VPOS(3)
        VPOSR1=VPOS(1)
        VPOSR2=VPOS(2)
        VDIREJE=VDIR(3)
        VDIRR1=VDIR(1)
        VDIRR2=VDIR(2)
c     PREGUNTAR A JAVIER QUE HACER EN CASO DE NINGUN OTRO
      else
        VPOSEJE=VPOS(1)
        VPOSR1=VPOS(2)
        VPOSR2=VPOS(3)
        VDIREJE=VDIR(1)
        VDIRR1=VDIR(2)
        VDIRR2=VDIR(3)
      END IF

C Radio posicion
      RP2=VPOSR1**2+VPOSR2**2
      RP=DSQRT(RP2)

C Caso (muy particular) en el que el neutron viaja en direccion del eje

      IF (DABS(VDIREJE).EQ.1.D0) THEN
        DO I=2,4
         IF (RP.GT.RADIO(I)) GOTO 100
        END DO
100     IF (DABS(VPOSEJE).LE.H/2.D0) THEN
            DIS(10-I)=H/2.D0-VDIREJE*VPOSEJE
            DIS(8)=DH
        ELSE
            IF(VDIREJE*VPOSEJE)101,103,102
101               DIS(8)=DH
                  DIS(2)=DABS(-VDIREJE*VPOSEJE-H/2.D0)
                  DIS(10-I)=H
                  GOTO 103
102               DIS(8)=H/2.D0+DH-VDIREJE*VPOSEJE
103               CONTINUE
        END IF
      RETURN
      END IF

C (Coseno del) angulo posicion-direccion, en el plano perpendicular al eje
      COSALF=1.D0
      IF (RP.GT.0.D0)
     # COSALF=(VDIRR1*VPOSR1+VDIRR2*VPOSR2)/RP/
     #DSQRT(1.D0-VDIREJE**2)

C Investiga las intersecciones con cada radio
      DO 2 I=1,4
      RM2=RADIO(I)**2
      IF (RM2.GT.RP2) THEN
        X=-RP*COSALF+DSQRT(RM2-RP2*(1.D0-COSALF**2))
        DIS(9-I)=X/DSQRT(1.D0-VDIREJE**2)
      END IF
      IF (RM2.LE.RP2) THEN
        IF(COSALF.GE.0.D0) GOTO 2
        IF(COSALF.LT.0.D0) THEN
          IF(RM2.LE.RP2*(1.D0-COSALF**2)) THEN
            GOTO 2
          ELSE
            X=-RP*COSALF+DSQRT(RM2-RP2*(1.D0-COSALF**2))
            DIS(9-I)=X/DSQRT(1.D0-VDIREJE**2)
            X=-RP*COSALF-DSQRT(RM2-RP2*(1.D0-COSALF**2))
            DIS(I)=X/DSQRT(1.D0-VDIREJE**2)
          END IF
        END IF
      END IF
2     CONTINUE

C Investiga las intersecciones con la tapa del cilindro
C (si no hay interseccion, saltea el siguiente paso)

      X1=0.D0
      X2=0.D0
      X3=0.D0
      IF (VDIREJE) 3,6,5
3     X1=-(H/2.D0+VPOSEJE)/VDIREJE
      X2=-(H/2.D0+DH+VPOSEJE)/VDIREJE
      X3=0.D0
      IF (VPOSEJE.GT.H/2.D0) X3=(H/2.D0-VPOSEJE)/VDIREJE
      GOTO 6
5     X1=(H/2.D0-VPOSEJE)/VDIREJE
      X2=(H/2.D0+DH-VPOSEJE)/VDIREJE
      X3=0.D0
      IF (VPOSEJE.LT.-H/2.D0) X3=(-H/2.D0-VPOSEJE)/VDIREJE

C Si X1 es negativo, el neutron esta en la tapa y saliendo hacia afuera,
C por lo que recorre un solo medio que se pone en DIS(8)

6      IF (X1.LT.0.D0) THEN
       IF (DIS(8).NE.0.D0 .AND. X2.GT.0.D0 .AND. X2.LT.DIS(8)) DIS(8)=X2
        DO I=1,7
          DIS(I)=0.D0
        END DO
        RETURN
       END IF

C Si X3 es no nula porque el neutron esta volando desde una tapa,
C resta X3 a todas las intersecciones con los cilindros

      DO I=1,8
         IF (DIS(I).GT.0.D0) DIS(I)=DIS(I)-X3
      END DO

C Compara las intersecciones con los cilindros y con la tapa
C a fin de determinar por donde sale antes

      DO 7 I=1,8
      IF(DIS(I).GT.X1-X3.AND.X1.NE.0.D0) GOTO 8      ! CONDICION MODIFICADA 4/10/01
7     CONTINUE
      GOTO 4                              ! LINEA AGREGADA EL 20/3/01
8     DIS(I)=X1-X3
       IF (I.LT.8) THEN
         DO 9 J=I+1,8
9        DIS(J)=0.D0
       END IF

C Hace las diferencias de distancias sucesivas a fin de ver
C la distancia recorrida en cada material

4     ZZ(1)=0.D0
      DO 10 I=2,8
      ZZ(I)=0.D0
10    IF (DIS(I).GT.0.D0) ZZ(I)=DIS(I)-DIS(I-1)
      DO 11 I=1,8
11    DIS(I)=ZZ(I)

C Agrega las distancias recorridas en las tapas de los cilindros

      IF ((X3.NE.0.D0.AND.DIS(2).EQ.0.D0) .OR. X3.LT.DIS(2)) DIS(2)=X3
      DX=X2-X1
      IF (DIS(8).NE.0.D0 .AND. DX.NE.0.D0 .AND. DX.LT.DIS(8)) DIS(8)=DX
      IF (DIS(8).EQ.0.D0) DIS(8)=DX

      RETURN
      END SUBROUTINE DISTAN


C--------------------------------------------------------------------
C     Calcula el camino libre medio del neutron y la distancia a la
C     siguiente interaccion. Modifica el peso del neutr-n.
C--------------------------------------------------------------------
      SUBROUTINE LAMBDA(XL,RANDOM,WGT,ICAN,IMED,TRANS)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /DATOSN/EIN,U,THETA,ALFA,DIS(8)
      COMMON /DATOSSAM/SE,SE0,SM,PS
      COMMON /DATOSCAN/SEC,SE0C,SMC,PSC

C Caso de corrida general (unico habilitado)

      TRANS=DEXP(-SMC*(DIS(2)+DIS(4)+DIS(6)+DIS(8))-SM*(DIS(3)+DIS(7)))
      if (wgt.lt.0) then
            write(15,*) wgt
      end if
      if (trans.GE.1.d0) then
	 wgt=0.d0
	 return
      end if
      A=1.D0/(1.D0-TRANS)
      X1=(1.D0-DEXP(-SMC*DIS(2)))*A
      X2=(1.D0-DEXP(-SMC*DIS(2)-SM*DIS(3)))*A
      X3=(1.D0-DEXP(-SMC*(DIS(2)+DIS(4))-SM*DIS(3)))*A
      X4=(1.D0-DEXP(-SMC*(DIS(2)+DIS(4)+DIS(6))-SM*DIS(3)))*A
      X5=(1.D0-DEXP(-SMC*(DIS(2)+DIS(4)+DIS(6))-SM*(DIS(3)+DIS(7))))*A

      IF (RANDOM.LT.X1) THEN
            XL=-DLOG(1.D0-RANDOM/A)/SMC
            ICAN=1
C     SCATT. EN CAN
            IMED=1
      END IF
      IF (RANDOM.GE.X1.AND.RANDOM.LT.X2) THEN
            XL=DIS(2)-(DLOG(1.D0-RANDOM/A)+DIS(2)*SMC)/SM
C	SCATT. EN AGUA
            IMED=0
      END IF
      IF (RANDOM.GE.X2.AND.RANDOM.LT.X3) THEN
            XL=DIS(2)+DIS(3)-
     #(DLOG(1.D0-(RANDOM/A))+DIS(2)*SMC+DIS(3)*SM)/SMC
            ICAN=1
C        SCATT. EN CAN
            IMED=1
      END IF
      IF (RANDOM.GE.X3.AND.RANDOM.LT.X4) THEN
            XL=DIS(2)+DIS(3)+DIS(4)+DIS(5)-
     #(DLOG(1.D0-(RANDOM/A))+(DIS(2)+DIS(4))*SMC+DIS(3)*SM)/SMC
            ICAN=1
C        SCATT. EN CAN
            IMED=1
      END IF
      IF (RANDOM.GE.X4.AND.RANDOM.LT.X5) THEN
            XL=DIS(2)+DIS(3)+DIS(4)+DIS(5)+DIS(6)-
     #(DLOG(1.D0-(RANDOM/A))+(DIS(2)+DIS(4)+DIS(6))*SMC+
     #DIS(3)*SM)/SM
C      SCATT. EN AGUA
            IMED=0
      END IF
      IF (RANDOM.GE.X5.AND.RANDOM.LT.1.D0) THEN
            XL=DIS(2)+DIS(3)+DIS(4)+DIS(5)+DIS(6)+DIS(7)-
     #(DLOG(1.D0-(RANDOM/A))+(DIS(2)+DIS(4)+DIS(6))*SMC+
     #(DIS(3)+DIS(7))*SM)/SMC
C      SCATT. EN CAN
            ICAN=1
            IMED=1
      END IF

      IF (WGT.LT.0) THEN
            write(13,*) PS, PC, WGT
      END IF

      WGT=WGT*(1.D0-TRANS)


      IF (IMED.EQ.0) THEN
        WGT=WGT*PS
      ELSE
        WGT=WGT*PSC
      END IF

      RETURN
      END SUBROUTINE LAMBDA

C--------------------------------------------------------------------
C     Calcula la nueva ubicaci-n de neutr-n a partir de la anterior y
C     de los -ngulos correspondientes.
C--------------------------------------------------------------------
      SUBROUTINE POSNUE(R,XL,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(2,3),D(2,3)
      DO 10,I=1,3
10     R(2,I)=R(1,I)+XL*D(2,I)
      RETURN
      END

C--------------------------------------------------------------------
C     Genera un numero al azar entre 0 y 1.
C--------------------------------------------------------------------
      SUBROUTINE GGUBFS(DSEED,RANDOM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RANDOM(1000)
      DATA D2P31M/2147483647.D0/
      DATA D2P31/2147483711.D0/
      DO 1 I=1,1000
      DSEED=DMOD(16807.D0*DSEED,D2P31M)
      R=DSEED/D2P31
      IF (R.EQ.0.D0) R=.0000001D0
      RANDOM(I)=R
1     CONTINUE
      RETURN
      END SUBROUTINE GGUBFS

C--------------------------------------------------------------------
C     Calcula los -ngulos de dispersi-n cuando ocurre una interacci-n.
C     Cambia los valores de posici-n y -ngulo, conservando los que
C     surgen de la -ltima iteraci-n.
C     Calcula el nuevo versor direcci-n.
C--------------------------------------------------------------------
      SUBROUTINE ANGULOS(A,R,D,U0,U,IMED,RAND1,RAND2)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER :: II, DELTAT
      COMMON /MUESTRA/BKT
      DIMENSION A(2,2),R(2,3),D(2,3),V(90)
      DATA PI/3.1415926535897932D0/
      E0=10.D0**U0
      E=10.D0**U

      STTT=0.D0
      S1=0.D0
      I=1
      DELTAT=2
      DO 100 II=2,180,DELTAT
      T=FLOAT(II)
      THETA=T/180.D0*PI
      IF (IMED.EQ.1) THEN
        CALL SDDSCAN(E0,E,THETA,S2)
      ELSE
        CALL MODEL(E0,E,THETA,BKT,3,S2)
      END IF
      S2=S2*E
      V(I)=S2
      S2=S2*2.D0*PI*DSIN(THETA)*PI/180.D0
      STTT=.5D0*(S1+S2)*DELTAT+STTT
      I=I+1
100   S1=S2

      P=0.D0
      S1=0.D0
      I=1
      DO 1 II=2,180,DELTAT
      T=FLOAT(II)
      THETA=T/180.D0*PI
      S2=V(I)
      S2=S2*2.D0*PI*DSIN(THETA)*PI/180.D0/STTT
      DS=.5D0*(S1+S2)*DELTAT
      IF (P+DS-RAND1)2,3,3
2     P=P+DS
      S1=S2
      I=I+1
1     CONTINUE
      TH=180.D0
      GOTO 4

3     AA=(S2-S1)/DELTAT
      B=(S1*T-S2*(T-DELTAT))/DELTAT
      BA=B/AA
      RAD=BA*BA+2.D0*(RAND1-P)/AA+(T-DELTAT)**2.D0+2*BA*(T-DELTAT)
      TH=-BA+DSQRT(RAD)
      IF (.NOT.(TH.GT.T-DELTAT.AND.TH.LT.T)) TH=-BA-DSQRT(RAD)

4     A(2,1)=TH
      A(2,2)=360.D0*RAND2

      DO 20,J=1,3
20    R(1,J)=R(2,J)
      A(1,2)=0.D0
      IF (DABS(D(2,3)).LT.1.D0) THEN
      A(1,2)=DACOS(D(2,1)/DSQRT(1.D0-D(2,3)*D(2,3)))*180.D0/PI
      END IF
      A(1,1)=DACOS(D(2,3))*180.D0/PI

      T1=A(1,1)/180.D0*PI
      T2=A(2,1)/180.D0*PI
      F1=A(1,2)/180.D0*PI
      F2=A(2,2)/180.D0*PI
      D(2,1)=DCOS(T1)*DCOS(F1)*DSIN(T2)*DCOS(F2)-
     #     DSIN(F1)*DSIN(T2)*DSIN(F2)+DSIN(T1)*DCOS(F1)*DCOS(T2)
      D(2,2)=DCOS(T1)*DSIN(F1)*DSIN(T2)*DCOS(F2)+
     #     DCOS(F1)*DSIN(T2)*DSIN(F2)+DSIN(T1)*DSIN(F1)*DCOS(T2)
      D(2,3)=DCOS(T1)*DCOS(T2)-DSIN(T1)*DSIN(T2)*DCOS(F2)

      RETURN
      END SUBROUTINE ANGULOS
C--------------------------------------------------------------------
C     Genera la distribuci-n angular para scattering en el can y la
C     normaliza.
C--------------------------------------------------------------------
      SUBROUTINE GENDANG(U,QEMP,DQ,A0,A1,A2,DANG)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer :: II
      PARAMETER(NDQE=400)
      DIMENSION QEMP(NDQE),DQ(NDQE),DANG(181)
      DATA PI/3.1415926535897932D0/
C Genera distribuci-n angular
      I=0
      DO 3 II=0,180,1
       T = float(II)
       I=I+1
       E=10.D0**U
       Q=2.D0*DSQRT(E/0.002072002D0)*DSIN(T/180.D0*PI/2.D0)
       IF (Q.GT.QEMP(NDQE)) THEN
       DANG(I)=A0+A1*Q+A2*Q**2
       GOTO 3
       END IF
       IF (QEMP(1).GT.Q) THEN
           DANG(I)=DQ(1)
           GOTO 3
       END IF
       DO 4 J=2,NDQE
         IF(QEMP(J).GT.Q) THEN
           Z=(Q-QEMP(J-1))/(QEMP(J)-QEMP(J-1))
           DANG(I)=DQ(J-1)+Z*(DQ(J)-DQ(J-1))
           GOTO 3
         END IF
4     CONTINUE
3     CONTINUE
C Calcula la integral de la distribucion angular
      STTT=0.D0
c      S1=DANG(1)*2.D0*PI*DSIN(THETA)*PI/180.D0
      s1=0.d0
      I=1
      DO 100 II=1,180,1
      T = float(II)
      THETA=T/180.D0*PI
      S2=DANG(I)
      S2=S2*2.D0*PI*DSIN(THETA)*PI/180.D0
      STTT=.5D0*(S1+S2)*1.D0+STTT
      I=I+1
100   S1=S2

C Normaliza
      DO 200 I=1,181
200   DANG(I)=DANG(I)/STTT
c      open(7,file='x')
c      do i=1,181
c      write(7,*)i-1,2*pi*dsin((i-1.d0)*pi/180.d0)*dang(i)*pi/180
c      end do
c      write(6,*)'escribi- dang para u=',u
c      close(7)
c      read(5,*)
      RETURN
      END SUBROUTINE GENDANG
C--------------------------------------------------------------------
C     Integra la distribuci-n en Q.
C     Da la constante de normalizaci-n.
C--------------------------------------------------------------------
      SUBROUTINE NORMDQ(U,QEMP,DQ,A0,A1,A2,DQNOR)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NDQE=400)
      DIMENSION QEMP(NDQE),DQ(NDQE)
      DATA PI/3.1415926535897932D0/
C CALCULA EL LIMITE SUPERIOR DE INTEGRACION QUE ES 2*K0
      E=10.D0**U
      XK0=DSQRT(E/0.002072002D0)
      QMAX=2.D0*XK0
C BUSCA EN QUE POSICION DEL VECTOR ESTA EL MAXIMO Q
      CALL LOCATE(QEMP,NDQE,QMAX,JMAX)
C HACE QUE JMAX SEA IMPAR PARA PODER INTEGRAR BIEN POR SIMPSON
      IF (2*(JMAX/2).EQ.JMAX) JMAX=JMAX+1
      IF (JMAX.EQ.NDQE) JMAX=NDQE-1
C INTEGRA EN EL POSIBLE RANGO NO TABULADO SUPERIOR A QEMP(NDQE)
      Q1=2.D0*DSQRT(E/0.002072002D0)
      Q0=QEMP(NDQE)
      SEXT=0.D0
      IF (Q1.GT.Q0) THEN
        SEXT=A0*(Q1-Q0)+A1/2.D0*(Q1**2-Q0**2)+A2/3.D0*(Q1**3-Q0**3)
      END IF

C INTEGRA

      DELTAQ=QEMP(2)-QEMP(1)
      S2=0.D0
      S4=0.D0
      DO I=1,JMAX/2
         J=2*I
         S4=S4+2.D0*PI/XK0**2*QEMP(J)*DQ(J)
         J=2*I+1
         S2=S2+2.D0*PI/XK0**2*QEMP(J)*DQ(J)
      END DO
         SIG1=2.D0*PI/XK0**2*QEMP(1)*DQ(1)
         SIG2=2.D0*PI/XK0**2*QEMP(JMAX)*DQ(JMAX)
         DQNOR=DELTAQ/3.D0*(SIG1+4.D0*S4+2.D0*S2-SIG2)+SEXT
      RETURN
      END SUBROUTINE NORMDQ
C--------------------------------------------------------------------
C     Calcula los -ngulos de dispersi-n cuando ocurre una interacci-n.
C     Cambia los valores de posicion y -angulo, conservando los que
C     surgen de la ultima iteracion.
C     Calcula el nuevo versor direcci-n.
C     Caso de distribucion angular empirica.
C--------------------------------------------------------------------
      SUBROUTINE ANGUEMP(A,R,D,DANG,RAND1,RAND2)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER::II,DELTAT
      DIMENSION A(2,2),R(2,3),D(2,3),DANG(181)
      DATA PI/3.1415926535897932D0/

      P=0.D0
      S1=0.D0
      I=1
      DELTAT=1
      DO 1 II=1,180,DELTAT
      T=FLOAT(II)
      THETA=T/180.D0*PI
      S2=DANG(I)
      S2=S2*2.D0*PI*DSIN(THETA)*PI/180.D0
      DS=.5D0*(S1+S2)*DELTAT
      IF (P+DS-RAND1)2,3,3
2     P=P+DS
      S1=S2
      I=I+1
1     CONTINUE
      TH=180.D0
      GOTO 4

3     AA=(S2-S1)/DELTAT
      B=(S1*T-S2*(T-DELTAT))/DELTAT
      BA=B/AA
      RAD=BA*BA+2.D0*(RAND1-P)/AA+(T-DELTAT)**2.D0+2.D0*BA*(T-DELTAT)
      TH=-BA+DSQRT(RAD)
      IF (.NOT.(TH.GT.T-DELTAT.AND.TH.LT.T)) TH=-BA-DSQRT(RAD)


4     A(2,1)=TH
      A(2,2)=360.D0*RAND2
      DO 20,J=1,3
20    R(1,J)=R(2,J)
      A(1,2)=0.D0
      IF (DABS(D(2,3)).LT.1.D0) THEN
      A(1,2)=DACOS(D(2,1)/DSQRT(1.D0-D(2,3)*D(2,3)))*180.D0/PI
      END IF
      A(1,1)=DACOS(D(2,3))*180.D0/PI

      T1=A(1,1)/180.D0*PI
      T2=A(2,1)/180.D0*PI
      F1=A(1,2)/180.D0*PI
      F2=A(2,2)/180.D0*PI
      D(2,1)=DCOS(T1)*DCOS(F1)*DSIN(T2)*DCOS(F2)-
     #     DSIN(F1)*DSIN(T2)*DSIN(F2)+DSIN(T1)*DCOS(F1)*DCOS(T2)
      D(2,2)=DCOS(T1)*DSIN(F1)*DSIN(T2)*DCOS(F2)+
     #     DCOS(F1)*DSIN(T2)*DSIN(F2)+DSIN(T1)*DSIN(F1)*DCOS(T2)
      D(2,3)=DCOS(T1)*DCOS(T2)-DSIN(T1)*DSIN(T2)*DCOS(F2)

      RETURN
      END SUBROUTINE ANGUEMP


C-----------------------------------------------------------------------
C     Ruleta rusa
C-----------------------------------------------------------------------
      SUBROUTINE RULRUS(RANDOM,WGT,WCO,IFLAG)
      IMPLICIT REAL * 8 (A-H,O-Z)

      IFLAG=1
      if (wgt.LE.0.d0) then
	iflag=0
	return
      end if
      IF (WGT.LT.WCO) THEN
        IF (RANDOM.LT..5D0) THEN
           WGT=WGT*2.D0
          ELSE
           IFLAG=0
        END IF
      END IF
c      write(14,*) WGT c THIS IS OK WITH WGT
      RETURN
      END SUBROUTINE RULRUS
C----------------------------------------------------------------------
C
C                       SMPSNU.FOR
C                             U2
C  Evalua                    /
C                            |  I(E) EdU
C                            |
C                            /
C                             U1
C-----------------------------------------------------------------------

      SUBROUTINE SMPSNU(U1,U2,VAL,IFLAG,IMED)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      PARAMETER (N=50)
      DU=(U2-U1)/DBLE(N)
      S2=0.D0
      S4=0.D0
      DO I=1,N/2
         U=U1+(2*I-1)*DU
         CALL SIGMUE(IFLAG,U,IMED,SIG)
         S4=S4+SIG
         U=U1+2*I*DU
         CALL SIGMUE(IFLAG,U,IMED,SIG)
         S2=S2+SIG
      END DO
         CALL SIGMUE(IFLAG,U1,IMED,SIG1)
         CALL SIGMUE(IFLAG,U2,IMED,SIG2)
         VAL=DU/3.D0*(SIG1+4.D0*S4+2.D0*S2-SIG2)
      RETURN
      END SUBROUTINE SMPSNU

C-----------------------------------------------------------------------
C
C  Calcula el valor de la funci-n
C                  2
C            1    d s
C           ---  ---- (E,E') * exp( -s (E')*D ) Ef(E')
C            s   dOdE                 T
C             T
C
C  donde E' es la energia de salida.
C
C  IFLAG=0 : Calcula sin atenuaci-n ni correcci-n por eficiencia del
C            detector
C  IFLAG distinto de 0: Calculo completo.
C-----------------------------------------------------------------------

      SUBROUTINE SIGMUE(IFLAG,X,IMED,SIG)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /DATOSN/EIN,U,THETA,ALFA,DIS(8)
      COMMON /DATOSSAM/SE,SE0,SM,PS
      COMMON /DATOSCAN/SEC,SE0C,SMC
      COMMON /TEMPERATURA/XKBT
      DATA HQS2M/2.07219093D-3/         ! eV
      E0=10**U
      E=10.D0**X
      Q=DSQRT((EIN+E-2.D0*DSQRT(EIN*E)*DCOS(ALFA))/HQS2M)
      ! CALL VALINTERP(X)

      IF (IMED.EQ.1) THEN
         CALL SDDSCAN(E0,E,THETA,S2DIF)
      ELSE
CC         CALL SDDSSAM(E0,E,THETA,S2DIF)
           CALL MODEL(E0,E,THETA,XKBT,3,S2DIF)   !NUEVO
      END IF

      IF (IFLAG.EQ.0) THEN
         IF (IMED.EQ.0) SIG=DLOG(10.D0)*E*S2DIF/SE0
         IF (IMED.EQ.1) SIG=DLOG(10.D0)*E*S2DIF/SE0C
         RETURN
      END IF

      TRANS=DEXP(-SMC*(DIS(2)+DIS(4)+DIS(6)+DIS(8))-
     #          SM*(DIS(3)+DIS(7)))
      IF (IMED.EQ.0) SIG=DLOG(10.D0)*E*S2DIF*TRANS*EFI(X)/SE0
      IF (IMED.EQ.1) SIG=DLOG(10.D0)*E*S2DIF*TRANS*EFI(X)/SE0C
      RETURN
      END SUBROUTINE SIGMUE

C----------------------------------------------------------------------
C     Funcion eficiencia aproximada para un tubo de 3He de una pulgada
C     de diametro, a incidencia normal
C----------------------------------------------------------------------
      REAL*8 FUNCTION EFI(X)
      IMPLICIT REAL * 8 (A-H,O-Z)
C      IF (X.LE.-0.5D0) THEN
C        EFI=1.D0/(1.D0+DEXP(2.D0*(X+0.5D0)))
C      ELSE
C        EFI=1.D0/(1.D0+DEXP(1.3D0*(X+0.5D0)))
C      END IF
c     PREGUNTAR A JAVIER
      X = X
      EFI=1.D0
      RETURN
      END FUNCTION EFI

C----------------------------------------------------------------------
C     Dado un vector ordenado, y un numero, encuentra entre cuales
C     valores del vector se halla ese numero
C----------------------------------------------------------------------
      SUBROUTINE LOCATE(XX,N,X,J)
      IMPLICIT REAL * 8 (A-H,O-Z)
      DIMENSION XX(N)
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GE.XX(1)).EQV.(X.GE.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GOTO 10
      ENDIF
      IF(X.EQ.XX(1))THEN
        J=1
      ELSE IF(X.EQ.XX(N))THEN
        J=N-1
      ELSE
        J=JL
      ENDIF
      RETURN
      END SUBROUTINE LOCATE

C----------------------------------------------------------------------------
      SUBROUTINE SBRSIGCAN(E0,E,S0TOT)
      IMPLICIT REAL * 8 (A-K,M-Z)
      COMMON /DATOSC/AC,KBTAUC,SIGBC
      S0TOT=0.D0
      X=DSQRT(E/KBTAUC)
      X0=DSQRT(E0/KBTAUC)
      ETA=(AC+1.D0)/2.D0/DSQRT(AC)
      RHO=(AC-1.D0)/2.D0/DSQRT(AC)
      Z=X**2.D0-X0**2.D0
      IF (E.GE.E0) THEN
        S0TOT=0.D0
        IF (DABS(Z).LT.5000)
     #    S0TOT=AC*SIGBC/8.D0/E0*(DFER(ETA*X-RHO*X0)-
     #    DFER(ETA*X+RHO*X0)+
     #    DEXP(-Z)*(DFER(RHO*X+ETA*X0)-
     #    DFER(RHO*X-ETA*X0)))
      ELSE
        S0TOT=0.D0
        IF (DABS(Z).LT.5000)
     #    S0TOT=AC*SIGBC/8.D0/E0*(DFER(ETA*X+RHO*X0)+
     #    DFER(ETA*X-RHO*X0)+
     #    DEXP(-Z)*(DFER(-RHO*X+ETA*X0)-
     #    DFER(RHO*X+ETA*X0)))
      END IF

      S0TOT=S0TOT*DLOG(10.D0)*E

      RETURN
      END SUBROUTINE SBRSIGCAN
C--------------------------------------------------------------------
      SUBROUTINE SDDSCAN(E0,E,ZETA,TT)
      IMPLICIT REAL * 8 (A-K,M-Z)
      COMMON /DATOSC/AC,KBTAUC,SIGBC
      DATA PI/3.1415926535897932D0/

      XHW=E0-E
      UM=(E0+E-2.*DSQRT(E0*E)*DCOS(ZETA))/AC
      UT=4.D0*KBTAUC
      IF(UM.NE.0.D0) THEN
       S=DEXP(-(XHW-UM)**2/(UM*UT))/DSQRT(UM*PI*UT)
      ELSE
       S=0
      END IF


      TT=SIGBC*S
      TT=TT*DSQRT(E/E0)/4.D0/PI

      RETURN
      END SUBROUTINE SDDSCAN

C----------------------------------------------------------------------------
      SUBROUTINE MODEL(E0,EF,ZZR,KBT,LINDIC,SALIDA)
C============================================================================
C SUBRUTINA MODEL:                                                 23/09/93
C     INPUT:
C       E0 = ENERGIA INCIDENTE (OBLIGATORIO)
C       EF = ENERGIA FINAL (LINDIC=2,3)
C       ZZR = ANGULO DE DISPERSION EN RADIANES (LINDIC=3)
C       KBT= TEMPERATURA EN eV (OBLIGATORIO)
C       LINDIC = 0 : CALCULA SECCION EFICAZ TOTAL (CON ABSORCION)
C       LINDIC = 1 : CALCULA SECCION EFICAZ TOTAL DE SCATTERING
C       LINDIC = 2 : CALCULA NUCLEO DE TRANSFERENCIA 0 DE E0 A EF
C       LINDIC = 3 : CALCULA SECCION EFICAZ DOBLE DIFERENCIAL
C NO!!  LINDIC = 4 : CALCULA EL COSENO MEDIO DE DISPERSION (TERMINO PPAL)
C     OUTPUT:
C       SALIDA = RESULTADO DEL CALCULO PEDIDO
C..................................................................
      IMPLICIT REAL * 8 (A-K,M-Z)
      PARAMETER (LAZ=5, LMZ=5)

      COMMON /DATOS/SABS,
     1         HW(LAZ,LMZ),HSIG(LAZ,LMZ),ML(LAZ,LMZ),
     1         SIGB(LAZ),MAT(LAZ),MMOL(LAZ),PP100(LAZ),
     1         LELEM,LATS,LMOD(LAZ)


      COMMON /INICIA/NWPA(LAZ,LMZ),NWPC(LAZ,LMZ),PE0P(LAZ,LMZ),
     1         GAMP(LAZ,LMZ),KBTP(LAZ),KBTPP(LAZ,LMZ),FPRIN(LAZ)

      COMMON /DE_E0/E0AN,PE0(LAZ,LMZ),NULA(LAZ),ETALA(LAZ),ROLA(LAZ),
     1         MULA(LAZ),TAULA(LAZ),PEPF(LAZ,LMZ)

      COMMON /PIES/PI,RPI,RADIAN

c      write(6,*)'e0,e0an',e0,e0an
c      read(5,*)

      IF(E0.NE.E0AN)CALL SBRE0(E0,KBT)

      LABSOR=0
      IF(LINDIC.EQ.0)THEN
        LINDIC=1
        LABSOR=1
      ENDIF

      SALIDA=0.D0
      RES   =0.D0
      KERAT =0.D0
      KERCT =0.D0
      SS0T  =0.D0
      SS1T  =0.D0
      ZZ=ZZR*180.D0/PI
c el angulo venia en radianes y lo pasa a grados
      DO 100 LA=1,LATS
        MU=MULA(LA)
        TAU=TAULA(LA)
        NU=NULA(LA)
        ETA=ETALA(LA)
        RO=ROLA(LA)

C TERMINO PRINCIPAL:
        SELECT CASE (LINDIC)
C -------------------------------------------------------
        CASE(1)
c      write(6,*)'1 1'
C   SECCION EFICAZ TOTAL:
          Z=MU*E0/TAU
          IF(NU.NE.1.D0)THEN
            C=(1.D0-NU*NU)/(1.D0+2.D0*MU*NU+MU*MU)
            RES=PI/Z*MU*MU*4.D0/(1.D0-NU*NU)*
     +        (DFER(DSQRT(Z))-
     +         DSQRT(1.D0-C)*DEXP(-C*Z)*DFER(DSQRT(Z*(1.D0-C))))
          ELSE
C     PREGUNTAR A JAVIER C=0
            C=0
            RES=4.D0*PI/(1.D0+2.D0/MU+1.D0/MU/MU)*
     +        ((1.D0+.5D0/Z)*DFER(DSQRT(Z))+
     +         DEXP(-Z)/DSQRT(PI*Z))
          ENDIF
C -------------------------------------------------------
        CASE(2)
C   NUCLEO DE TRANSFERENCIA DE ENERGIA DE ORDEN CERO DE E0 A EF:
          X0=DSQRT(E0/TAU)
          XF=DSQRT(EF/TAU)
          LUPDO=-1
          IF(E0.GT.EF)LUPDO=1
          RES=PI*MU/2.D0/E0*
     +      (DEXP(-(1.D0-NU)/2.D0*(E0-EF)/TAU)*
     +       (DFER(ETA*XF-RO*X0)+LUPDO*DFER(ETA*XF+RO*X0))-
     +       DEXP((1.D0+NU)/2.D0*(E0-EF)/TAU)*
     +       (DFER(RO*XF-ETA*X0)+LUPDO*DFER(RO*XF+ETA*X0)))
C -------------------------------------------------------
        CASE(3)
C   SECCION EFICAZ DOBLE DIFERENCIAL:
          W=E0-EF
          Q=E0+EF-2.D0*DSQRT(E0*EF)*DCOS(ZZ*RADIAN)
          RES=DSQRT(EF/E0)*DEXP(-(1.D0-NU)/2.D0*W/TAU)*
     +      DSQRT(MU/4.D0/PI/Q/TAU)*
     +      DEXP(-(W-Q/MU)*(W-Q/MU)*MU/4.D0/Q/TAU)
C -------------------------------------------------------
        CASE(4)
C   COSENO MEDIO DE DISPERSION (TERMINO PRINCIPAL EXCLUSIVAMENTE)
C          LF=0
C          CALL SBRMU(E0,MU,TAU,NU,ETA,RO,SS0,SS1,LF)
C          SS0T=SS0T+SS0*FPRIN(LA)
C          SS1T=SS1T+SS1*FPRIN(LA)
C -------------------------------------------------------
        END SELECT
        RES=RES*FPRIN(LA)

C TERMINOS FONONICOS:
        IF((NU*NU-1.D0)*(LINDIC-4)*LMOD(LA).EQ.0)GOTO 111
        DO 110 L=1,LMOD(LA)
          KERA=0.D0
          KERC=0.D0
          P0=PE0(LA,L)
          IF(P0.EQ.0.D0)GOTO 110
          PF=PEPF(LA,L)
C
C @@@@@@@@@@@@@@@
C ANIQUILICION:
          E0A=E0+HW(LA,L)-DABS(HSIG(LA,L))
          SELECT CASE (LINDIC)
C -------------------------------------------------------
          CASE (1)
c      write(6,*)'1 2'
C SECCION EFICAZ TOTAL
            Z=E0A/TAU*MU
            KERA=PI/E0A*(4.D0*MU*TAU/(1.D0-NU*NU))**2*
     +        (DFER(DSQRT(Z))-C*(1.D0-C)*DSQRT(Z/PI)*DEXP(-Z)-
     +         DSQRT(1.D0-C)*(1.D0+C/2.D0+C*(1.D0-C)*Z)*DEXP(-C*Z)*
     +         DFER(DSQRT((1.D0-C)*Z)))
            KERA=KERA*NWPA(LA,L)*P0*DSQRT(E0A/E0)
C -------------------------------------------------------
          CASE (2)
C   NUCLEO DE TRANSFERENCIA DE ENERGIA DE ORDEN CERO DE E0 A EF:
            X0=DSQRT(E0A/TAU)
            LUPDO=-1
            IF(E0A.GT.EF)LUPDO=1
            KERA=2.D0*MU*TAU*PI*MU/2.D0/E0A*
     +        ((1.D0+(E0A-EF)/2.D0/TAU)*
     +          DEXP(-(1.D0-NU)/2.D0*(E0A-EF)/TAU)*
     +          (DFER(ETA*XF-RO*X0)+LUPDO*DFER(ETA*XF+RO*X0))-
     +         (1.D0-(E0A-EF)/2.D0/TAU)*
     +          DEXP((1.D0+NU)/2.D0*(E0A-EF)/TAU)*
     +          (DFER(RO*XF-ETA*X0)+LUPDO*DFER(RO*XF+ETA*X0))-
     +         2.D0/DSQRT(PI*MU)*DEXP(-(1.D0-NU)/2.D0*(E0A-EF)/TAU)*
     +          ((XF+X0)*DEXP(-(ETA*XF-RO*X0)**2)+
     +           LUPDO*(XF-X0)*DEXP(-(ETA*XF+RO*X0)**2)))
            KERA=KERA*NWPA(LA,L)*P0*DSQRT(E0A/E0)
C -------------------------------------------------------
          CASE(3)
C   SECCION EFICAZ DOBLE DIFERENCIAL:
          W=E0A-EF
          Q=E0A+EF-2.D0*DSQRT(E0A*EF)*DCOS(ZZ*RADIAN)
          KERA=DSQRT(EF/E0A)*Q*DEXP(-(1.D0-NU)/2.D0*W/TAU)*
     +      DSQRT(MU/4.D0/PI/Q/TAU)*
     +      DEXP(-(W-Q/MU)*(W-Q/MU)*MU/4.D0/Q/TAU)
          KERA=KERA*NWPA(LA,L)*P0*DSQRT(E0A/E0)
C -------------------------------------------------------
          END SELECT
C
C @@@@@@@@@@@@@
C CREACION:
          IF(HSIG(LA,L).GT.0)THEN
C   MODO "CLASICO" DEL TERMINO DE CREACION:
            E0C=E0-HW(LA,L)+HSIG(LA,L)
            IF(E0C.GT.0.)THEN
              SELECT CASE (LINDIC)
C -------------------------------------------------------
              CASE (1)
c      write(6,*)'1 3'
C SECCION EFICAZ TOTAL
                Z=E0C/TAU*MU
                KERC=PI/E0C*(4.D0*MU*TAU/(1.-NU*NU))**2*
     +            (DFER(DSQRT(Z))-C*(1.D0-C)*DSQRT(Z/PI)*DEXP(-Z)-
     +             DSQRT(1.D0-C)*(1.D0+C/2.D0+C*(1.D0-C)*Z)*DEXP(-C*Z)*
     +             DFER(DSQRT((1.D0-C)*Z)))
                KERC=KERC*NWPC(LA,L)*P0*DSQRT(E0C/E0)
C -------------------------------------------------------
              CASE (2)
C   NUCLEO DE TRANSFERENCIA DE ENERGIA DE ORDEN CERO DE E0 A EF:
                X0=DSQRT(E0C/TAU)
                LUPDO=-1
                IF(E0C.GT.EF)LUPDO=1
                KERC=2.D0*MU*TAU*PI*MU/2.D0/E0C*
     +            ((1.D0+(E0C-EF)/2.D0/TAU)*
     +             DEXP(-(1.D0-NU)/2.D0*(E0C-EF)/TAU)*
     +              (DFER(ETA*XF-RO*X0)+LUPDO*DFER(ETA*XF+RO*X0))-
     +             (1.D0-(E0C-EF)/2.D0/TAU)*
     +              DEXP((1.D0+NU)/2.D0*(E0C-EF)/TAU)*
     +              (DFER(RO*XF-ETA*X0)+LUPDO*DFER(RO*XF+ETA*X0))-
     +             2.D0/DSQRT(PI*MU)*DEXP(-(1.D0-NU)/2.D0*(E0C-EF)/TAU)*
     +              ((XF+X0)*DEXP(-(ETA*XF-RO*X0)**2)+
     +               LUPDO*(XF-X0)*DEXP(-(ETA*XF+RO*X0)**2)))
                KERC=KERC*NWPC(LA,L)*P0*DSQRT(E0C/E0)
C -------------------------------------------------------
              CASE(3)
C   SECCION EFICAZ DOBLE DIFERENCIAL:
                W=E0C-EF
                Q=E0C+EF-2.D0*DSQRT(E0C*EF)*DCOS(ZZ*RADIAN)
                KERC=DSQRT(EF/E0C)*Q*DEXP(-(1.D0-NU)/2.D0*W/TAU)*
     +            DSQRT(MU/4.D0/PI/Q/TAU)*
     +            DEXP(-(W-Q/MU)*(W-Q/MU)*MU/4.D0/Q/TAU)
                KERC=KERC*NWPC(LA,L)*P0*DSQRT(E0C/E0)
C -------------------------------------------------------
              END SELECT
            ELSE
              KERC=0.D0
            ENDIF
          ELSE
C   MODO FONON "ANCHO" PARA EL TERMINO DE CREACION:
            SELECT CASE (LINDIC)
C -------------------------------------------------------
            CASE (1)
c      write(6,*)'1 4'
C SECCION EFICAZ TOTAL
              KERC=(4.D0*MU*TAU/(1.D0+2.D0*MU*NU+MU*MU))**2*
     +          DSQRT(MU/4.D0/PI/TAU)*4.D0*PI
              KERC=KERC*NWPC(LA,L)*PF*DSQRT(1.D0/E0)
C -------------------------------------------------------
            CASE (2)
C   NUCLEO DE TRANSFERENCIA DE ENERGIA DE ORDEN CERO DE E0 A EF:
              KERC=DEXP(-EF/4.D0/MU/TAU*(1.D0+2.D0*MU*NU+MU*MU))*
     +          DSQRT(MU/4.D0/PI/TAU)*EF*4.D0*PI
              KERC=KERC*NWPC(LA,L)*PF*DSQRT(1.D0/E0)
C -------------------------------------------------------
            CASE(3)
C   SECCION EFICAZ DOBLE DIFERENCIAL:
              KERC=DEXP(-EF/4.D0/MU/TAU*(1.D0+2.D0*MU*NU+MU*MU))*
     +          DSQRT(MU/4.D0/PI/TAU)*EF
              KERC=KERC*NWPC(LA,L)*PF*DSQRT(1.D0/E0)
C -------------------------------------------------------
            END SELECT
          ENDIF
          KERAT=KERAT+KERA
          KERCT=KERCT+KERC
  110   CONTINUE
  111   CONTINUE
        SALIDA=SALIDA+KERAT+KERCT+RES
  100 CONTINUE
      IF(LINDIC.EQ.4)SALIDA=SS1T/SS0T
      IF(LABSOR.EQ.1)THEN
        SALIDA=SALIDA+SABS/DSQRT(E0)
        LINDIC=0
      ENDIF
      RETURN
      END SUBROUTINE MODEL


C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C --------------------------------------------------------------------------
C --------------------------------------------------------------------------
C SUBRUTINA SBRE0(E0,KBT)
C
C CALCULA TODOS LOS PARAMETROS DEL MODELO PARA UNA DADA ENERGIA INICIAL
C Y TEMPERATURA.
C
C INVOCADA POR SBRSETMU, SBRKER05
C PUEDE NECESITAR INICIA, INDAT Y DFER
C
C LOS DATOS DE INPUT VIENEN POR COMMON ("DATOS","INICIA")
C LOS PARAMETROS DE OUTPUT VAN POR COMMON ("DE_E0")
C --------------------------------------------------------------------------
      SUBROUTINE SBRE0(E0,KBT)
      IMPLICIT REAL * 8 (A-K,M-Z)
      PARAMETER (LAZ=5, LMZ=5)

      COMMON /DATOS/SABS,
     1         HW(LAZ,LMZ),HSIG(LAZ,LMZ),ML(LAZ,LMZ),
     1         SIGB(LAZ),MAT(LAZ),MMOL(LAZ),PP100(LAZ),
     1         LELEM,LATS,LMOD(LAZ)

      COMMON /INICIA/NWPA(LAZ,LMZ),NWPC(LAZ,LMZ),PE0P(LAZ,LMZ),
     1         GAMP(LAZ,LMZ),KBTP(LAZ),KBTPP(LAZ,LMZ),FPRIN(LAZ)

      COMMON /DE_E0/E0AN,PE0(LAZ,LMZ),NULA(LAZ),ETALA(LAZ),ROLA(LAZ),
     1         MULA(LAZ),TAULA(LAZ),PEPF(LAZ,LMZ)
      COMMON /PIES/PI,RPI,RADIAN

      DATA KBTAN/-1.D0/

      E0AN=E0

      IF (KBT.NE.KBTAN) CALL INICIAL(KBT,0)
      KBTAN=KBT

C PARA MODOS ANCHOS (HSIG(LA,L)<=0):
C    SE TIENE FUNCION "DENSIDAD DE ESTADOS ANCHOS":
C                 F(E0)=3*RPI*(E0/HW)^2*(1-DFER(E0/HW))
C    QUE DEFINE PE0(X1=E0/HW)=(1.D0-2.D0*X1*X1)*(1.D0-DFER(X1))+
C                    2.D0/RPI*X1*DEXP(-X1*X1)
C    PARA TERMINOS DE CREACION SE UTILIZA
C         SIGMA=SIGB*F(E0)*NWPC*SIGMA(E0=0,E,THETA)

      DO 2000 LA=1,LATS

        MM=MAT(LA)
        TAU=0.D0
        MU=0.D0
        GAM=0.D0

        DO 2100 L=1,LMOD(LA)

          IF(HSIG(LA,L).LE.0)THEN
            X1=E0/HW(LA,L)
            PE0(LA,L)=(1.D0-2.D0*X1*X1)*(1.D0-DFER(X1))+
     1        2.D0/RPI*X1*DEXP(-X1*X1)
            PEPF(LA,L)=3.D0*RPI*X1*X1*(1.D0-DFER(X1))
            IF(X1.LT..5D0)THEN
              PE0(LA,L)=1.D0
              PEPF(LA,L)=0.D0
            ENDIF
            IF (X1.GT.1.5D0)THEN
              PE0(LA,L)=0.D0
              PEPF(LA,L)=0.D0
            ENDIF
          ELSE
            X1=-(E0-HW(LA,L))/HSIG(LA,L)
            X2=HW(LA,L)/HSIG(LA,L)
            PE0(LA,L)=(RPI*X2*(1.D0+DFER(X1))+DEXP(-X1*X1))/PE0P(LA,L)
            PEPF(LA,L)=0.D0
          ENDIF
          TAU=TAU+PE0(LA,L)*KBTPP(LA,L)
          MU=MU+PE0(LA,L)/ML(LA,L)
          GAM=GAM+PE0(LA,L)*GAMP(LA,L)
 2100   CONTINUE
        MU=1.D0/(1.D0/MM-MU)
        TAU=MU*(KBTP(LA)-TAU)

        NU=1.D0/DSQRT(1.D0+4.D0*GAM*TAU*MU)

        MU=NU*MU
        TAU=TAU*NU

        AUX=DSQRT(MU)
        ETA=.5D0*(AUX+1.D0/AUX)
        RO=.5D0*(AUX-1.D0/AUX)

        NULA(LA)=NU
        ETALA(LA)=ETA
        ROLA(LA)=RO
        MULA(LA)=MU
        TAULA(LA)=TAU

 2000 CONTINUE

      RETURN
      END SUBROUTINE SBRE0

C --------------------------------------------------------------------------
C --------------------------------------------------------------------------
C --------------------------------------------------------------------------
      SUBROUTINE INICIAL(KBT,LPRINT)

      IMPLICIT REAL *8 (A-K,M-Z)

      PARAMETER (LAZ=5 , LMZ=5)

      COMMON /DATOS/SABS,
     1         HW(LAZ,LMZ),HSIG(LAZ,LMZ),ML(LAZ,LMZ),
     1         SIGB(LAZ),MAT(LAZ),MMOL(LAZ),PP100(LAZ),
     1         LELEM,LATS,LMOD(LAZ)

      COMMON /INICIA/NWPA(LAZ,LMZ),NWPC(LAZ,LMZ),PE0P(LAZ,LMZ),
     1         GAMP(LAZ,LMZ),KBTP(LAZ),KBTPP(LAZ,LMZ),FPRIN(LAZ)

c      COMMON /TITULO/TITULO(18)
      COMMON /PIES/PI,RPI,RADIAN

      COMMON /DE_E0/E0AN,PE0(LAZ,LMZ),NULA(LAZ),ETALA(LAZ),ROLA(LAZ),
     1         MULA(LAZ),TAULA(LAZ),PEPF(LAZ,LMZ)


      IF (LATS.EQ.0)      CALL INDAT
      E0AN=-1.D0
C      VALORES CONSTANTES
      DO 888 LL=1,LATS
        KBTP(LL)=0.D0
        FPRIN(LL)=SIGB(LL)*PP100(LL)/4.D0/PI
        DO 2000 L=1,LMOD(LL)
          NW=1.D0/(DEXP(HW(LL,L)/KBT)-1.D0)
          EWL=HW(LL,L)*(NW+.5D0)
          X2=HW(LL,L)/HSIG(LL,L)
          PE0P(LL,L)=RPI*X2*(1.D0+DFER(X2))+DEXP(-X2*X2)
          KBTP(LL)=KBTP(LL)+EWL/ML(LL,L)
          KBTPP(LL,L)=EWL/ML(LL,L)
          GAMP(LL,L)=(2.D0*NW+1.D0)/HW(LL,L)/ML(LL,L)
          NWPA(LL,L)=NW/HW(LL,L)/ML(LL,L)*FPRIN(LL)
          NWPC(LL,L)=(1.D0+NW)/HW(LL,L)/ML(LL,L)*FPRIN(LL)
 2000   CONTINUE

        KBTP(LL)=KBT/MMOL(LL)+KBTP(LL)
  888 CONTINUE

      IF (LPRINT.EQ.0) RETURN

c      WRITE (LPRINT,601) TITULO
C601   FORMAT (1X,//,10X,18A4,//)

      TEMPK=KBT*11604.5D0

      WRITE (LPRINT,610) LELEM,LATS,SABS,TEMPK-273.15D0,TEMPK,KBT
  610 FORMAT (          10X,'SUSTANCIA               : ',I2,/,
     1      10X,'NUMERO DE ATOMOS        : ',I2,/,
     3      10X,'SIGMA DE ABSORCION TOTAL: ',1PE12.5,/,
     6      10X,'TEMPERATURA ( C, K, eV) : ',3(1PE12.5,3X)/)

      DO 889 LL=1,LATS
      WRITE (LPRINT,620)LL,LMOD(LL),MMOL(LL)
     1      ,MAT(LL),SIGB(LL)
     1      ,PP100(LL),(ML(LL,L),L=1,LMOD(LL))
  620 FORMAT (10X,'PARA EL ATOMO NUMERO ',I2,' : ',/,
     1      16X,'NUMERO DE MODOS     :',I2,/,
     1      16X,'MASA MOLECULAR      :',1PE12.5,/,
     2      16X,'MASA ATOMICA        :',1PE12.5,/,
     3      16X,'SIGMA bound         :',1PE12.5,/,
     1      16X,'FACTOR DE CONTRIBUC.:',1PE12.5,/,
     4  16X,'MASA DE CADA MODO   :',5(1PE12.5,3X),//)
      WRITE (LPRINT,*)
      DO L=1,LMOD(LL)
      WRITE (LPRINT,630) L,HW(LL,L),HSIG(LL,L)
      ENDDO
  630 FORMAT (10X,'ENERGIA DEL MODO ',I2,' = ',1PE12.5,
     1      ' CON SIGMA = ',1PE12.5)
      WRITE (LPRINT,*)
  889 CONTINUE

      WRITE (LPRINT,590)
 590  FORMAT (10X,'ATOMO  MODO '
     1,'    MASA          Kb TAU          GAMA            NU')

      DO 200 LL=1,LATS
      MM=MAT(LL)
      DO 150 LKK=LMOD(LL)+1,1,-1
      GAM=0.D0
      KBTAU=0.D0
      MU=0.D0
      DO 100 L=LMOD(LL),LKK,-1
      KBTAU=KBTAU+KBTPP(LL,L)
      MU=MU+1.D0/ML(LL,L)
      GAM=GAM+GAMP(LL,L)
100   CONTINUE
      MUKK=1.D0/(1.D0/MM-MU)
      KBTAUKK=MUKK*(KBTP(LL)-KBTAU)
      NUKK=1.D0/DSQRT(1.D0+4.D0*GAM*MUKK*KBTAUKK)
      WRITE (LPRINT,600) LL,L,MUKK,KBTAUKK,GAM,NUKK
150   CONTINUE
      WRITE (LPRINT,*)
200   CONTINUE
600   FORMAT (10X,I3,4X,I3,2X,5(1PE12.5,3X))

      RETURN
      END SUBROUTINE INICIAL
C========================================================================
C========================================================================

      SUBROUTINE INDAT
      IMPLICIT REAL * 8 (A-K,M-Z)
      PARAMETER (LAZ=5 , LMZ=5)

      COMMON /DATOS/SABS,
     1         HW(LAZ,LMZ),HSIG(LAZ,LMZ),ML(LAZ,LMZ),
     1         SIGB(LAZ),MAT(LAZ),MMOL(LAZ),PP100(LAZ),
     1         LELEM,LATS,LMOD(LAZ)

c      COMMON /TITULO/TITULO(18)
      COMMON /PIES/PI,RPI,RADIAN

      DIMENSION SAB(LAZ)

      RPI=DSQRT(PI)
      RADIAN=PI/180.D0


c      READ(1,1) TITULO
C1     FORMAT (18A4)
      READ(1,*) LELEM,LATS
      DO 100 LL=1,LATS
        READ(1,*)LMOD(LL)
        IF (LMOD(LL).GT.0)THEN
          READ(1,*) (HW(LL,L),L=1,LMOD(LL))
          READ(1,*) (HSIG(LL,L),L=1,LMOD(LL))
          READ(1,*) (ML(LL,L),L=1,LMOD(LL))
        ENDIF
        READ(1,*) SIGB(LL),SAB(LL),MAT(LL),MMOL(LL),PP100(LL)
        SABS=SABS+SAB(LL)*PP100(LL)
  100 CONTINUE
      RETURN
      END SUBROUTINE INDAT



C--------------------------------------------------------------------
C  Contribucion angular. Calcula la contibucion de la historia actual
C  a un bin angular dado. Para ello interpola en la distribucion en escala
C  Q que viene como input (DQMUE y DQCAN en el programa ppal.), le aplica
c  la constante de normalizaci-n correspondiente a esa energ-a, y le
C  aplica atenuacion y eficiencia de los detectores.
c  Modificada el 2/7/07
C--------------------------------------------------------------------

      SUBROUTINE CONTRANG(IFLAG,QEMP,DQ,DQNOR,A0,A1,A2,VAL)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NDQE=400)
      DIMENSION QEMP(NDQE),DQ(NDQE)
      COMMON /DATOSN/EIN,U,THETA,ALFA,DIS(8)
      COMMON /DATOSSAM/SE,SE0,SM,PS
      COMMON /DATOSCAN/SEC,SE0C,SMC,PSC

      E=10.D0**U
      Q=2.D0*DSQRT(E/0.002072002D0)*DSIN(THETA/2.D0)
      IF (Q.GT.QEMP(NDQE)) THEN
      VAL=(A0+A1*Q+A2*Q**2)/DQNOR
      GOTO 3
      END IF
      IF (QEMP(1).GT.Q) THEN
          VAL=DQ(1)/DQNOR
          GOTO 3
      END IF
      DO 4 J=2,NDQE
        IF(QEMP(J).GT.Q) THEN
          Z=(Q-QEMP(J-1))/(QEMP(J)-QEMP(J-1))
          VAL=(DQ(J-1)+Z*(DQ(J)-DQ(J-1)))/DQNOR
          GOTO 3
        END IF
4     CONTINUE
3     CONTINUE

      TRANS=DEXP(-SMC*(DIS(2)+DIS(4)+DIS(6)+DIS(8))-
     #          SM*(DIS(3)+DIS(7)))
      IF (IFLAG.EQ.1) VAL=VAL*TRANS !*EFI(U)

      RETURN
      END SUBROUTINE CONTRANG

      SUBROUTINE BOX_MULLER(U1, U2, X, Y)
      real(8) :: R, theta
      real(8), intent(in) :: U1, U2
      real(8), intent(out) :: X, Y
      DATA PI/3.1415926535897932D0/
      DIMENSION RANDOM(1000)
c     usaremos la transformación de Box-Muller para samplear una gaussiana
c     https://bjlkeng.io/posts/sampling-from-a-normal-distribution/

      R = sqrt(-2 * log(U1))
      theta = 2 * PI * U2
      X = R * cos(theta)
      Y = R * sin(theta)
      END SUBROUTINE BOX_MULLER

      SUBROUTINE INELASTICO_GENERALIZADO(Er, sigma, RAND1, RAND2, Uout)
      real(8), intent(in) :: Er, sigma, RAND1, RAND2
      real(8), intent(out) :: Uout
      real(8) :: X, Y
C     P(w) = 1/(sqrt(2*pi*sigma**2)) * exp(-(w - wr)**2/(2*simga**2))

      WRITE(*,*) "working", RAND1, RAND2
      call BOX_MULLER(RAND1, RAND2, X, Y)
      write(*,*) "Random numbers: ", X, Y
      Eout = X * sigma + Er
      write(*,*) "Eout: ", Eout
      Uout = log10(Eout)
      write(*,*) "Uout: ", Uout

      END SUBROUTINE INELASTICO_GENERALIZADO

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   IMSL ROUTINE NAME   - MERFD  =DFER
C
C----------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - EVALUATE THE ERROR FUNCTION OF A DOUBLE
C                           PRECISION ARGUMENT
C
C   USAGE               - RESULT = DFER(Y)
C
C   ARGUMENTS    Y      - INPUT DOUBLE PRECISION ARGUMENT OF THE ERROR
C                           FUNCTION.
C                DFER   - OUTPUT DOUBLE PRECISION VALUE OF THE ERROR
C                           FUNCTION. DFER MUST BE TYPED DOUBLE
C                           PRECISON IN THE CALLING PROGRAM.
C
C   PRECISION/HARDWARE  - DOUBLE/H32
C                       - NOT AVAILABLE/H36,H48,H60
C                         NOTE - DFER MAY NOT BE SUPPLIED BY IMSL IF
C                           IT RESIDES IN THE MATHEMATICAL SUBPROGRAM
C                           LIBRARY SUPPLIED BY THE MANUFACTURER.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C----------------------------------------------------------------------
C
      FUNCTION DFER(Y)
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DFER,Y
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DIMENSION          P(5),Q(4),P1(9),Q1(8),P2(6),Q2(5)
      DOUBLE PRECISION   P,Q,P1,Q1,P2,Q2,XMIN,XLARGE,SQRPI,X,
     *                   RES,XSQ,XNUM,XDEN,XI
      INTEGER            ISW,I
C                                  COEFFICIENTS FOR 0.0 .LE. Y .LT.
C                                  .477
      DATA               P(1)/113.8641541510502D0/,
     *                   P(2)/377.4852376853020D0/,
     *                   P(3)/3209.377589138469D0/,
     *                   P(4)/.1857777061846032D0/,
     *                   P(5)/3.161123743870566D0/
      DATA               Q(1)/244.0246379344442D0/,
     *                   Q(2)/1282.616526077372D0/,
     *                   Q(3)/2844.236833439171D0/,
     *                   Q(4)/23.60129095234412D0/
C                                  COEFFICIENTS FOR .477 .LE. Y
C                                  .LE. 4.0
      DATA               P1(1)/8.883149794388376D0/,
     *                   P1(2)/66.11919063714163D0/,
     *                   P1(3)/298.6351381974001D0/,
     *                   P1(4)/881.9522212417691D0/,
     *                   P1(5)/1712.047612634071D0/,
     *                   P1(6)/2051.078377826071D0/,
     *                   P1(7)/1230.339354797997D0/,
     *                   P1(8)/2.153115354744038D-8/,
     *                   P1(9)/.5641884969886701D0/
      DATA               Q1(1)/117.6939508913125D0/,
     *                   Q1(2)/537.1811018620099D0/,
     *                   Q1(3)/1621.389574566690D0/,
     *                   Q1(4)/3290.799235733460D0/,
     *                   Q1(5)/4362.619090143247D0/,
     *                   Q1(6)/3439.367674143722D0/,
     *                   Q1(7)/1230.339354803749D0/,
     *                   Q1(8)/15.74492611070983D0/
C                                  COEFFICIENTS FOR 4.0 .LT. Y
      DATA               P2(1)/-3.603448999498044D-01/,
     *                   P2(2)/-1.257817261112292D-01/,
     *                   P2(3)/-1.608378514874228D-02/,
     *                   P2(4)/-6.587491615298378D-04/,
     *                   P2(5)/-1.631538713730210D-02/,
     *                   P2(6)/-3.053266349612323D-01/
      DATA               Q2(1)/1.872952849923460D0/,
     *                   Q2(2)/5.279051029514284D-01/,
     *                   Q2(3)/6.051834131244132D-02/,
     *                   Q2(4)/2.335204976268692D-03/,
     *                   Q2(5)/2.568520192289822D0/
C                                  CONSTANTS
      DATA               XMIN/1.0D-10/,XLARGE/6.375D0/
      DATA               SQRPI/.5641895835477563D0/
C                                  FIRST EXECUTABLE STATEMENT
      X = Y
      ISW = 1
      IF (X.GE.0.0D0) GO TO 5
      ISW = -1
      X = -X
    5 IF (X.LT..477D0) GO TO 10
      IF (X.LE.4.0D0) GO TO 25
      IF (X.LT.XLARGE) GO TO 35
      RES = 1.0D0
      GO TO 50
C                                  ABS(Y) .LT. .477, EVALUATE
C                                  APPROXIMATION FOR DFER
   10 IF (X.LT.XMIN) GO TO 20
      XSQ = X*X
      XNUM = P(4)*XSQ+P(5)
      XDEN = XSQ+Q(4)
      DO 15 I=1,3
         XNUM = XNUM*XSQ+P(I)
         XDEN = XDEN*XSQ+Q(I)
   15 CONTINUE
      RES = X*XNUM/XDEN
      GO TO 50
   20 RES = X*P(3)/Q(3)
      GO TO 50
C                                  .477 .LE. ABS(Y) .LE. 4.0
C                                  EVALUATE APPROXIMATION FOR DFER
   25 XSQ = X*X
      XNUM = P1(8)*X+P1(9)
      XDEN = X+Q1(8)
      DO 30 I=1,7
         XNUM = XNUM*X+P1(I)
         XDEN = XDEN*X+Q1(I)
   30 CONTINUE
      RES = XNUM/XDEN
      GO TO 45
C                                  4.0 .LT. ABS(Y), EVALUATE
C                                  MINIMAX APPROXIMATION FOR DFER
   35 XSQ = X*X
      XI = 1.0D0/XSQ
      XNUM = P2(5)*XI+P2(6)
      XDEN = XI+Q2(5)
      DO 40 I=1,4
         XNUM = XNUM*XI+P2(I)
         XDEN = XDEN*XI+Q2(I)
   40 CONTINUE
      RES = (SQRPI+XI*XNUM/XDEN)/X
   45 RES = RES*DEXP(-XSQ)
      RES = 1.0D0-RES
   50 IF (ISW.EQ.-1) RES = -RES
      DFER = RES
      RETURN
      END FUNCTION DFER


