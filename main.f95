Module MathLib

contains

    SUBROUTINE LEGS (A_,Nn,B_,Xx)
      INTEGER, INTENT (IN) :: Nn
      INTEGER :: iI,iJ
      INTEGER, DIMENSION (Nn) :: INDX
      REAL, INTENT (IN), DIMENSION (Nn,Nn) :: A_
      REAL, INTENT (IN), DIMENSION (Nn) :: B_
      REAL, DIMENSION (Nn,Nn) :: Aa
      REAL, DIMENSION (Nn) :: Bb
      REAL, INTENT (OUT), DIMENSION (Nn) :: Xx

      Aa  = A_
      Bb  = B_
    !
      CALL ELGS (Aa,Nn,INDX)
    !
      DO iI = 1, Nn-1
        DO iJ = iI+1, Nn
          Bb(INDX(iJ)) = Bb(INDX(iJ))-Aa(INDX(iJ),iI)*Bb(INDX(iI))
        END DO
      END DO
    !
      Xx(Nn) = Bb(INDX(Nn))/Aa(INDX(Nn),Nn)
      DO iI = Nn-1, 1, -1
        Xx(iI) = Bb(INDX(iI))
        DO iJ = iI+1, Nn
          Xx(iI) = Xx(iI)-Aa(INDX(iI),iJ)*Xx(iJ)
        END DO
        Xx(iI) =  Xx(iI)/Aa(INDX(iI),iI)
      END DO
    !
    END SUBROUTINE LEGS
    !
    SUBROUTINE ELGS (Aa,Nn,INDX)
      INTEGER, INTENT (IN) :: Nn
      INTEGER :: iI, iJ, iK = 0, ITMP
      INTEGER, INTENT (OUT), DIMENSION (Nn) :: INDX
      REAL :: C1, PI, PI1,PJ
      REAL, INTENT (INOUT), DIMENSION (Nn,Nn) :: Aa
      REAL, DIMENSION (Nn) :: C
    !
    ! Initialize the index
      DO iI = 1, Nn
        INDX(iI) = iI
      END DO
    !
    ! Find the rescaling factors, one from each row
      DO iI = 1, Nn
        C1= 0.0
        DO iJ = 1, Nn
          C1 = AMAX1(C1,ABS(Aa(iI,iJ)))
        END DO
        C(iI) = C1
      END DO
    !
    ! Search the pivoting (largest) element from each column
      DO iJ = 1, Nn-1
        PI1 = 0.0
        DO iI = iJ, Nn
          PI = ABS(Aa(INDX(iI),iJ))/C(INDX(iI))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            iK   = iI
          ENDIF
        END DO
    !
    ! Interchange the rows via INDX(N) to record pivoting order
        ITMP     = INDX(iJ)
        INDX(iJ) = INDX(iK)
        INDX(iK) = ITMP

        DO iI = iJ+1, Nn
          PJ  = Aa(INDX(iI),iJ)/Aa(INDX(iJ),iJ)
    !
    ! Record pivoting ratios below the diagonal
    !
          Aa(INDX(iI),iJ) = PJ
    !
    ! Modify other elements accordingly
          DO iK = iJ+1, Nn
            Aa(INDX(iI),iK) = Aa(INDX(iI),iK)-PJ*Aa(INDX(iJ),iK)
          END DO
        END DO
      END DO
    !
    END SUBROUTINE ELGS
    !
    character(len=20) function rStr(rk)
    implicit none

    ! Convertir un intier en chaine de caractères.
    real, intent(in) :: rk

    write(rStr, *) rk
    rStr = adjustl(rStr)
    end function rStr

    character(len=20) function Str(ik)
    implicit none

    ! Convertir un intier en chaine de caractères.
    integer, intent(in) :: ik

    write(Str, *) ik
    Str = adjustl(str)
    end function Str

    Integer function Str2int(str_)
       implicit none

       character(len=*), intent(in) :: str_
       integer :: int_
       read(str_, *)  int_
       Str2int = int_
    end function Str2int

    Real function Str2float(str_)
       implicit none

       character(len=*), intent(in) :: str_
       Real :: f_

       read(str_, "(F9.3)")  f_
       Str2float = f_
    end function Str2float

    character(len=20) function Formatte(s_, nc_)
    ! Convertir un intier en chaine de caractères.
    integer, intent(in) :: nc_
    integer, intent(in) :: s_
    character(len=20) :: sfmt

        if (s_ == 0) then
            sfmt = REPEAT(" ", int(nc_/2)) // Str(s_) // REPEAT(" ", int(nc_/2))
        else
            sfmt = REPEAT(" ", int(nc_/2)) // "-" // REPEAT(" ", int(nc_/2))
        endif
        write(Formatte, *) s_
        Formatte = adjustl(Formatte)
    end function Formatte

    ! Processus d'importation à partir de fichier IEEE CDF

End module MathLib

program newtonrap

Use MathLib

Implicit none

!> Declarations des variables
real,parameter :: pi=3.14159
integer :: N, I, J, W, NB, k, JN, Np, Ng, Kmax = 100, Cond, narg, slack = 1, io
character(len=30) :: NName

!> Variables principales
complex, allocatable, Dimension(:,:) :: Ya, Z, Spq, Ys
complex, allocatable, Dimension(:) :: V, Sg, Sd
real, allocatable, dimension(:,:) :: YB, YG, Ja, c !, Vs
Integer, allocatable, dimension(:) :: Nfrom, Nto, Indx
real, allocatable, dimension(:) :: P, Q, Fd, X, DX, ModV, Theta
character(1), allocatable, dimension(:) :: Typen
real :: eps = 0.001 , MaxDx, t1, t2, PL, QL, Sbase
character(100) :: A2Fmt, AFmt, A3Fmt, fnamein, fnameout, name, hlpfilename = "NR.hlp"
character(len=100) :: hlpstr
logical :: lookForIn = .false., fileExist, LookForInOut = .FALSE., Yes
logical :: lookForOut= .false., lookForCDF = .FALSE.
Character(len=5) :: CalcWith = "INV"
Character(1) :: Y_
! variables génériques
real :: a, b, d, e

! ici commence mon programme

! Passage des argument, correspondance au nom de fichier. narg représente le nombre d'arguments
narg=command_argument_count()
! Vérification des arguments
! fnamein et fnameout sont respectivement les noms de fichiers d'entrée et de sortie (noms complets)
if (narg > 0) then
   do i= 1, narg
      call get_command_argument(i, name)
      select case(adjustl(name))
            case("--in")
                lookForIn  = .TRUE. !change logical value
            case("--out")
                lookForOut = .TRUE.
            case("--cdfin")
                lookForIn  = .TRUE.
                lookForCDF = .TRUE.
            case("--INV")
                calcWith = "INV"
            case("--GAUSS")
                calcWith = "GAUSS"
            case("--LU")
                calcWith = "LU"
            case("--help")
                inquire(file = hlpfilename, exist = fileExist) ! Vérifier si le fichier d'aide existe
                    if (.not. fileExist) then
                        write(*, *) "fichier d'aide ", hlpfilename, ' inéxistant'
                        write(*, *) "Utilisation de NRPower :"
                        write(*, *) "NRPower [Option] [Attribut]"
                        write(*, *) "Option : --help, --in, --out"
                        write(*, *) "Attribut : Nom du fichier (d'Aide, d'entrée de données, de sortie des résultats"
                        stop
                    else
                        open(10, file=hlpfilename, action='write')
                        io = 0
                        do while (io == 0)
                           read(10, *, IOSTAT=io) hlpstr
                           write(*,*) adjustl(trim(hlpstr))
                        end do
                    endif
            case default
                if (LookForIn) then
                    ! Nom du fichier d'entrée des données réseau
                    fnamein = adjustl(name)
                    ! Vérifier son existance dans le chemin spécifié
                    inquire(file = fnamein, exist = fileExist)
                    if (.not. fileExist) then
                        write(*,*) "fichier de données ", fnamein, " est inéxistant dans l'emplacement désigné ..."
                        stop
                    else
                        write(*,*) "Accès au fichier de données ", fnamein, "..."
                        open(50, file = fnamein, action = "read")
                        if (LookForCDF) Then
                            Call CDFImport1(50)
                            ! Dimensionnement des matrices et vecteurs
                            print*, "Allocation mémoire des matrices du réseau"
                            allocate (Ya(N,N), Z(N,N), Ys(N, N), YB(N,N), YG(N,N), Spq(N,N))
                            allocate (V(N), P(N), Q(N), Typen(N), ModV(N), Theta(N), Sg(N), Sd(N), Indx(N))
                            ! Initialisation des données
                            ! Matrice admittance
                            Ya = cmplx(0., 0.)
                            Call CDFImport2(50)
                            allocate (Nfrom(NB), Nto(NB))
                            ! Vecteurs Nfrom et Nto à 0
                            Nfrom = 0
                            Nto   = 0
                            Call CDFImport3(50)
                            Call Mef
                            Call SaveImport(fnamein)
                            print*, "Dé-allocation dynamique de la mémoire des données CDF ..."
                            deallocate (Ya, Z, Ys, YB, YG, Spq)
                            deallocate (V, P, Q, Typen, ModV, Theta, Sg, Sd, Indx)
                            deallocate (Nfrom, Nto)
                            print*, "Ouverture du fichier de données conforme NMSS ... traitement."
                            LookForCDF = .FALSE.
                            close(50)
                            open(50, file = fnamein, action = "read")
                        endif
                    endif
                    LookForIn = .FALSE.
                    LookForInOut = .TRUE.
                elseif (LookForOut) then
                    ! Nom du fichier d'entrée des données réseau
                    fnameout = adjustl(name)
                    ! Vérifier son existance dans le chemin spécifié
                    inquire(file = fnameout, exist = fileExist)
                    if (fileExist) then
                        write(*,*) "fichier de sortie ", fnameout, " éxiste dans l'emplacement désigné ..."
                        Yes = .FALSE.
                        Do While (.not. Yes)
                           Print*, "Comfirmer l'ecrasement (O/N) ? "
                           read(*, "(A1)") Y_
                           if (Y_ == "O") Then
                               Yes = .TRUE.
                               write(*,*) "Accès au fichier de sortie ", fnameout, "..."
                               open(60, file = fnameout, action='write')
                           else
                               print*,"Impossibilité de continuer ... Merci de réessayer."
                               stop
                           endif
                        EndDo
                    else
                        write(*,*) "Accès au fichier de sortie ", fnameout, "..."
                        open(60, file = fnameout, action='write')
                    endif
                    LookForOut = .FALSE.
                    LookForInOut = .TRUE.
                else
                    write(*,*)"Option [",adjustl(name),"] inconnue"
                endif
      end select
   end do
endif

! saisie Dimension réseau

if (.not. LookForIn .and. .not. LookForInOut) Then
    write(*,*) "Execution des données TEST ... "
    open(50, file='/home/mtamali/NR/NR/obj/Debug/data.dat')
end if

if (.not. LookForOut .and. .not. LookForInOut) Then
    write(*,*) "Ecriture des résultats selon données TEST ... "
    open(60, file='/home/mtamali/NR/NR/obj/Debug/result.res', action='write')
end if

! Choix de laméthode d'entrée des données (CDF ou standard NMSS)
If (.not. lookForCDF) Then
    ! Lecture des nombres de noeuds et de branches
    read(50, *) NName, N, NB, Sbase
    ! Ecriture des données du réseau.
    print*, "Lecture de l'identification du réseau"
    write(60, *) 'Identification : ', NName
    write(60, *) "le nombre de noeuds N= ", trim(Str(N)), ", le nombre de branches Nb= ", trim(Str(NB)),&
   &" et la puissance de base Sb= ", Sbase, " MVA."

    print*, "Allocation mémoire des matrices du réseau"
    allocate (Ya(N,N), Z(N,N), Ys(N, N), YB(N,N), YG(N,N), Spq(N,N), Indx(N))
    allocate (V(N), P(N), Q(N), Typen(N), ModV(N), Theta(N), Sg(N), Sd(N), Nfrom(NB), Nto(NB))
    ! Initialisation des données
    ! Matrice admittance
    Ya = cmplx(0., 0.)
    ! Vecteurs Nfrom et Nto à 0
    Nfrom = 0
    Nto   = 0
    ! Vecteur tension nodal
    print*, "Initialisation du vecteur tension nodal du réseau ..."
    V=cmplx(1.,0.)
    ! Dimensionnement des matrices et vecteurs
    ! Initialisation des puissances
    Sg=(0.0, 0.0)
    Sd=(0.0, 0.0)
    ! Lecture des données des lignes de transport
    print*, "Lecture des donées lignes de transport du réseau"
    do w=1,NB
        ! Bfrom BTo R X Gsh Bsh
        Read(50,*) I, J, a, b, d, e
        Nfrom(w) = I
        Nto(w) = J
        Z(I,J)=cmplx(a,b)
        Ys(I,J)=cmplx(d,e)
    end do
    ! Initialisation des variables Np et Ng : Nombre de noeuds consommateurs et générateurs respectivement
    Np    = 0
    Ng    = 0
    ! Initialisation des vecteur ModV et Theta, Module de tension nodal et Argument
    ModV  = 1.0
    Theta = 0.0

    ! Saisie type des noeud, valeurs spécifiées.
    ! Noeud Générateur G (Lire V, Pg), Consommateur C (Lire Pd, Qd) et Bilan B (Affecter V)
    print*, "Lecture des données caractéristiques des noeuds du réseau (pu)"
    do I = 1, N
       Read(50,*) typen(I), a, b
       if (typen(I) == "G") then
           Ng      = Ng + 1
           ModV(i) = a
           V(I)    = ModV(I)*(V(i)/ModCmplx(V(I)))
           Sg(I)=Cmplx(b, 0.0)
           P(I)    = Real(Sg(I) - Sd(I))
           write(60,"(A,A,A,F9.3,A,A,A,F9.3)")  "Noeud générateur  |V|(", trim(Str(I)), ")=", &
           &ModCmplx(V(I)), " | Pg(", trim(Str(I)), ")=", P(I)
       elseif (typen(I)=="C") then
           Np      = Np + 1
           Sd(I)=Cmplx(a, b)
           P(I)    = Real(Sg(I) - Sd(I))
           Q(I)    = Imag(Sg(I) - Sd(I))
           write(60,"(A,A,A,F9.3,A6,A,A,F9.3)") "Noeud consommateur Pd(", trim(Str(I)), ")=", &
           &P(I), " | Qd(", trim(Str(I)), ")=", Q(I)
       else
           slack    = I
           ModV(I)  = a
           Theta(I) = b
           V(I)=cmplx(a*cos(b), a*sin(b))
           write(60, "(A,A,A,F9.3,A,A,A,F9.3)") "Noeud balancier   |V|(", trim(Str(I)), ")=", ModCmplx(V(I)), &
           &" | Theta(", trim(Str(I)), ")=", Argum(V(I))
       end if
    enddo
    !c:saisie de epsilon et de Kmax
    print*, "Lecture des éléments de tolérance/limite des itérations des traitements ..."
    read(50,*) a, Kmax
    eps=a
endif

! Mise en forme des indices de noeuds.
Call Mef
!Call SaveImport
!stop

write(60,*) '------------------------------------------------------------------'
write(60,"(A)") 'le vecteur tension nodal après initialisation ...'
write(60,"(2F6.3,A5)") (V(Indx(I)),"*i pu", I = 1, N)
write(60,*) '------------------------------------------------------------------'

! Subroutine de calcul de la matrice admittance Ya
Call CalcYa

! Formation des matrices G et B
print*, "Initialisation de la matrice admitance du réseau"
YG=real(Ya)
YB=aimag(Ya)
! Initialisation des vecteurs P et Q
P=0.
Q=0.

A2Fmt= "(" // Trim(Str(N)) // "(A1,(F9.3,1X,F9.3),A3))"

! Impression de la matrice admittance (ligne/ligne)
!write(60,*) 'la matrice admittance Y ...'
!do I=1,N
!   write(60, A2Fmt) ("(", Ya(I,J), "*i)", J=1,N)
!end do
!write(60,*) '------------------------------------------------------------------'

write(60,*) 'la matrice conductance G=Real(Y)'
AFmt= "(" // Trim(Str(N)) // "F9.3,1X)"
do I=1,N
   write(60, AFmt) (YG(I,J), J=1,N)
end do

write(60,*) '------------------------------------------------------------------'
write(60,*) 'la matrice susceptance B=Imag(Y)'
do I=1,N
   write(60, AFmt) (YB(I,J), J=1,N)
end do
write(60,*) '------------------------------------------------------------------'

! Alloue dynamiquement la vraie dimension de Ja, Fd et DX
JN = Ng+2*Np
print*, "Initialisation de la matrice Jacobienne du réseau ... "
Allocate (Ja(JN, JN), Fd(JN), DX(JN), X(JN), c(JN, JN)) !, INDX(JN)

!-----------------------------------------------------------------------------------------------------
write(60,"(A28, A, A23, F9.3)") "Nombre max des itérations : ", trim(Str(Kmax)), ", Précision demandée : ", eps

! Affectation et formation de la valeur ESTIME du veteur des inconnues [Teta, |V|] à l'itération k=0
X = 0.
print*, "Initialisation du vecteur INCONNUES du réseau ..."
write(60,"(A39)") "Vecteur inconnues après initialisation ..."
write(60,*) '------------------------------------------------------------------'
do I = 1, N-1
   X(I) = Argum(V(Indx(I+1)))
end do
do I = N, JN
   X(I) = ModCmplx(V(Indx(I - Np + 1)))
end do

do I=1, JN
   if (I .le. N-1) then
       write(60, "(A2,A,A4,F9.3,A)") "X(", trim(Str(I)), ") = ", X(I), " radian(s)"
   else
       write(60, "(A2,A,A4,F9.3,A)") "X(", trim(Str(I)), ") = ", X(I), " pu"
   endif
end do
write(60,*) '=================================================================='

! Début des itérations
Cond = 0 ! Condition d'arrêt relative à Epsilloon
k=1 ! Indice de l'itération
print*, "Début du traitement du réseau ", NName
! Evaluation du moment du début des traitements
call cpu_time ( t1 )
write(*, "(A,f9.3,A2)") "Temps début du traitement évalué à : ", t1, " s"

Do while (Cond==0 .and. k<=Kmax)
    write(60, "(A26, A)") "Début de l'itération, k = ", trim(Str(k))

    !Calcul des composantes de la matrice Jacobienne
    Ja = 0.0
    ! Element de la sous-matrice J1
    do i = 2, N
       do j = 2, N
          if (Indx(i) == Indx(j)) then
              Ja(i-1, j-1) = DP_Dtetaii(Indx(i), N)
          else
              Ja(i-1, j-1) = DP_Dtetaij(Indx(i), Indx(j))
          end if
       end do
    enddo
    ! Element de la sous-matrice J2
    do i = 1, N-1
       do j = N, JN
          if (Indx(i+1) == Indx(j-Np+1)) then
              Ja(i, j) = DP_DModVii(Indx(i+1), N)
          else
              Ja(i, j) = DP_DModVij(Indx(i+1), Indx(j-Np+1))
          end if
       end do
    enddo
    ! Element de la sous-matrice J3
    do i = N, JN
       do j = 1, N-1
          if (Indx(i-Np+1) == Indx(j+1)) then
              Ja(i, j) = DQ_Dtetaii(Indx(i-Np+1), N)
          else
              Ja(i, j) = DQ_Dtetaij(Indx(i-Np+1), Indx(j+1))
          endif
       end do
    enddo
    ! Element de la sous-matrice J4
    do i = N, JN
       do j = N, JN
          if (Indx(i-Np+1) == Indx(j-Np+1)) then
              Ja(i, j) = DQ_DModVii(Indx(i-Np+1), N)
          else
              Ja(i, j) = DQ_DModVij(Indx(i-Np+1), Indx(j-Np+1))
          endif
       end do
    enddo

    write(60,*) '------------------------------------------------------------------'
    ! Impression des elements de la matrice Jacobienne
    A3Fmt= "(" // Trim(Str(JN)) // "F9.3,1X)"
    write(60, *) "Matrice jacobienne à l'itération k=", trim(Str(k))
    do i = 1, JN
        write(60, A3Fmt) (Ja(i, j), j = 1, JN)
    end do

    write(60,*) '------------------------------------------------------------------'
    ! Calcul des fonctions DeltaF
    write(60, *) "Vecteur second membre DF ( itération k :",trim(Str(k)),")"
    Fd = 0
    ! Les DeltaP de tous le noeuds sauf balancier (Ng + Np)
    do I = 1, N-1
       ! Correspond à P (Noeuds de type PV et PQ)
       Fd(I)   = - P(Indx(I+1)) + calcP(Indx(I+1))
    enddo
    ! DeltaQ des neux consommateurs seulement (Np)
    do I = N, JN
           Fd(I) = - Q(Indx(I-Np+1)) + CalcQ(Indx(I-Np+1))
    enddo

    ! Impression des DeltaP, DeltaQ résultants
    write(60,*) "Impression des DeltaP, DeltaQ résultants"
    do I=1, JN
       write(60, "(A3,A,A2,F9.3)") "DF(",trim(Str(I)),")= ",Fd(I)
    end do

    if (CalcWith == "INV") Then
    !    ! Inversion de la matrice Jacobienne, en utilisant La méthode de Gauss-Jordan
        Call InvMat(Ja, c, JN)
        ! Calcul de la solution du système linéaire J.DX = F
        DX = MatMul(c, Fd)
    elseif (CalcWith == "GAUSS") then
        Call GaussMat(Ja, Fd, DX, JN)
    elseif (CalcWith == "LU") then
        Call LEGS(Ja, JN, Fd, DX)
    endif

    ! Calculer les nouvelles valeurs de Theta et ModV
    X = X - DX
!
    write(60,*) '------------------------------------------------------------------'
    write(60,*) "Vecteur corrections à l'itération k= ", trim(Str(k))
    write(60, "(E9.3)") (DX(I), I = 1, JN)
!
    write(60,*) '------------------------------------------------------------------'
    write(60,*) "Vecteur solution intermédiaire à l'itération k= ", trim(Str(k))
    write(60, "(E9.3)") (X(I), I = 1, JN)
!
    ! Ré-établissement des valeurs des vecteurs ModV et de Theta
    do i=1, N-1 ! les N premières valeurs sont pour les arguments Theta
       Theta(i+1) = X(i)
    end do
    do i = N, JN ! Les suivantes sont pour les modules de V
       ModV(i-Np+1) = X(i)
    end do
    ! Ré-évaluation du vecteur tension nodal, à partir des corrections obtenues
    do i = 1, N
       V(i) = cmplx(ModV(i)*cos(Theta(i)), ModV(i)*sin(Theta(i)))
    enddo
!
    ! Déterminer le maximum des DX
    MaxDX=Maxval(DX)
    write(60,*) '------------------------------------------------------------------'
    write (60,"(A8,E9.3,A12,E9.3)") "MaxDX = ",MaxDX, " Epsillon = ",eps

    ! Comparer avec Eps
    if (MaxDX <= eps) then
        ! si petit terminer, imprimer resultats
        Cond=1
        write(60,*) '------------------------------------------------------------------'
        write(60,*) " Vecteur tension nodal à l'itération k= ", trim(Str(k))
        write(60,*) '------------------------------------------------------------------'
        write(60,*) '         V(I) (pu)               |V| (pu)     Theta (Deg)'
        write(60,*) '------------------------------------------------------------------'
        do i = 1, N
           write(60, "(A1,2F12.8,A6, F12.8, 2X, F12.8)") "(",V(i),"*i) - ", ModCmplx(V(i)), Argum(V(i))*180/pi
        enddo

        write(60,*) '------------------------------------------------------------------'
        ! Calcul du transit des puissance P et Q dans le réseau et des pertes avec les dernières corrections
        PL=0.
        QL=0.
        do i=1, N
           if (Typen(i) == "B") Then
                P(i) = CalcP(i)
                Q(i) = calcQ(i)
           elseif (Typen(i) == "G") Then
                Q(i) = CalcQ(i)
           end if
           PL = PL + P(i)
           QL = QL + Q(i)
        enddo

        ! Calcul du transit de puissances Spq
        Spq = 0.
        do i=1, N
           do j=1, N
              if (i /= j) Then
                  Spq(i, j) = (V(i) - V(j))*conjg((V(i) - V(j))*Ya(i, j) - Ys(i, j)*V(i))
                  Spq(j, i) = (V(j) - V(i))*conjg((V(j) - V(i))*Ya(i, j) - Ys(j, i)*V(j))
              end if
           end do
        end do

        ! Impression du transit de puissance et des pertes par transmission
        write(60, "(A,F9.3,A5)") "            Tableau des puissances nodales  (pu), Sbase = ", Sbase, " MVA."
        write(60, *) "------------------------------------------------------------------------------------"
        write(60, *) "|       Si         |      Pgi     |      Qgi      |      Pdi     |       Qdi       |"
        write(60, *) "------------------------------------------------------------------------------------"
        AFmt = "(A2,2F9.3,A1,F14.8,A1,F15.8,A1,F14.8,A1,F17.8,A1)"
        do i=1, N
            write(60, AFmt) "|", Cmplx(P(i),Q(i)), "|", Real(Sg(i)), "|", Imag(Sg(i)), "|", &
            &Real(Sd(i)), "|", Imag(Sd(i)), "|"
        end do
        write(60, "(A)") "          Tableau des puissances transmises  (pu)"
        write(60, *) "------------------------------------------------------------------------------------"
        write(60, *) "| i | j |        Sij (pu)        |        Sij (pu)        |        SLij (pu)       |"
        write(60, *) "------------------------------------------------------------------------------------"
        AFmt = "(A2,I3,A1,I3,A1,2F12.8,A1,2F12.8,A1,2F12.8,A1)"
        do i=1, N
            do j=i+1, N
                write(60, AFmt) "|", i,"|",j,"|",Spq(i,j),"|",Spq(j,i),"|",(Spq(i,j)-Spq(j,i)),"|"
            end do
        end do
        write(60, *) "------------------------------------------------------------------------------------"
        write(60, "(A)") "              Pertes par transmission  (pu)"
        write(60, *) "------------------------------------------------------------------------------------"
        write(60,*) "           ", Cmplx(PL, QL)
        write(60, *) "------------------------------------------------------------------------------------"
        call cpu_time ( t2 )
        write(*, "(A,f6.3)") "Arrêt des traitement à ", t2
        print*, '=================================================================='
        write(60, "(A,f6.3,A)") "Arrêt des traitement à ", t2, " s"
    else    ! sinon refaire le calcul Jacob Ja
        k=k+1
    end if
    !Cond=1
enddo ! Fin de la boucle

! Fin du programme Newton raphson pour le calcul du Load Flow
write(*,"(A,f6.3,A2)") "Durée du traitement est évaluée à ", (t2 - t1), " s"
write(60, "(A,f6.3,A2)") "Durée du traitement est évaluée à ", (t2 - t1), " s"
write(60,*) '=================================================================='
! Impression des résultats ans le fichier de sortie (60)
!
Close(50)
Close(60)

contains

    ! Fonction de calcul de la amatrice admittance Ya
    Subroutine CalcYa()
        Integer :: w, i_, j_

        do w=1,NB
            i_ = Nfrom(w)
            j_ = Nto(w)
            Ya(i_,j_)=-1/Z(i_,j_)
            ! ce ne sont pas les meme valeurs, Yii et Yjj
            Ya(i_,i_)=Ya(i_,i_)+(-Ya(i_,j_)+Ys(i_,j_)/2)
            Ya(j_,j_)=Ya(j_,j_)+(-Ya(i_,j_)+Ys(i_,j_)/2)
            Ya(j_,i_)=Ya(i_,j_)
            Z(j_,i_)=Z(i_,j_)
        end do
    end Subroutine

    ! Fonction du calcul de la puissance active
    Real function calcP(iI)
        integer :: ik, k_
        integer, intent(IN) :: iI
        Real :: S

        S = 0.
        Do ik = 1, N
           k_ = Indx(ik)
           S = S + ModCmplx(V(k_))*(YG(iI,k_)*cos(Teta(iI,k_)) + YB(iI,k_)*sin(Teta(iI,k_)))
        enddo
        calcP = ModCmplx(V(iI))*S
    end function

    ! Fonction du calcul de la puissance réactive
    Real function calcQ(iI)
        integer :: ik, k_
        integer, intent(IN) :: iI
        Real :: S

        S = 0.
        Do ik = 1, N
           k_ = Indx(ik)
           S = S + ModCmplx(V(k_))*(YG(iI,k_)*sin(Teta(iI,k_)) - YB(iI,k_)*cos(Teta(iI,k_)))
        enddo
        calcQ = ModCmplx(V(iI))*S
    end function

    ! Fonction du calcul du module de la tension nodale en un noeud
    Real function ModCmplx(z)
        complex, intent(IN) :: z

        !ModCmplx = (sqrt(Real(z)**2 + imag(z)**2))
        ModCmplx = (cabs(z))
    end function ModCmplx

    ! Fonction du calcul de l'argument de la tension nodale en un noeud
    Real function Argum(z)
        real :: theta

        complex, intent(IN) :: z

        theta = atan(imag(z)/real(z))
        !Argum = (theta*180)/pi
        Argum = theta
    end function
    !
    ! Fonction du calcul de la solution du système AX=B selon GAUSS
    Subroutine GaussMat(Jaa, Bb, Xx, Nn)

        integer, intent(IN) :: Nn
        Real, intent(IN) :: Jaa(:,:), Bb(:)
        Real, allocatable, intent(OUT) :: Xx(:)

        Real, allocatable, Dimension(:,:) :: Aa
        integer :: ii, ij, ik
        Real :: Xsum, xmult
        !character(len=40) :: A3Fmt

        ! Formation de la matrice a = [Jaa | I]
        allocate(Aa(Nn,Nn+1), Xx(Nn))

        Xx=0.
        Aa=0.
        do ii=1, Nn
            do ij=1, Nn+1
               Aa(ii,ij) = Jaa(ii,ij)
            end do
            Aa(ii, Nn+1) = Bb(ii)
        end do

        ! Elimination de Gauss
        do ik = 1, Nn-1
           do ii = ik+1, Nn
              xmult = Aa(ii, ik) / Aa(ik, ik)
              Aa(ii, ik) = 0.
              do ij = ik+1, Nn+1
                 Aa(ii, ij) = Aa(ii, ij) - xmult * Aa(ik, ij)
              end do
           end do
        end do

        Xx(Nn) = Aa(Nn, Nn+1) / Aa(Nn, Nn)
        do ii = Nn-1, 1, -1
           Xsum = 0.
           do ij= ii+1, Nn
              Xsum = Xsum + Aa(ii, ij) * Xx(ij)
           end do
           Xx(ii) = (Aa(ii, Nn+1) - Xsum) / Aa(ii, ii)
        end do

        deallocate(Aa)
    End Subroutine

    ! Procédure d'inversion matrice selon GAUSS-JORDAN
    Subroutine InvMat(Aa, Cc, Nn)
        !implicit none

        integer, Intent(IN) :: Nn
        Real, Intent(INOUT) :: Aa(Nn,Nn)

        Real, Intent(OUT) :: Cc(Nn,Nn)
        Real :: L(Nn,Nn), U(Nn,Nn), b(Nn), d(Nn), Xx(Nn)
        Real :: coeff
        integer ii, ij, ik

        L = 0.0
        U = 0.0
        b = 0.0

        do ik = 1, Nn-1
           do ii = ik+1, Nn
              coeff = Aa(ii,ik) / Aa(ik,ik)
              L(ii,ik) = coeff
              do ij = ik + 1, Nn
                 Aa(ii,ij) = Aa(ii,ij) - coeff*Aa(ik,ij)
              end do
           end do
        end do

        do ii = 1, Nn
          L(ii,ii) = 1.0
        end do

        do ij = 1, Nn
          do ii = 1, ij
            U(ii, ij) = Aa(ii,ij)
          end do
        end do

        do ik = 1, Nn
          b(ik) = 1.0
          d(1) = b(1)
          do ii = 2, Nn
            d(ii) = b(ii)
            do ij = 1, ii-1
              d(ii) = d(ii) - L(ii,ij)*d(ij)
            end do
          end do
          Xx(Nn)=d(Nn)/U(Nn,Nn)
          do ii = Nn-1, 1, -1
            Xx(ii) = d(ii)
            do ij = n, ii+1, -1
              Xx(ii) = Xx(ii) - U(ii, ij)*Xx(ij)
            end do
            Xx(ii) = Xx(ii)/U(ii,ii)
          end do
          do ii = 1, Nn
            Cc(ii, ik) = Xx(ii)
          end do
          b(ik)=0.0
        end do
    end subroutine InvMat
    !
    ! Function de calcul de dP/dTeta ii
    Real function DP_Dtetaii(ii, Nn)
        Integer, intent(IN) :: ii, Nn
        integer :: ij, j_
        Real :: S

        S = 0.
        do ij = 1, Nn
           j_ = Indx(ij)
           if (ii /= j_) then
               S = S + ModCmplx(V(j_))*(YB(ii,j_)*cos(Teta(ii,j_)) - YG(ii,j_)*sin(teta(ii,j_)))
           endif
        enddo
        DP_Dtetaii = S*ModCmplx(V(ii))
    end function

    ! Function de calcul de dP/dTeta ij
    Real function DP_Dtetaij(ii, ij)
        Integer, intent(IN) :: ii, ij

        DP_Dtetaij = ModCmplx(V(ii))*ModCmplx(V(ij))*(YG(ii,ij)*sin(teta(ii,ij)) - YB(ii,ij)*cos(teta(ii,ij)))
    end function

    ! Function de calcul de dQ/dTeta ii
    Real function DQ_Dtetaii(ii, Nn)
        Integer, intent(IN) :: ii, Nn

        Integer :: ij,j_
        Real :: S

        S=0.
        do ij=1,Nn
           j_ = Indx(ij)
           if (ii/=j_) then
               S = S + ModCmplx(V(j_))*(YG(ii,j_)*cos(teta(ii,j_)) + YB(ii,j_)*sin(teta(ii,j_)))
           endif
        enddo
        DQ_Dtetaii = CalcP(ii) - YG(ii,ii)*(ModCmplx(V(ii))**2)
    end function

    ! Function de calcul de dQ/dTeta ij
    Real function DQ_Dtetaij(ii, ij)
        Integer, intent(IN) :: ii, ij

        DQ_Dtetaij = - ModCmplx(V(ii))*ModCmplx(V(ij))*(YG(ii,ij)*cos(teta(ii,ij)) + YB(ii,ij)*sin(teta(ii,ij)))
    end function

    ! Function de calcul de dP/dMod ii
    Real function DP_DModVii(ii, Nn)
        Integer, intent(IN) :: ii, Nn

        Integer :: ij, j_
        Real :: S

        S=0.
        do ij=1,Nn
           j_ = Indx(ij)
           if (ii/=j_) then
               S = S + ModCmplx(V(j_))*(YG(ii,j_)*cos(teta(ii,j_)) + YB(ii,j_)*sin(teta(ii,j_)))
           endif
        enddo
        DP_DModVii = ModCmplx(V(ii))*S + 2*YG(ii,ii)*(ModCmplx(V(ii)))**2
    end function

    ! Function de calcul de dP/dMod ij
    Real function DP_DModVij(ii, ij)
        Integer, intent(IN) :: ii, ij

        DP_DModVij = ModCmplx(V(ii))*ModCmplx(V(ij))*(YG(ii,ij)*cos(teta(ii,ij)) + YB(ii,ij)*sin(teta(ii,ij)))
    end function

    ! Function de calcul de dQ/dMod ii
    Real function DQ_DModVii(ii, Nn)
        Integer, intent(IN) :: ii, Nn

        Integer :: ij, j_
        Real :: S

        S = 0.
        do ij=1, Nn
           j_ = Indx(ij)
           if (ii /= j_) then
               S = S + ModCmplx(V(j_))*(YG(ii,j_)*sin(teta(ii,j_)) - YB(ii,j_)*cos(teta(ii,j_)))
           endif
        enddo
        DQ_DModVii = ModCmplx(V(ii))*S - 2*YB(ii,ii)*(ModCmplx(V(ii)))**2
    end function

    ! Function de calcul de dP/dMod ij
    Real function DQ_DModVij(Ii, Ji)
        Integer, intent(IN) :: Ii, Ji

        DQ_DModVij = ModCmplx(V(Ii))*ModCmplx(V(Ji))*(YG(Ii,Ji)*sin(teta(Ii,Ji)) - YB(Ii,Ji)*cos(teta(Ii,Ji)))
    end function

    ! Fonction de calcul de la diférence des arguments de deux tensions nodale
    Real function Teta(Ii,Ji)
        Integer, intent(IN) :: Ii,Ji

        Teta = Argum(V(Ii)) - Argum(V(Ji))
    end function

    ! Fonction Import de données à partir de fichier CDF (Common Data Format)
    subroutine CDFImport1(f_)
        !implicit none

        Integer, Intent(In) :: f_
        character(len=200) :: a_

        ! Lecture des informations de la première section (Identification du réseau)
        read(f_,'(A)') a_
        !NName = adjustl(a_(11:30))
        NName = adjustl(a_(46:73))
        SBase = Str2float(a_(32: 37))

        ! Lecture de la deuxième section (Pas importante) mention du nombre de neuds dans le réseau.
        read(f_, '(A)') a_
        N = Str2int(a_(30: 46))
    end subroutine CDFImport1

    subroutine CDFImport2(f_)

        Integer, Intent(In) :: f_
        Integer :: ii, t_, inn = 0
        character(len=200) :: a_
        Real :: Pd_, Qd_, Pg_, Qg_

        ! Lecture de la troisième section, Données des noeuds
        loop1 : Do
           read(f_, "(A)") a_
           ii  = Str2Int(a_(1: 4))
           If (ii /= -999) then
               inn = inn +1
               ii  = Str2Int(a_(1: 4))
               t_ = Str2Int(a_(25: 26))
               If (t_ .le. 1) then
                  Typen(ii) = "C"
               elseif (t_ == 2) then
                  Typen(ii) = "G"
               else
                  Typen(ii) = "B"
               endif
               Pd_   = Str2Float(a_(41: 49))
               Qd_   = Str2Float(a_(50: 59))
               Pg_   = Str2Float(a_(60: 67))
               Qg_   = Str2Float(a_(68: 75))
               V(ii) = Str2Float(a_(28: 33))
               Sd(ii) = Cmplx(Pd_/Sbase, Qd_/Sbase)
               Sg(ii) = Cmplx(Pg_/Sbase, Qg_/Sbase)
           else
               exit loop1
           endif
           ! Affectation du nombre de noeuds, exact.
           N = inn
        enddo loop1

        ! Après detection de la fin de la section 3, lecture de la section quatrième : Données es lignes
        read(f_, '(A)') a_
        Nb = Str2Int(a_(30: 46))
    end subroutine CDFImport2

    subroutine CDFImport3(f_)
        Integer, Intent(In) :: f_
        Integer :: ii, jj, inb = 0
        character(len=200) :: a_
        Real :: r_ , x_

        loop2: Do
           read(f_, "(A)") a_
           ii  = Str2Int(a_(1: 4))
           If (ii /= -999) then
               inb = inb + 1
               Nfrom(inb) = ii
               jj = Str2Int(a_(6: 9))
               Nto(inb) = jj

               r_ = Str2Float(a_(20: 29))
               x_ = Str2Float(a_(30: 40))
               Z(ii, jj) = Cmplx(r_, x_)
               Ys(ii, jj) = Cmplx(0., Str2Float(a_(41: 50))/2)
           else
               exit loop2
           endif
        enddo loop2
    end subroutine CDFImport3

    Subroutine Mef()

    Integer :: i_
    !character(30) :: fm

        Indx = 0
        If (LookForCDF) Then
            ! Domaine du noeud Bilan B, un seul noeud.
            do i_= 1, N
               If (Typen(i_) == "B") Then
                   Indx(1) = i_
                   !Indx(i_) = 1
               endif
            enddo
            ! Domaine des noeuds générateurs, à partir de 2 (1 réservé au noeud Bilan).
            Ng=0
            do i_= 1, N
               if (Typen(i_) == "G") then
                   Ng = Ng +1
                   Indx(Ng + 1) = i_
                   !Indx(i_) = Ng + 1
               endif
            end do
            ! Domaine des neuds consommateurs, à partir de (Ng+1).
            Np=0
            do i_= 1, N
               If (Typen(i_) == "C") Then
                   Np = Np +1
                   Indx(Np + (Ng + 1)) = i_
                   !Indx(i_) = Np + (Ng + 1)
               end if
            end do
        else
           do i_ = 1, N
              Indx(i_) = i_
           enddo
        end if

!        fm = "(" // Str(Ng+1) // "(I3,2X))"
!        write(*, fm) (Indx(i_), i_=1,N)

    end Subroutine Mef

    Subroutine SaveImport(f_)
        implicit none

        character(len=100), intent(out) :: f_
        integer :: i_, j_, k_

        f_ = adjustl(trim(f_(1: len(f_)-4)) // trim(GenerateRndCar()) // ".dat")
        write(*,*) "Ecriture temporaire dans fichier conforme standard NMSS ..."

        open (unit=70, file=f_, action="write")
        ! Retrait de tous les espaces de Nname
        Call removesp(Nname)

        write(70, "(A,1X,I4,1X,I4,1X,F6.1)") Nname, N, NB, Sbase
        do i_=1,NB
!           j_ = Indx(Nfrom(i_))
!           k_ = Indx(Nto(i_))
           j_ = Nfrom(i_)
           k_ = Nto(i_)
           write(70, "(I4,1X,I4,1X,2F10.4,1X,2F10.4, I3,1X,I3)") Indx(j_), Indx(k_), Z(j_, k_), Ys(j_, k_)
        enddo
        do j_ = 1, N
           i_ = Indx(j_)
           if (Typen(i_) == "B") Then
                write(70, "(A1,1X,F9.3, 1X, F9.3, 1X, I4)") Typen(i_), ModCmplx(V(i_)), Argum(V(i_)), i_
           elseif (Typen(i_) == "G") Then
                write(70, "(A1,1X,F9.3, 1X, F9.3, 1X, I4)") Typen(i_), ModCmplx(V(i_)), Real(Sg(i_)), i_
           elseif (Typen(i_) == "C") Then
                write(70, "(A1,1X,F9.3, 1X, F9.3, 1X, I4)") Typen(i_), Real(Sd(i_)), Imag(Sd(i_)), i_
           endif
        enddo
        write(70, "(F10.6,1X,I4)") eps, Kmax
        print*, "Fermeture du fichier temporaire ..."
        close(70)
    End Subroutine

    subroutine Removesp(str)
    implicit none

    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str))::outstr
    integer :: i_, k_, ich, lenstr

    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    k_=0

    do i_=1,lenstr
      ch=str(i_:i_)
      ich=iachar(ch)
      select case(ich)
        case(0:32)  ! space, tab, or control character
             cycle
        case(33:)
          k_=k_+1
          outstr(k_:k_)=ch
      end select
    end do

    str=adjustl(outstr)
    end subroutine Removesp

    character(len=30) Function GenerateRndCar()
        implicit none

        real :: r_ = 2008.1960

        call RANDOM_NUMBER(r_)
        GenerateRndCar = Adjustl(Str(floor(r_*10000)))

    end function

end program newtonrap
