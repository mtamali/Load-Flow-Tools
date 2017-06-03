Module ToolsLib

contains

    ! Fonction du calcul du module de la tension nodale en un noeud
    Real function ModCmplx(z)
        complex, intent(IN) :: z

        !ModCmplx = (sqrt(Real(z)**2 + Aimag(z)**2))
        ModCmplx = (cabs(z))
    end function ModCmplx

    ! Fonction du calcul de l'argument de la tension nodale en un noeud
    Real function Argum(z)
        real :: theta

        complex, intent(IN) :: z

        theta = Atan(Aimag(z)/real(z))
        !Argum = (theta*180)/pi
        Argum = theta
    end function

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
            do ij = Nn, ii+1, -1
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

    ! Gauss-Jordan avec pivotation de lignes
    SUBROUTINE GaussJordanPivation (A_,Nn,B_,Xx)
      Implicit none

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
      CALL PivotationLU (Aa,Nn,INDX)
      DO iI = 1, Nn-1
        DO iJ = iI+1, Nn
          Bb(INDX(iJ)) = Bb(INDX(iJ))-Aa(INDX(iJ),iI)*Bb(INDX(iI))
        END DO
      END DO

      Xx(Nn) = Bb(INDX(Nn))/Aa(INDX(Nn),Nn)
      DO iI = Nn-1, 1, -1
        Xx(iI) = Bb(INDX(iI))
        DO iJ = iI+1, Nn
          Xx(iI) = Xx(iI)-Aa(INDX(iI),iJ)*Xx(iJ)
        END DO
        Xx(iI) =  Xx(iI)/Aa(INDX(iI),iI)
      END DO
    END SUBROUTINE GaussJordanPivation

    ! Subroutine de pivotation de lignes
    SUBROUTINE PivotationLU (Aa,Nn,INDX)
      Implicit none
      INTEGER, INTENT (IN) :: Nn
      INTEGER :: iI, iJ, iK = 0, ITMP
      INTEGER, INTENT (OUT), DIMENSION (Nn) :: INDX
      REAL :: C1, PI, PI1,PJ
      REAL, INTENT (INOUT), DIMENSION (Nn,Nn) :: Aa
      REAL, DIMENSION (Nn) :: C

      DO iI = 1, Nn
        INDX(iI) = iI
      END DO

      DO iI = 1, Nn
        C1= 0.0
        DO iJ = 1, Nn
          C1 = AMAX1(C1,ABS(Aa(iI,iJ)))
        END DO
        C(iI) = C1
      END DO

      DO iJ = 1, Nn-1
        PI1 = 0.0
        DO iI = iJ, Nn
          PI = ABS(Aa(INDX(iI),iJ))/C(INDX(iI))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            iK   = iI
          ENDIF
        END DO

        ITMP     = INDX(iJ)
        INDX(iJ) = INDX(iK)
        INDX(iK) = ITMP

        DO iI = iJ+1, Nn
          PJ  = Aa(INDX(iI),iJ)/Aa(INDX(iJ),iJ)
          Aa(INDX(iI),iJ) = PJ
          DO iK = iJ+1, Nn
             Aa(INDX(iI),iK) = Aa(INDX(iI),iK)-PJ*Aa(INDX(iJ),iK)
          END DO
        END DO
      END DO
    END SUBROUTINE PivotationLU

    ! Convertir un intier en chaine de caractères.
    character(len=20) function rStr(rk)
    implicit none

    real, intent(in) :: rk

    write(rStr, *) rk
    rStr = adjustl(rStr)
    end function rStr

    ! Convertir un intier en chaine de caractères.
    character(len=20) function Str(ik)
    implicit none

    integer, intent(in) :: ik

    write(Str, *) ik
    Str = adjustl(str)
    end function Str

    ! Convertit une chaine de caractères (numerique) en entier
    Integer function Str2int(str_)
       implicit none

       character(len=*), intent(in) :: str_
       integer :: int_
       read(str_, *)  int_
       Str2int = int_
    end function Str2int

    ! Convertit une chaine de caractères (Réelle) en réel
    Real function Str2float(str_)
       implicit none

       character(len=*), intent(in) :: str_
       Real :: f_

       read(str_, "(F9.3)")  f_
       Str2float = f_
    end function Str2float

    ! Subroutine de nettoyage de caractère spéciaux d'une chaine de caractères.
    Subroutine Removesp(str)
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

    ! Génére une chaine de caractère (à partir de nombres aléatoires)
    character(len=30) Function GenerateRndCar()
        implicit none

        real :: r_ = 2008.1960

        call RANDOM_NUMBER(r_)
        GenerateRndCar = Adjustl(Str(floor(r_*10000)))

    end function GenerateRndCar
End module ToolsLib
