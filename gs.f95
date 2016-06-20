      program GaussSeidel
        implicit none

      ! Diclarations
        character :: FileImport, LINE, NetDataDate,NetOrigName,NetSeason,&
                &NetCaseId,FData, CharNbBus, CharNbBranch
        character(len=106) :: FMT1
        character(len=124) :: FMT2
        character,dimension(500) :: BusName
        character(len=32)::NetDataFileIn, NetDataFileOut, arg, NetFileFormat
        character(len=100) :: OutMessage
        integer,dimension(500) :: NbNode,LFAreaN,LZoneN,Ntype,Ib,Jb,LFAreaB,&
                &LZoneB,Circ,BType,LineRate1,LineRate2,LineRate3,BCntl,BSide
        integer :: N,I,J,Kmax,S,K,M,W,NB,NetYear,ii,jj
        real,dimension(500) :: FV,FDelta,Pd,Qd,Pg,Qg,Vb,Vspec,Qmax,Qmin,&
                &ChargLine,TrRatio,TrPhase,Tapmin,Tapmax,SteSize,Vmin,Vmax
        real::a,b,d,e,NetMVABase
        real,Dimension(500):: P,Q
        complex,dimension(500) :: Ybs
        complex,dimension(500,500) :: Z,Ya,Ys
        complex,Dimension(500)::V
        complex::SS,ST,Va,eps,DVamax

      ! arg represente les arguments transmits à l'application
        NetDataFileOut=""
      if (iargc()/=0 .and. iargc()>0) then
         CALL getarg(1, arg)
         NetDataFileIn=arg
         if (iargc()>1) then
             CALL getarg(2, arg)
             NetDataFileOut=arg
         endif
         if (iargc()>2) then
             CALL getarg(3, arg)
             NetFileFormat=arg
             NetFileFormat=trim(NetFileFormat)
         endif
         DO i = 4, iargc()
            CALL getarg(i, arg)
            selectcase (arg)
                Case ("GS")
                Case ("NR")
                Case ("RL")
                Case ("FD")
            endselect
         enddo
      ! NetDataFile : Nom du fichier de donnees du reseau
        open(unit=100,file=NetDataFileIn,status='old',action='read')
      ! Data.res : Nom du fichier de resultats du reseau
        if (NetDataFileOut=="") Then
            open(unit=200,file='data.html')
        else
            NetDataFileOut=trim(NetDataFileOut) // ".html"
            open(unit=200,file=NetDataFileOut)
        endif
        ! Preparation du fichier à soumettre : RESULTATS et insertion des BALISES HTML necessaire
        write(200,*) "<!DOCTYPE html>"
        write(200,*) "<html lang='fr'>"
        write(200,*) "<head>"
        write(200,*) "<meta charset='utf-8'>"
        write(200,*) "<title>Net Platform Version 1.0 | UTMB, SimulIA Team ENERGARID Lab., 2015</title>"
        write(200,*) "<meta name='description' content='Electrical Network Platform 1.0.'>"
        write(200,*) "<meta name='keywords' content='Electrical, Network, GS, NR, OPF, Load Flow'>"
        write(200,*) "<meta name='author' content='Etudiants Master ETT 2Annee, 2015 | TAMALI Mohammed'>"
        write(200,*) "<meta name='viewport' content='width=device-width, initial-scale=1.0'>"
        write(200,*) "<link rel='stylesheet' href='style.css' type='text/css' />"
        write(200,*) "<script type='text/javascript' src='jquery.min.js'></script>"
        write(200,*) "</head><body class='blurBg-true' style='background-color:#ebdea9'>"
        write(200,*) "<div style='margin-left:50px; width:80%;'>"
        write(200,*) "<H1>DONNEES DU RESEAU [",trim(NetDataFileIn),"]</H1><hr><br >"

        if(NetFileFormat=='CDF')then
            call CDFImport(100)
        else
            Call SimulIANetRead(100,200)
        endif

        ! Debut du traitement selon la methode GAUSS-SEIDEL
        K=0
        write(200,*) '<H3>RESULTATS INTERMEDIAIRES</H3><br />'
        100 write(200,'(A,I4,A)') 'Iteration : ',K,'<br />'
        ! Boucle principale du calcul des element du vecteur tensions nodales
        do I=1,N
           ! On saute le noeud BILAN numero s
           if (I/=S)then
               ! Calcul de la somme des courant de contributions  SS
               SS=(0.,0.)
               do J=1,N
                  if(I/=J)then
                     SS=SS+Ya(I,J)*V(J)
                  endif
               enddo
               ! Si TYPE du noeud courant I est Producteur (1), on calcul la puissance reactive Q(i)
               if (Ntype(I)==2) then
                  ST=cmplx(0.,0.)
                  do J=1,N
                     ST=ST+Ya(I,J)*V(J)
                  enddo
                  Q(I)=aimag(V(I)*ST)
               endif
               ! On sauvegarde l'ancienne valeur de V(i) dans Va
               Va=V(I)
               ! On calcul la nouvelle valeur de V(i)
               V(I)=(1/Ya(I,I))*(cmplx(P(I),-Q(I))/V(I)-SS)
               ! On compare la difference des deux valeur Va et V(i) à Delta Vmax, si superieur on ecrase Delta Vmax
               if (real(Va-V(I))>real(DVamax) .and. aimag(Va-V(I))>aimag(DVamax)) then
                  DVamax=(Va-V(I))
               endif
           endif
        enddo
        ! Fin de la boucle principale
        ! On imprime les valeurs intermediaires de V
        write(200,*) '<H3>TENSIONS NODALES (TMP)</H3><br />'
        do M=1,N
            write(200,'(2F12.6,A)') V(M),'<br />'
        enddo
        if (cabs(DVamax)>cabs(eps)) then
           K=K+1
           if (K<=Kmax) then
              Go to 100
           else
              write(200,*) '<H3>Le calcul a diverge, le maximum des iterations ',&
                    &Kmax,' a ete atteint</H3><br />'
              GoTo 2
           endif
        else
           ! Calcul de P(s) et Q(s), du noeud BILAN
           ST=cmplx(0.,0.)
           do I=1,N
              ST=ST+Ya(S,I)*V(S)
           enddo
           P(S)=Real(V(S)*ST)
           Q(S)=aimag(V(S)*ST)
           write(200,*) '<H3>PUISSANCE APPARENTE NODALES DU NOEUD BILAN</H3><br />'
           write(200,'(2F14.6,A)') P(S),Q(S),'<br />'
           write(200,*) '<H2>RESULTATS FINAUX (en p.u.)</H2><br />'
           write(200,'(A,I4,A)') "<H3>Pour un nombre d'iterations gal à ",K,"</H3><br />"
           write(200,*) '<hr>'
           do I=1,N
              write(200,'(2F12.6,A)') V(I),'<br />'
           enddo
           ! calcul des puissance des pertes par transmission PL et QL
           SS=cmplx(0.,0.)
           ST=cmplx(0.,0.)
           do I=1,N
              SS=SS+P(S)
              ST=ST+Q(S)
           enddo
           write(200,*) '<hr>'
           write(200,*) '<H3>PERTES PAR TRANSMISSION PL et QL (en p.u.)</H3><br />'
           write(200,'(2F12.6,A)') SS, ST,'<br />'
        endif
 2      write(200,*) "</div></body>"
        write(200,*) "</html>"
      else
         OutMessage="Donner au moins un argument, NOM DU FICHIER\n"
         open(unit=300,file='Error')
         write(300,*) OutMessage
      endif

      CONTAINS

      Subroutine SimulIANetRead(InFile, OutFile)
            implicit none
            integer,intent(In) :: InFile, OutFile
            ! Lecture des nombres de noeuds et de lignes
          ! saisie de la dimension reseau
            Read(InFile,'(2I2)') N,NB
            ! Eciture de N et de Nb dans le fichier à soumettre
            write(OutFile,*)'<H3>le nombre =',N,' le nombre de branches=',NB,'</H3><br />'
            ! Initialisation des elements de la matrice admittance
            do I=1,N
               do J=1,N
                  Ya(I,J)=cmplx(0.,0.)
               enddo
            enddo
            ! Lecture et ecriture des donnees des lignes
            write(OutFile,*) '<H3>La matrice admittance, elements extra-diagonaux</H3><br />'
            write(OutFile,*) '  Ni  Nj        Ya[i,j]<br />'
            do w=1,NB
               Read(InFile,*) I,J,a,b,d,e
               Z(I,J) = cmplx(a,b)
               Ys(I,J)= cmplx(d,e)
               if(real(Z(I,J))>1.0e-6 .and. aimag(Z(I,J))>1.0e-6)then
                  Ya(I,J)=-1/Z(I,J)
               else
                  Ya(I,J)=cmplx(1.0e-6, 1.0e-6)
               endif
               Ya(J,I)=Ya(I,J)
               Z(J,I)=Z(I,J)
               write(OutFile, '(2I4,2F12.6,A)') I,J,Ya(I,J),'<br />'
            enddo
            ! Calcul les element de diagonale
            write(OutFile,*)'<H3>La matrice admittance, les elements diagonaux</H3><br />'
            do I=1,N
                Ya(I,I)=(0.,0.)
                 do J=1,N
                   if(Z(I,J)/=0)then
                      Ya(I,I)=Ya(I,I)+(-Ya(I,J)+Ys(I,J)/2)
                   endif
                 enddo
                 write(OutFile, '(2I4,2F12.6,A)') I,I,Ya(I,I),'<br />'
            enddo
            ! Initialisation de Delta_Vmax
            DVamax=cmplx(-100000, -100000)
            ! Initialisation du vecteur tensions nodales
            do I=1,N
               V(I)=cmplx(1.,0.)
            enddo
            ! Ici commence la saisie des tableaux de planification
            ! Saisie numero nud balancier

            Read(InFile,*),S
            write(OutFile,*) '<H3>Numero du noeud BILAN ',S,'</H3><br />'
            Ntype(S)=0
            write(OutFile,*) '<H3>TABLEAU DE PLANIFICATION</H3><br />'
            write(OutFile,*) ' TYPE      P         Q/|V|'
            ! Saisie du vecteur type des nud
            do I = 1,N
               if (I/=S) then
                  Read(InFile,*) Ntype(I)
                  ! Si le nud est de type consommateur : 2 sinon Producteur : 1, 0 pour le noeud Bilan
                  if (Ntype(I) <= 1) then
                      Read(InFile,*) P(I),Q(I)
                      write(OutFile,'(I4,2F12.6,A)') Ntype(I),P(I),Q(I),'<br />'
                  endif
                  if (Ntype(I) == 2) then
                      Read(InFile,*) a,b,P(I)
                      V(I)=cmplx(a,b)
                      write(OutFile,'(I4,2F12.6,A)') Ntype(I),P(I),cabs(V(I)),'<br />'
                  endif
               endif
            enddo
            ! Saisie de Epsilon,Kmax
            Read(InFile,*)a,b,Kmax
            eps=cmplx(a,b)
            write(OutFile,*) '<H3>   Epsillon             Kmax</H3><br />'
            write(OutFile,'(A,2F10.6,A,I10)') '(',eps,')', Kmax
        end subroutine SimulIANetRead


      ! Processus d'importation à partir de fichier IEEE CDF
      subroutine CDFImport(f)
         implicit none
         integer,Intent(In):: f

         !FileImport=trim(f)
         !open(unit=10,file=FileImport,status='old',action='read')
         ! Lecture des informations de la première section (Identification du réseau)
         read(f,'(1X,A8,1X,A19,2X,F6.1,1X,I4,1X,A1,1X,A)') NetDataDate,NetOrigName,&
                &NetMVABase,NetYear, NetSeason, NetCaseId
         ! Lecture de la deuxième section (Pas importante) mention du nombre de neuds dans le réseau.
         read(f, '(A16,1X,A)') FData, CharNbBus
         CharNbBus=trim(CharNbBus)
         ! Lecture de la troisième section, Données des noeuds
         FMT1='(2X,A10,1X,I2,I3,1X,I2,1X,F6.4,F5.2,F8.2,F9.2,F7.2,F8.2,1X,&
                &F7.2,1X,F6.4,F7.2,F7.2,2F8.4,1X,I4,1X,I4)'
         i=0
         loop1: Do
            read(f, '(I4)'), w
            If (w/=-999) then
                i=i+1
                NbNode(i)=w
                read(f, FMT1) BusName(i),LFAreaN(i),LZoneN(i),Ntype(i),&
                        &FV(i),FDelta(i),Pd(i),Qd(i),Pg(i),Qg(i),Vb(i),&
                        &Vspec(i),Qmax(i),Qmin(i),Ybs(i),w
            else
                exit loop1
            endif
         enddo loop1
         ! Après detection de la fin de la section 3, lecture de la section quatrième : Données es lignes
         read(f, '(A19,1X,A)') FData, CharNbBranch
         CharNbBranch=trim(CharNbBranch)
         FMT2='(1X,I2,I2,2X,I1,1X,F10.6,F11.6,F9.5,F6.0,1X,F6.0,1X,F6.0,&
                &1X,I4,1X,I1,2X,F6.4,1X,F7.2,F7.4,F6.4,1X,F6.5,1X,2F7.4,2X,I4)'
         i=0
         loop2: Do
            read(f, '(I4)'), w
            If (w/=-999) then
                i=i+1
                Ib(i)=w
                read(f, '(1X,I4)'),Jb(i)
                ii=Ib(i)
                jj=Jb(i)
                read(f, FMT2) LFAreaB(i),LZoneB(i),Circ(i),BType(i),Z(ii,jj),&
                &Ys(ii,jj),LineRate1(i),LineRate2(i),LineRate3(i),BCntl(i),&
                &BSide(i),TrRatio(i),TrPhase(i),Tapmin(i),Tapmax(i),SteSize(i),Vmin(i),Vmax(i)
            else
                exit loop2
            endif
         enddo loop2
      end subroutine CDFImport
      end program

