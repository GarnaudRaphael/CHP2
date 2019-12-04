module gradientconjugue
  use matrices
  implicit none
  contains


        subroutine gc(nnz,colonnes,ic,b,NxNy,Nx,Ny,xk)
          implicit none
          double precision , dimension(:),allocatable ::b,xk, p, Ap, rk, produitmatriciel,nnz,colonnes,ic,pk,pkA


          double precision :: alpha, rr, rr1, pAp,produitvectoriel,alphak,betak, produitvectoriel2,produitvectoriel3,rkrk,pkAp
          integer :: j,i,k,statinfo,Np,me, NxNy,Nx,Ny,i1,iN,he,e,d
          allocate (rk(NxNy))
          allocate (pkA(NxNy))
          allocate (pk(NxNy))
          allocate(Ap(NxNy))
          e=0
          d=4
          xk=0!x0=0
          rk = b
          pk = b
          k=0
          do while(norm2(rk)>10**(-e))
            call produitvec(Nx,Ny,rk,rk,rkrk,statinfo,me,i1,iN)
            call produitmat(nnz,colonnes,ic,pk,Nx,Ny,pkA,i1,iN,statinfo)
            call produitvec(Nx,Ny,pkA,pk,pkAp,statinfo,me,i1,iN)
            alphak=rkrk/pkAp
            xk=xk+alphak*pk
            !call produitmat(nnz,colonnes,ic,xk,Nx,Ny,Ap,i1,iN,statinfo)
            rk=rk-alphak*pkA
            call produitvec(Nx,Ny,rk,rk,produitvectoriel,statinfo,me,i1,iN)
            betak=produitvectoriel/rkrk
            !print*," pkAp = ",pkAp
            pk=rk+betak*pk
            k=k+1
          !  print*,rk
            if(k>10**d)then
              print*,"erreur"
              exit
            end if
          end do
          !print*," k = ",k," x = ",xk
          deallocate(pkA)
          deallocate(Ap)
          deallocate(rk)
          deallocate(pk)


      end subroutine gc


    subroutine produitmat(nnz,colonnes,ic,xx,Nx,Ny,produitmatriciel,i1,iN,statinfo)
      implicit none
      integer::Nx,Ny,i,i1,iN,he,me,statinfo,Np,j1,jN,k1,k2,k3,j
      double precision,dimension(:,:),allocatable::A
      double precision,dimension(:),allocatable::xx,pm,produitmatriciel,nnz,colonnes,ic

      !allocate(produitmatriciel(1:(Nx+2)*(Ny+2)))
      produitmatriciel=0
      !print*,"xx = ",xx
      !open(unit=7, file="bm1.txt",form="formatted",access="sequential")
      do i=1,(Nx+2)*(Ny+2)
        produitmatriciel(i)=0
        k1=ic(i)!indice du première éléments non nul de la ligne
        k2=ic(i+1)-1!indice du dernière éléement de la lgine non nul
        do j=k1,k2
          k3=(colonnes(j))!colonnes de la matrice qui correspond à l'indice dans x
          produitmatriciel(i)=produitmatriciel(i) +nnz(j)*xx(k3)!j correponds à l'indice de la colonnes considérée
        !  write(7,*)produitmatriciel(i),nnz(j),k3,xx(k3)
        end do
      end do
      !close(7)
    end subroutine produitmat

    subroutine produitvec(Nx,Ny,v,vv,produitvectoriel,statinfo,me,i1,iN)
      implicit none
      integer::Nx,Ny,	i1,iN,i,me,statinfo,Np
      double precision,dimension((nx+2)*(ny+2))::v,vv
      double precision:: pv,produitvectoriel

      produitvectoriel=0
      do i=1,(Nx+2)*(Ny+2)
        produitvectoriel= produitvectoriel + v(i)*vv(i)
      end do

    end subroutine produitvec




end module gradientconjugue
