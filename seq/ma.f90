module matrices
		use fonctions
		implicit none



		contains

			subroutine matrice(Nx,Ny,hx,hy,nnz,colonnes,ic,D,i1,iN,dt)
			use mpi
			implicit none
			integer,intent(in)::Nx,Ny
			integer::i,j,cont,diim,Np,me,statinfo,i1,iN,n1,n2,n3
			double precision,intent(in) :: hx,hy,dt,D
			double precision,dimension(:),allocatable::nnz,colonnes,ic
			double precision,dimension(:,:),allocatable::A,AA

			double precision,dimension(Nx+2)::x
			double precision, dimension(Ny+2)::y




				call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
				call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)

				allocate (AA(5,(Nx+2)*(Ny+2)))
				allocate (A((Nx+2)*(Ny+2),(Nx+2)*(Ny+2)))




!print*,D
! construction des 5 diagonales
				AA=0.d0
				AA(3,:)=1.d0
				do j=1,Ny
					do i=1,Nx
						AA(2,i+j*(Nx+2))=-dt*D*(1.d0/(hx**2))
					end do
					do i=3,Nx+2
						AA(4,i+j*(Nx+2))=-dt*D*(1.d0/(hx**2))
					end do
					do i=2,Nx+1
						AA(1,i+(j-1)*(Nx+2))=-dt*D*(1.d0/(hy**2))
						AA(3,i+j*(Nx+2))=1.d0-dt*D*(-2.d0/(hx**2) -2.d0/(hy**2))
						AA(5,i+(j+1)*(Nx+2))=-dt*D*(1.d0/(hy**2))
					end do
				end do


! stockage de ces 5 diagonales dans la matrice de taille (Nx+2)*(Ny+2)
				A=0.d0
				A(1,1)=AA(3,1)
				A((Nx+2)*(Ny+2),(Nx+2)*(Ny+2))=AA(3,(Nx+2)*(Ny+2))
				DO i=2,(Nx+2)*(Ny+2)-1
					A(i,i)=AA(3,i)	;A(i,i-1)=AA(2,i-1) ; A(i,i+1)=AA(4,i+1)
					if (i>=(Nx+2+1)) then
						A(i,i-Nx-2)=AA(1,i-Nx-2)
					end if
					if (i<=(Nx+2)*(Ny+2)-Nx-2-1) then
						A(i,i+Nx+2)=AA(5,i+Nx+2)
					end if

				ENDDO



				deallocate(AA)
				call laplacien(Nx,Ny,hx,hy,dt,D,A)


! nombre d'éléments non nuls dans la matrice


				diim=0
				do i=i1,iN
					do j=1,(Nx+2)*(Ny+2)
						if (A(i,j)/=0) then
							diim=diim+1
						end if
					end do
				end do
! construction de la matrice csr
				allocate (nnz(diim))
				allocate (colonnes(diim))
				cont=0
				allocate (ic( iN-i1+2))

				do i=i1,iN
					ic(i)=(cont+1)
					do j=1,(Nx+2)*(Ny+2)
						if (A(i,j)/=0) then
							cont=cont+1
							nnz(cont)=A(i,j)
							colonnes(cont)=j
						end if
					end do
				end do
				ic(iN-i1+2)=cont+1

!				n3=size(nnz)
!				open(unit=5, file="ic.txt",form="formatted",access="sequential")
!					do i=1,n3
!							write(5,*)ic(i)
!					enddo
!					close(5)
				deallocate(A)


			end subroutine matrice

			subroutine laplacien(Nx,Ny,dx,dy,dt,D,M)
			  implicit none
			  integer :: i,Nx,Ny,n
			  double precision :: h,dx,dy,dt,D
			  double precision, dimension(:,:), allocatable :: A,B
			  double precision, dimension(:,:), allocatable :: M,K
			  allocate(A(1:Nx+2,1:Nx+2));allocate(B(1:Nx+2,1:Nx+2));allocate(K(1:Nx+2,1:3*(Nx+2)));!allocate(M(1:(Nx+2)*(Ny+2),1:(Nx+2)*(Ny+2)));
			  M=0.d0
			  A=0.d0
			  K=0.d0
				B=0.d0
			  A(1,1)=1.0
			  do i=2,Nx+1
			    A(i,i-1)=-dt*(1.0/(dx*dy))
			    A(i,i)=1+2*dt*(1.0/(dx**2))+2*dt*(1.0/(dy**2))
			    A(i,i+1)=-dt*(1.0/(dx*dy))
			  enddo
			  A(Nx+2,Nx+2)=1.0
			  do i=2,Nx+1
			    B(i,i)=-dt*1.0/(dx*dy)
			  enddo
			  !open(unit=2, file="A.txt",form="formatted",access="sequential")
			  !do i=1,n+1
			  !  write(2,*)A(i,:)
			  !enddo
			  !close(2)0
			  K(:,1:Nx+2)=B
			  K(:,Nx+3:2*(Nx+2))=A
			  K(:,2*(Nx+2)+1:3*(Nx+2))=B

			  do i=1,Ny
			    !M(i*(Nx+1)+1:(i+1)*(Nx+1),(i-1)*(Nx+1)+1:(i+2)*(Nx+1))=K
					M(i*(Nx+2)+1:(i+1)*(Nx+2),(i-1)*(Nx+2)+1:(i+2)*(Nx+2))=K
			  enddo
				do i=1,Nx+2
						M(i,i)=1.0
						M(i+(Ny+1)*(Nx+2),i+(Ny+1)*(Nx+2))=1.0
				enddo
			  !M(:(n+1)**2,n*(n+1):(n+1)**2)=1
			  !print*,M
!				open(unit=6, file="M.txt",form="formatted",access="sequential")
!					do i=1,(Nx+2)*(Ny+2)
!							write(6,*)M(i,:)
!					enddo
!					close(6)
			  deallocate(A)
			  deallocate(B)
			  deallocate(K)
			return
			end subroutine


			subroutine laplaciencsr(hx,hy,dt,ic,colonnes,nnz,Nx,Ny)
			  implicit none
			  integer :: i,n,Nx,Ny,j,Njic,Nj,Njc,i1,iN,d,n1,n2,n3
			  double precision :: dx,dy,dt,hx,hy
				double precision, dimension(5) ::fg
			  double precision, dimension(:), allocatable :: colonnes,nnz,ic,ic1
				allocate(ic(1:nx*ny+1))
				allocate(nnz(1:(2*(2*3+4*(Nx-2))+(Ny-2)*(2*4+5*(Nx-2)))))
				allocate(colonnes(1:(2*(2*3+4*(Nx-2))+(Ny-2)*(2*4+5*(Nx-2)))))
				nnz=0.d0
				ic=0.d0
				colonnes=0.d0

				!disparition de g10 et de g01
				ic(1)=1
				nnz(1:3)=[1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2)),-1*dt*(1/(hy**2))]
				colonnes(1:3)=[1,2,2+nx-1]
				ic(2)=ic(1)+3
				!k=i*nx+j
				Nj=4
				Njic=2
				do i=1,nx-2
					!on se balade de l'indice 2 à l'indice nx-1, donc on est en i+1 et non en i
					nnz((i-1)*4+Nj:i*4+nj)=[-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2)),-1*dt*(1/(hy**2))]
					ic(1+Njic)=ic(Njic)+4
					colonnes((i-1)*4+Nj:i*4+nj)=[i+1-1,i+1,i+1+1,i+1+nx]
					Njic=Njic+1
				end do
				Nj=Nj+4*(nx-2)
				nnz(Nj:Nj+2)=[-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-1*dt*(1/(hy**2))]
				colonnes(Nj:Nj+2)=[nx-1,nx,nx+nx]
				ic(Njic+1)=ic(Njic)+3 !i=i+1, à la fin i=nx-1
				Nj=Nj+3
				Njic=Njic+1
				!print*,"njic = ",njic
				fg=[-1*dt*(1/(hy**2)),-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2)),-1*dt*(1/(hy**2))]
				do j=1,ny-2
						nnz(Nj:Nj+3)=[-1*dt*(1/(hy**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2)),-1*dt*(1/(hy**2))]
						colonnes(Nj:Nj+3)=[j*nx+1-nx,j*nx+1,j*nx+1+1,j*nx+1+nx]
						ic(Njic+1)=ic(Njic)+4
						!k=i*nx+j
						Nj=Nj+4
						Njic=Njic+1
						do i=1,nx-2
							!on se balade de l'indice 2 à l'indice nx-1, donc on est en i+1 et non en i
							nnz((i-1)*5+Nj:i*5+nj)=fg
							ic(Njic+1)=ic(Njic)+5
							colonnes((i-1)*5+Nj:i*5+nj)=[j*nx+i+1-nx,j*nx+i+1-1,j*nx+i+1,j*nx+i+1+1,j*nx+i+1+nx]
							Njic=Njic+1
						end do
						Nj=Nj+5*(nx-2)
						nnz(Nj:Nj+3)=[-1*dt*(1/(hy**2)),-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-1*dt*(1/(hy**2))]
						colonnes(Nj:Nj+3)=[j*nx+1+nx-1-nx,j*nx+1+nx-1-1,j*nx+1+nx-1,j*nx+1+nx-1+nx]
						ic(Njic+1)=ic(Njic)+4 !i=i+1, à la fin i=nx-1

						Njic=Njic+1 !décalage pour i
						Nj=Nj+4 !
						Njc=Njc+2 ! valeurs aux bords
				end do

				nnz(Nj:Nj+2)=[-1*dt*(1/(hy**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2))]
				colonnes(Nj:Nj+2)=[(ny-1)*nx+1-nx,(ny-1)*nx+1,(ny-1)*nx+1+1]
				ic(Njic+1)=ic(Njic)+3
				!k=i*nx+j
				Nj=Nj+3
				Njic=Njic+1

				do i=1,nx-2
					!on se balade de l'indice 2 à l'indice nx-1, donc on est en i+1 et non en i
					nnz((i-1)*4+Nj:i*4+nj)=[-1*dt*(1/(hy**2)),-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2))]
					ic(Njic+1)=ic(njic)+4
					colonnes((i-1)*4+Nj:i*4+nj)=[(ny-1)*nx+1+i-nx,(ny-1)*nx+1+i-1,(ny-1)*nx+1+i,(ny-1)*nx+1+i+1]
					Njic=Njic+1
				end do
				Nj=Nj+4*(nx-2)
				nnz(Nj:Nj+2)=[-1*dt*(1/(hy**2)),-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2))]
				colonnes(Nj:Nj+2)=[ny*nx-nx,ny*nx-1,ny*nx]
				ic(Njic+1)=ic(Njic)+3 !i=i+1, à la fin i=nx-1
				!ic(Njic+2)=ic(Njic+1)+1

			return
			end subroutine

			subroutine sm1(x,y,Nx,Ny,ssm,dt,u,dx,dy)  !second memebre
				implicit none
				integer::Nx,Ny,i,j,k,me,i1,iN,Np,he,statinfo
				double precision :: dt,dx,dy
				double precision, dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension((Ny+2)*(Nx+2))::u
				!call charge(me, Np, Ny,j1,jN)
				double precision, dimension(:),allocatable::ssm
				!allocate(ssm(1:(Nx+2)*(Ny+2)))
				ssm=0.d0
				do j=1,Ny+2 ! dimension Nx+2
					do i=1,Nx+2
						k=i+(j-1)*(Nx+2)
						if (i==1) then!gauche
							ssm(k)=0.d0
						else if (i==Nx+2) then!droite
							ssm(k)=0.d0
						else if ((j==1)) then!bas
							ssm(k)=0.d0
						else if (j==Ny+2) then!haut
							ssm(k)=0.d0
						else
							ssm(k)=u(k)+f1(x(i),y(j))*dt
						end if
					end do
				end do
!				open(unit=1, file="sm1.txt",form="formatted",access="sequential")
!					do i=1,Nx+2
!						do j=1,Ny+2
!							write(1,*)ssm(i+(j-1)*(Nx+2)),i,j
!						enddo
!					enddo
!					close(1)

			end subroutine sm1

			subroutine sm3(x,y,Nx,Ny,ssm,dt,u,Nt,dx,dy)  !second memebre
				implicit none
				integer::Nx,Ny,i,j,k,me,i1,iN,Np,he,statinfo,n,Nt
				double precision :: dt,Lx,Ly,dx,dy
				double precision, dimension(Nx)::x
				double precision, dimension(Ny)::y
				double precision,dimension(Nt)::t
				double precision, dimension(Ny*Nx)::u
				!call charge(me, Np, Ny,j1,jN)
				double precision, dimension(:),allocatable::ssm
				Lx=1
				Ly=1
				ssm=0
				do j=1,Ny+2
					do i=1,Nx+2
						k=i+(j-1)*(Nx+2)
						if ((j==1)) then
							ssm(k)=0.d0
						else if (j==Ny) then
							ssm(k)=0.d0
						else if (i==1) then!gouche
							ssm(k)=1.d0*dt
						else if (i==Nx) then!droite
							ssm(k)=1.d0*dt
						else
							ssm(k)=ff(x(i+1),y(j+1),t(n),Ly,Lx)*dt+u(k)
						end if
					end do
				end do

			end subroutine sm3



			subroutine sm2(x,y,Nx,Ny,ssm,dt,u,dx,dy)  !second memebre
				implicit none
				integer::Nx,Ny,i,j,k,me,i1,iN,Np,he,statinfo
				double precision :: dt,dx,dy
				double precision, dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension((Nx+2)*(Ny+2))::u
				!call charge(me, Np, Ny,j1,jN)
				double precision, dimension(:),allocatable::ssm

				!allocate(ssm((Nx+2)*(Ny+2)))
				ssm=0
				do j=1,Ny+2
					do i=1,Nx+2
						k=i+(j-1)*(Nx+2)
						if ((j==1)) then
							ssm(k)=f(x(i),y(j))*dt
						else if (j==Ny+2) then
							ssm(k)=f(x(i),y(j))*dt
						else if (i==1) then!gouche
							ssm(k)=f(x(i),y(j))*dt
						else if (i==Nx+2) then!droite
							ssm(k)=f(x(i),y(j))*dt
						else
							ssm(k)=f(x(i),y(j))*dt+u(k)
						end if
					end do
				end do

			end subroutine sm2

			subroutine charge(me,Np,N,i1,iN)
				implicit none
				integer ::i1,iN,Me,Np,N,r,q
				q=N/np
				r=N-q*np
				if (me<r) then
					i1=me*(q+1) +1
					iN=(me+1)*(q+1)
				else
					i1= 1 + r + me*q
					iN= i1 +q -1
				endif
			end subroutine charge


end module  matrices
