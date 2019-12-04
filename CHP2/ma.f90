module matrices
		use mpi
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

							allocate (A((Nx+2)*(Ny+2),(Nx+2)*(Ny+2)))


							call laplacien(Nx,Ny,hx,hy,dt,D,A)


			! nombre d'éléments non nuls dans la matrice


							diim=0
							do i=1,(Nx+2)*(Ny+2)
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

							do i=1,(Nx+2)*(Ny+2)
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
!										open(unit=6, file="M.txt",form="formatted",access="sequential")
!											do i=1,(Nx+2)*(Ny+2)
!													write(6,*)M(i,:)
!											enddo
!											close(6)
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

			subroutine sm11(x,y,Nx,Ny,ssm,dt,u,dx,dy)  !proc 0
				implicit none
				integer::Nx,Ny,i,j,k,me,i1,iN,Np,he,statinfo,tag0,tag1,request,r
				double precision :: dt,dx,dy
				integer, dimension(Nx+2) :: haut, bas
				double precision, dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension((Ny+2)*(Nx+2))::u
				  integer,dimension(MPI_STATUS_SIZE) :: status
				!call charge(me, Np, Ny,j1,jN)
				double precision, dimension(:),allocatable::ssm
				!allocate(ssm(1:(Nx+2)*(Ny+2)))
				!call MPI_INIT(statinfo)
				call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
				call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
				ssm=0.d0
				!do j=1,Ny+2 ! dimension Nx+2
				!	do i=1,Nx+2

				do j=1,(Ny+2)
					do i=1,(Ny+2)
						k=i+(j-1)*(Nx+2)
						if (i==1) then!bas
							ssm(k)=0
						else if (i==(Nx+2)) then!haut
							ssm(k)=0
						else if ((j==1)) then!gauche
							ssm(k)=0
						else if (j==(Ny+2)) then!droite
							ssm(k)=0
						else
							ssm(k)=f1(x(i),y(j))
						end if
					end do
				end do
				!call MPI_SEND(ssm((bas(2:Nx-1))),Nx,MPI_DOUBLE_PRECISION,1,tag1,MPI_COMM_WORLD,statinfo )
				!call MPI_RECV(ssm(haut(2:Nx-1)),Nx,MPI_DOUBLE_PRECISION,1,tag0,MPI_COMM_WORLD,status,statinfo)
!				open(unit=1, file="sm1.txt",form="formatted",access="sequential")
!					do i=1,Nx+2
!									if(me==0) then	do j=1,Ny+2
!							write(1,*)ssm(i+(j-1)*(Nx+2)),i,j
!						enddo
!					enddo
!					close(1)

end subroutine sm11

			subroutine sm12(x,y,Nx,Ny,ssm,dt,u,dx,dy)  !proc1
				implicit none
				integer:: Nx,Ny,i,j,k,me,i1,iN,Np,he,statinfo,tag0,tag1,request,r
				double precision :: dt,dx,dy
				integer, dimension(Nx+2) :: haut, bas
				double precision, dimension(:),allocatable::x
				double precision, dimension(:),allocatable::y
				double precision, dimension(:),allocatable::u
				double precision, dimension(:),allocatable::ssm
				  integer,dimension(MPI_STATUS_SIZE) :: status
				!call charge(me, Np, Ny,j1,jN)
				!call MPI_INIT(statinfo)
				call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
				call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
				!allocate(ssm(1:(Nx+2)*(Ny+2)))
				ssm=0.d0
				bas=0
				haut=0
				do j=1,(Ny+2) !
					do i=1,(Nx+2)
						k=i+(j-1)*(Nx+2)!
						if (i==1) then!bas normalement ça devrait être i=(Nx+2)/2 mais on le rajoute dans l'expression de x
							ssm(k)=0.d0
						else if (i==(Nx+2)) then!haut
							ssm(k)=0.d0
						else if ((j==1)) then!gauche
							ssm(k)=0.d0
						else if (j==(Ny+2)) then!droite
							ssm(k)=0.d0
						else
							!print*,"i = ",i," j = ",j
							ssm(k)=f1(x(i-r+(Nx+2)/2),y(j-r))!on fait le recouvrement
						end if
					end do
				end do
				!print*,"haut = ",haut(2:Nx+1)
				!print*,"bas = ",bas(2:Nx+1)
				!call MPI_SEND(ssm(bas(2:Nx-1)),Nx,MPI_DOUBLE_PRECISION,0,tag0,MPI_COMM_WORLD,statinfo )
				!call MPI_RECV(ssm(haut(2:Nx-1)),Nx,MPI_DOUBLE_PRECISION,0,tag1,MPI_COMM_WORLD,status,statinfo)
!				open(unit=1, file="sm1.txt",form="formatted",access="sequential")
!					do i=1,Nx+2
!						do j=1,Ny+2
!							write(1,*)ssm(i+(j-1)*(Nx+2)),i,j
!						enddo
!					enddo
!					close(1)

			end subroutine sm12

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
