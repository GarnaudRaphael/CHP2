program programme
	use mpi
	use fonctions
	use matrices
	use gradientconjugue
	implicit none
	integer,parameter :: Nx=100,Ny=100,NxNy=Nx*Ny
	integer::j1,jN,k1,k2,k3,j,i,i1,iN,statinfo,Np,me,cont,diim,Nt,n2,n5,n1=4,k
	double precision,parameter :: hx=1.d0/(Nx+1),hy=1.d0/(Ny+1),Lx=1.d0,Ly=1.d0,dt=0.01d0,D=1.d0
	double precision, dimension(:),allocatable :: b,xn,ssm,ssme,produitmatriciel,v,vv,nnz,colonnes,ic,xx
	double precision, dimension(:),allocatable :: solutionexacte,smf,U,Uo,Uinitial,xk
	double precision, dimension((Ny+2),(Nx+2))::fg,hg,gg1,Ufir
	integer, dimension( MPI_STATUS_SIZE) :: status
	double precision,dimension(:,:),allocatable::A,AA
	double precision:: produitscalaire,temps_debut,temps_fin,temps_fin_max
	double precision,dimension(Nx+2)::x
	double precision, dimension(Ny+2)::y
	double precision, dimension(100)::temps

	x=0.d0
	y=0.d0

	call MPI_INIT(statinfo)

	temps_debut= MPI_WTIME()

	call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
	call charge(me,Np,(Nx+2)*(Ny+2),i1,iN)
	!print*,i1,iN
	allocate (b(1:(Nx+2)*(Ny+2)))
	allocate (U(1:(Nx+2)*(Ny+2)))
	allocate (Uo(1:(Nx+2)*(Ny+2)))
	allocate(xn((Nx+2)*(Ny+2)))
	allocate(ssm((Nx+2)*(Ny+2)))
	Nt=5
	! axe des abscisses
	do i=1,Nx+2
		x(i)=(i-1)*hx
	end do
	! axe des ordonnées
	do i=1,Ny+2
		y(i)=(i-1)*hy
	end do
	!print*," x = ",x
	! échelle de temps
	do i=1,100
		temps(i)=(i-1)*dt
	end do
	b=0.d0
	U=0.d0
	do j=1,Ny+2 ! dimension Nx+2
		do i=1,Nx+2
			k=i+(j-1)*(Nx+2)
			if (i==1) then!gauche
				Uo(k)=1.d0
			else if (i==Nx+2) then!droite
				Uo(k)=1.d0
			else if ((j==1)) then!bas
				Uo(k)=0.d0
			else if (j==Ny+2) then!haut
				Uo(k)=0.d0
			else
				Uo(k)=ff(x(i),y(j),temps(1),Ly,Lx)
			end if
		end do
	end do

	!print*,Uo
	! appel du vecteur initial U contenant au bord les fonctions g et h
	call matrice(Nx,Ny,hx,hy,nnz,colonnes,ic,D,i1,iN,dt)
	!call sm1(x,y,Nx,Ny,ssm,dt,Uo,hx,hy)
	b=Uo

	!print*,"hey"
	do i=1,Nt

		call gc(nnz,colonnes,ic,b,(Nx+2)*(Ny+2),Nx,Ny,xn)

		U=xn
		!deallocate(ssm)
		call sm3(x,y,Nx,Ny,ssm,dt,U,Nt,hx,hy)
		!call sm3(x,y,Nx,Ny,ssm,dt,u,n,Nt,dx,dy)
		b=ssm
	end do
	!print*,x

	n2=size(U)
open(unit=2, file="sol3.dat",form="formatted",access="sequential")
	do i=1,Nx+2
		do j=1,Ny+2
			write(2,*)x(i),x(j),U(i+(j-1)*(Nx+2))
		enddo
	enddo
	close(2)
!print*,"what"
	temps_fin= (MPI_WTIME()-temps_debut)
	deallocate (b)
	deallocate (U)
	deallocate (Uo)
	deallocate(xn)
	deallocate(ssm)

	call MPI_REDUCE (temps_fin,temps_fin_max,1, MPI_DOUBLE_PRECISION , MPI_MAX ,0,MPI_COMM_WORLD ,statinfo)
	if (me == 0) then
		print*,"Temps : ",temps_fin_max," secondes"
	end if
			!print*,me,U
		call MPI_FINALIZE(statinfo)


end program programme
