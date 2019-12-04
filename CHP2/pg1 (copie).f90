program programme
	use mpi
	use fonctions
	use matrices
	use gradientconjugue
	implicit none
	integer,parameter :: Nx=30,Ny=30,NxNy=Nx*Ny
	integer::j1,jN,k1,k2,k3,j,i,i1,iN,statinfo,Np,me,cont,diim,Nt,n2,n5,n1=4,k
	double precision,parameter :: hx=1.d0/(Nx+1),hy=1.d0/(Ny+1),Lx=1.d0,Ly=1.d0,dt=0.01d0,D=1.d0
	double precision, dimension(:),allocatable :: b,xn,ssm,ssme,produitmatriciel,v,vv,nnz,colonnes,ic,xx,xn1,xn2,ssm1,ssm2
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
	!allocate (b((Nx+2)*(Ny+2)))
	!allocate (U((Nx+2)*(Ny+2)))
	!allocate (Uo((Nx+2)*(Ny+2)))
	!allocate (Uo2((Nx+2)*(Ny+2)))
if(me==0) then
	allocate(xn1((Nx+2)*(Ny+2)))
	allocate(ssm1((Nx+2)*(Ny+2)))
	do i=1,Nx+2
		x(i)=(i-1)*hx
	end do
	! axe des ordonnées
	do i=1,Ny+2
		y(i)=(i-1)*hy
	end do
else if(me==1) then
	allocate(xn2((Nx+2)*(Ny+2)))
	allocate(ssm2((Nx+2)*(Ny+2)))
	do i=1,Nx+2
		x(i)=(i-1+Nx+2)*hx
	end do
	! axe des ordonnées
	do i=1,Ny+2
		y(i)=(i-1+Nx+2)*hy
	end do
end if

	Nt=10
	! axe des abscisses

	!print*," x = ",x
	! échelle de temps
	do i=1,100
		temps(i)=(i-1)*dt
	end do
	b=0.d0
	U=0.d0
	if(me==0) then
	do j=1,Ny+2 ! dimension Nx+2
		do i=1,Nx+2
			k=i+(j-1)*(Nx+2)
			if (i==1) then!bas
				ssm1(k)=0
			else if (i==Nx+2) then!haut
				ssm1(k)=0
			else if ((j==1)) then!gauche
				ssm1(k)=0
			else if (j==Ny+2) then!droite
				ssm1(k)=0
			else
				ssm1(k)=f1(x(i),y(j))
			end if
		end do
	end do

else if(me==1) then
	do j=1,Ny+2 ! dimension Nx+2
		do i=1,Nx+2
			k=i+(j-1)*(Nx+2)!on passe au se
			if (i==1) then!bas
				ssm2(k)=0
			else if (i==Nx+2) then!haut
				ssm2(k)=0
			else if ((j==1)) then!gauche
				ssm2(k)=0
			else if (j==Ny+2) then!droite
				ssm2(k)=0
			else
				ssm2(k)=f1(x(i),y(j))
			end if
		end do
	end do
end if

	call matrice(Nx,Ny,hx,hy,nnz,colonnes,ic,D,i1,iN,dt)
	!call sm1(x,y,Nx,Ny,ssm,dt,Uo,hx,hy)


	!print*,"hey"
	do i=1,Nt
		if(me==0) then
		call gc(nnz,colonnes,ic,ssm2,(Nx+2)*(Ny+2),Nx,Ny,xn1)
		call sm11(x,y,Nx,Ny,ssm1,dt,xn1,hx,hy)
	else if(me==1) then
		call gc(nnz,colonnes,ic,ssm1,(Nx+2)*(Ny+2),Nx,Ny,xn2)
		call sm12(x,y,Nx,Ny,ssm2,dt,xn2,hx,hy)
	end if
	end do
	!print*,x

	n2=size(U)
open(unit=2, file="sol1.dat",form="formatted",access="sequential")
	do i=1,Nx+2
		do j=1,Ny+2
			write(2,*)x(i),x(j),U(i+(j-1)*(Nx+2))
		enddo
	enddo
	close(2)
	!print*,maxval(U)
!print*,"what"
	temps_fin= (MPI_WTIME()-temps_debut)
	deallocate (xn1)
	deallocate(xn2)
	deallocate (ssm1)
	deallocate (ssm2)


	call MPI_REDUCE (temps_fin,temps_fin_max,1, MPI_DOUBLE_PRECISION , MPI_MAX ,0,MPI_COMM_WORLD ,statinfo)
	if (me == 0) then
		print*,"Temps : ",temps_fin_max," secondes"
	end if
			!print*,me,U
		call MPI_FINALIZE(statinfo)


end program programme
