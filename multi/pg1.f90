program programme
	use mpi
	use fonctions
	use matrices
	use gradientconjugue
	implicit none
		character*13 ::name
	integer,parameter :: Nx=50,Ny=50,NxNy=Nx*Ny
	integer::j1,jN,k1,k2,k3,j,i,i1,iN,statinfo,Np,me,cont,diim,Nt,n2,n5,n1=4,k,r,Nx2,Ny2
	double precision,parameter :: hx=1.d0/(Nx+1),hy=1.d0/(Ny+1),Lx=1.d0,Ly=1.d0,dt=0.01d0,D=1.d0
	double precision, dimension(:),allocatable :: b,xn,ssm,ssme,produitmatriciel,v,vv,nnz,colonnes,ic,xx,xn1,xn2,ssm1,ssm2
	double precision, dimension(:),allocatable :: solutionexacte,smf,U,Uo,Uinitial,xk
	double precision, dimension((Ny+2),(Nx+2))::fg,hg,gg1,Ufir
	integer, dimension( MPI_STATUS_SIZE) :: status
	double precision,dimension(:,:),allocatable::A,AA
	double precision:: produitscalaire,temps_debut,temps_fin,temps_fin_max
	double precision,dimension(:),allocatable::x
	double precision, dimension(:),allocatable::y
	double precision, dimension(100)::temps


	call MPI_INIT(statinfo)

	temps_debut= MPI_WTIME()

	call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
	call charge(me,Np,Nx+2,i1,iN)! i1
	!print*,"i1 = ",i1,"i2 = ",iN
	!allocate (b((Nx+2)*(Ny+2)))
	!allocate (U((Nx+2)*(Ny+2)))
	!allocate (Uo((Nx+2)*(Ny+2)))
	!allocate (Uo2((Nx+2)*(Ny+2)))
	allocate(y(Ny/Np+2))
	allocate(x(Nx/Np+2))
	allocate(xn((Nx/Np+2)*(Ny/Np+2)))
	allocate(ssm((Nx/Np+2)*(Ny/Np+2)))
	Nx2=Nx/Np;
	Ny2=Ny/Np;
	Nt=1
	r=0
	x=0.d0
	y=0.d0
	ssm=0.d0
	xn=0.d0
	do i=1,Nx2+2!fois le numéro du proc
		x(i)=(i-1+me*(Nx2+2))*hx
	end do
	do i=1,Ny2+2
		y(i)=(i-1)*hy
	end do

	! axe des abscisses

	!print*," x = ",x
	! échelle de temps
	do i=1,100
		temps(i)=(i-1)*dt
	end do
	ssm=0.d0
	xn=0.d0
	do j=1,Ny2+2 !
		do i=1,Nx2+2
			k=tradinv(i,j,Nx,Np)!
			if (i==1) then!bas normalement ça devrait être i=(Nx+2)/2 mais on le rajoute dans l'expression de x
				ssm(k)=0.d0
			else if (i==(Nx2+2)) then!haut
				ssm(k)=0.d0
			else if ((j==1)) then!gauche
				ssm(k)=0.d0
			else if (j==(Ny2+2)) then!droite
				ssm(k)=0.d0
			else
				!print*,"i = ",i," j = ",j
				ssm(k)=f1(x(i-r),y(j))!on fait le recouvrement
			end if
		end do
	end do



	call matrice(Nx2,Ny2,hx,hy,nnz,colonnes,ic,D,i1,iN,dt)
	!call sm1(x,y,Nx,Ny,ssm,dt,Uo,hx,hy)


	!print*,"hey"
	do i=1,Nt
		call gc(nnz,colonnes,ic,ssm,(Nx2+2)*(Ny2+2),Nx2,Ny2,xn)! comme il n'y pas pas deux point extérieur mais un seul on inclu dans le Nx le point extérieur
		call sm1(x,y,Nx2,Ny2,ssm,dt,xn,hx,hy)
		!xn=ssm
		!call sm11(x,y,Nx,Ny,ssm,dt,xn,hx,hy)
	end do
	!print*,x
	!print*,maxval(U)
!print*,"what"

!xn=ssm
	call Rename(Me,name)
	open(unit=2, file=name,access="sequential",form="formatted")
		do i=1,(Nx2+2)
			do j=1,(Ny2+2)
				write(2,*)x(i),y(j),xn(i+(j-1)*(Nx+2)/Np)
			enddo
		enddo
		close(2)


	temps_fin= (MPI_WTIME()-temps_debut)
	deallocate (x)
	deallocate (y)
	deallocate (xn)
	deallocate (ssm)
	deallocate(ic)
	deallocate(nnz)
	deallocate(colonnes)


	call MPI_REDUCE (temps_fin,temps_fin_max,1, MPI_DOUBLE_PRECISION , MPI_MAX ,0,MPI_COMM_WORLD ,statinfo)
	if (me == 0) then
		print*,"Temps : ",temps_fin_max," secondes"
	end if
			!print*,me,U
		call MPI_FINALIZE(statinfo)

contains

	subroutine Rename(Me,name)
		implicit none
		integer :: Me
		character*13 ::name
		character*3 :: tn
		integer :: i1,i2,i3
		i1 = Me/100
		i2 =( Me - 100*i1)/10
		i3 = Me - 100*i1 -10*i2
		tn = char(i1+48)//char(i2+48)//char(i3+48)
		name='sol'//tn//'.dat'
	end subroutine Rename
end program programme
