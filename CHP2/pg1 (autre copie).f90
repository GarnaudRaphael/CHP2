program programme
	use mpi
	use fonctions
	use matrices
	use gradientconjugue
	implicit none
		character*13 ::name
	integer,parameter :: Nx=40,Ny=40,NxNy=Nx*Ny
	integer::j1,jN,k1,k2,k3,j,i,i1,iN,statinfo,Np,me,cont,diim,Nt,n2,n5,n1=4,k,r
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
	call charge(me,Np,2*(Nx+2)*(Ny+2),i1,iN)
	!print*,"i1 = ",i1,"i2 = ",iN
	!allocate (b((Nx+2)*(Ny+2)))
	!allocate (U((Nx+2)*(Ny+2)))
	!allocate (Uo((Nx+2)*(Ny+2)))
	!allocate (Uo2((Nx+2)*(Ny+2)))
	allocate(xn((Nx+2)*(Ny+2)))
	allocate(ssm((Nx+2)*(Ny+2)))
if(me==0) then
	do i=1,Nx+2
		x(i)=(i-1)*hx
	end do
	! axe des ordonnées
	do i=1,Ny+2
		y(i)=(i-1)*hy
	end do
else if(me==1) then
	do i=(Nx+2)/2,Nx+2
		x(i)=(i-1+Nx+2)*hx
	end do
	! axe des ordonnées
	do i=(Ny+2)/2,Ny+2
		y(i)=(i-1+Nx+2)*hy
	end do
end if

	Nt=1
	r=0
	! axe des abscisses

	!print*," x = ",x
	! échelle de temps
	do i=1,100
		temps(i)=(i-1)*dt
	end do
	ssm=0.d0
	xn=0.d0
	if(me==0) then
	do j=1,Ny+2 ! dimension Nx+2
		do i=1,Nx+2
			k=i+(j-1)*(Nx+2)
			if (i==1) then!bas
				ssm(k)=0
			else if (i==Nx+2) then!haut
				ssm(k)=0
			else if ((j==1)) then!gauche
				ssm(k)=0
			else if (j==Ny+2) then!droite
				ssm(k)=0
			else
				ssm(k)=f1(x(i),y(j))
			end if
		end do
	end do

else if(me==1) then
	do j=(Ny+2)/2,Ny+2 ! dimension Nx+2
		do i=1,(Nx+2)/2,Nx+2
			k=i+(j-1)*(Nx+2)!on passe au se
			if (i==1) then!bas
				ssm(k)=0
			else if (i==Nx+2) then!haut
				ssm(k)=0
			else if ((j==1)) then!gauche
				ssm(k)=0
			else if (j==Ny+2) then!droite
				ssm(k)=0
			else
				ssm(k)=f1(x(i),y(j))
			end if
		end do
	end do
end if

	call matrice(Nx/2,Ny/2,hx,hy,nnz,colonnes,ic,D,i1,iN,dt)
	!call sm1(x,y,Nx,Ny,ssm,dt,Uo,hx,hy)


	!print*,"hey"
	do i=1,Nt
		if(me==0) then
		call gc(nnz,colonnes,ic,ssm,(Nx+2)*(Ny+2),Nx/2,Ny/2,xn)
		call sm11(x,y,Nx,Ny,ssm,dt,xn,hx,hy)
	else if(me==1) then
		call gc(nnz,colonnes,ic,ssm,(Nx+2)*(Ny+2),Nx/2,Ny/2,xn)
		call sm12(x,y,Nx/2,Ny/2,ssm,dt,xn,hx,hy)
	end if
	end do
	!print*,x

	!print*,maxval(U)
!print*,"what"
call Rename(Me,name)
open(unit=2, file="sol.dat",Access = 'append',form="formatted")
	do i=1,Nx+2
		do j=1,Ny+2
			write(2,*)x(i),x(j),xn(i+(j-1)*(Nx+2))
		enddo
	enddo
	close(2)
	temps_fin= (MPI_WTIME()-temps_debut)
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
