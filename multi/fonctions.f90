		module fonctions
			implicit none
			contains
				function f1(x,y)
					implicit none
					double precision, intent(in)::x,y
					double precision :: f1
					f1=2*(y-y**2 + x - x**2)
				end function f1

				function f(x,y)
					implicit none
					double precision, intent(in)::x,y
					double precision :: f
					f=sin(x)+cos(y)
				end function f

				function g(x,y)
					implicit none
					double precision, intent(in)::x,y
					double precision :: g
					g=0
				end function g

				function h(x,y)
					implicit none
					double precision,intent(in)::x,y
					double precision :: h
					h=0
				end function h

				function ff(x,y,t,Ly,Lx)
					implicit none
					double precision, parameter :: pi=3.14159265
					double precision, intent(in)::x,y,t,Ly,Lx
					double precision::ff
					ff=exp(-((x-Lx/2.d0)**2))*(exp(-((y-Ly/2.d0)**2)))*cos((pi/2.d0)*t)
				end function ff

				function uex2(x,y)
					implicit none
					double precision, parameter :: pi=3.14159265
					double precision, intent(in)::x,y
					double precision::ff
					double precision::uex2
					uex2=sin(x)+cos(y)
				end function uex2

				function uex1(x,y)
					implicit none
					double precision, parameter :: pi=3.14159265
					double precision, intent(in)::x,y
					double precision::ff
					double precision::uex1
					uex1=sin(x)+cos(y)
				end function uex1

				function gg(x,y,t)
					implicit none
					double precision, intent(in)::x,y,t
					double precision :: gg
					gg=0
				end function gg


				function hh(x,y,t)
					implicit none
					double precision, intent(in)::x,y,t
					double precision :: hh
					hh=1.
				end function hh

				function trad(k,Nx,i1,iN)! k est la numérotation globale et i numérotation locale
					implicit none
					integer:: trad,k,i1,Nx,Ny,iN
					trad=k
				end function


				function tradinv(i,j,Nx,Np)! k est la numérotation globale et i numérotation locale
					implicit none
					integer:: tradinv,i,j,Nx,Np
					tradinv=i+(j-1)*(Nx+2)/Np
				end function


		end module
