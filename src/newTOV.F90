PROGRAM newTOV
  
! Creates a 1D profile for a star with given EoS table.
! This 1D profile then can be used for e.g. SPH simulations.
! It can be switched between Newtonian stellar equations and TOV equations.
!
! The EOS driver routine of Evan O'Connor et al. is used.
! source: https://github.com/evanoconnor/EOSdriver
!
! The form of the TOV equations used with the baryon density can be found in
! https://arxiv.org/abs/1306.4034 by Kaplan et al.
! EOS tables are available on https://stellarcollapse.org/equationofstate
!
! The code is build upon better understanding the equations and process
! of building neutron star models. For sure it can be made much faster 
! and efficient. 

	USE eosmodule
	IMPLICIT NONE

	DOUBLE PRECISION,PARAMETER :: msun=1.989d33,rsun=6.96d10
	CHARACTER*24 :: EOStab
	CHARACTER*50 :: filename

	! variables for EOSDriver
	INTEGER keytemp,keyerr,TOV,points,i
	DOUBLE PRECISION xrho,xye,xtemp,xtemp2
	DOUBLE PRECISION xenr,xprs,xent,xcs2,xdedt,xmunu
	DOUBLE PRECISION xdpderho, xdpdrhoe
	DOUBLE PRECISION xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
	DOUBLE PRECISION xxa,xxh,xxn,xxp,iprs

	DOUBLE PRECISION rhoc,dmutot,h
	DOUBLE PRECISION,DIMENSION(4) :: Y
	DOUBLE PRECISION idens,imass,irho,irad,ipress,iye,itemp
	DOUBLE PRECISION rhoc_min,rhoc_max,step

!----INPUT VALUES ------------------------------------------------------

	EOStab="LS220.h5"  		! path to EOS table
	
	filename="1D_star_profile"	! output filename
	
	TOV=0	! 0 - use Newtonian stellar equations
		    ! 1 - use TOV equations (GR)

	itemp=0.5d0      ! temperature T in MeV
	xent=1.0d0       ! entropy S in kB
	keytemp=1        ! 1 - for temperature T-slice, input: rho,T,Ye
	                 ! 2 - for entropy S-slice, input: rho,S,Ye
	rhoc=0.4893269d15	 ! central density rho_c in cgs
	!rho_c=0.4685d15 for entropy slice with S=1.0d0 kbar -> M=1.4 M_sun
	h=75.d0		     ! stepsize for RK in cm		

!-----------------------------------------------------------------------    

	! initial conditions
	call readtable(EOStab)
	call betaeqYe(rhoc,itemp,iye,xprs,dmutot,xenr,xent,keytemp)

	irad=0.d0	! radius		r
	irho=rhoc	! central density	rhoc
	Y(1)=0.d0	! baryon mass 		mb 
	Y(2)=xprs	! central pressure 	Pc
	Y(3)=0.d0	! gravitational mass	mg	[only for full TOV]
	Y(4)=0.d0	! gravitational pot.	Phi [only for full TOV]

	! Solving ODE system with RK4 while P >= 0
	open(1,file=filename)
	do while((Y(2).ge.0.0d0).and.(irho.ge.1.5d3)) 
		write(1,'(7(1x,es23.16))') irad,irho,Y(1),Y(2),itemp,iye,xenr
		call RK4(irad,irho,itemp,iye,xenr,dmutot,Y,h,TOV,xent,keytemp)
		irad=irad+h
	enddo
	close(1)
	
	print '(a,e15.5)', 'central density: ',rhoc
	print '(a,f15.5)', 'Star Mass [M_sun]: ',Y(1)/msun
	print '(a,i15)', '# of iterations: ',int(irad/h+1)

END PROGRAM newTOV

!************************************************************************

SUBROUTINE betaeqYe(rhoe,temp,yee,xprs,dmutot,xenr,xent,keytemp)

	! find Ye for beta-equilibrium: dmutot=0=mu_e+mu_p-mu_n
	! for low densities finds Ye where dmutot is minimal
	
	USE eosmodule
	IMPLICIT NONE

	INTEGER keytemp,keyerr,i,itermax
	DOUBLE PRECISION xprs,rhoe,temp,dmutot
	DOUBLE PRECISION yel,yeg,yee,tol
	DOUBLE PRECISION xtemp2,xenr,xent,xcs2,xdedt,xmunu
	DOUBLE PRECISION xdpderho, xdpdrhoe
	DOUBLE PRECISION xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
	DOUBLE PRECISION xxa,xxh,xxn,xxp
	
	keyerr=0
	tol=1.0d-9
	i=0
	itermax=50
	yel=eos_yemin		
	yeg=eos_yemax
	yee=(yel+yeg)*0.5d0
	call nuc_eos_full(rhoe,temp,yee,xenr,xprs,xent,xcs2,xdedt,&
		xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
		xmuhat,keytemp,keyerr,precision)
	dmutot=xmu_p+xmu_e-xmu_n  
	do while(abs(dmutot).ge.tol.and.(i.le.itermax))
		if(dmutot.le.tol)then
			yel=yee
			yee=(yee+yeg)*0.5d0
		else	
			yeg=yee
			yee=(yee+yel)*0.5d0
		endif
		call nuc_eos_full(rhoe,temp,yee,xenr,xprs,xent,xcs2,xdedt,&
			xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,&
			xmu_p,xmuhat,keytemp,keyerr,precision)
		dmutot=xmu_p+xmu_e-xmu_n
		i=i+1
	enddo
		
end SUBROUTINE betaeqYe

SUBROUTINE findrho(press,rho,xtemp,iye,xenr,dmutot,xent,keytemp)

	! find rho for given p, T, Ye with bisection

	USE eosmodule
	IMPLICIT NONE

	INTEGER keytemp,keyerr,j,itermax
	DOUBLE PRECISION xtemp,xprs,dmutot
	DOUBLE PRECISION yel,yeg,iye,rhol,rhog,rho,press,tol
	DOUBLE PRECISION xtemp2,xenr,xent,xcs2,xdedt,xmunu
	DOUBLE PRECISION xdpderho,xdpdrhoe
	DOUBLE PRECISION xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
	DOUBLE PRECISION xxa,xxh,xxn,xxp
	keyerr=0
	j=0
	itermax=50
	rhog=rho*1.2d0
	rhol=eos_rhomin
	rho=(rhol+rhog)*0.5d0
	tol=1.0d-9
	call nuc_eos_full(rho,xtemp,iye,xenr,xprs,xent,xcs2,xdedt,&
		xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
		xmuhat,keytemp,keyerr,precision)
	do while(abs((xprs-press)/press).ge.tol.and.(j.le.itermax))     
		call betaeqYe(rho,xtemp,iye,xprs,dmutot,xenr,xent,keytemp)
		if(((xprs-press)/press).le.tol)then
			rhol=rho
			rho=(rho+rhog)*0.5d0
		else
			rhog=rho
			rho=(rho+rhol)*0.5d0 
		endif
		call nuc_eos_full(rho,xtemp,iye,xenr,xprs,xent,xcs2,xdedt,&
			xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,&
			xmu_p,xmuhat,keytemp,keyerr,precision)
		j=j+1
	enddo
	call betaeqYe(rho,xtemp,iye,xprs,dmutot,xenr,xent,keytemp)
	
end SUBROUTINE findrho

SUBROUTINE RK4(irad,irho,itemp,iye,xenr,dmutot,Y,h,TOV,xent,keytemp)

	! Iteration for solving ODE system: dY/dr = f(r,Y) 
	! with 4th order Runge-Kutta method
	!
	! Y(1) 	mass		m   
	!		baryon mass	mb  [if TOV]
	! Y(2) 	pressure 	p
	! Y(3) 	grav. mass	mg  [if TOV]
	! Y(4)	potential   Phi [if TOV]
	! 
	! variables rho, T and Ye are calculated with findrho & EOSDriver
	
	USE eosmodule
	IMPLICIT NONE
	INTEGER TOV,keytemp
	DOUBLE PRECISION, DIMENSION(4) :: Y,iY
	DOUBLE PRECISION, DIMENSION(4) :: k1,k2,k3,k4
	DOUBLE PRECISION :: irad,irho,itemp,iye,h,dmutot,xenr,xent

	iY=Y
	call stellarEq(irad,irho,itemp,iY,Y,TOV,xent,keytemp)
	k1=h*Y
	call findrho(iY(2)+k1(2)/2,irho,itemp,iye,xenr,dmutot,xent,keytemp)
	call stellarEq(irad+h/2,irho,itemp,iY+k1/2,Y,TOV,xent,keytemp)
	k2=h*Y
	call findrho(iY(2)+k2(2)/2,irho,itemp,iye,xenr,dmutot,xent,keytemp)
	call stellarEq(irad+h/2,irho,itemp,iY+k2/2,Y,TOV,xent,keytemp)
	k3=h*Y
	call findrho(iY(2)+k3(2)/2,irho,itemp,iye,xenr,dmutot,xent,keytemp)
	call stellarEq(irad+h,irho,itemp,iY+k3,Y,TOV,xent,keytemp)
	k4=h*Y
	Y=iY+(k1+2*k2+2*k3+k4)/6
	call findrho(Y(2),irho,itemp,iye,xenr,dmutot,xent,keytemp)

end SUBROUTINE RK4

SUBROUTINE stellarEq(irad,irho,itemp,iY,Y,TOV,xent,keytemp)

	! ODE system for stellar structure: dY/dr = f(r,Y)
	! with Y(1) mass	 m   
	!	        baryon mass mb  [if TOV]
	!      Y(2) pressure 	 P
	!      Y(3) grav. mass	 mg  [if TOV]
	!      Y(4) potential   Phi [if TOV]
	
	USE eosmodule
	IMPLICIT NONE
	INTEGER TOV,keytemp
	DOUBLE PRECISION,DIMENSION(4) :: iY,Y
	DOUBLE PRECISION irad,irho,c2,slant
	DOUBLE PRECISION itemp,xenr,iye,dmutot,xent

    c2=clight*clight
    
    if(TOV.eq.0)then
		if(iY(1).eq.0.0d0)then
			Y(1) = 4.d0/3*pi*irho*irad**2
			Y(2) = -2.d0/3*ggrav*irho**2*pi*irad
		else
			Y(1) = 4.d0*pi*irho*irad**2
			Y(2) = -ggrav*irho*iY(1)/irad**2
		endif
    else
		call betaeqYe(irho,itemp,iye,iY(2),dmutot,xenr,xent,keytemp)
		slant=(1.d0+xenr/c2+iY(2)/(irho*c2))
		if(iY(1).eq.0.0d0)then
			Y(1) = 4.d0*pi*irho*irad**2
			Y(2) = -ggrav*irho*slant*4.d0*pi*iY(2)*irad/c2
			Y(3) = 4.0d0*pi*irad**2*(irho-iY(2)/c2)
			Y(4) = 4.0d0*ggrav*pi*iY(2)*irad/(c2*c2)
		else
			Y(1) = 4.d0*pi*irho*irad**2*(1.0d0-2.0d0*ggrav*iY(3)/ &
				&(irad*c2))**(-0.5d0)
			Y(2) = -ggrav*irho*slant*iY(3)/irad**2*(1.0d0+4.0d0*pi*iY(2)/c2* &
				&irad**3/iY(3))*(1.0d0-2.0d0*ggrav*iY(3)/c2/irad)**(-1)
			Y(3) = 4.0d0*pi*irad**2*(irho*slant-iY(2)/c2)
			Y(4) = (ggrav/c2)*iY(3)/irad**2*(1.0d0+4.0d0*pi*iY(2)/c2* &
				&irad**3/iY(3))*(1.0d0-2.0d0*ggrav*iY(3)/c2/irad)**(-1)
		endif
	endif
end SUBROUTINE stellarEq

