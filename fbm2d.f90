program fbm2d
	implicit none
	
	integer :: l
	real :: h
	integer :: seed
	character(100) :: filename
	character(100) :: chbuf
	
	real,allocatable :: p(:,:)
	
	integer :: i,j,k
	
	call getarg(1,chbuf)
	read(chbuf,*)l
	call getarg(2,chbuf)
	read(chbuf,*)h
	call getarg(3,chbuf)
	read(chbuf,*)seed
	seed=-seed
	call getarg(4,chbuf)
	read(chbuf,*)filename
	
	allocate(p(l,l))
	
	call ffm2d(l,seed,h,p)
	
	open(1,file='fbm2d_'//trim(filename)//'.dat')
	write(1,'(2i8)')l,l
	do j=1,l
		do i=1,l
			write(1,'(2i8,f16.7)')i-1,j-1,p(i,j)
		end do
	end do
	close(1)
	
end program fbm2d

subroutine ffm2d(l,seed,h,z)
	implicit none
	
	integer :: l
	integer :: seed
	real :: h
	real,dimension(l,l) :: z
	
	complex,allocatable :: a(:,:)
	real,parameter :: pi=acos(-1.)
	real :: phase
	real :: rad
	integer :: i0,j0
	
	integer,allocatable :: nn(:)
	real, allocatable :: b(:)
	
	real :: ran2, gasdev
	integer :: i,j,k
	
	allocate(a(0:l-1,0:l-1))
	allocate(nn(2))
	allocate(b(2*l*l))
	a=0
	nn=0
	b=0
	
	do i=0,l/2
		do j=0,l/2
			phase=2*pi*ran2(seed)
			if(i/=0.or.j/=0)then
				rad=(i*i+j*j)**(-(h+1)/2)*gasdev(seed)
			else
				rad=0.
			end if
			a(i,j)=cmplx(rad*cos(phase),rad*sin(phase))
			
			if(i==0)then
				i0=0
			else
				i0=l-i
			end if

			if(j==0)then
				j0=0
			else
				j0=l-j
			end if
			a(i0,j0)=cmplx(rad*cos(phase),-rad*sin(phase))
		end do
	end do
	a(l/2,0)=cmplx(real(a(l/2,0)),0.0)
	a(0,l/2)=cmplx(real(a(0,l/2)),0.0)
	a(l/2,l/2)=cmplx(real(a(l/2,l/2)),0.0)
	do i=1,l/2-1
		do j=1,l/2-1
			phase=2*pi*ran2(seed)
			rad=(i*i+j*j)**(-(h+1)/2)*gasdev(seed)
			a(i,l-j)=cmplx(rad*cos(phase),rad*sin(phase))
			a(l-i,j)=cmplx(rad*cos(phase),-rad*sin(phase))
		end do
	end do
	
	nn(1)=l
	nn(2)=l
	k=0
	do i=0,l-1
		do j=0,l-1
			k=k+1
			
			b(2*k-1)=real(a(i,j))
			b(2*k)=aimag(a(i,j))
		end do
	end do
	
	call fourn(b,nn,2,-1)
	
	k=0
	do i=1,l
		do j=1,l
			k=k+1
			z(i,j)=b(2*k-1)
		end do
	end do
	
	deallocate(a)
	deallocate(nn)
	deallocate(b)
end subroutine ffm2d

subroutine fourn(data,nn,ndim,isign)
	integer isign,ndim,nn(ndim)
	real data(*)
	integer i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3
	integer k1,k2,n,nprev,nrem,ntot
	real tempi,tempr
	real theta,wi,wpi,wpr,wr,wtemp
	
	ntot=1
	do idim=1,ndim
		ntot=ntot*nn(idim)
	enddo
	
	nprev=1
	do idim=1,ndim
		n=nn(idim)
		nrem=ntot/(n*nprev)
		ip1=2*nprev
		ip2=ip1*n
		ip3=ip2*nrem
		i2rev=1
		do i2=1,ip2,ip1
			if(i2.lt.i2rev)then
				do i1=i2,i2+ip1-2,2
					do i3=i1,ip3,ip2
						i3rev=i2rev+i3-i2
						tempr=data(i3)
						tempi=data(i3+1)
						data(i3)=data(i3rev)
						data(i3+1)=data(i3rev+1)
						data(i3rev)=tempr
						data(i3rev+1)=tempi
					enddo
				enddo
			endif
			ibit=ip2/2
1			if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
				i2rev=i2rev-ibit
				ibit=ibit/2
				goto 1
			endif
			i2rev=i2rev+ibit
		enddo
		ifp1=ip1
2		if(ifp1.lt.ip2)then
			ifp2=2*ifp1
			theta=isign*6.28318530717959d0/(ifp2/ip1)
			wpr=-2.d0*sin(0.5d0*theta)**2
			wpi=sin(theta)
			wr=1.d0
			wi=0.d0
			do i3=1,ifp1,ip1
				do i1=i3,i3+ip1-2,2
					do i2=i1,ip3,ifp2
						k1=i2
						k2=k1+ifp1
						tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
						tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
						data(k2)=data(k1)-tempr
						data(k2+1)=data(k1+1)-tempi
						data(k1)=data(k1)+tempr
						data(k1+1)=data(k1+1)+tempi
					enddo
				enddo
				wtemp=wr
				wr=wr*wpr-wi*wpi+wr
				wi=wi*wpr+wtemp*wpi+wi
			enddo
			ifp1=ifp2
			goto 2
		endif
		nprev=n*nprev
	enddo
end subroutine fourn

function gasdev(idum)
	integer idum
	real gasdev
	
	integer iset
	real fac,gset,rsq,v1,v2,ran2
	save iset,gset
	data iset/0/
	
	if (idum.lt.0)iset=0
	if (iset.eq.0)then
1		v1=2.*ran2(idum)-1.
		v2=2.*ran2(idum)-1.
		rsq=v1**2+v2**2
		if(rsq.ge.1..or.rsq.eq.0.)goto 1
		fac=sqrt(-2.*log(rsq)/rsq)
		gset=v1*fac
		gasdev=v2*fac
		iset=1
	else
		gasdev=gset
		iset=0
	endif
end function gasdev

function ran2(idum)
	integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
	real ran2,am,eps,rnmx
	parameter (im1=2147483563,im2=2147483399,am=1.d0/im1,imm1=im1-1,&
		ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,&
		ntab=32,ndiv=1+imm1/ntab,eps=3.d-16,rnmx=1.d0-eps)
	integer idum2,j,k,iv(ntab),iy
	save iv,iy,idum2
	data idum2/123456789/, iv/ntab*0/, iy/0/
	if (idum.le.0) then
		idum=max(-idum,1)
		idum2=idum
		do j=ntab+8,1,-1
			k=idum/iq1
			idum=ia1*(idum-k*iq1)-k*ir1
			if (idum.lt.0) idum=idum+im1
			if (j.le.ntab) iv(j)=idum
		end do
		iy=iv(1)
	endif
	k=idum/iq1
	idum=ia1*(idum-k*iq1)-k*ir1
	if (idum.lt.0) idum=idum+im1
	k=idum2/iq2
	idum2=ia2*(idum2-k*iq2)-k*ir2
	if (idum2.lt.0) idum2=idum2+im2
	j=1+iy/ndiv
	iy=iv(j)-idum2
	iv(j)=idum
	if(iy.lt.1)iy=iy+imm1
	ran2=min(am*iy,rnmx)
end function ran2
