! This module contains the algebra routines
!
!   20 eigenvalues
!   40 determinant
!   70 inverse
!  100 printmatrix
!  120 exponential of a diagonal matrix
!
!
!
!
module matrixmod

implicit none
integer, private, parameter :: i4=selected_int_kind(9)
integer, private, parameter :: r8=selected_real_kind(15,9)
complex(kind=r8), private, parameter :: czero = (0.0_r8,0.0_r8)
complex(kind=r8), private, parameter :: cone  = (1.0_r8,0.0_r8)
complex(kind=r8), private, parameter :: ci    = (0.0_r8,1.0_r8)
!real(kind=r8),private,save,allocatable :: A_N_zeros(:,:),A_M_zeros(:,:)
integer(kind=i4),private,save :: N_s


contains

subroutine bubblesort(vec,N)
real(kind=r8)::vec(N)
real(kind=r8) temp
integer(kind=i4) :: bubble, N, j,lsup
lsup=N
!lsup is the size of the array to be used

do while (lsup > 1)
bubble = 0 !bubble in the greatest element out of order
do j = 1, (lsup-1)
if (vec(j) > vec(j+1)) then
temp = vec(j)
vec(j) = vec(j+1)
vec(j+1) = temp
bubble = j
endif
enddo
lsup = bubble
enddo

endsubroutine bubblesort

!DIAGONALIZE A SYMMETRIC REAL MATRIX
subroutine eigenrs(eigvec,eigval,n)
!
! compute eigenvalues and eigenvectors of a real symmetric matrix
! input matrix in eigvec (only lower triangle used), output eigenvectors
! in eigvec, output eigenvalues in eigval
!
integer(kind=i4) :: n,info
real(kind=r8) :: eigvec(n,n),eigval(n),work(130*n)
!
! lwork >= (nb+2)*n so assume block size 128 is more than enough
!
if (n.lt.1) return
call dsyev('v','u',n,eigvec,n,eigval,work,130*n,info)
if (info.ne.0) then
write (6,'(/,''Error in dsyev'',i10)') info
stop
endif
return
end subroutine eigenrs







subroutine eigengenright(A,wr,wi,vr,n)
integer(kind=i4) :: n,info
real(kind=r8) :: A(n,n),wr(n),wi(n),work(4*n),vl(N,N),vr(N,N),A_temp(N,N)
A_temp=A
call DGEEV ('N','V',n,A_temp,n,wr,wi,vl,n,vr,n,WORK,4*n,info)
if(info.ne.0) then
    write(*,*) "prob"
endif
end subroutine eigengenright

subroutine eigengenleft(A,wr,wi,vl,n)
integer(kind=i4) :: n,info
real(kind=r8) :: A(n,n),wr(n),wi(n),work(4*n),vl(N,N),vr(N,N),A_temp(N,N)
A_temp=A
call DGEEV ('V','N',n,A_temp,n,wr,wi,vl,n,vr,n,WORK,4*n,info)
if(info.ne.0) then
write(*,*) "prob",info
endif
endsubroutine eigengenleft





subroutine inv(a,n)
!
! This calculate the inverse of a matrix can be optimized
!
integer(kind=i4) :: n,ipiv(n),info,i
real (kind=r8) :: a(n,n),cwork(n,n)
!
! lapack routine for lu factorization
!
call dgetrf(n,n,a,n,ipiv,info)
if (info.ne.0) then
write (6,'(/,''Error in dgetrf'',i10)') info
stop
endif
!
! lapack routine to calculate inverse from factorization
!
call dgetri(n,a,n,ipiv,cwork,n*n,info)
if (info.ne.0) then
write (6,'(/,''Error in dgetri'',i10)') info
stop
endif

return
end subroutine inv

subroutine mprint(A,N)
!
! This routine print a square matrix
!
real(kind=r8)::A(N,N)
integer(kind=i4)::i,j,N

do i = 1,N
write (*,*) (A(i,j), j=1,N)
enddo

end subroutine mprint



!fermi_distribution of diagonal matrix
subroutine fermi_diagonal(F_E,E_mat,beta,mu,N)
real(kind=r8) :: F_E(N,N),E_mat(N,N)
real(kind=r8) :: beta,mu
integer(kind=i4)::N,i
F_E=E_Mat
do i=1,N
F_E(i,i)=1./(exp(beta * (E_mat(i,i)-mu)     )+1)
enddo

endsubroutine fermi_diagonal

subroutine solve_linear_system(M,inhom,solu,N)
integer(kind=i4)::N,ipiv(N),info
real(kind=r8)::M(N,N),inhom(N),solu(N)
solu=inhom
call dgesv(N,1,M,N,ipiv,solu,N,info)
write(*,*) "info", info
endsubroutine solve_linear_system



subroutine deter(a,det,n)
!
!compute the determinat of a matrix via LU decomposition
!det(L*U)=det(L)det(U)  product of the element on the diagonal LoL
!
integer(kind=i4) :: n,ipiv(n),info,i
real (kind=r8) :: a(n,n),det,cwork(n,n)

!
! lapack routine for L.U factorization
!
call DGETRF(n,n,a,n,ipiv,info)
if (info.ne.0) then
write (6,'(/,''Error in zgetrf'',i10)') info
stop
endif
!
! calculate determinant
!
det=1
do i=1,n
det=det*a(i,i)
if (ipiv(i).ne.i) det=-det
enddo
return
end subroutine deter


subroutine gen_mprint(A,N_r,N_col)
!
! This routine print a square matrix
!
real(kind=r8)::A(N_r,N_col)
integer(kind=i4)::i,j,N_r,N_col

do i = 1,N_r
write (*,*) (A(i,j), j=1,N_col)
enddo
end subroutine gen_mprint

subroutine generate_Q(Q,Kc,E,beta,mu,N)

real(kind=r8)::Q(N,N),Kc(N,N),beta,mu
real(kind=r8)::F_E(N),E(N)  
integer(kind=i4)::N,i,j

if(beta.gt.50) then
	do i=1,N
		if(E(i).gt.mu) then
			F_E(i)=0 
		else
			F_E(i)=1
		endif 
	enddo
else
	do i=1,N
		F_E(i)=1./( exp(beta*(E(i)-mu) ) + 1. )
	enddo 
endif 

do i=1,N
	do j=1,N
		if(Kc(i,j).eq.0) then
			Q(i,j)=0
		else 
			Q(i,j)=( ( F_E(i)-F_E(j) )/(E(i)-E(j)) )*Kc(i,j)
			if(abs(E(i)-E(j)).lt.1E-7.and.Kc(i,j).gt.1E-7) then
				write(*,*) "strange",i,j,E(i),E(j),Kc(i,j),Q(i,j)   
			endif
		endif       
	enddo
enddo 
endsubroutine generate_Q

subroutine generate_FR(FR,R,E,beta,mu,N) 
real(kind=r8)::FR(N,N),R(N,N),E(N),F_E(N) 
real(kind=r8)::beta,mu
integer (kind=i4)::N,i,j
if(beta.gt.50) then
	do i=1,N
		if(E(i).gt.mu) then
			F_E(i)=0
		else
			F_E(i)=1
		endif
	enddo
else
	do i=1,N
		F_E(i)=1./( exp(beta*(E(i)-mu) ) + 1. ) 
	enddo
endif


do i=1,N
	do j=1,N
		FR(i,j)=F_E(i)*R(i,j)  
	enddo
enddo 
endsubroutine generate_FR











end module
