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
!
!I need to build the K_corr in the psi_basis blocking the stuff
!this should be r
!
subroutine generate_Q(E,K_corr,Q,beta,mu,N)
integer(kind=i4)::N,i,j,s
real(kind=r8)::E(N),K_corr(N,N),Q(N,N)
real(kind=r8)::delta_E,pdf_is,delta_sign,beta,mu

do i=1,N
    do j=1,N
        if(abs(K_corr(i,j)).lt.1.E-13) then
            K_corr(i,j)=0
        endif
    enddo
enddo

do i=1,N
    do s=1,N
        if(K_corr(i,s).eq.0)  then
            !simple cases
            Q(i,s)=0
        else
            !double pole cases or similar one
            delta_E=abs( E(i)-E(s) )
            delta_sign=E(i)-E(s)
            if(delta_E.lt.1.E-6) then
                write(*,*) "r4ok"
                if(beta.gt.600) then
                    if(mu.ne.E(i).and.mu.ne.E(s)) then
                        pdf_is=0
                    else
                    write(*,*) "whowowo"
                    pdf_is= -beta/( (1+exp(-beta*(E(i)-mu)))**2 )
                    endif
                else
                !if the temperature is high then we can use the formal derivative of the
                !fermi distributin maybe this expression can be regularized
                    pdf_is=-beta* (exp(beta*(E(i)-mu))) / (exp(beta*(E(i)-mu))+1)**2
                endif
                Q(i,s)=pdf_is*K_corr(i,s)
            else
            !no double pole case maybe we should do somenthing also here
                if(beta.gt.600)
                    if(E(i).lt.mu.and.E(s).lt.mu) then
                        pdf_is=0
                    endif
                    if(E(i).gt.mu.and.E(s).gt.mu)then
                        pdf_is=0
                    endif
                    if(E(i).gt.mu.and.E(s).lt.mu) then
                        pdf_is=-1
                    endif
                    if(E(i).lt.mu.and.E(s).gt.mu) then
                        pdf_id=1
                    endif
                else
                    pdf_is=( 1./(exp(beta*(E(i)-mu))+1) - 1./(exp(beta*E(j)-mu)+1) )/(E(i)-E(s))
                    Q(i,s)= pdf_is*K_corr(i,s)
                endif
                    Q(i,s)=pdf_is*K_corr    
            endif
        endif

    enddo
enddo

endsubroutine generate_Q
!
!
!
subroutine generate_MP1(MP1,R,N)
integer(kind=i4)::N
real(kind=r8)::MP1(N,N),R(N,N)
MP1=transpose(R)
endsubroutine generate_MP1

subroutine generate_MP2(MP2,E,M,R,beta,mu,N)
integer(kind=i4)::i,j,N,q
real(kind=r8):: f_e,beta,mu
real(kind=r8)::MP2(N,N),M(N,N),R(N,N),E(N)
do i=1,N
    do j=1,N
        MP2(i,j)=0
    enddo
enddo

do i=1,N
    f_e=1./(exp(beta*(E(i)-mu)  )+1)
    do j=1,N
        do q=1,N
            MP2(i,j)=MP2(i,j)+R(q,i)*f_e*M(q,j)
        enddo
    enddo
enddo

endsubroutine generate_MP2


subroutine generate_MP3(MP3,Q,R,M,N)
integer(kind=i4)::i,j,p,s,N
real(kind=r8)::MP3(N,N),Q(N,N),R(N,N),M(N,N)

do i=1,N
    do j=1,N
        MP3(i,j)=0
    enddo
enddo



do i=1,N
    do j=1,N

        do s=1,N
            do p=1,N
                MP3(i,j)=MP3(i,j)+Q(i,s)*R(p,s)*M(p,j)
            enddo
        enddo


    enddo
enddo

endsubroutine generate_MP3

subroutine generate_v1(v1,kn_term,E,R,beta,mu,N)
real(kind=r8)::kn_term(N),E(N),f_e,v1(N),R(N,N)
real(kind=r8)::mu,beta
integer(kind=i4)::i,q,N

do i=1,N
    v1(i)=0
enddo


do i=1,N
    f_e=1./(exp(beta*(E(i)-mu)  )+1)
    do q=1,N
        v1(i)=v1(i)+f_e*R(q,i)*kn_term(q)
    enddo
enddo



endsubroutine generate_v1

subroutine generate_v2(v2,kn_term,Q,R,beta,mu,N)
real(kind=r8)::Q(N,N),R(N,N),beta,mu,kn_term(N),v2(N)
integer(kind=i4)::i,j,p,s,N
do i=1,N
    v2(i)=0
enddo

do i=1,N
    do p=1,N
        do s=1,N
            v2(i)=v2(i)+Q(i,s)*R(p,s)*kn_term(p)
        enddo
    enddo
enddo

endsubroutine generate_v2



end module
