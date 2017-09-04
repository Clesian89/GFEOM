
! Two operator c_k (d'dc)_k Correlation function functional thoery.
!
! Created by Francesco Catalano on 12/09/16.
!







program FC_CFFT
use matrixmod
!variables type
!---------------------------------------------------------------------------------
implicit none
integer, parameter :: i1=selected_int_kind(1)
integer, parameter :: i4=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(15)
integer, parameter :: r8=selected_real_kind(15,9)
!---------------------------------------------------------------------------------

!Physical System variables
!---------------------------------------------------------------------------------
integer(kind=i4)::info,N_site_x,N_site_y,N_tot,M_site,cor
real(kind=r8):: U,beta,t_x,t_y,mu,t_kb,pi,tot_weight
real(kind=r8):: bar_n_old,e_tot,muinit,E_init,E_old,temp,det
!---------------------------------------------------------------------------------


!Matrix Variables
!--------------------------------------------------------------------------------
real(kind=r8),allocatable::Q(:,:),kn_term(:)
real(kind=r8),allocatable::K_mat(:,:),K12(:,:),K21(:,:),K11(:,:),K22(:,:)
real(kind=r8),allocatable::R1(:,:),R2(:,:),omega1(:,:),omega1old(:,:)
real(kind=r8),allocatable::omega8(:,:),exppomega8(:,:),expmomega8(:,:),ID8(:,:)
real(kind=r8),allocatable::wr(:),wi(:),vl(:,:),wrord(:),invR2(:,:),K_teff(:,:)
real(kind=r8),allocatable::KeffP(:,:),L_mat(:,:),R_mat(:,:),E_mat(:,:),M_mat(:,:)
real(kind=r8),allocatable::F_E(:,:),ID20(:,:),RFL(:,:),inhom(:),L_syst(:,:),solu(:)
real(kind=r8),allocatable::anticom_av(:),spectral(:),reig_test(:),ieig_test(:),rot_C_mat(:,:)
real(kind=r8),allocatable::psi_K_corr(:,:),E_n(:),M1(:,:),M2(:,:)
real(kind=r8),allocatable::mat_A(:,:),mat_B(:,:),mat_C(:,:),omega_A(:,:),omega_B(:,:),temp_BB(:,:)
real(kind=r8),allocatable::wr_A(:),wr_B(:),wi_A(:),wi_B(:),temp_AA(:,:),omegatrA(:,:),omegatrB(:,:)
real(kind=r8),allocatable::MP1(:,:),MP2(:,:),MP3(:,:),v1(:),v2(:),M_tilde(:,:),v_tilde(:)
real(kind=r8),allocatable::doub_pol_corr(:)
!----------------------------------------------------------------------------------

!Dummy_variables or cycle varibles
!-----------------------------------------------------------------------------------
integer(kind=i4):: i,j,k,l,k_i,iter,i_check,k_check,q_i,s,p
real(kind=r8):: Es

!File variables
!---------------------------------------------------------------------------------
character(len=100)::fn,fn2,fn3,fn4,fn5,fn6,fn7
integer(kind=i4)::dy,dy2,dy3,dy4,dy5,dy6,dy7
logical:: useit
!---------------------------------------------------------------------------------
!
! Set the matrix The physical constant H=t_x*Hop+U*Hint
!
!
!---------------------------------------------------------------------------------
t_x=-1
U=10
mu=7
beta=1000

!---------------------------------------------------------
!
!FILE STUFF
!
!---------------------------------------------------------

dy=8
dy2=9
dy3=10
dy4=11
dy5=12
dy6=13

write(fn,'(a)') '1_eig_ite.dat'
write(fn2,'(a)') '2_eig_ite.dat'
write(fn3,'(a)')'3_eig_ite.dat'
write(fn4,'(a)')'4_eig_ite.dat'
write(fn5,'(a)')'excitation_ite.dat'
write(fn6,'(a)')'lefteigenvectors.dat'
!write(fn7,'(a)')'comparisonofeigenvaluesEQKQ'


!open(dy,file=fn,action="write",status="new")
!open(dy2,file=fn2,action="write",status="new")
!open(dy3,file=fn3,action="write",status="new")
!open(dy4,file=fn4,action="write",status="new")
!open(dy5,file=fn5,action="write",status="new")
!open(dy6,file=fn6,action="write",status="new")

!---------------------------------------------------------
!
!ALLOCATION
!
!---------------------------------------------------------
allocate(K_mat(20,20))
allocate(ID20(20,20))
allocate(Q(20,20))
allocate(omega8(8,8))
allocate(expmomega8(8,8))
allocate(exppomega8(8,8))
allocate(ID8(8,8))
allocate(K_teff(8,8))


allocate(K12(4,4))
allocate(K21(4,4))


allocate(R2(4,4))
allocate(KeffP(4,4))
allocate(invR2(4,4))
allocate(omega1old(4,4))
allocate(omega1(4,4))

!TEST
allocate(R1(2,2))
allocate(K11(2,2))
allocate(reig_test(2))
allocate(ieig_test(2))


allocate(wr(20))
allocate(wi(20))
allocate(kn_term(20))
allocate(inhom(20))
allocate(solu(20))
allocate(vl(20,20))
allocate(wrord(20))
allocate(L_mat(20,20))
allocate(R_mat(20,20))
allocate(E_mat(20,20))
allocate(M_mat(20,20))
allocate(F_E(20,20))
allocate(RFL(20,20))
allocate(L_syst(20,20))
allocate(anticom_av(20))
allocate(spectral(20))

allocate(K22(20,20))

!ALLOCATION 
allocate(mat_A(12,12))
allocate(mat_B(8,8))
allocate(mat_C(12,8))
allocate(omega_A(12,12))
allocate(omega_B(8,8))
allocate(wr_A(12))
allocate(wr_B(8))
allocate(wi_A(12))
allocate(wi_B(8))
allocate(temp_AA(12,12))
allocate(omegatrA(12,12))
allocate(omegatrB(8,8))
allocate(temp_BB(8,8))
allocate(rot_C_mat(12,8))
allocate(psi_K_corr(20,20))
allocate(M_tilde(20,20))
!THE LINEAR SYSTEM
allocate(E_n(20))
allocate(MP1(20,20))
allocate(MP2(20,20))
allocate(MP3(20,20))
allocate(v1(20))
allocate(v2(20))

allocate(v_tilde(20))
!---------------------------------------------------------
!
!MATRIX  DEFINITION
!
!---------------------------------------------------------

do i=1,20
        kn_term(i)=0
        spectral(i)=0
    do j=1,20
        K_mat(i,j)=0
        M_mat(i,j)=0
        F_E(i,j)=0
        ID20(i,j)=0
        psi_K_corr(i,j)=0
    enddo
enddo
!can be optimized
do i=1,20
    ID20(i,i)=1
enddo
!DEFINITION OF K
!---------------------------
K_mat(1,7)=t_x
K_mat(1,2)=U
K_mat(2,4)=t_x
K_mat(2,8)=t_x
K_mat(2,3)=-t_x
K_mat(2,2)=U
K_mat(3,2)=-t_x
K_mat(3,5)=t_x
K_mat(3,9)=t_x
K_mat(3,3)=U
K_mat(3,13)=U
K_mat(4,2)=t_x
K_mat(4,5)=-t_x
K_mat(4,10)=t_x
K_mat(4,14)=U
K_mat(5,4)=-t_x
K_mat(5,11)=t_x
K_mat(5,3)=t_x
K_mat(5,6)=-U
K_mat(6,12)=t_x
K_mat(6,15)=-U
K_mat(7,1)=t_x
K_mat(7,11)=U
K_mat(8,2)=t_x
K_mat(8,9)=-t_x
K_mat(8,10)=t_x
K_mat(8,12)=-U
K_mat(9,8)=-t_x
K_mat(9,11)=t_x
K_mat(9,3)=t_x
K_mat(9,15)=U
K_mat(10,4)=t_x
K_mat(10,8)=t_x
K_mat(10,11)=-t_x
K_mat(10,10)=U
K_mat(10,16)=U
K_mat(11,10)=-t_x
K_mat(11,9)=t_x
K_mat(11,5)=t_x
K_mat(11,11)=U
K_mat(12,6)=t_x
K_mat(12,12)=U
K_mat(13,15)=-t_x
K_mat(13,17)=t_x
K_mat(14,16)=-t_x
K_mat(14,17)=t_x
K_mat(14,14)=U
K_mat(15,18)=t_x
K_mat(15,13)=-t_x
K_mat(15,15)=U
K_mat(16,14)=-t_x
K_mat(16,18)=t_x
K_mat(17,18)=-t_x
K_mat(17,14)=2*t_x
K_mat(17,13)=2*t_x
K_mat(17,19)=U
K_mat(18,16)=2*t_x
K_mat(18,17)=-t_x
K_mat(18,15)=2*t_x
K_mat(18,20)=U
K_mat(19,20)=t_x
K_mat(19,18)=-t_x
K_mat(19,13)=t_x
K_mat(19,14)=t_x
K_mat(19,19)=U
K_mat(20,19)=t_x
K_mat(20,16)=t_x
K_mat(20,15)=t_x
K_mat(20,17)=-t_x
K_mat(20,20)=U
!---------------------------

! THE UPPER MATRIX
do i=1,12
    do j=1,12
        mat_A(i,j)=K_mat(i,j)
    enddo
enddo


!THE DOWN MATRIX
do i=1,8
    do j=1,8
        mat_B(i,j)=K_mat(i+12,j+12)
    enddo
enddo

!THE CORRELATION MATRIX
do i=1,12
    do j=1,8
        mat_C(i,j)=K_mat(i,j+12)
    enddo
enddo


call eigengenleft(mat_B,wr_B,wi_B,omega_B,8)
call eigengenleft(mat_A,wr_A,wi_A,omega_A,12)

write(*,*) "The eigenvalues of A"
do i=1,12
    write(*,*) i,wr_A(i),wi_A(i)
enddo
write(*,*) " "
write(*,*) "The eigenvalues of B"
do i=1,8
    write(*,*) i,wr_B(i),wi_B(i)
enddo
write(*,*) " "
write(*,*) " "
!do i=1,12
!    do j=1,12
!        omegatrA(i,j)=omega_A(j,i)
!    enddo
!enddo
omegatrB=transpose(omega_B)
omegatrA=transpose(omega_A)

call inv(omegatrA,12)
call inv(omegatrB,8)

temp_BB=matmul(transpose(omega_B),matmul(mat_B,omegatrB))
temp_AA=matmul(transpose(omega_A),matmul(mat_A,omegatrA))
write(*,*) "Check Simmetry transform B "
write(*,*) " "
write(*,*) " "
call mprint(temp_BB,8)
write(*,*) " "
write(*,*) " "
write(*,*) "Check Simmetry transform A "
write(*,*) " "
write(*,*) " "
call mprint(temp_AA,12)
write(*,*) "THE CORRELATED BLOCK"
write(*,*) " "
write(*,*) " "
call gen_mprint(mat_C,12,8)
write(*,*) " "
write(*,*) " "
write(*,*) "THE ROTATED CORRELATED BLOCK"
write(*,*) " "
write(*,*) " "
rot_C_mat=matmul(transpose(omega_A),matmul(mat_C,omegatrB))
call gen_mprint(rot_C_mat,12,8)

write(*,*) "PROBLEMATIC EIGENVECTOR OF A"

do i=1,12
write(*,*) i, omega_A(i,3),omega_A(i,6),omega_A(i,7),omega_A(i,8)
enddo
write(*,*) " "
write(*,*) " "
write(*,*) "PROBLEMATIC EIGENVECTOR OF B"
do i=1,8
write(*,*) i+12, omega_B(i,1),omega_B(i,3),omega_B(i,5),omega_B(i,6)
enddo
write(*,*) " "
write(*,*) " "
write(*,*) "NON PROBLEMATIC EIGENVECTOR OF A"
do i=1,12
write(*,*) i, omega_A(i,2),omega_A(i,1),omega_A(i,4),omega_A(i,5)
enddo
write(*,*) "NON PROBLEMATIC EIGENVECTOR OF B"
do i=1,8
write(*,*) i+12, omega_B(i,2),omega_B(i,4),omega_B(i,7),omega_B(i,8)
enddo
write(*,*) " debug 1"
!matrix of the correlation
do i=1,12
    do j=1,8
        if(abs(rot_C_mat(i,j)).lt.1.E-12) then
            psi_K_corr(i,j+12)=0
        else
            psi_K_corr(i,j+12)=rot_C_mat(i,j)
        endif
    enddo
enddo
write(*,*) " debug2"
!LER US BUILD THE EIGENVALUES OF THE BASIS
do i=1,12
    E_n(i)=wr_A(i)
enddo
do i=1,8
    E_n(i+12)=wr_B(i)
enddo
write(*,*) "debug 3"
!LET US BUILD THE FULL MATRIX R
do i=1,12
    do j=1,12
        R_mat(i,j)=omega_A(i,j)
    enddo
enddo

do i=1,8
    do j=1,8
        R_mat(i+12,j+12)=omega_B(i,j)
    enddo
enddo
write(*,*) "chjeijqi"
!Now I have to write the matrix M 
kn_term(1)=1
!the known term
M_mat(2,1)=1
M_mat(3,7)=1
M_mat(4,7)=1
M_mat(5,1)=1
M_mat(6,5)=-2
M_mat(6,9)=2
M_mat(13,4)=-1
M_mat(14,3)=1
M_mat(15,9)=-1
M_mat(16,10)=1
M_mat(17,2)=-1
M_mat(17,5)=1
M_mat(18,8)=1
M_mat(18,11)=-1
M_mat(19,5)=1
M_mat(19,6)=1
M_mat(20,11)=-1
M_mat(20,12)=-1
!At this point I have to write 
write(*,*) "feheoefhdp"
call generate_Q(E_n,psi_K_corr,Q,beta,mu,20)
call mprint(Q,20)
!At this point the MP1
call generate_MP1(MP1,R_mat,20)
call generate_MP2(MP2,E_n,M_mat,R_mat,beta,mu,20)
call generate_MP3(MP3,Q,R_mat,M_mat,20)
call generate_v1(v1,kn_term,E_n,R_mat,beta,mu,20)
call generate_v2(v2,kn_term,Q,R_mat,beta,mu,20)
M_tilde=MP1-MP2-MP3
v_tilde=v1+v2
write(*,*) "THE MATRIX Q"
call mprint(Q,20)

!let us solve the thing.
call solve_linear_system(M_tilde,v_tilde,solu,20)
write(*,*) "the solution to the linear system is"
do i=1,20
    write(*,*) i,solu(i)
enddo
!Let us evaluate the spectral weight of 
write(*,*) "The average of anticomm"
anticom_av=matmul(M_mat,solu)+kn_term
do i=1,20
write(*,*) i,anticom_av(i)
enddo
write(*,*) "spectrum of the theory"
do i=1,20
    write(*,*) i,E_n(i)
enddo
write(*,*) "double pole contribution",psi_K_corr(8,15)

Es=0
    do p=1,20
        Es=Es+R_mat(p,15)*anticom_av(p)
    enddo
write(*,*) "correction due to double pole",Es*psi_K_corr(8,15)
!Let us write the inverse of R
L_mat=R_mat
call inv(L_mat,20)
write(*,*) "inverse of the L matrix"

end program



















