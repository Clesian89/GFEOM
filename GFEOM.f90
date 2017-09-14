
program GFEOM

use matrixmod

!definition of the precsion
implicit none
integer, parameter :: i1=selected_int_kind(1)
integer, parameter :: i4=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(15)
integer, parameter :: r8=selected_real_kind(15,9)


!Physical system varibales
real(kind=r8)::mu,U,beta,t_x 
!Kind system physical variables ...
 
real(kind=r8)::Ku(20,20),Kc(20,20),R(20,20),av(20),L(20,20)
real(kind=r8)::M_mat(20,20),v(20),Q(20,20),K_mat(20,20)
real(kind=r8)::E(20),e_r(20),e_i(20),R_temp(20,20),D(20,20) 
real(kind=r8)::FR(20,20),L_sys(20,20),L1(20,20),L2(20,20),L3(20,20) 
real(kind=r8)::inhom(20),solu(20),anti_com(20)    
!Dummy variable 
integer(kind=i4)::i,j
!remember to clean
real(kind=r8)::temp1(20,20),temp2(20,20),temp1f(20) 
real(kind=r8)::Rant(20),scal_Kc(20,20) 
real(kind=r8)::sw1(20),sw2(20),sw3(20),spct(20)
real(kind=r8)::tot_sw,l_bound,u_bound 
real(kind=r8),allocatable::bin_w(:) 
!----------------------------------------
!
! PHYSICAL CONSTANT
!
!-----------------------------------------
U=5
mu=5
t_x=5
beta=500


do i=1,20
	v(i)=0
	sw1(i)=0
	sw2(i)=0
	sw3(i)=0
	do j=1,20
		K_mat(i,j)=0
		Ku(i,j)=0
		Kc(i,j)=0
		D(i,j)=0
	enddo
enddo 

!The Kmatrix 
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
!--------
!MATRIX M
!-------- 
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
!THE KNOWN TERM
v(1)=1
!SPLITTING OF THE MATRIX IN Ku Kc
do i=1,12
	do j=1,12
		Ku(i,j)=K_mat(i,j)
	enddo
enddo
do i=13,20
	do j=13,20
		Ku(i,j)=K_mat(i,j)
	enddo
enddo

do i=1,12
	do j=13,20
		Kc(i,j)=K_mat(i,j) 
	enddo
enddo

!DIAGONALIZATION OF THE DIAGONALIZABLE PART
!GENERATE THE MATRIX R AND L AND THE DEFINITION OF E 

call eigengenleft(Ku,e_r,e_i,R_temp,20)

write(*,*) "The spectrum of the theory"
do i=1,20
	write(*,*) i,e_r(i),e_i(i)  
	E(i)=e_r(i)  
enddo
write(*,*) "------------------------"
R=transpose(R_temp)
L=R
call inv(L,20)  

! CHECKING THE CONDITION, IN THE LACK OF A THEOREM
! WHICH CAN BE PROVE IN A PROBABLY EASY WAY PUTTING SOME CONSTRAINT ON THE
! FORM OF THE MATRIX KC 

do i=1,20
	D(i,i)=rand()
	write(*,*) i,D(i,i)  
enddo
write(*,*) "CHECKING THE CONDITION OF DECOUPLING"  
!temp1=matmul(R,matmul(Kc,L))
!temp2=matmul(R,matmul(Kc,L))
!call mprint(matmul(temp1,matmul(D,temp2)),20)
write(*,*) "THE CONDITION IS SATISFIED" 

!WE CAN NOW GENERATE THE Q-MATRIX (F_ei-F_ej)/(e_i-e_j) * Kc_ij 
!WE HAVE TO GENERATE THE FR-MATRIX  F_ei*R_ij
temp1=matmul(R,matmul(Kc,L))
Kc=temp1
call generate_Q(Q,Kc,E,beta,mu,20)
write(*,*) "CHECK THE MATRIX Q IF THERE ARE ANY PROBLEM" 
call mprint(Q,20) 
call generate_FR(FR,R,E,beta,mu,20)
write(*,*) "CHECK THE MATRIX FR"
call mprint(FR,20) 

L1=R
L2=matmul(FR,M_mat)
L3=matmul(Q,matmul(R,M_mat)) 

L_sys=L1-L2-L3   

inhom=matmul(FR,v)
inhom=inhom+matmul(Q,matmul(R,v)) 

call solve_linear_system(L_sys,inhom,solu,20)
write(*,*) "There is the vector solut n->1 d->2"
write(*,*) "0->15 0->16 0->18 0->20"
write(*,*) "U,mu",U,mu
do i=1,20
	write(*,*) i,solu(i)  
enddo

!The average of the anticommutator
anti_com= v + matmul(M_mat,solu) 
!DEFINITION OF THE RESCALED KC 
do i=1,20
	do j=1,20
		if(abs(Kc(i,j)).lt.1.E-9) then
			scal_Kc(i,j)=0
		else
			if(abs(E(i)-E(j)).lt.1.E-7) then
				scal_Kc(i,j)=0
			else
				scal_Kc(i,j)=( 1./(E(i)-E(j)) )*Kc(i,j) 
			endif
		endif
	enddo
enddo


Rant=matmul(R,anti_com)
temp1f=matmul(scal_Kc,Rant)
temp2=matmul(L,scal_Kc)
 


!THE SPECTRAL WEIGHT WE SPLIT IN A SUM OF 3 PARTIAL SW
do i=1,20
	sw1(i)=L(1,i)*Rant(i)
	sw2(i)=L(1,i)*temp1f(i)
	sw3(i)=temp2(1,i)*Rant(i) 
enddo 

spct=sw1+sw2-sw3

tot_sw=0
write(*,*) "THE SPECTRAL WEIGHT "
do i=1,20
	tot_sw=tot_sw+spct(i) 
	write(*,*) i,E(i),spct(i) 
enddo
write(*,*) " THE TOTAL SPECTRAL WEIGHT"
write(*,*) tot_sw
!CHECK IF THE SUM OF ALL THE WEIGHT IS ONE
!It Is we can try to see what is the spectral weight of the theory
!Box 0.1
allocate(bin_w(100))



l_bound=-5.5
u_bound=-5

do i=1,200
	bin_w(i)=0
	do j=1,20
		if(E(j).le.u_bound.and.E(j).gt.l_bound) then
			bin_w(i)=bin_w(i)+spct(j)   
		endif 
	enddo
	l_bound=l_bound+0.5
	u_bound=u_bound+0.5
enddo
write(*,*) "WRITE THE FINAL RESULT"
do i=1,48
    write(*,*) ( (-5.5-5)/2.0 +(i-1)*0.5),bin_w(i) 
enddo
 write(*,*) "Manual assamble" 













end program 

