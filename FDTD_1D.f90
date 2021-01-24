Program FDTD_1D
Implicit none

Integer :: i    !!X coordinate
Integer :: j    !!Y coordinate
Integer :: n    !!Time!!

Real, Parameter :: C_S=5. !!Cell size
Real, Parameter :: T_D=0.5 !!Time Difference

!E-FIELD!!!!!!!!!!!!
Real :: E_z(10,10,10)
Real :: E_x(10,10,10)
Real :: J_z(10,10,10)

!H-FIELD!!!!!!!!!!!!
Real :: H_y(10,10,10)
Real :: M_y(10,10,10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Real, Parameter :: pi=Acos(-1.0)

!Material Propertiy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!I assume the lunar environment. 

Real, Parameter :: Conduct = 7*10**(-4)!! Conductivity  
Real, Parameter :: Magloss = 0         !! Magnetic loss  
Real, Parameter :: Permit  = 2.87      !! Permitivity  Free space=1/(36*pi))*10**(-9) 
Real, Parameter :: Permeat = 1.012     !! Permeativity Free space=(4*pi)*10**(-7) 

!E-Filed Component!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Real :: C_a(10,10,10) 
C_a=(1-(Conduct*T_D)/(2.*Permit)) / (1.+(Conduct*T_D)/(2.*Permit))

Real :: C_b(10,10,10)  
C_b=(T_D/(Permit*C_S)) / (1.+(Conduct*T_D)/(2.*Permit))


!H-Field Component!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Real :: D_a(10,10,10) 
D_a=(1-(Magloss*T_D)/(2.*Permeat)) / (1.+(Magloss*T_D)/(2.*Permeat))
    
Real :: D_b(10,10,10) 
D_b=(T_D/(Permeat*C_S)) / (1.+(Magloss*T_D)/(2.*Permeat))



!Initial Values!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
n=1
j=1
i=1

E_z(i/2-1/2,j/2+1/2,n/2-1/2)=10
H_y(i/2,j/2+1,n/2)=10
H_y(i/2-1,j/2+1/2,n/2)=8

!Path!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Open (unit=10,file='FDTD_1D.txt', status='replace')

!Time Step
 Do n=1,10
!Space Step
  Do j=1,10
   Do i=1,10
 
!Finite Difference!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
E_z(i/2-1/2,j/2+1/2,n/2+1/2)=C_a(i,j,k)*E_z(i/2-1/2,j/2+1/2,n/2-1/2)&
                            &+C_b(i,j,k)&
                            &*H_y(i/2,j/2+1,n/2)-H_y(i/2-1,j/2+1/2,n/2)&
                            &-J_z(i,j,k)*C_S

H_y(i/2,j/2+1/2,n/2+1) = D_a(i,j,k)*H_y(i/2,j/2+1/2,n/2)&
                       &+D_b(i,j,k),&
                       &*E_x(i/2,j/2+1/2,n/2+1/2)-E_x(i/2,j/2+1/2,n/2+1/2)&
                       &-M_y(i,j,k)*C_S
   End do 
  End do
 End do

End program FDTD_1D
