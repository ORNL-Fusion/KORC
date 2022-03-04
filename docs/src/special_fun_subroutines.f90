FUNCTION chebev_s(a,b,c,x)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b,x
	REAL(SP), DIMENSION(:), INTENT(IN) :: c
	REAL(SP) :: chebev_s
	INTEGER(I4B) :: j,m
	REAL(SP) :: d,dd,sv,y,y2
	if ((x-a)*(x-b) > 0.0) call nrerror('x not in range in chebev_s')
	m=size(c)
	d=0.0
	dd=0.0
	y=(2.0_sp*x-a-b)/(b-a)
	y2=2.0_sp*y
	do j=m,2,-1
		sv=d
		d=y2*d-dd+c(j)
		dd=sv
	end do
	chebev_s=y*d-dd+0.5_sp*c(1)
END FUNCTION chebev_s


FUNCTION chebev_v(a,b,c,x)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(:), INTENT(IN) :: c,x
	REAL(SP), DIMENSION(size(x)) :: chebev_v
	INTEGER(I4B) :: j,m
	REAL(SP), DIMENSION(size(x)) :: d,dd,sv,y,y2
	if (any((x-a)*(x-b) > 0.0)) call nrerror('x not in range in chebev_v')
	m=size(c)
	d=0.0
	dd=0.0
	y=(2.0_sp*x-a-b)/(b-a)
	y2=2.0_sp*y
	do j=m,2,-1
		sv=d
		d=y2*d-dd+c(j)
		dd=sv
	end do
	chebev_v=y*d-dd+0.5_sp*c(1)
END FUNCTION chebev_v


SUBROUTINE beschb_s(x,gam1,gam2,gampl,gammi)
	USE nrtype
	USE nr, ONLY : chebev
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x
	REAL(DP), INTENT(OUT) :: gam1,gam2,gampl,gammi
	INTEGER(I4B), PARAMETER :: NUSE1=5,NUSE2=5
	REAL(SP) :: xx
	REAL(SP), DIMENSION(7) :: c1=(/-1.142022680371168_sp,&
		6.5165112670737e-3_sp,3.087090173086e-4_sp,-3.4706269649e-6_sp,&
		6.9437664e-9_sp,3.67795e-11_sp,-1.356e-13_sp/)
	REAL(SP), DIMENSION(8) :: c2=(/1.843740587300905_sp,&
		-7.68528408447867e-2_sp,1.2719271366546e-3_sp,&
		-4.9717367042e-6_sp, -3.31261198e-8_sp,2.423096e-10_sp,&
		-1.702e-13_sp,-1.49e-15_sp/)
	xx=8.0_dp*x*x-1.0_dp
	gam1=chebev(-1.0_sp,1.0_sp,c1(1:NUSE1),xx)
	gam2=chebev(-1.0_sp,1.0_sp,c2(1:NUSE2),xx)
	gampl=gam2-x*gam1
	gammi=gam2+x*gam1
END SUBROUTINE beschb_s


SUBROUTINE beschb_v(x,gam1,gam2,gampl,gammi)
	USE nrtype
	USE nr, ONLY : chebev
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(:), INTENT(OUT) :: gam1,gam2,gampl,gammi
	INTEGER(I4B), PARAMETER :: NUSE1=5,NUSE2=5
	REAL(SP), DIMENSION(size(x)) :: xx
	REAL(SP), DIMENSION(7) :: c1=(/-1.142022680371168_sp,&
		6.5165112670737e-3_sp,3.087090173086e-4_sp,-3.4706269649e-6_sp,&
		6.9437664e-9_sp,3.67795e-11_sp,-1.356e-13_sp/)
	REAL(SP), DIMENSION(8) :: c2=(/1.843740587300905_sp,&
		-7.68528408447867e-2_sp,1.2719271366546e-3_sp,&
		-4.9717367042e-6_sp, -3.31261198e-8_sp,2.423096e-10_sp,&
		-1.702e-13_sp,-1.49e-15_sp/)
	xx=8.0_dp*x*x-1.0_dp
	gam1=chebev(-1.0_sp,1.0_sp,c1(1:NUSE1),xx)
	gam2=chebev(-1.0_sp,1.0_sp,c2(1:NUSE2),xx)
	gampl=gam2-x*gam1
	gammi=gam2+x*gam1
END SUBROUTINE beschb_v
