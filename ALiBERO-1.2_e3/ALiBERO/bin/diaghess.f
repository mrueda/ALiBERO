c----------------------------------------------------------

c     Manuel Rueda, PhD 2006

      program diaghess

c     Diagonalizes a Hessian derived from EN-NMA
c     matrix, real, symmetrical

c     Adapted from diagstd.f included in elNemo software

c     Diagonalization routine: TQLI (EISPACK).
c    -simple, public domain, but slow (all eigenvectors are computed).

      implicit none

c----------maximum size of the system-----------------------
      integer maxnatom, maxdim
      parameter(maxnatom=2500,maxdim=3*maxnatom)

c----------name of file-------------------------------------
      character hessian_file*20
      character evec_file*20
      parameter(hessian_file='hessian.dat')
      parameter(evec_file='eigenvec.dat')

c----------variables for LAPACK-----------------------------
      character jobz, uplo
      integer lwork
      parameter(lwork=3*maxdim)
      real*8 work(lwork)
      integer dim
      integer info
c-----------------------------------------------------------
      real*8 a(maxdim,maxdim)
      real*8 w(maxdim)
      integer i, j, ii, jj
      real*8 f
      
c------------initialising the array--------------------------
      do i=1,maxdim
         do j=1,maxdim
            a(i,j)=0
         enddo
      enddo
      
c-----the Hessian does not necessarily contains all the elements
      open(file=hessian_file,form='FORMATTED',unit=10,status='OLD',
     $     err=999)
      dim=0
      do
         read(10,*,end=100) i, j, f
         a(i,j) = f
c-----symmetric---------------------------------------------------
         a(j,i) = f
         if( dim < i ) then
            dim = i
         endif
         if( dim < j ) then
            dim = j
         endif
      enddo
      
 100  continue

      write(*,'(A,I6)') 'The dimension: ',dim
      
      jobz = 'V'
      uplo = 'L'

c-----Diagonalization---------------------------------------------------
      call dsyev( jobz, uplo, dim, a, maxdim, w, work, lwork, info)

c-----Printing evecs-evals in ptraj format
      write(*,'(A,I6)') 'return info from dsyev: ',info
      open(file=evec_file,unit=12)      
      write(12,*) 'Eigenvector file: NMA'
      write(12,*) 'Contains ',dim, ' eigenvectors'
c     do j=7,dim
      do j=7,106
         jj=j-6
         write(12,*) '****'
         write(12,'(I5,F12.5)') jj,w(j)
         write(12,110) (a(i,j),i=1,dim)
      end do
      close(12)
 110  format(7f11.5)
      write(*,'(A,A)') 'ptraj based evectors in: ',evec_file
      
      stop

 999  continue
      write(*,'(A,A)') 'failed to open hessian file: ',hessian_file
      end
