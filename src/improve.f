      subroutine imprpara(maxcomp,maxorder,order,mix,andir,orient,lev,
     1                    vext,par,npar)
      implicit logical (a-z)
      integer maxcomp,maxorder,order(3,3,3),npar
      real*8 mix(maxorder,3,3,3),andir(3,maxorder,3,3,3),
     1       orient(2,maxorder,3,3,3),lev(2,3,3,3),par(npar),vext(3)
      integer i1,j1,k1,i2,j2,k2,o,o2,ord,ord2,icomp,besti,bestj,
     1        bestk,besto
      real*8 krit(6,3,3,3),dir(3),bestor(2),th,z,bestkrit
C   restricts maximum order to 6
      DO i1=1,3
         i2=4-i1
         DO j1=1,3
            j2=4-j1
            DO k1=1,3
               k2=4-k1
               if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
               dir(1)=(i1-2)*vext(1)
               dir(2)=(j1-2)*vext(2)
               dir(3)=(k1-2)*vext(3)
               ord = order(i1,j1,k1)
               DO o=1,ord
                  z = abs(dir(1)*andir(1,o,i1,j1,k1)+
     1                    dir(2)*andir(2,o,i1,j1,k1)+
     1                    dir(3)*andir(3,o,i1,j1,k1))
C  this is large if the component shows in direction of the voxel
                  ord2 = order(i2,j2,k2)
                  DO o2=1,ord2
                     z=z+mix(o,i2,j2,k2)*
     1               abs(andir(1,o,i1,j1,k1)*andir(1,o2,i2,j2,k2)+
     1                   andir(2,o,i1,j1,k1)*andir(2,o2,i2,j2,k2)+
     1                   andir(3,o,i1,j1,k1)*andir(3,o2,i2,j2,k2))
C this adds something if the same information is present on the opposite side
                  END DO
                  krit(o,i1,j1,k1)=mix(o,i1,j1,k1)*z
               END DO
            END DO
         END DO
      END DO
C  now prepare initial parameters
      icomp=maxcomp
      th=lev(1,2,2,2)-lev(2,2,2,2)
      DO while (icomp.gt.0) 
C
C   search for highest ranked component (direction)
C
         bestkrit=krit(1,1,1,1)
         bestor(1) = orient(1,1,1,1,1)
         bestor(2) = orient(2,1,1,1,1)
         besti=1
         bestj=1
         bestk=1
         besto=1
         DO i1=1,3
            DO j1=1,3
               DO k1=1,3
                  if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
                  ord = order(i1,j1,k1)
                  DO o=1,ord
                     if(krit(o,i1,j1,k1).gt.bestkrit) THEN
                        bestkrit=krit(o,i1,j1,k1)
                        bestor(1) = orient(1,o,i1,j1,k1)
                        bestor(2) = orient(2,o,i1,j1,k1)
                        besti=i1
                        bestj=j1
                        bestk=k1
                        besto=o
                     END IF
                  END DO
               END DO
            END DO
         END DO
         th=th+lev(1,besti,bestj,bestk)-lev(2,besti,bestj,bestk)
         par(2*(maxcomp-icomp+1))=bestor(1)
         par(2*(maxcomp-icomp+1)+1)=bestor(2)
C  mark this entry and all that are to close as uninteresting
         dir(1)=andir(1,besto,besti,bestj,bestk)
         dir(2)=andir(2,besto,besti,bestj,bestk)
         dir(3)=andir(3,besto,besti,bestj,bestk)
         DO i1=1,3
            DO j1=1,3
               DO k1=1,3
                  if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
                  ord = order(i1,j1,k1)
                  DO o=1,ord
                     if(dir(1)*andir(1,o,i1,j1,k1)+
     1                  dir(2)*andir(2,o,i1,j1,k1)+
     2                  dir(3)*andir(3,o,i1,j1,k1).gt..9) THEN
                        krit(o,i1,j1,k1)=0.d0
                     END IF
                  END DO
               END DO
            END DO
         END DO
         icomp=icomp-1
      END DO
      par(1)=th/(maxcomp+1)
      RETURN
      END
      subroutine imprparb(maxcomp,maxorder,order,mix,andir,orient,lev,
     1                    vext,par,npar,npar1)
      implicit logical (a-z)
      integer maxcomp,maxorder,order(3,3,3),npar,npar1
      real*8 mix(maxorder,3,3,3),andir(3,maxorder,3,3,3),
     1       orient(2,maxorder,3,3,3),lev(2,3,3,3),par(npar),vext(3)
      integer i1,j1,k1,i2,j2,k2,o,o2,ord,ord2,icomp,besti,bestj,
     1        bestk,besto,ord0,i,ord1
      real*8 krit(6,3,3,3),dir(3),bestor(2),th,z,bestkrit
C   restricts maximum order to 6
      ord0=order(2,2,2)
C  keep directions that are found to be informative
      th=2.5d0
      IF(ord0.gt.0) THEN
         th=(lev(1,2,2,2)-lev(2,2,2,2))
         DO i=1,ord0
            par(2*(maxcomp-i+1))=orient(1,i,2,2,2)
            par(2*(maxcomp-i+1)+1)=orient(2,i,2,2,2)            
         END DO
      END IF
      DO i1=1,3
         i2=4-i1
         DO j1=1,3
            j2=4-j1
            DO k1=1,3
               k2=4-k1
               if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
               dir(1)=(i1-2)*vext(1)
               dir(2)=(j1-2)*vext(2)
               dir(3)=(k1-2)*vext(3)
               ord1 = order(i1,j1,k1)
               DO o=1,ord1
                  z = abs(dir(1)*andir(1,o,i1,j1,k1)+
     1                    dir(2)*andir(2,o,i1,j1,k1)+
     1                    dir(3)*andir(3,o,i1,j1,k1))
C  this is large if the component shows in direction of the voxel
                  ord2 = order(i2,j2,k2)
                  DO o2=1,ord2
                     z=z+mix(o,i2,j2,k2)*
     1               abs(andir(1,o,i1,j1,k1)*andir(1,o2,i2,j2,k2)+
     1                   andir(2,o,i1,j1,k1)*andir(2,o2,i2,j2,k2)+
     1                   andir(3,o,i1,j1,k1)*andir(3,o2,i2,j2,k2))
C this adds something if the same information is present on the opposite side
                  END DO
                  krit(o,i1,j1,k1)=mix(o,i1,j1,k1)*z
               END DO
            END DO
         END DO
      END DO
C  make directions that are close to existing ones as uninteresting 
      DO i1=1,3
         DO j1=1,3
            DO k1=1,3
               if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
               ord = order(i1,j1,k1)
               DO o=1,ord
                  z=0.d0
                  if(ord0.gt.0) THEN
                     DO i=1,ord0
                        z=max(z,abs(
     1                     andir(1,i,2,2,2)*andir(1,o,i1,j1,k1)+
     1                     andir(2,i,2,2,2)*andir(2,o,i1,j1,k1)+
     2                     andir(3,i,2,2,2)*andir(3,o,i1,j1,k1)))
                     END DO
                  END IF
                  krit(o,i1,j1,k1)=krit(o,i1,j1,k1)*(1-z*z)
               END DO
            END DO
         END DO
      END DO
C  now prepare initial parameters
      icomp=maxcomp
      npar1=2*ord0+1
      DO while (icomp.gt.ord0) 
C
C   search for highest ranked component (direction)
C
         bestkrit=krit(1,1,1,1)
         bestor(1) = orient(1,1,1,1,1)
         bestor(2) = orient(2,1,1,1,1)
         besti=1
         bestj=1
         bestk=1
         besto=1
         DO i1=1,3
            DO j1=1,3
               DO k1=1,3
                  if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
                  ord = order(i1,j1,k1)
                  DO o=1,ord
                     if(krit(o,i1,j1,k1).gt.bestkrit) THEN
                        bestkrit=krit(o,i1,j1,k1)
                        bestor(1) = orient(1,o,i1,j1,k1)
                        bestor(2) = orient(2,o,i1,j1,k1)
                        besti=i1
                        bestj=j1
                        bestk=k1
                        besto=o
                     END IF
                  END DO
               END DO
            END DO
         END DO
         if(bestkrit.gt.0.2d0) THEN
            th=max(th,lev(1,besti,bestj,bestk)-lev(2,besti,bestj,bestk))
            par(2*(maxcomp-icomp+1))=bestor(1)
            par(2*(maxcomp-icomp+1)+1)=bestor(2)
            npar1=npar1+2
C  mark this entry and all that are to close as uninteresting
            dir(1)=andir(1,besto,besti,bestj,bestk)
            dir(2)=andir(2,besto,besti,bestj,bestk)
            dir(3)=andir(3,besto,besti,bestj,bestk)
            DO i1=1,3
               DO j1=1,3
                  DO k1=1,3
                     if(i1.eq.2.and.j1.eq.2.and.k1.eq.2) CYCLE
                     ord = order(i1,j1,k1)
                     DO o=1,ord
                        if(dir(1)*andir(1,o,i1,j1,k1)+
     1                     dir(2)*andir(2,o,i1,j1,k1)+
     2                     dir(3)*andir(3,o,i1,j1,k1).gt..9) THEN
                           krit(o,i1,j1,k1)=0.d0
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         ELSE
C  no more intersting direction reduce order 
            DO i=1,icomp
               par(2*(maxcomp-i+1))=0.d0
               par(2*(maxcomp-i+1)+1)=0.d0
            END DO
C            npar=2*(maxcomp-icomp)-1
            CONTINUE
            RETURN
         END IF
         icomp=icomp-1
      END DO
      par(1)=th
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine sweepimp(si,s0,n,ng0,ng1,siq,ms0)
C
C   for voxel in mask:
C   calculate mean s0 value
C   sweep s0 from si to generate  siq
C   calculate variance of siq
C
      integer n,ng0,ng1,si(n,ng1),s0(n,ng0)
      real*8 siq(n,ng1),ms0(n)
      integer i,k
      real*8 s0mean
      DO i=1,n
         z=0.d0
         DO k=1,ng0
            z=z+s0(i,k)
         END DO
         s0mean = z/ng0
         ms0(i) = s0mean
         DO k=1,ng1
            siq(i,k)=min(si(i,k)/s0mean,0.99d0)
         END DO
      call rchkusr()
      END DO
      RETURN
      END
      subroutine outlier1(si,mask,n,nb,s0ind,ls0,sinew,ind,lind)
C
C   same as subroutine outlier but restricted to a mask
C   replace physically meaningless Si values by mean S0
C
      implicit logical(a-z)
      integer n,nb,ls0,sinew(n,nb),ind(n),lind
      real*8 si(n,nb)
      logical s0ind(nb),mask(n)
      integer i,j,ls0m1
      real*8 s0
      logical changed
      ls0m1=ls0-1
      lind=0
      DO i=1,n
         if(mask(i)) THEN
            s0=0
            DO j=1,nb
               if(s0ind(j)) THEN
                  s0=s0+si(i,j)
                  sinew(i,j)=si(i,j)
               END IF
            END DO
            s0=(s0+ls0m1)/ls0
            changed=.FALSE.
            DO j=1,nb
               if(.not.s0ind(j)) THEN 
                  if(si(i,j).gt.s0) THEN
                     sinew(i,j)=s0
                     changed=.TRUE.
                  ELSE 
                     sinew(i,j)=si(i,j)
                  END IF
               END IF
            END DO
            if(changed) THEN
               lind=lind+1
               ind(lind)=i
            END IF
         END IF
      END DO
      call rchkusr()
      RETURN
      END
