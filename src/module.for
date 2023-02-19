      module codat  ! rr
c-----------------------------------------------------------------------
c.... Purpose: codat.h
c-----------------------------------------------------------------------
      logical coflg
      end module

      module conval ! rr
c-----------------------------------------------------------------------
c.... Purpose: conval.h
c        vvv(i, 1-26) - Character a-z 
c        vvv(i,27-36) - Character 0-9 
c        vvv(i,    0) - Character ' ' 
c        vvv(1-26,k)  - Character a-z
c        www(26)      - Upper case letters with values assigned
c     update ww KIT 11/14
c-----------------------------------------------------------------------
      real*8 vvv(26,0:36),www(26)
      end module

      module errchk ! rr
c-----------------------------------------------------------------------
c.... Purpose: errchk.h
c-----------------------------------------------------------------------
      logical         errck
      end module

      module iofile ! rr
c-----------------------------------------------------------------------
c.... Purpose: iofile.h
c-----------------------------------------------------------------------
      integer ior,iow
      end module

      module iosave
c-----------------------------------------------------------------------
c.... Purpose: iosave.h
c-----------------------------------------------------------------------
      integer lfile
      logical lread,lsave
      end module

