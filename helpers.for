      subroutine pintio(y,ifld)
c----------------------------------------------------------------------
c      Purpose: Character string input routine for data
c               N.B. This routine has largely been superceded by
c                    dinput functions.

c      Inputs:
c         lfld   - Field width to separate input data items into.

c      Outputs:
c         y      - Character string input, in field widths of ifld
c
c----------------------------------------------------------------------
      USE iofile
      USE iosave
      CHARACTER*80 X,Y
      CHARACTER*1 XX(80),YY(80)
      EQUIVALENCE (X,XX)
c.... read a record from the input file
100   if (ior.gt.0)  then
        read (ior,1000,err=901,end=902) x
      else
        read (  *,1000,err=901,end=902) x
      end if
      if(lsave) write(lfile,1000) x
c.... adjust and move the record to the y-array
c.kneb       call acheck(x,y,ifld,80,80)
      call acheck(xx,yy,ifld,80,80)
      DO 50 I = 1,80
  50  Y(I:I) = YY(I)
      return
c.... read error encountered
901   call  errclr ('PINTIO')
      stop 'ERROR in PINTIO at data input' !ww 
      goto  100
c.... eof encountered
c.kne 902   call  endclr ('PINTIO',x(1))
902   call  endclr ('PINTIO',x)
      goto  100
c.kne 1000  format(80a1)
1000  format(a)
      end

      subroutine acheck(x,y,n0,nl,nlc)
c----------------------------------------------------------------------
c      Purpose:   Parse a string to find fields separated by commas

c      Inputs:
c         x(*) -  Character string of data to parse
c         n0   -  Field width for parsed data
c         nl   -  Total number of characters in x-array
c         nlc  -  Total number of characters in y-array

c      Outputs:
c         y(*) -  Parsed data in field widths of n0
c----------------------------------------------------------------------
      character*1 x(*),y(*),macr*4

      do 100 ii = nl,1,-1
       if(x(ii).ne.' ') go to 110
100   continue

110   do 150 i = 1,nlc
       y(i) = ' '
150   continue

      k = 0
      il= 0
      do 200 i = 1,ii
        if(x(i).eq.',') then
          k  = k + n0
          if(k.gt.nlc-n0) go to 210
          il = k - i
        else
          y(i+il) = x(i)
        end if
200   continue
      k  = k + n0
cww210   call just(y,k,n0)
210   continue
c.... justify alphanumeric data in a string:- numbers right - alphanumeric remain left
c>>w
c.... do not in case of macros REST,1234 and SAVE,1234  (1234 is then a file extension)
      do im=1,4
        macr(im:im)=y(im)
      end do
      if(macr.eq.'save'.or.macr.eq.'rest'.or.
     +   macr.eq.'SAVE'.or.macr.eq.'REST') return
cww<<
      call just(y,k,n0)
      return
      end

      subroutine just(y,k,n0)
c----------------------------------------------------------------------
c
c      Purpose: Justify alphanumeric data in a string:
c               - Numbers are right justified
c               - Alphanumerics remain left justified
c
c      Inputs:
c         y*(*) - Unjustified string of data
c         k     - Length of string
c         n0    - Field width for justification
c
c      Outputs:
c         y*(*) - Justified string of data
c----------------------------------------------------------------------
      character*1 y(*)
      n1 = n0 - 1
      do i = 1,k,n0
c.... find last character in string
        do j = i,i+n1
          if(y(j).ne.' ') go to 100
        end do   
        y(i+n1) = '0'
100     if(y(i+n1).eq.' ') then
c....     identify a number in the field and right justify
          if((y(i).ge.'0'.and.y(i).le.'9') 
     +                   .or. (y(i).eq.'-')
     +                   .or.(y(i).eq.'+') 
     +                   .or. (y(i).eq.'.')) then
            do j = i+n1-1,i,-1
              if(y(j).ne.' ') go to 110
            end do  
110         nl = n1 + i - j
            do l = j,i,-1
              y(l+nl) = y(l)
              y(l) = ' '
            end do
          end if
        end if
      end do
      return

      end

      subroutine endclr (subnam,chr)
c----------------------------------------------------------------------
c
c      Purpose: end-of-file clearing routine
c
c      Inputs:
c         subnam   - name of Subroutine with error
c         chr      - error data
c
c      Outputs:    - error message 
c
c----------------------------------------------------------------------
      USE iofile
      character     subnam*(*),chr*(*)
c
      if (ior.gt.0)  then
         write(iow,2000) subnam
cww      stop
       return
      else
         backspace  5
         chr = ' '
      end if
      return
c.... error message, only for batch mode
 2000 format (' ** ERROR in ',a6,' ** end of file encountered')
      end

      subroutine    errclr (subnam)
c----------------------------------------------------------------------
c
c      Purpose: input error clearing routine
c
c      Inputs:
c         subnam   - name of Subroutine with error
c
c      Outputs:    - error message 
c
c----------------------------------------------------------------------
      USE iofile
      character     subnam*(*), yyy*80
c
      write(iow,2002) subnam
      write(  *,2002) subnam
 
cww   write(yyy,2002) subnam
cww   call drawmess(yyy,1,-2)
      return
2002  format (' ** ERROR in ',a6,' ** reinput last record')
      end

      logical function pcomp(a,b,n)
c-----------------------------------------------------------------
c      Purpose: Compare character strings for match
c               Ignores upper/lower case differences.
c
c      Inputs:
c         a(*)   - Character string 1
c         b(*)   - Character string 2
c         n      - Number of characters to compare
c
c      Outputs:
c         pcomp  - Flag, true if a = b
c-----------------------------------------------------------------

      character*(*) a,b
c      character*4 a,b
      pcomp = .false.
c.... compute the increment between an upper and a lower case letter
      inc = ichar('A') - ichar('a')
c.... compare for a match
      do 100 i = 1,n
         ia = ichar(a(i:i))
         ib = ichar(b(i:i))
c.... test all permutations of characters for a match
         if(ia.ne.ib .and. ia+inc.ne.ib .and. ia.ne.ib+inc ) return
100   continue
      pcomp = .true.
      return
      end
Cccc
      subroutine pconst(prt)
c-----------------------------------------------------------------
c
c     Purpose: Input parameter expressions:  let = expression
c              arithmetic free input routine
c              start from PMESH 
c
c      Inputs:
c         prt    - Print input values if true
c
c      Outputs:
c        Values of parameters a-z are stored in array vvv(i,j)
c        vvv(i, 1-26) Character a-z 
c        vvv(i,27-36) Character 0-9 
c        vvv(i,    0) Character ' ' 
c        vvv(1-26,k)  Character a-z
c        character position:
c        0=048 ... 9=057   ... i-47    =0
c        A=065 ... Z=090   ... i-64    =A 
c        a=097 ... z=122   ... i-64-32 =a
        
c        modified WW 11/05 to two characters
c-----------------------------------------------------------------
      USE conval
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      character*1 let*2,eql,x(75),xx*4
      logical pcomp,prt
cww   integer n,i,j,ial,izl,iau,izu,id,iq, i0,i9

c.... character positions
      i0  = ichar('0') ! 48 
      i9  = ichar('9') ! 57 
      iq  = ichar('=') ! 61
      ial = ichar('a') ! 97
      izl = ichar('z') !122
      iau = ichar('A') ! 65
      izu = ichar('Z') ! 90
      id  = ial - iau  ! 32
c
      if (prt) then
                     write(iow,2001)
        if(ior.lt.0) write(*  ,2001)
      end if
c.... input a record from file ior or the keyboard *
1     if(ior.gt.0) then
          read (ior,1000,err=901,end=902) let,eql,x
      else
          write(  *,3000)
          read (  *,1000,err=901,end=902) let,eql,x
      end if

c.... check length of string
      num=ipos(x,75) 
      if(num.eq.75) then 
         write(*,*) 'Warning: string too long'
         write(*,*) x      
         stop   
      end if

c     check for macro with one character 
      if(let(2:2).eq.'=') then
        do i = 75,2,-1
          x(i) = x((i-1))
        end do
        x(1)     = eql   
        eql      = '='
        let(2:2) = ' '
      end if

              
c...  check for a blank character or the null character = blank line
10    if(let(1:1).eq.' '.or.ichar(let(1:1)).eq.0) then
          x(1) = ' '
          let  = ' '
          return
      end if
c...  compare 'xx' for a match to 'list' = list values to screen
      xx(1:1) = let(1:1)
      xx(2:2) = let(2:2)
      xx(3:3) = eql
      xx(4:4) = x(1)
      if(pcomp(xx,'list',4)) then
        if(ior.lt.0) then
          write(*,2001)
          do i = 1,26
            if(vvv(i,0).ne.0.0d0) then  ! character,-
              write(*,2000) char(i+96),' ',vvv(i,0)
            end if
            do j = 1,26                 ! character,character
              if(vvv(i,j).ne.0.0d0) then
                write(*,2000) char(i+96),char(j+96),vvv(i,j)
              end if
            end do
            do j = 27,36
              if(vvv(i,j).ne.0.0d0) then ! character,digit
                write(*,2000) char(i+96),char(j+i0-27),vvv(i,j)
              end if
            end do
          end do
        end if
        go to 1
      end if
c.... Locate correct location for the addition in lower case letters
      i = ichar(let(1:1))
      if(ial.le.i .and. i.le.izl) i = i-96
      if(iau.le.i .and. i.le.izu) i = i-64
      if(i.gt.26) go to 901

      if(let(2:2).eq.' ') then
        j        = 0
      else
        j = ichar(let(2:2))
        if(ial.le.j .and. j.le.izl) then
          j = j-96
          if(j.gt.26) go to 901
        end if
        if(iau.le.j .and. j.le.izu) then
          j = j-64
          if(j.gt.26) go to 901
        end if
        if( i0.le.j .and. j.le.i9 ) then 
          j = j-21  !-47+26
          if(j.lt.27.or.j.gt.36) go to 901
        end if
      end if

      errck = .false.
      call setval(x,75, val)


c.... store value at calculated position
      vvv(i,j) = val

      if (prt) then
                     write(iow,2000) let(1:1),let(2:2),vvv(i,j)
        if(ior.lt.0) write(  *,2000) let(1:1),let(2:2),vvv(i,j)
      end if
      go to 1
c.... error on a read
901   call  errclr ('PCONST')
      if (ior.lt.0)  goto 1
      return
c.... eof encountered
902   call  endclr ('PCONST',let)
      goto  10
c.... formats
 1000 format(a2,76a1)
 2000 format(5x,'Parameter ',a1,a1,' = ',e15.8)
 2001 format(/'  p a r a m e t e r   v a l u e s')
 3000 format(' Use "list" to give current values - <CR> to exit'/
     1       ' Input: letter=expression (no blanks)'/'  -->',$)
      end
c
      subroutine setval(xi,num, val)
c----------------------------------------------------------------------
c
c      Purpose: Represent character constants by free inputs in strings
c      Inputs:
c        xi - input string - num-characters long
c        xs - string to search
c        xt - temporary string
c
c      Outputs:
c         val       - Value of string
c
c----------------------------------------------------------------------
      USE codat
      USE conval
      USE errchk
      USE iofile
      implicit double precision (a-h,o-z)
      logical errco
      character*1 xt(75),xs(75),xi(num)
      dimension v(25)
c.... check length of string
      if(num.gt.75) then
         write(*,*) 'Warning: string too long'
         write(*,*) xi
         stop
      end if
c.... read the value if no constants have been set
cwd   changed i=1,15 to i=1,21 to allow double precision input
      if(.not.coflg) then
        do 40 i = 1,21
cwd        do 40 i = 1,15
          xs(i) = ' '
40      continue
cwd     replaced 15 by 21->see above
        nex = 21 - num
cwd        nex = 15 - num
        do 50 i = 1,num
          xs(i+nex) = xi(i)
50      continue
        errco = .false.
        call pinval(xs,val,errco)
        if(errco) go to 60
        return
      end if
c.... find the letter number for this parameter
60    do 100 i = 1,75
        xs(i) = ' '
100   continue
      do 110 i = 1,num
        xs(i) = xi(i)
110   continue
c.... evaluate expression
1     nex = 0
      call parexp(xs,xt,v,nex,errck)
      if(errck) go to 150
      call pfuncs(xs,v,val,nex,errck)
      if(errck) go to 150
      return
c.... an error has been detected in the statement respecify
150   if(ior.lt.0) then
        write(*,2001) (xi(i),i=1,num)
151     read (*,1000,err=152,end=153) xt
        go to  154
152     call  errclr ('SETVAL')
        go to  151
153     call  endclr ('SETVAL',xt)
154     write(iow,2002) xt
        call pcheck(1,xt,errck)
        do 160 i = 1,74
          xs(i) = xt(i)
160     continue
        xs(75) = ' '
        errck = .false.
        go to 1
c.... error on a read
      else
        call  errclr ('SETVAL')
      end if
      return
c.... formats
 1000 format(75a1)
 2001 format(2x,a1,' = ',74a1/'   >',$)
 2002 format('  Correction:>',75a1)
      end
c
      subroutine pinval(xss,val,errco)
c----------------------------------------------------------------------
c
c      Purpose:  Moves character string into real value
c
c      Inputs:
c         xss(*)  - Character string
c
c      Outputs:
c         val     - Value extracted from character string
c         error   - Flag, true if error occurs
c
cwd  15 replaced by 21
c----------------------------------------------------------------------
      logical errco
      character*21 xs
      character*1  xss(21)
      real*8 val
      do 10 i=1,21
10    xs(i:i) = xss(i)
      read(xs,1000,err=100) val
      return
100   errco = .true.
      return
1000  format(f21.0)
      end
c
      subroutine parexp(x,xs,v,nex,error)
c----------------------------------------------------------------------
c      Purpose: Identify parenthetical expressions and evaluate
c
c      Inputs:
c         x(*)     - String containing expression to evaluate

c      Scratch:
c         xs(*)    - Array used to temporarily store expression
c         v(*)     - Array to hold values

c      Outputs:
c         x(*)     - Expression replaced by upper case letter
c         nex      - Number of upper case letters used
c         error    - Flag, true if error occurs
c      Common returns:
c         www(*)   - Upper case letters with values assigned
c----------------------------------------------------------------------

      USE conval
      implicit double precision (a-h,o-z)
      logical error
      character*1 x(*),xs(*)
      dimension v(*)
c.... find parenthetical expressions and remove
      do 130 i = 1,75
        if(x(i).eq.'(') then
          i1 = i + 1
          do 120 j = i1,75
            if(x(j).eq.'(') then
              call errclr('PAREXP')
cww           stop
            return
            else if(x(j).eq.')') then
              do 50 l = 1,j-i+1
                    xs(l) = ' '
50            continue
              i2 = j - 1
              if(i2.lt.i1) then
                call errclr('PAREXP')
cww                 stop
              return
              else
                k = 0
                    do 100 l = i1,i2
                      k = k + 1
                      xs(k) = x(l)
                      x(l)      = ' '
100                 continue
                    x(i2+1)  = ' '
c.... evaluate the expression in the parenthesis
                    call evalex(xs,v,val,k,error)
                    if(error) return
                    nex = nex + 1
                    www(nex) = val
c.... put an upper case letter in the expression and close up remainder
                    x(i) = char(nex +64)
                    i2 = i2 -i1 + 2
                    do 110 l = i1,75
                      x(l) = ' '
                      if(l+i2.le.75) then
                        x(l) = x(l+i2)
                      end if
110                 continue
              end if
              go to 125
            end if
120       continue
125       continue
        end if
130   continue
      return
      end
c
      subroutine evalex(xs,v,val,ns,error)
c----------------------------------------------------------------------
c     Purpose: Identify expression in character string and evaluate

c     Inputs:
c         xs(*) - Character string of input data
c         v(*)  - Array of real values
c         ns    - Length of character string

c     Outputs:
c        val   - Value of expression
c        error - Error indicator if true
c
c.... evaluate an expression: (+) add,   (-) subtract, (*) multiply,
c                             (/) divide,(^) power    
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical error
      character*1 xs(ns),x(80),y,op(25)
      dimension v(*)
c.... first pack the expression by removing any ' ' characters
      call pevpak(xs,ns, x,num)
      k  = 0
      do 100 i = 1,num
        if(k.eq.0 .and. x(i).eq.' ') go to 100
        if((x(i).eq.'+') .or. (x(i).eq.'-')) then
          if(i.gt.2) then
            y = x(i-1)
            if(y.eq.'e'.or.y.eq.'d'.or.y.eq.'E'.or.y.eq.'D') then
              if((x(i-2).ge.'0'.and.x(i-2).le.'9').or.
     1            x(i-2).eq.'.') go to 100
            end if
          end if
          k     = k + 1
          op(k) = x(i)
          x(i)  = ','
        end if
        if((x(i).eq.'*') .or. (x(i).eq.'/') .or. (x(i).eq.'^')) then
          k     = k + 1
          op(k) = x(i)
          x(i)  = ','
        end if
100   continue
      call dcheck(x,v,num,error)
      if(error) return
c.... compute the value of the expression
      val = v(1)
      if(k.ge.1) then
c..... 1. evaluate all exponentiations
        i = 1
120     continue
       if(op(i).eq.'^') then
          v(i) = v(i) ** v(i+1)
          k = k - 1
          do 110 j = i,k
            v(j+1) = v(j+2)
            op(j) = op(j+1)
110         continue
       end if
          i = i + 1
        if(i.le.k) go to 120
c..... 2. evaluate all multiplications and divisions
        i = 1
140     continue
          if(op(i).eq.'*') v(i) = v(i) * v(i+1)
          if(op(i).eq.'/') v(i) = v(i) / v(i+1)
          if(op(i).eq.'*' .or. op(i).eq.'/') then
            k = k - 1
            do 130 j = i,k
              v(j+1) = v(j+2)
              op(j) = op(j+1)
130         continue
          else
            i = i + 1
          end if
        if(i.le.k) go to 140
c..... 3. evaluate all additions and subtractions
        val = v(1)
        if(k.gt.0) then
          do 160 i = 1,k
            if(op(i).eq.'+') val = val + v(i+1)
            if(op(i).eq.'-') val = val - v(i+1)
160       continue
        end if
      end if
      return
      end
c
      subroutine pevpak(xs,ns, x,n)
c----------------------------------------------------------------------
c
c      Purpose: Remove unwanted blanks from character string
c      Inputs:
c         xs(*)   - Unpacked character string
c         ns      - length of unpacked string

c      Outputs:
c         x(*)    - Packed character string
c         n       - length of packed string
c
c----------------------------------------------------------------------
      character*1 xs(ns), x(ns)
      n = 0
      do 100 i = 1,ns
        if(xs(i).ne.' ') then
          n = n + 1
          x(n) = xs(i)
        end if
100   continue
      return
      end
c
      subroutine dcheck(x,vd,nt,error)
c----------------------------------------------------------------------
c      Purpose: Internal input of values from string data

c      Inputs:
c         x(*)  - Character array from inputs
c         nt    - Length of character string

c      Outputs:
c         vd(*) - Numerical values from string inputs
c         error - True of error occurs during inputs
c----------------------------------------------------------------------
      USE conval
      USE iofile
      implicit double precision (a-h,o-z)
      logical error
      character*1 x(nt)
      character*255 y
      character*1 yy(255)
      dimension ivd(2,16),vd(*)
c
c.... expression substitution using predefined constants
c
      n0 = 15
      nn = 1
      do 50 i = 1,nt
        if(x(i).eq.',') nn = nn + 1
50    continue
      call acheck(x,yy,n0,nt,255)
      call pcharr(yy,ivd,n0,nn*n0)
      error = .false.
      do 12 i=1,255
12    y(i:i) = yy(i)
      read(y,'(17f15.0)',err=200) (vd(i),i=1,nn)
      do 60 i = 1,nn
        if(ivd(1,i).gt.0) then
          vd(i) = vvv(ivd(1,i),ivd(2,i))
        else if(ivd(1,i).lt.0) then
          vd(i) = www(-ivd(1,i))
        end if
 60   continue

      return
c.... error
 200    write(iow,2000) nn,y(1:nn*15)
      if(ior.lt.0) then
        write(  *,2000) nn,y(1:nn*15)
      end if
      call errclr('DCHECK')
      error = .true.
      return
2000  format(/' *ERROR* attempting to input',i3,' value(s).',
     &        '  Input is:'/a)
      end
c
      subroutine pcharr(y,ivd,n0,nt)
c-----------------------------------------------------------------
c      Purpose: Identify 'character' variables in the string
c               N.B. Lower case variables have been input,
c                    upper case computed

c      Inputs:
c         y(*)     - String to search
c         n0       - Field width
c         nt       - Number of characters in y-array

c      Outputs:
c         ivd(2,*) - Number of characters found
c-----------------------------------------------------------------
      character*1 y(*)
      integer   n0,nt,n1,k,i,j,n,ivd(2,*),kk

      n1 = n0 - 1
      k = 0
      do i = 1,nt,n0
        k = k + 1
        ivd(1,k) = 0
        ivd(2,k) = 0
        n = ichar( y(i) ) - 64
        if(n.gt.0) then
          if(n.gt.58) go to 200
          if(n.gt.26) then
            ivd(1,k) = n - 32
            j = ichar(y(i+1))
            if(j.eq.32) then
              ivd(2,k) = 0
            elseif(j.ge.ichar('a') .and. j.le.ichar('z')) then
              ivd(2,k) = j - ichar('a') + 1
            elseif(j.ge.ichar('0') .and. j.le.ichar('9')) then
              ivd(2,k) = j - ichar('0') + 27
            endif
          else
            ivd(1,k) = - n
          end if


          do j = i,i+n1
            y(j) = ' '
c           1 arbitrary statement necessary for SALFORD without DEBUG Option!!?? 
            kk = 42           
          end do  

c....  set value
          y(i+n1) = '0'

        end if

      end do   


      return
c.... error
 200  call errclr('PCHARR')
      write(*,*) 'error character no: ',n,'field ',i, 'value',y(i) 
 
      return
      end
c
      subroutine pfuncs(x,v,val,nex,error)
c----------------------------------------------------------------------
c
c      Purpose: Evaluate expressions with functions in input records for
c               inc dec  int  abs   sqrt
c               sin sind asin asind sinh
c               cos cosd acos acosd cosh
c               tan tand atan atand tanh
c               log exp
c      extension d: input in degrees 
c
c      Inputs:
c         x(*)   - String of input
c         v(*)

c      Outputs:
c         nex    - Number of www(*) used to hold function value
c         error  - Flag, true if error occurs
c         val    - Expression value
c----------------------------------------------------------------------
      USE conval
      implicit double precision (a-h,o-z)
      character*1 x(*),yy(75),xx*5
      logical pcomp,error
      dimension v(*)

c.... rad->deg
      rtod = datan(1.d0)/45.d0 

c.... evaluate functions in expressions

      yy = ' '
      k  = 0
      i  = 1
140   continue
       k = k + 1
        xx(1:1) = x(i)
        xx(2:2) = x(i+1)
        xx(3:3) = x(i+2)
        xx(4:4) = x(i+3)
        xx(5:5) = x(i+4)
        if(     pcomp(xx,'atand',5))then
          kk = 1
        else if(pcomp(xx,'asind',5))then
          kk = 2
        else if(pcomp(xx,'acosd',5))then
          kk = 3
        else if(pcomp(xx,'atan',4))then
          kk = 4
        else if(pcomp(xx,'asin',4))then
          kk = 5
        else if(pcomp(xx,'acos',4))then
          kk = 6
        else if(pcomp(xx,'cosh',4))then
          kk = 7
        else if(pcomp(xx,'sinh',4))then
          kk = 8
        else if(pcomp(xx,'tanh',4))then
          kk = 9
        elseif(pcomp(xx,'cosd',4))then
          kk = 10
        else if(pcomp(xx,'sind',4))then
          kk = 11
        else if(pcomp(xx,'tand',4))then
          kk = 12
        else if(pcomp(xx,'sqrt',4))then
          kk = 13
        else if(pcomp(xx,'exp',3)) then
          kk = 14
        else if(pcomp(xx,'sin',3)) then
          kk = 15
        else if(pcomp(xx,'cos',3)) then
          kk = 16
        else if(pcomp(xx,'tan',3)) then
          kk = 17
        else if(pcomp(xx,'abs',3)) then
          kk = 18
        else if(pcomp(xx,'int',3)) then
          kk = 19
        else if(pcomp(xx,'log',3)) then
          kk = 20
        else if(pcomp(xx,'inc',3)) then
          kk = 21
        else if(pcomp(xx,'dec',3)) then
          kk = 22
        else
          kk = 0
        endif

c       Evaluate functions
        if(kk.ne.0 ) then
c         Functions 1 to 3
          if(kk.le.3) then
            j = 5
c         Functions 4 to 13
          elseif(kk.le.13) then
            j = 4
c         Functions 14 to 22
          else
            j = 3
          endif

          nn = ichar(x(i+j)) - 64
          if(nn.gt.26) then
            jj = nn + ichar(x(i+j+1))
            if    (jj.ge.ichar('a') .and. jj.le.ichar('z')) then
              ii = jj - ichar('a') + 1
            elseif(jj.ge.ichar('0') .and. jj.le.ichar('9')) then
              ii = jj - ichar('0') + 27
            else
              ii = 0
            endif

            val = vvv(nn-32,ii)
          else
            val = www(nn)
          endif

          nex = nex + 1
          if(     kk.eq.1) then
            www(nex) = atan(val)/rtod
          else if(kk.eq.2) then
            www(nex) = asin(val)/rtod
          else if(kk.eq.3) then
            www(nex) = acos(val)/rtod
          else if(kk.eq.4) then
            www(nex) = atan(val)
          else if(kk.eq.5) then
            www(nex) = asin(val)
          else if(kk.eq.6) then
            www(nex) = acos(val)
          else if(kk.eq.7) then
            www(nex) = cosh(val)
          else if(kk.eq.8) then
            www(nex) = sinh(val)
          else if(kk.eq.9) then
            www(nex) = tanh(val)
          else if(kk.eq.10) then
            www(nex) = cos(val*rtod)
          else if(kk.eq.11) then
            www(nex) = sin(val*rtod)
          else if(kk.eq.12) then
            www(nex) = tan(val*rtod)
          else if(kk.eq.13) then
            www(nex) = sqrt(val)
          else if(kk.eq.14) then
            www(nex) = exp(val)
          else if(kk.eq.15) then
            www(nex) = sin(val)
          else if(kk.eq.16) then
            www(nex) = cos(val)
          else if(kk.eq.17) then
            www(nex) = tan(val)
          else if(kk.eq.18) then
            www(nex) = abs(val)
          else if(kk.eq.19) then
            www(nex) = int(val)
          else if(kk.eq.20) then
            www(nex) = log(val)
          else if(kk.eq.21) then
            www(nex) = padd(val)
          else if(kk.eq.22) then
            www(nex) = psub(val)
          endif
          yy(k) = char(nex+64)
          i = i + j
        else
          yy(k) = x(i)
        endif
        i = i + 1
        if(i.lt.75) go to 140   ! fehler mit fullcheck 
c       if(i.lt.73) go to 140   

c.... final evaluation of expression
      call evalex(yy,v,val,k,error)
      return
      end
c
      subroutine dinput(d,nn)
c----------------------------------------------------------------------
c
c      Purpose: Data input subprogram for real values
c
c      Inputs:
c         nn    - Number of data items to input
c
c      Outputs:
c         d(nn) - Values for nn items input
c
c      Comments:
c         possible characters in line: 200
c         max number of items:          16
c         length(print format) of item: 15
c         for input with variables
c           max. length in line:        75
c
c----------------------------------------------------------------------
      USE errchk
      USE iofile
      USE iosave
      implicit double precision(a-h,o-z)
cww      character xxx(80)*1
      character xxx(201)*1
      dimension d(nn)
c.... check on number of items
      errck = .false.
      if(nn.gt.16) then
          if(ior.lt.0) then
            write(*,2000) nn
            errck = .true.
            return
          else
            write(iow,2000) nn
cww         stop
          return
          end if
      end if
      do 50 n = 1,nn
50      d(n) = 0.0d0
51    if(ior.gt.0) then
          read (ior,1000,err=901,end=902) xxx
      else
          read (*  ,1000,err=901,end=902) xxx
      end if
c.... check length of string
      if(xxx(201).ne.' ') then
         write(*,*) 'Warning: string too long'
         write(*,*) xxx
         stop
      end if
52    if(lsave) write(lfile,1000) xxx
c.... if no characters in list return
      if (xxx(1).eq.char(0))  return
cww      do 60 nl = 80,1,-1
      do 60 nl = 201,1,-1
          if(xxx(nl).ne.' ') go to 70
60    continue
c.... if blank characters in list return
      if (xxx(1).eq.' ')  return
70    no = 1
      nv = 1
c.... skip leading blank characters
      do 80 n = 1,nl
          if(xxx(n).ne.' ') go to 90
80    continue
c.... format separators are blanks or commas
90    if(xxx(n).eq.' '.or.xxx(n).eq.',') then
          if(n.gt.no) call setval(xxx(no),n-no,d(nv))
92        n = n + 1
          if(n.lt.nl.and.xxx(n).eq.' ') go to 92
          no = n
          nv = nv + 1
      else
          n  = n + 1
      end if
      if(n.le.nl.and.nv.le.nn) go to 90
c.... fill in last value if needed
      if(n.gt.no.and.nv.le.nn) call setval(xxx(no),n-no,d(nv))
      return
c.... read error encoutered
901   call  errclr ('DINPUT')
      goto  51
c.... eof encountered
902   call  endclr ('DINPUT',xxx(1))
      goto  52
c
cww1000  format(80a1)
1000  format(201a1)
2000  format(' ** ERROR ** too many items requested, limit = 16')
      end

      integer function ipos(file,nn)
c----------------------------------------------------------------------
c
c      Purpose: determine length of a character string
c
c      Inputs:
c         file(nn) - name
c
c      Outputs:
c         ipos     - length of character string
c    
c----------------------------------------------------------------------
      character*1 file(nn)
      do 100 n = nn,1,-1
        if(file(n).ne.' ') go to 200
100   continue
      ipos = 0
      return
200   ipos = n
      return
      end
c
      function padd(val)
c----------------------------------------------------------------------
c     for pfuncs: calculate increase of value      
c----------------------------------------------------------------------
c
      implicit  none
      real*8    padd, val, xval
      data      xval /0.0d0/
c     Look at parameter
      if(val.eq.0.0d0) then
        xval = 0.0d0
      else
        xval = xval + val
      endif
      padd = xval
      end
c
      subroutine pcheck(nc,xs,error)
c-----------------------------------------------------------------
c      Purpose: Check that input string contains admissible data
c               and parentheses match.  Convert all input letters
c               to lower case for further processing

c      Inputs:
c         nc     - Number of characters to check
c         xs(*)  - Character array

c      Outputs:
c         error  - Flag, true if error occurs
c-----------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical error
      character*1 x(75),xs(75)
c.... make sure that all of x is in lower case and blanks are removed
      i = 0
      do 100 j = 1,75
        x(j) = ' '
        if(xs(j).ne.' ' .and. xs(j).ne.'=' .and. xs(j).ne.',') then
          i = i + 1
          x(i)  = xs(j)
          xs(j) = ' '
          n = ichar( x(i) )
          if(n.ge.65 .and. n.le.90) x(i) = char(n + 32)
        end if
100   continue
c.... move back and check the characters for incorrect parenthesis
      error = .false.
      n = 0
      do 200 j = 1,i
        xs(j) = x(j)
        if(xs(j).eq.'(') n = n+1
        if(xs(j).eq.')') n = n-1
        if(n.lt.0 .or. n.gt.1 ) error = .true.
200   continue
      if(n.ne.0) error = .true.
      n = ichar(xs(1))
      if(n.lt.97 .or. n.gt.122) error = .true.
c.... check the characters for incorrect parameters
      if(.not.error) then
        do 210 j = 2,i
          n = ichar(xs(j))
          if(.not.(n.ge.97 .and. n.le.122) .and.
     1       .not.(n.ge.40 .and. n.le.57) ) then
            error = .true.
          end if
210     continue
      end if
      if(error) then
        write(*,2000)
      else
        write(*,2001) nc,(xs(j),j=1,i)
      end if
      return
 2000 format(' Incorrect statement - reinput ')
 2001 format('   No.',i3,'>',a1,' = ',74a1)
      end
c
      function psub(val)
c----------------------------------------------------------------------
c     for pfuncs: calculate decrease of value      
c----------------------------------------------------------------------
c
      implicit  none
      real*8    psub, val, xval
      data      xval / 0.0d0 /
c     Look at parameter
      if(val.eq.0.0d0) then
        xval = 0.0d0
      else
        xval = xval - val
      endif
      psub = xval
      end
c
