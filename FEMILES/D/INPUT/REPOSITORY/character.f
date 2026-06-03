    character*80 h,hh
    character*3 x
    h='hola'
    hh=h
    x='aa_'
c    hh=x//h
        write(hh(lnblnk(h)+1:lnblnk(hh)+1+5),'(A5)')'_LINE'
    write(*,*)hh
    end
