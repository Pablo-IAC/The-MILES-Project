c*******SUBRUTINA QUE LOCALIZA EL TIEMPO MAS PROXIMO DE LA ISOCRONA********
	SUBROUTINE loca(xx,n,x,k)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION xx(n)
	jl=0
	k=0
	ju=n+1
10	if (ju-jl.gt.1) then
		jm=(ju+jl)/2.d0
		if ((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
			jl=jm
		else
			ju=jm
		endif	
		go to 10
	endif
	j=jl
	if (abs(xx(j+1)-x).lt.abs(xx(j)-x)) then
		k=j+1
	else
		k=j
	endif
	return
	end
