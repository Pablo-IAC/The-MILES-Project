        vktcs=0.059
	vjtcs=0.043
	jktcs=vktcs-vjtcs
	dkjtcs=0.043-(0.013/0.907)*(jktcs-0.01)
	vkj=(vktcs-0.05)/0.994
	write(*,*)'Kj-Ktcs =',dkjtcs,'V-Kj =',vkj
	end
