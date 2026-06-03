#!bin/csh
################################################################################ 
# python: it might be required to install numpy package: sudo pip3 install numpy
################################################################################ 
# compiler:
F_COMPILER = gfortran
#FFLAGS = -e -O5
#FFLAGS = -m64 -mcmodel=medium -pedantic -Wall -Wsurprising -ffixed-line-length-132 -g
FFLAGS = -Wsurprising -mcmodel=large -ffixed-line-length-132 -O3
#FFLAGS = -pedantic -Wall -Wsurprising -ffixed-line-length-132 -g
# directorio de instalacion del ejecutable:
#dir_inst = ../
#dir_prog = ./lib
#
FSOURCE = aaa.f \
          a_alpha.f \
          a_IRTF.f \
          a_XSL_276A.f \
          a_NGSL.f \
	  dirs.f\
          hrsl_alpha.f \
          hrsl_IRTF.f \
          hrsl_XSL_276A.f \
          hrsl_NGSL.f \
          STU_alpha.f \
          STU_IRTF.f \
          STU_XSL_276A.f \
          STU_NGSL.f \
	  sigmam_alpha.f \
	  sigmam_IRTF.f \
	  sigmam_XSL_276A.f \
	  sigmam_NGSL.f \
          e_alpha.f \
          e_IRTF.f \
          e_XSL_276A.f \
          e_NGSL.f \
          fra_u_usdss.f\
          fraU.f\
          respu.f\
          respk.f\
          frausdss.f\
          frausdss004.f\
          respusdss.f\
          xflabs.f\
          fracat.f\
          respij.f\
          respic.f\
          FNORM.f \
          NICO.f \
          loca.f \
          matinx.f \
          limf.f \
          lbusc.f \
          sigmac.f \
          d4000.f \
          nIR.f \
          cattri.f \
          lcn1.f \
          lcn2.f \
          ca4227.f \
          lgband.f \
          fe4383.f \
          ca4455.f \
          fe4531.f \
          fe4668.f \
          hbeta.f \
          fe5015.f \
          mg1.f \
          mg2.f \
          mgb.f \
          fe5270.f \
          fe5335.f \
          fe5406.f \
          fe5709.f \
          fe5782.f \
          lnad.f \
          ltio1.f \
          ltio2.f \
          hgamaf.f \
          hdeltf.f \
          hgamaa.f \
          hdelta.f \
	  dpythag.f \
	  dsvbksb.f \
	  dsvdcmp.f \
	  respv.f \
	  respb.f \
	  hunt.f \
	  polint.f \
	  indexx.f \
	  remn.f \
	  CO.f
#          a_XSL.f \
#          hrsl_XSL.f \
#          STU_XSL.f \
#          sigmam_XSL.f \
#          e_XSL.f \

FOBJECT = aaa.o \
          a_alpha.o \
          a_IRTF.o \
          a_XSL_276A.o \
          a_NGSL.o \
	  dirs.o\
          hrsl_alpha.o \
          hrsl_IRTF.o \
          hrsl_XSL_276A.o \
          hrsl_NGSL.o \
          STU_alpha.o \
          STU_IRTF.o \
          STU_XSL_276A.o \
          STU_NGSL.o \
	  sigmam_alpha.o \
	  sigmam_IRTF.o \
	  sigmam_XSL_276A.o \
	  sigmam_NGSL.o \
          e_alpha.o \
          e_IRTF.o \
          e_XSL_276A.o \
          e_NGSL.o \
          fra_u_usdss.o\
          fraU.o\
          respu.o\
          respk.o\
          frausdss.o\
          frausdss004.o\
          respusdss.o\
          xflabs.o\
          fracat.o\
          respij.o\
          respic.o\
          FNORM.o \
          NICO.o \
          loca.o \
          matinx.o \
          limf.o \
          lbusc.o \
          sigmac.o \
          d4000.o \
          nIR.o \
          cattri.o \
          lcn1.o \
          lcn2.o \
          ca4227.o \
          lgband.o \
          fe4383.o \
          ca4455.o \
          fe4531.o \
          fe4668.o \
          hbeta.o \
          fe5015.o \
          mg1.o \
          mg2.o \
          mgb.o \
          fe5270.o \
          fe5335.o \
          fe5406.o \
          fe5709.o \
          fe5782.o \
          lnad.o \
          ltio1.o \
          ltio2.o \
          hgamaf.o \
          hdeltf.o \
          hgamaa.o \
          hdelta.o \
	  dpythag.o \
	  dsvbksb.o \
	  dsvdcmp.o \
	  respv.o \
	  respb.o \
	  hunt.o \
	  polint.o \
	  indexx.o \
	  remn.o \
	  CO.o
#          a_XSL.o \
#          hrsl_XSL.o \
#          STU_XSL.o \
#          sigmam_XSL.o \
#          e_XSL.o \
	  
# default rule to create program
a: $(FOBJECT)
	$(F_COMPILER) -o $@ $(FOBJECT) $(FFLAGS)
	mv $@ ../aaa

.f.o: $(FSOURCE)
	$(F_COMPILER) -c $(FSOURCE) $(FFLAGS)

