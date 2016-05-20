##	Makefile
#
CC      = gcc#
CFLAGS	=-Wall -pedantic#

rm	= rm -f#				delete file

all:	_
clean::		;@ $(MAKE) T='$T' _clean
_clean:	_	;  $(rm) *.o $T a.out core *.tmp *.ps *.bak
run::	_
_:		;@ echo -------------------- $D --------------------


D =	Network analysis

MA=nw_reactant.o nw_massaction.o nw_massaction_system.o nw_cons.o

T = $(MA) JonesMann_JBC1994$x michaelismenten$x HockinMann_JBC2002$x HockinMann_JBC2002.o AplusBisC$x AplusBisC.o nw_cons$x reducedHM$x nw_data.o nw_data$x 
#nw_model.o nw_inhibit.o nw_enzyme.o nw_enzyme_ode.o nw_enzyme_ode$x
# nw_testMCMC$x
#compareHM_JBC2002$x compareIndividuals_JBC2002$x

all:	$T

GSLINC= `gsl-config --cflags`
#MCMC= $(HOME)c/misspiggy/mp_parameter.o $(HOME)c/misspiggy/mp_mcmc.o $(HOME)c/misspiggy/likelihood.o $(HOME)c/misspiggy/mp_proposal.o
GSL = `gsl-config --cflags` `gsl-config --libs` -lm

# nw_model.o:	nw_model.c nw_model.h
# 		$(CC) $(CFLAGS) -c -o $@ nw_model.c -I$(INCLUDE) -I$(INC) -I$(INC2)

# nw_enzyme.o:	nw_inhibit.h nw_inhibit.c nw_enzyme.c nw_enzyme.h
# 		$(CC) $(CFLAGS) -c -o $@ nw_enzyme.c -I$(INCLUDE) -I$(INC) -I$(INC2)

# nw_inhibit.o:	nw_inhibit.c nw_inhibit.h
# 		$(CC) $(CFLAGS) -c -o $@ nw_inhibit.c -I$(INCLUDE) -I$(INC) -I$(INC2)

# nw_enzyme_ode.o:	nw_enzyme_ode.c nw_enzyme_ode.h
# 		$(CC) $(CFLAGS) -c -o $@ nw_enzyme_ode.c -I$(INCLUDE) -I$(INC) -I$(INC2)

# nw_enzyme_ode$x:	nw_model.o nw_inhibit.o nw_enzyme.o nw_enzyme_ode.h nw_enzyme_ode.c
# 		$(CC) $(CFLAGS) -DTEST -o $@ nw_model.o nw_inhibit.o nw_enzyme.o nw_enzyme_ode.c -I$(INCLUDE) -I$(INC) -I$(INC2) $(GSL)

# nw_testMCMC$x:	nw_model.o nw_inhibit.o nw_enzyme.o nw_enzyme_ode.h nw_enzyme_ode.o $(MCMC) nw_testMCMC.c
# 		$(CC) $(CFLAGS) -DTEST -o $@ nw_model.o nw_inhibit.o nw_enzyme.o nw_enzyme_ode.o $(MCMC) nw_testMCMC.c -I$(INCLUDE) -I$(INC) -I$(INC2) $(GSL)

nw_massaction.o:	nw_massaction.c nw_massaction.h
			$(CC) $(CFLAGS) -c -o $@ nw_massaction.c $(GSLINC)

nw_massaction_system.o:	nw_massaction_system.c nw_massaction_system.h
			$(CC) $(CFLAGS) -c -o $@ nw_massaction_system.c $(GSLINC)

JonesMann_JBC1994$x:	JonesMann_JBC1994.h JonesMann_JBC1994.c $(MA)
			$(CC) $(CFLAGS) -o $@ $(MA) JonesMann_JBC1994.c $(GSL)

HockinMann_JBC2002.o:	HockinMann_JBC2002.h HockinMann_JBC2002.c
			$(CC) $(CFLAGS) -c -o $@ HockinMann_JBC2002.c $(GSLINC)

HockinMann_JBC2002$x:	HockinMann_JBC2002.h HockinMann_JBC2002.c $(MA)
			$(CC) $(CFLAGS) -DTEST -o $@ $(MA) HockinMann_JBC2002.c $(GSL)

michaelismenten$x:	michaelismenten.c $(MA)
			$(CC) $(CFLAGS) -o $@ $(MA) michaelismenten.c $(GSL)

nw_reactant.o:	nw_reactant.c nw_reactant.h
		$(CC) $(CFLAGS) -c -o $@ nw_reactant.c $(GSLINC)

AplusBisC$x:	AplusBisC.c $(MA)
			$(CC) $(CFLAGS) -DTEST -o $@ $(MA) AplusBisC.c $(GSL)

AplusBisC.o:	AplusBisC.h AplusBisC.c
		$(CC) $(CFLAGS) -c -o $@ AplusBisC.c $(GSLINC)

nw_cons.o:	nw_cons.c #nw_cons.h
		$(CC) $(CFLAGS) -c -o $@ nw_cons.c $(GSLINC)

nw_cons$x:	nw_cons.c nw_reactant.o nw_massaction.o nw_massaction_system.o
		$(CC) $(CFLAGS) -DTEST -o $@ nw_reactant.o nw_massaction.o nw_massaction_system.o nw_cons.c $(GSL)

reducedHM$x:	reducedHM.c HockinMann_JBC2002.o $(MA)
		$(CC) $(CFLAGS) -o $@ HockinMann_JBC2002.o $(MA) reducedHM.c $(GSL)

nw_data.o:	nw_data.h nw_data.c
		$(CC) $(CFLAGS) -c -o $@ nw_data.c $(GSLINC) 

nw_data$x:	nw_data.h nw_data.c
		$(CC) $(CFLAGS) -DTEST -o $@ nw_data.c $(GSLINC)

# compareHM_JBC2002$x:	compareHM_JBC2002.c $(MA) HockinMann_JBC2002.o nw_data.o
# 			$(CC) $(CFLAGS) -o $@ $(MA) HockinMann_JBC2002.o nw_data.o compareHM_JBC2002.c $(GSL) 

# compareIndividuals_JBC2002$x:	compareIndividuals_JBC2002.c $(MA) HockinMann_JBC2002.o nw_data.o
# 			$(CC) $(CFLAGS) -o $@ $(MA) HockinMann_JBC2002.o nw_data.o compareIndividuals_JBC2002.c $(GSL)
