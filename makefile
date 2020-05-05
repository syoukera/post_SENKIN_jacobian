COMPILE = gfortran

COMLIB1 =

COMOTT1 = -O0 -Wall -Wextra

COMCONV = 

FILES = \
	calc_rop.o \
	cklib.o \



.SUFFIXES: .f .o

.f.o:
	$(COMPILE) $(COMOTT1) $(COMCONV) -c $<

calc_rope: $(FILES)
	$(COMPILE) $(COMOTT1) $(COMCONV) $(FILES) -o calc_rope $(COMLIB1)

calc_rop.o: calc_rop.f
cklib.o: cklib.f