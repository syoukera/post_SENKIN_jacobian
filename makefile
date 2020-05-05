COMPILE = gfortran

COMLIB1 =

COMOTT1 = -O0 -Wall -Wextra

COMCONV = 

FILES = \
	src/calc_jac.o \
	src/cklib.o \



.SUFFIXES: .f .o .f90

.f.o:
	$(COMPILE) $(COMOTT1) $(COMCONV) -c $< -o $@

.f90.o:
	$(COMPILE) $(COMOTT1) $(COMCONV) -c $< -o $@ -J src

calc_jac: $(FILES)
	$(COMPILE) $(COMOTT1) $(COMCONV) $(FILES) -o calc_jac $(COMLIB1) -I src