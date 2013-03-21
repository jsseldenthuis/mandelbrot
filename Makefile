FC = gfortran
FFLAGS = -Wall -Wextra -march=native -O3 -ffast-math -fopenmp
LDFLAGS = -fopenmp
LIBS =

FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LIBS += $(shell pkg-config --libs plplotd-f95)

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

OBJS = \
	fractal.o

all: fractal

fractal: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMPILE) -o $@ -c $<

clean:
	$(RM) fractal $(OBJS) *.mod

