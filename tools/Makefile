CC=gcc
DESTDIR=/staff/edp/dwagon/bin.ultrix

all: texproc ldiv demo.dvi

texproc: proc.c
	$(CC) -o texproc proc.c

install:
	install -s texproc $(DESTDIR)

demo.dvi: tst.tex
	latex demo

tst.tex: ldiv
	ldiv tst > tst.tex
