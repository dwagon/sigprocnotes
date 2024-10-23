DESTDIR=$(HOME)/bin
all: texproc

texproc: proc.c
	$(CC) -o texproc proc.c

install:
	install -s texproc $(DESTDIR)
