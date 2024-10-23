DESTDIR=$(HOME)/bin
all: install

texproc: proc.c
	$(CC) -o texproc proc.c

install: texproc
	install -s texproc $(DESTDIR)
