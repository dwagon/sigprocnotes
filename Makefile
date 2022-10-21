# $Id: Makefile,v 1.5 1993/06/09 00:58:49 dwagon Exp $
# $Header: /ccstaff2/edp/dwagon/signotes/RCS/Makefile,v 1.5 1993/06/09 00:58:49 dwagon Exp $
# $Log: Makefile,v $
# Revision 1.5  1993/06/09  00:58:49  dwagon
# Force removal of old tex files before creating new ones to prevent
# overwrite failure due to read only chmoding
#
# Revision 1.4  1993/06/09  00:49:40  dwagon
# Added outline.tex to TEXFILES
#
# Revision 1.3  1993/06/07  04:35:55  dwagon
# TeX files chmoded read only after creation to prevent accidental editing
#
# Revision 1.2  1993/06/07  03:23:36  dwagon
# Changed file names to correspond with chapter names not numbers
#
# Revision 1.1  1993/03/24  03:54:10  dwagon
# Initial revision
#
TEXFILES= notes.tex intro.tex fourseries.tex fourtrans.tex laptrans.tex \
		  anasigproc.tex digitanalog.tex discsig.tex filterdes.tex DFT.tex \
		  compDFT.tex outline.tex advDSP.tex

TEX=/usr/local/texlive/2007/bin/i386-darwin/latex
#TEX=amslatex
CAT=/bin/cat
PROC=tools/texproc
CHMOD=/bin/chmod 0400
RM=/bin/rm

all: notes.dvi #notes.ps

.SUFFIXES:.zd .tex .dvi .ps

.zd.tex:
	$(RM) -f $*.tex
	$(CAT) $< | $(PROC) > $*.tex
	$(CHMOD) $*.tex

.tex.dvi:
	$(TEX) $<

.dvi.ps:
	dvips $< 

notes.dvi: $(TEXFILES)
	$(TEX) notes.tex
	makeindex notes.idx

index: $(TEXFILES)
	$(TEX) notes.tex
	$(TEX) notes.tex
	makeindex notes.idx
	$(TEX) notes.tex

it: $(TEXFILES)
	$(TEX) notes.tex

notes.ps: notes.dvi
	dvips notes.dvi

clean: 
	rm -f *.dvi *.log *.aux *.toc *.ps *.err *.lot *.lof *.idx *.ilg *.tex \
	*.ind
