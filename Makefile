TEXFILES= notes.tex intro.tex fourseries.tex fourtrans.tex laptrans.tex \
		  anasigproc.tex digitanalog.tex discsig.tex filterdes.tex DFT.tex \
		  compDFT.tex outline.tex advDSP.tex

TEX=latex
CAT=/bin/cat
PROC=texproc
CHMOD=/bin/chmod 0400
RM=/bin/rm

all: notes.pdf

.SUFFIXES:.zd .tex .dvi .pdf

.zd.tex:
	$(RM) -f $*.tex
	$(CAT) $< | $(PROC) > $*.tex
	$(CHMOD) $*.tex

.tex.dvi:
	$(TEX) $<

.dvi.pdf:
	dvipdf $< 

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

notes.pdf: notes.dvi
	dvipdf notes.dvi

clean: 
	rm -f *.dvi *.log *.aux *.toc *.pdf *.err *.lot *.lof *.idx *.ilg *.tex *.ind
