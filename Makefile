HOST	=	nalab.mind.meiji.ac.jp
DIR	=	labo/text
DEST	=	$(HOST):Sites/$(DIR)
DEST2	=	$(HOST):Sites/program/fdm
PUSH	=	bin/pushtosakura
FILES	=	heat-fdm-1.dvi heat-fdm-1.pdf
FILES2	=	adi.pdf theta.pdf

default: all
all: $(PROGS) $(FILES)

heat-fdm-1.dvi: heat-fdm-1.tex heat1g-v2.c ../reference/reference.tex heat2d/*.c
	uplatex heat-fdm-1.tex
	upbibtex heat-fdm-1
	uplatex heat-fdm-1.tex
	uplatex heat-fdm-1.tex

heat-fdm-1.pdf: heat-fdm-1.dvi
	dvipdfmx -d 5 -O 2 heat-fdm-1.dvi

heat-fdm-1.ps.gz: heat-fdm-1.dvi
	dvips -f heat-fdm-1.dvi | gzip -9 > $@

heat2d.pdf: heat2d.tex
	uplatex heat2d.tex
	upbibtex heat2d
	uplatex heat2d.tex
	uplatex heat2d.tex
	dvipdfmx -d 5 -O 2 heat2d.dvi

copy: $(FILES)
	scp -p $(FILES) $(DEST)
	ssh $(HOST) $(PUSH) $(DIR)

copyprog: heat1g-v2.c
	scp -p heat1g-v2.c $(DEST2)
	ssh $(HOST) $(PUSH) program/fdm

adi.dvi: adi.tex
	uplatex adi.tex

theta.dvi: theta.tex
	uplatex theta.tex

adi.pdf: adi.dvi
	dvipdfmx -d 5 -O 2 adi.dvi

theta.pdf: theta.dvi
	dvipdfmx -d 5 -O 2 theta.dvi

copy2: $(FILES2)
	scp -p $(FILES2) $(DEST)
