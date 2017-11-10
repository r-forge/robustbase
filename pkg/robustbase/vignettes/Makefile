###--- Makefile (adapted from pcalg package)

DRIVERS = lmrob_simulation psi_functions fastMcd-kmini

Extras = lmrob_simulation_gg2

#R_PROG = R-devel
R_PROG = R

## "Writing R Extensions" requires: 1st target must also clean up
doc: all
	${MAKE} clean

all: $(DRIVERS:=.tex) $(DRIVERS:=.pdf)

All: all $(Extras:=.pdf)

## Recreate the figures by running ggplot2 functions
## {MK gave up updating the ggplot2 function calls every few months
##  when Hadley decided to "improve" ggplot2 in incompatible ways}
lmrob_simulation_gg2.Rnw: lmrob_simulation.Rnw
	sed '/^<</s/,eval=FALSE,/,eval=TRUE,/' $< > $@

%.tex: %.Rnw
	${R_PROG} CMD Sweave $<

%.R: %.Rnw
	${R_PROG} CMD Stangle $<

## render pdf and compress it using qpdf
%.pdf: %.tex *.bib
	texi2pdf $<
	qpdf --stream-data=compress $@ temp$@
	mv temp$@ $@

purge:
	$(RM) Rplots.pdf $(DRIVERS:=.out) *.aux *.log *.blg
clean: purge
	$(RM) psi-*.pdf $(DRIVERS:=.tex) $(DRIVERS:=.bbl) $(DRIVERS:=-*.pdf)

realclean: clean
	$(RM) $(DRIVERS:=.pdf)
