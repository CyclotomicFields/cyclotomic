paper.pdf: paper.tex
	pdflatex $<

paper.tex: ../src/cyclotomic.rs
	./src2latex.py < $< > $@

.PHONY: clean
clean:
	rm -f *.pdf *.aux *.log *.tex
