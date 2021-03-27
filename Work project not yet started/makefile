all:
	make compile
	make clean
	make github

compile:
	latex Ha19_difference_full;
	bibtex Ha19_difference_full;
	latex Ha19_difference_full;
	pdflatex Ha19_difference_full;

clean:
	rm *.aux *.bbl *.blg *.log *.out *.svn *.snm
	rm -f *~ *.bak *.dvi *.ps *.bbl *.blg *.aux *.log *.thm  *.snm *.toc *.out *.nav

github:	
# Enter the command "sudo su" before running github	
	sudo git add --all
	git commit --all -m "I am still working on these topic, do I?"
	git push --all


