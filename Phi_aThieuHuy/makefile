 all:
	make compile
	make clean
	make github

compile:
	latex Huy_Phi_2019;
	bibtex Huy_Phi_2019;
	latex Huy_Phi_2019;
	pdflatex Huy_Phi_2019;

clean:
	rm *.aux *.bbl *.blg *.log *.out *.svn *.snm
	rm *~ *.bak *.dvi *.ps *.bbl *.blg *.aux *.log *.thm  *.snm *.toc *.out *.nav

github:	
# Enter the command "sudo su" before running github	
	sudo git add --all
	git commit --all -m "I am writing"
	git push --all

