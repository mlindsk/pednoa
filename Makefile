doc:
	Rscript -e "devtools::document()"

check: 
	Rscript -e "devtools::check()"

rhub_windows: 
	Rscript -e "rhub::check_on_windows()"

rhub_linux: 
	Rscript -e "rhub::check_on_linux()"

rhub_solaris: 
	Rscript -e "rhub::check_on_solaris()"

rhub_cran:
	Rscript -e "rhub::check_for_cran()"

install:
	cd ..; R CMD build pednoa/; \
	R CMD INSTALL pednoa_0.1.0.9999.tar.gz

readme:
	Rscript -e "rmarkdown::render('README.Rmd')"

clean:
	rm -f README.html
