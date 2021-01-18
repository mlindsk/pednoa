.PHONY: compile doc check install build_site run_main test clean

all: doc check

doc:
	Rscript -e "devtools::document()"

check: 
	Rscript -e "devtools::check()"

install:
	sudo Rscript -e "devtools::install()"

readme:
	Rscript -e "rmarkdown::render('README.Rmd')"

clean:
	rm -f README.html
