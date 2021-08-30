PANDOC_FLAGS = --toc --standalone

all: README-PYPI.html README-test.html

README-PYPI.html: README.md
	python -m readme_renderer $< -o $@	

README-test.html: README.md
	pandoc -o $@ $< 

%.html: %.md
	pandoc $(PANDOC_FLAGS) $< -o $@  && open -g -a Safari $@
	fswatch -e ".*\.html" -o . | while read num ; do pandoc $(PANDOC_FLAGS) $< -o $@ && open -g -a Safari $@; done

.PHONY: all clean

clean:
	-rm -rf .nox
	-rm README*.html
	-find . -type d -name "__pycache__" -exec rm -rf {} \;
	-find . -type f -name "*.pyc" -delete
