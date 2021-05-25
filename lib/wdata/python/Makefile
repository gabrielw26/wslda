all: README-PYPI.html README-test.html

README-PYPI.html: README.md
	python -m readme_renderer $< -o $@	

README-test.html: README.md
	pandoc -o $@ $< 

.PHONY: all clean

clean:
	rm README*.html
