SAGES=$(wildcard *.sage)
PYS=$(SAGES:%.sage=%.py)

test: $(PYS)
	@echo "Testing files"
	sage -c 'import sage_coding_project' 1>/dev/null

%.py: %.sage
	-rm $*.py $*.pyc
	sage -preparse $*.sage

clean:
	-rm -rf $(PYS) $(SAGES:%sage=%pyc)

README.html: README.md
	markdown README.md > README.html

