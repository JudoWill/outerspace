testrun: a.txt

c.txt: 
	touch c.txt

b.txt: c.txt
	touch b.txt

a.txt: b.txt
	touch a.txt

cleantest:
	rm -f ?.txt
	#? means any one character
