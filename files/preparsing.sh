for FILE in FFTpreproc FFTa field_iso testtwostep twostepRMFE twostepinstance polyintereval generatormatrix matrices RMFEmapsfrommatrices

do
	sage -preparse $FILE.sage
	mv $FILE.sage.py $FILE.py
done

