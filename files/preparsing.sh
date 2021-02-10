for FILE in FFTpreproc FFTa RMFE_23 RMFE_24 RMFE_35 field_iso tests 

do
	sage -preparse $FILE.sage
	mv $FILE.sage.py $FILE.py
done

