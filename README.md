# RMFE 
Implementation in SageMath of polynomial-evaluation based Reverse Multiplication Friendly Embeddings (Cascudo/Cramer/Xing/Yuan, Crypto 18), see [article](https://eprint.iacr.org/2018/429.pdf)
SageMath is a free open-source math software based on Python.

**Requirements:**
- [SageMath](https://www.sagemath.org/download.html)
- Python 3


**Usage**
 - **Running the examples:**
 The file tests.py (folder files) can be run with the command:
`sage -python tests.py`

 - **Modifying it:**
 Easiest is to modify the .sage files, and then execute `./preparse.sh`(within the /files folder) . The .py files are automatically generated from the .sage files. 
This is except for lowlevel.py and pypolyfunctions.py which are directly written as Python files, and can be modified directly.
 

 




