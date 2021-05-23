# RMFE 
Implementation in SageMath of polynomial-evaluation based Reverse Multiplication Friendly Embeddings ([CCXY18] (https://eprint.iacr.org/2018/429.pdf)).

SageMath is a free open-source math software based on Python.

**Requirements:**
- [SageMath](https://www.sagemath.org/download.html) . See more info about usage [here](https://doc.sagemath.org/html/en/faq/faq-usage.html)
- Python 3

**Implementated notions:**
So far constructions of RMFE over F2 based on polynomial interpolation, both direct (as (k1,e1)2-RMFE, (k1,e1)2-RMFE where k1=2 or 3 and e1>=2+k1-1) and as concatenation of 2 direct interpolation-based construction, namely the concatenation of a (k1,e1)2-RMFE and a (k2,e2)2^e1-RMFE, where k1,e1 are as above, k2<=2^e1+1, and e2>=2*k2-1. See CCXY18.

**Main files:**
 - files/twostepinstance.sage:
 Parameters are specified by the class instance specified in file twostepinstance.sage. Defining an object in that class requires parameters k1,k2,e1,e2. It creates all parameters for a (k,e)-RMFE, stored as variables of the object, including k,e, the intermediate field F=F2^e1 and output field H=F2^e2 and their generator polynomials f, h (namely H=F_2[X]/(h), same with F and f). These polynomials are decided by SAGE, and can be modified as instance.h=instance.R(polynomial) where polynomial is written as X^i1+X^i2+...

 - files/twostepRMFE.sage:
 Contains functions computing the application of phi and psi of the RMFE (nevertheless, if the same RMFE is going to be used several times it is recommended to create the generator matrices with the file below, and compute results as a matrix vector multiplication)

 - files/matrices.sage:
 Creates generator matrices for a given instance of RMFE.



**Usage:**
 
 - **Modifying files:**
 Easiest is to modify the .sage files, and then execute `./preparse.sh`(within the /files folder) . The .py files are automatically generated from the .sage files. 
This is except for lowlevel.py and pypolyfunctions.py which are directly written as Python files, and can be modified directly.
 
 - **Examples:**
 A test file for the functions defined in files/twostepRMFE.sage can be found in files/
 Precomputed generator matrices and diverse data used for several articles can be found in folder files/output_data. They have been created with the file /files/matrices.sage. TODO: Explanation of the different files.
  
**Some References:**
 - [BMN17] (https://eprint.iacr.org/2018/395) Block, Maji, Nguyen: Secure computation with constant communication overhead using multiplication embeddings. Indocrypt 18.
 - [CCXY18] (https://eprint.iacr.org/2018/429.pdf) Cascudo, Cramer, Xing, Yuan: Amortized Complexity of Information-Theoretically Secure MPC Revisited. Crypto 18. 
 - [CG20] (https://eprint.iacr.org/2020/162.pdf) Cascudo, Gundersen: A Secret-Sharing Based MPC Protocol for Boolean Circuits with Good Amortized Complexity. TCC 20.
 - [PS21] (https://eprint.iacr.org/2020/1412.pdf) Polychroniadou, Song: Constant-Overhead Unconditionally Secure Multiparty Computation over Binary Fields. Eurocrypt 21.
 - [CG21] Cascudo, Giunta. Ongoing work. 




