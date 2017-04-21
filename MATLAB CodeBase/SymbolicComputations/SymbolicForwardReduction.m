syms a_1 a_2 a_3 a_4 a_5 a_6
assume(a_1, 'real')
assume(a_2, 'real')
assume(a_3, 'real')
assume(a_4, 'real')
assume(a_5, 'real')
assume(a_6, 'real')

b_1 = a_3 - a_1;
b_2 = a_5 - a_2;
b_3 = a_6 - a_3;
M = [
     0  b_1  0  0   b_2   0   0 b_3 0 ; 
     0  b_2  0  0   b_1   0   0 b_2 0 ; 
     0  b_3  0  0   b_2   0   0 b_1 0 ; 
     0  0    0  b_1 b_2  b_3  0 0   0 ; 
     0  0    0  b_2 b_1  b_2  0 0   0 ; 
     0  0    0  b_3 b_2  b_1  0 0   0 ; 
     0  0    0  b_2 b_3  b_2  0 0   0 ; 
     ]
 
 [L,U] = lu(M)
