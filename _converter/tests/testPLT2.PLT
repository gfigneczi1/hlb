

A nobit 1 key=A_1 unit=potatoe scmin=-1000 scmax=1000 bit 0
//First test comment below the first line.
// And another one, just for fun
A nobit 1 key=A_2 unit=tomatoe scmin=-1000 scmax=1000 bit 1 



+ exp=x/2


B nobit 1 key=B_3 unit=suffer scmin=-800 scmax=800 bit 0
+ exp(x1=A_2)=(x/2)+x1/x-1
C nobit 1 key=C_4 unit=pain scmin=-200 scmax=200 bit 1
// Hey, another comment
D nobit 1 key=D_5 unit=pony scmin=-150 scmax=150 bit 0
+ exp(x1=B_3,x2=A_2)=(x1+x2+x)/3
// Ho- ho another comment
// A comment, again
E nobit 2 key=E_6 unit=dragon scmin=-300 scmax=300 bit 1
B nobit 5 key=B_7 unit=the_winter_is_coming scmin=-250 scmax=250 bit 0
D nobit 8 key=D_8 unit=a_lannister_always_pays_his_debts scmin=-1000 scmax=1000 bit 1
+ exp=x/3
// Comment after comment
E nobit 100 key=E_9 unit=for_frodo scmin=-70 scmax=7 bit 0


