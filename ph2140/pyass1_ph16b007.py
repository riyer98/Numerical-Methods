#Narayana's Cows

def Ncow(n):
    (a,b,c)=(1,1,1)
    result=[1,1,1]
    for i in range(3,n):
        (a,b,c)=(b,c,c+a)
        result.append(c)
        i=i+1
    return result

#modulo program

def modfn(y,N):
    x=[y]
    a=0
    for i in range(0,N-1):
        a=(2*x[i])%1
        x.append(a)
        i+=1
    return x
#x(n>M) goes to zero because of the limit in the precision/number of decimal places in python,
#and the binary representation of the numbers. multiplying by 2 is multiplying by 10 in binary.
#by multiplying the remainder with 2 sufficient number of times, all the stored decimal places
#become integers in binary, thus giving a remainder of 0. this cannot happen in actuality in the
#unless the number is of the form a + (1/2^n) where a and n are integers(a,n>=0).
#the smallest possible number that python 3.6.2 recognises is 10^-323, which is close to 2^-1073.
#thus the max value of M  for which x(n>M)=0 is 1073.





#funtions for determining prime numbers less than a number

def prime1(n):
    if n<0:
        print('not applicable')
        return
    import math
    i=2;k=0;
    while i<n:
        j=2;
        while j<=math.sqrt(i):
            if i%j==0:
                break
            j+=1
        if j>math.sqrt(i):
            k+=1
        i+=1
    print(k)
    
#this program has an index which is set to 0 and increased by 1 for every prime no. < n
#the index j is for checking the factors for each i<n. we only need to check till j=greatest integer(sqrt(i))
#because suppose j*m=i, if j>sqrt(i) then m must be < sqrt(i), which has already been accounted for.
#if i%j==0 is false for 2<j<i-1 then i ust be prime, hence the index k is incremented by 1.

    
def prime2(n):
    if n<0:
        print("not applicable")
        return
    import math
    primes=[]
    i=2
    while i<n:
        j=2;
        while j<=math.sqrt(i):
            if i%j==0:
                break
            j+=1
        if j>(math.sqrt(i)):
            primes.append(i)
        i+=1
    return len(primes)
#this function uses the same logic as the previous one, except it appends all prime numbers to a list.
#then it returns the length of the list

#if pi(n) be no of primes < n, then:
#pi(10)=4
#pi(10**2)=25
#pi(10**3)=168
#pi(10**4)=1229
#pi(10**5)=9592

#for n=10**5 both programs take approximately 2 s
#for other n values both programs are nearly instantaneous.
