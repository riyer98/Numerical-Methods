def prime(n):
    primes=[2]
    i=2
    while i<n:
        j=0
        while j<len(primes):
            if primes[j]>i**0.5:
                primes.append(i)
                break
            if i%primes[j]==0:
                break
            j+=1
        i+=1
    return (len(primes)-1)
