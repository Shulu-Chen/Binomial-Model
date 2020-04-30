
```python
import math
from scipy.stats import norm
import time
import matplotlib.pyplot as plt
import numpy as np
```

### Binomial model function
```python
def Binomial(Option,K,T,S,sigma,r,q,N,Exercise):
    start=time.process_time()
    delta=T/N
    u=math.exp(sigma*math.sqrt(delta))
    d=math.exp(-sigma*math.sqrt(delta))
    p = (math.exp((r - q) * delta) - d) / (u - d)
    if Exercise=="E":
        if Option=="C":
            #European Call
            f = []
            for j in range(N + 1):
                S_Nj = S * pow(u, j) * pow(d, N - j)
                f_Nj = max(0, S_Nj - K)
                f.append(f_Nj)
            for n in reversed(range(N)):
                f_n = []
                for j in range(n+1):
                    f_nj = math.exp(-r*delta)*(p*f[j+1]+(1-p)*f[j])
                    f_n.append(f_nj)
                f = f_n
            option_price = f[0]
        elif Option=="P":
            #European Put#
            f=[]
            for j in range(N+1):
                S_Nj=S*pow(u,j)*pow(d,N-j)
                f_Nj=max(0,K-S_Nj)
                f.append(f_Nj)
            for n in reversed(range(N)):
                f_n=[]
                for j in range(n+1):
                    f_nj=math.exp(-r*delta)*(p*f[j+1]+(1-p)*f[j])
                    f_n.append(f_nj)
                f=f_n
            option_price=f[0]
        else:
            print("wrong Option input")
    elif Exercise=="A":
        if Option=="C":
            #American Call#
            f=[]
            for j in range(N+1):
                S_Nj=S*pow(u,j)*pow(d,N-j)
                f_Nj=max(0,S_Nj-K)
                f.append(f_Nj)
            for n in reversed(range(N)):
                f_n=[]
                for j in range(n+1):
                    f_nj=max(max(0,S*pow(u,j)*pow(d,n-j)-K),math.exp(-r*delta)*/
                             (p*f[j+1]+(1-p)*f[j]))
                    f_n.append(f_nj)
                f=f_n
            option_price=f[0]
        elif Option=="P":
            #American Put#
            f=[]
            for j in range(N+1):
                S_Nj=S*pow(u,j)*pow(d,N-j)
                f_Nj=max(0,K-S_Nj)
                f.append(f_Nj)
            for n in reversed(range(N)):
                f_n=[]
                for j in range(n+1):
                    f_nj=max(max(0,K-S*pow(u,j)*pow(d,n-j)),math.exp(-r*delta)*/
                             (p*f[j+1]+(1-p)*f[j]))
                    f_n.append(f_nj)
                f=f_n
            option_price=f[0]
        else:
            print("wrong Option input")
    else:
        print("wrong Exercise input")
    option_price=round(option_price,3)
    end=time.process_time()
    return option_price,start-end
```

### Function to Count step needed for $10^3$ accuracy
```python
def count_step(op,qt):
    T3=np.arange(1/12,13/12,1/12)
    step_list=[]
    p_tmp=0
    n=1
    for i in T3:
        for j in np.arange(n,5000,3):
            p=Binomial(op,100,i,100,0.2,0.05,qt,j,"A")
            if abs(p-(p+p_tmp)/2)<0.001:
                step_list.append(j)
                print(i,j)
                n=j
                break
            p_tmp=p
    return step_list
```

### Function to count early exercise boundary for puts
```python
def count_put_s(q,N):
    n=10000
    s_star=[]
    for j in range(0,12):
        for i in reversed(range(n)):
            intrinsic_value_tmp=max((100 - S[i]), 0)
            p_tmp=Binomial("P",100,T[j],S[i],0.2,0.05,q,100,"A")
            dif_tmp=abs(intrinsic_value_tmp-p_tmp)
            if dif_tmp<0.005:
                print(S[i],"i")
                st=i
                for k in reversed(range(st)):
                    p=Binomial("P",100,T[j],S[k],0.2,0.05,q,N,"A")
                    print(S[k],"k")
                    intrinsic_value=max((100-S[k]),0)
                    dif=abs(intrinsic_value-p)
                    print(dif)
                    if dif<0.005:
                        s_star.append(S[k])
                        print(S[k])
                        n=k
                        break
                break
    return s_star
```

### Function to count early exercise boundary for calls
```python
def count_call_s(q,N):
    n=10000
    s_star=[]
    for j in range(0,12):
        for i in range(n,20000):
            intrinsic_value_tmp=max((S[i]-100), 0)
            p_tmp=Binomial("C",100,T[j],S[i],0.2,0.05,q,100,"A")
            dif_tmp=abs(intrinsic_value_tmp-p_tmp)
            if dif_tmp<0.005:
                print(S[i],"i")
                st=i
                for k in range(st,20000):
                    p=Binomial("C",100,T[j],S[k],0.2,0.05,q,N,"A")
                    print(S[k],"k")
                    intrinsic_value=max((S[k]-100),0)
                    dif=abs(intrinsic_value-p)
                    print(dif)
                    if dif<0.005:
                        s_star.append(S[k])
                        print(S[k])
                        n=k
                        break
                break
    return s_star
```

### Compute Black-Scholes Option Price
```python
Option2="C"
K2=100
T2=1
S2=100
sigma2=0.2
r2=0.05
q2=0.04
N2=range(1,200)
Exercise2="E"
d1 = (math.log(S2/K2)+(r2-q2+0.5*pow(sigma2,2))*T2)/(sigma2*math.sqrt(T2))
d2 = d1-sigma2*math.sqrt(T2)
BS_price = S2*math.exp(-q2*T2)*norm.cdf(d1)-K2*math.exp(-r2*T2)*norm.cdf(d2)
```

### Compute Number of Time Steps needed for Accuracy
```python
step_list1=count_step("P",0)
step_list2=count_step("P",0.04)
step_list3=count_step("C",0.04)
step_list4=count_step("C",0.08)
```

### Compute Boundary $S^*$
```python
s_star1=count_put_s(0,940)
s_star2=count_put_s(0.04,3043)
s_star3=count_call_s(0.04,2305)
s_star4=count_call_s(0.08,1696)
```

### Plot Binomial Price compared with Black-Schole Price
```python
p_list=[]
error_list=[]
for i in N2:
    p=Binomial(Option2,K2,T2,S2,sigma2,r2,q2,i,Exercise2)
    p_list.append(p)
    error_list.append(abs(p-BS_price))
plt.subplot(2, 1, 1)
plt.plot(N2, p_list,label='Binomial price')
plt.axhline(y=BS_price,label='Black-Scholes Price',c="orange")
plt.annotate(s=round(BS_price,3) ,xy=(0,8) ,xytext=(190,7.7))
plt.ylabel('European Call Price')
plt.xlabel("Number of time steps")
plt.legend()
plt.subplot(2, 1, 2)
plt.plot(N2, error_list,label='Error',c="red")
plt.ylabel('Error')
plt.xlabel("Number of time steps")
plt.legend()
plt.show()
```

### Plot Number of time steps V.S running time figures
```python
Option2="P"
K2=100
T2=1
S2=100
sigma2=0.2
r2=0.05
q2=0.04
N2=range(1,2000,10)
time_e=[]
time_a=[]
time_ec=[]
time_ac=[]
for i in N2:
    p1,t1=Binomial("P",K2,T2,S2,sigma2,r2,q2,i,"E")
    p2,t2=Binomial("P",K2,T2,S2,sigma2,r2,q2,i,"A")
    p3,t3=Binomial("C",K2,T2,S2,sigma2,r2,q2,i,"E")
    p4, t4 = Binomial("C", K2, T2, S2, sigma2, r2, q2, i, "A")
    time_e.append(t1)
    time_a.append(t2)
    time_ec.append(t3)
    time_ac.append(t4)
plt.plot(N2, time_e,label='European Put Time')
plt.plot(N2, time_a,label='American Put Time')
plt.plot(N2, time_ec,label='European Call Time')
plt.plot(N2, time_ac,label='American Call Time')
plt.title('Compute time VS Step numbers')
plt.ylabel('Compute time')
plt.xlabel("Number of time steps")
plt.legend()
plt.show()
```

### Plot option price V.S. $S_0$
```python
k3=100
S3=np.arange(0,200,0.1)
intrinsic_value=[]
p_list3=[]
p_list4=[]
gap=[]
for i in range(2000):
    v=max((S3[i]-k3),0)
    intrinsic_value.append(v)
    p=Binomial("C",k3,1,S3[i],0.2,0.05,0.04,10,"A")
    p2=Binomial("C",k3,1,S3[i],0.2,0.05,0.08,10,"A")
    p_list3.append(p)
    p_list4.append(p2)
plt.plot(S3, intrinsic_value,label="Intrinsic Value")
plt.plot(S3,p_list3,label="American Call Price, q=0.04")
plt.plot(S3,p_list4,label="American Call Price, q=0.08")
plt.title('12-Month American Call Price VS S_0')
plt.ylabel('American Call Price')
plt.xlabel("S_0")
plt.legend()
plt.show()
```
