#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 17:52:29 2018

@author: Tyler
"""

from scipy.special import factorial as fact
import numpy as np
import matplotlib.pyplot as plt


def Lambda(l,k,n,m,Print=False):    
    p=choose(n-m,k-l)*choose(m,l)*fact(k)*fact(n-k)/(fact(n))
    if Print ==True:
        print(p,l,k,n,m)
    return p



def choose(n,k):
    return fact(n)/(fact(n-k)*fact(k))




def generate_distribution(n):
    P=[]
    #P.append(0.20)
    #P.append(0.20)
    #P.append(0.40)
    #P.append(0.20)
  
    for i in range(0,int(np.around(n-len(P)))):
        P.append(1/n)
    

    
    #for i in range(0,int(np.around(n))):
     #   P.append(1/n/2)

    

    return P




def He(l,k,n,m):
    Lam=0
        
    if l<=k and l<=m and m+k-2*l>=0 and n>=m+k-2*l:
        
        if l>=k+m-n:
        
            #Lam+=Lambda(l-s,k,n,m)*binom(k-l+s,p,s)
            Lam+=Lambda(l,k,n,m)
            
    return Lam




    

def Lambda_dom(l,k,n,m,p,other_parent,epsilon,mark):
    Lam=0
    

    H=(0,0,0)
    H=np.asarray(H,dtype=float)
    for s in range(0,l+1):
        
        if l-s<=k and l-s<=m and m+k-2*l+2*s>=s and n>=m+k-2*l+2*s:
            
            if l-s>=k+m-n:
            
                L=Lambda(l-s,k,n,m)*binom(m+k-2*l+2*s,p,s)
                Lam+=L
                
                if mark==True:
                    L1=L*(l-s)/(n)
                    L2=L*(n-k-m+l-s)/(n)
                    L3=L*(m+k-2*(l-s))/(n)
                    
                    if L1<0 or L2<0 or L3<0:
                        raise ValueError('Values less than Zero')
                    
                    H[0]+=L1
                    H[1]+=L2
                    H[2]+=L3
     
    LLam=Lam*(1+epsilon/(np.abs(other_parent-l)+0.01))

    return Lam,H



'''
D=[]
for i in range(0,6):
    D.append(Lambda_dom(i,5,10,5))
plt.plot(np.linspace(0,6,6),D)
'''
def binom(n,p,k):
    return choose(n,k)*p**k*(1-p)**(n-k)   
    
def generateP_distribution(H,n):
    P=[]
    for i in range(0,n):
        P.append(binom(n-1,H[0]+0.5*H[2],i))
        
    return P
        
    
    


def Season(P,H,n,b,p,x,epsilon,beta,alpha):
    P=generateP_distribution(H,n)
   
    
    d=0
    P_new=np.zeros(len(P))
    
    H=np.asarray((0,0,0),dtype=float)
    
    for i in range(0,n):     
        for j in range(0,n):           
            for K in range(0,n):                                   
                    for k2 in range(0,n): 
                        if j<=i:
                            Game=(game(K,b,p,player=1,beta=beta,alpha=alpha)+game(k2,b,p,player=2,beta=beta,alpha=alpha))/2
                            Lambda1=Lambda_dom(K,i,n-1,j,x,k2,epsilon,mark=True)
                            Lambda2=Lambda_dom(k2,i,n-1,j,x,K,epsilon,mark=False)
                            psi=Lambda1[0]*Lambda2[0]*P[i]*P[j]*Game
                                       
                            P_new[K]+=psi
                            d+=psi

                           
                            total_P=Lambda2[0]*P[i]*P[j]*Game
                            
                            Lambda1[1][0]=Lambda1[1][0]*total_P
                            Lambda1[1][1]=Lambda1[1][1]*total_P
                            Lambda1[1][2]=Lambda1[1][2]*total_P
                                   
                            H+=Lambda1[1]
                           


                        if j>i:
                            Lambda1=Lambda_dom(K,j,n-1,i,x,k2,epsilon,mark=True)
                            Lambda2=Lambda_dom(k2,j,n-1,i,x,K,epsilon,mark=False)
                            
                            psi=Lambda1[0]*Lambda2[0]*P[i]*P[j]*Game
                                       
                            P_new[K]+=psi
                            d+=psi
                            
                            
                            total_P=Lambda2[0]*P[i]*P[j]*Game
                            
                            Lambda1[1][0]=Lambda1[1][0]*total_P
                            Lambda1[1][1]=Lambda1[1][1]*total_P
                            Lambda1[1][2]=Lambda1[1][2]*total_P
                            H+=Lambda1[1]
                            
    
    
    return P_new/d,H/d




                    
def game(K,b,p,player,beta,alpha):
    #alpha=[alpha*K,K]
    
    #if alpha!= None:
     #   K+=np.random.randn()*alpha[0]+alpha[1]
      #  K=K/2
        
        
    if player ==2:
        return p+K*beta*b
    if player ==1:
        return p-K*b
    
    

def many_seasons_meanVar(H,n,b,p,x,num_seasons,beta,epsilon,alpha):
    
    P=generate_distribution(n)
    
    
    H=[0.33,0.34,0.33]
    P=generateP_distribution(H,n)
    
    mean_all=[]
    std_all=[]
    plt.figure(1)
    diff=0.1
    
    for i in range(0,num_seasons):
        
        
        beta1=np.random.randn()+beta
                             
        if i <8:
            beta1=beta
        
        old_mean=np.mean(P)

        P,H=Season(P=P,H=H,n=n,b=b,p=p,x=x,epsilon=epsilon,beta=beta,alpha=alpha)
        
        
        #plt.subplot(211)
        
        #plt.plot(np.linspace(0,b*len(P),len(P)),P,label=str(p))
        #plt.show()
        
        
        mean=np.sum(P*np.linspace(0,len(P),len(P)))
        
        mean_all.append(mean)
        std=np.std(P)
        std_all.append(std)
        
        if i>4:
            diff=mean-old_mean
        
        if np.around(diff,3)==0:
            i=num_seasons
       
    mean_all=np.asarray(mean_all)
    #plt.subplot(212)
    #plt.title(str(beta))
    plt.plot(np.linspace(0,len(mean_all),len(mean_all)),mean_all,label='x='+str(x)+'--beta:'+str(beta))
    #plt.show()
   
    #plt.show()
    #plt.plot(np.linspace(0,b*len(P),len(P)),P,label=str(x))
    #plt.show()
           
    return mean, std

for i in range(0,5):
    if i>0:
        plt.title('Varying the probability of gene recessivity')
        plt.legend()
        plt.show()
    for j in range(0,6):
        H=[0.5,0,0.5]
        mean, std=many_seasons_meanVar(H=H,n=6,b=0.05,p=0.4,x=0.5,num_seasons=1000,beta=1+j/15,epsilon=0,alpha=None)
        

''''
Mean=[]
Std=[]

for i in range(1,10):  
    
    beta=i
    mean, std=many_seasons_meanVar(P,n=len(P),b=0.05,p=0.3,x=0.48,num_seasons=50,beta=beta)
    Mean.append(mean)
    Std.append(std)
            
    
plt.subplot(211)
plt.title('Mean')
plt.plot(np.linspace(0,len(Mean),len(Mean)),Mean)

plt.subplot(212)
plt.title('Standard Deviation')
plt.plot(np.linspace(0,len(Std),len(Std)),Std)
            
    
    
plt.show()
    
    
    
    #plt.show()
    
    

   
    
#many_seasons(0.9,P)
'''      

    