#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 17:58:46 2018

@author: tylerclark
"""
from scipy.special import factorial as fact
import numpy as np
import matplotlib.pyplot as plt

def game(K,b,p,player,beta,alpha):

    '''
    This specifies the results of a two-player game
    '''
        
    if player ==2:
        return p+K*beta*b
    if player ==1:
        return p-K*b

def Lambda(l,k,n,m,Print=False):
    '''
    '''
    p=choose(n-m,k-l)*choose(m,l)*fact(k)*fact(n-k)/(fact(n))
    if Print ==True:
        print(p,l,k,n,m)
    return p


def choose(n,k):
    '''
    This function returns the combinations of choosing n out of k
    '''
    return fact(n)/(fact(n-k)*fact(k))

def binom(n,p,k):
    '''
    The probability calculation taken from the binomial distribution
    '''
    return choose(n,k)*p**k*(1-p)**(n-k)  




def generate_distribution(H,n):
    '''
    This function takes in the number of genes n and the frequencies of the 
    heterogeneous and homogeneous genes and outputs a distribution of 
    macrostates
    '''
    P=[]
    for i in range(0,n):
        P.append(binom(n-1,H[0]+0.5*H[2],i))
        
    return P


def Lambda_dom(l,k,n,m,p,other_parent,epsilon,mark):
    '''
    '''
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
     
    if epsilon != None:
        Lam=Lam*(1+epsilon/(np.abs(other_parent-l)+0.01))

    return Lam,H

def Season(P,H,n,b,p,x,epsilon,beta,alpha):
    '''
    '''
    P=generate_distribution(H,n)
   
    
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




def many_seasons_meanVar(parms):
    '''
    '''
    
    num_seasons=parms['num_seasons']
    beta=parms['beta']
    epsilon=parms['epsilon']
    alpha=parms['alpha']
    x=parms['x']
    p=parms['p']
    b=parms['b']
    n=parms['n']
    H=parms['H']
    
    plot_distribution_season=parms['plot_dist_season']
    plot_distribution_end=parms['plot_dist_end']
    
    plot_mean_var=parms['plot_mean_var']
     
    P=generate_distribution(H,n)
    
    mean_all=[]
    std_all=[]
    
    for i in range(0,num_seasons):
        
        if parms['rand_beta'][0] ==True:
            
            beta=np.random.randn()*parms['rand_beta'][1]+beta

        P,H=Season(P=P,H=H,n=n,b=b,p=p,x=x,epsilon=epsilon,beta=beta,alpha=alpha)
        
        if plot_distribution_season ==True:
        
            plt.title('Distribution of Macrostates')
            plt.plot(np.linspace(0,b*len(P),len(P)),P,label='Season= '+str(i))
            
        
        mean=np.sum(P*np.linspace(0,len(P),len(P)))
        
        mean_all.append(mean)
        std=np.std(P)
        std_all.append(std)
        
    if plot_distribution_season ==True:
        plt.show()
                
    if plot_mean_var ==True:
        
            mean_all=np.asarray(mean_all)
            plt.plot(np.linspace(0,len(mean_all),len(mean_all)),mean_all,label='x='+str(x)+'--beta:'+str(beta))
        
        
     
    if plot_distribution_end ==True:

        plt.title('Macrostate distribution after '+str(num_seasons))
        plt.plot(np.linspace(0,b*len(P),len(P)),P,label=str(x))
        
           
    return mean, std




def explore(parms,var1,var2):
    
    for i in range(0,var1['steps']):
        
        if i>0 and parms['plot_dist_end']==True or parms['plot_mean_var']==True:
            plt.legend()
            plt.show()
            
        for j in range(0,var2['steps']):
            parms[var1['var']]=i*var1['max']/var1['steps']+var1['min']
            parms[var2['var']]=j*var2['max']/var2['steps']+var2['min']
            mean, std=many_seasons_meanVar(parms)
            
            
            