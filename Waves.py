# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 23:33:59 2019

@author: Admin
"""


import  glob, os  
import numpy as np
import matplotlib.pyplot as plt
import math
import time
    
def solveRS(dx,dt,L,a,f,psi,fi,mu1,mu2,T, Action = None, solution = None):
    gamma = a*dt/dx
    M = int(round(T/dt))
    C2 = gamma**2
    N = int(L/dx)
    
    y = np.zeros(N + 1)
    y_1 = np.zeros(N + 1)
    y_2 = np.zeros(N + 1)
    x = np.linspace(0, L, N+1)
    t = np.linspace(0, M*dt, M+1)
    E = 0
    
    
    for i in range(N+1):    #нулевой слой
        y_1[i] = psi(x[i])
    
    if Action is not None:
        Action(y_1,x,t,0)
        
    if solution !=  None:
        u_e = solution(x, t[0])
        E  = max (np.abs(y_1 - u_e).max(), E)
        
    for i in range(1,N):    #первый слой
        y[i] = y_1[i] +dt*fi(x[i]) + 0.5*C2*(y_1[i+1] - 2*y_1[i] + y_1[i-1])+dt*dt/2*f(x[i],t[0])
    y[0] = y_1[0] +dt*fi(x[0]) + C2*(y_1[1] - y_1[0] - dx*mu1(t[0]))+dt*dt/2*f(x[0],t[0])
    y[N] = y_1[N] +dt*fi(x[N]) + C2*(y_1[N-1] - y_1[N] + dx*mu2(t[0]))+dt*dt/2*f(x[N],t[0])
    
    if Action is not None:
        Action(y,x,t,1)
     
    if solution != None:
        u_e = solution(x, t[1])
        E = max (np.abs(y - u_e).max(), E)   
    y_2[:], y_1[:] = y_1, y     #переносим слои
   
    for n in range(1,M):
	# Пересчитываем значения во внутренних узлах сетки на слое n+1
        for i in range(1, N):
            y[i] = 2*y_1[i] - y_2[i] + C2*(y_1[i+1] - 2*y_1[i] + y_1[i-1])+dt*dt*f(x[i],t[n])
        y[0] = 2*y_1[0] - y_2[0] + 2*C2*(y_1[1] - y_1[0] - dx*mu1(t[n]))+dt*dt*f(x[0],t[n])
        y[N] = 2*y_1[N] - y_2[N] + 2*C2*(y_1[N-1] - y_1[N] + dx*mu2(t[n]))+dt*dt*f(x[N],t[n])
        
        if Action is not None:
            if  Action(y,x,t,n+1):
                break
        if solution != None:
            u_e = solution(x, t[n+1])

            E = max (np.abs(y - u_e).max(), E)

    # Изменяем переменные перед переходом на следующий
    # временной слой
        y_2[:], y_1[:] = y_1, y
    
    

    return y, x, t, E


def viz(
    dx,dt,L,a,f,psi,fi,mu1,mu2,T,   # Параметры задачи
    umin, umax,              # Интервал для отображения u
    animate=True,             # Расчет с анимацией
    solution = None
    ):


    class PlotMatplotlib:
        def __call__(self, u, x, t, n):
            """Функция user_action для солвера."""
            if n == 0:
                plt.ion()
                self.lines = plt.plot(x, u, 'r-')
                plt.xlabel('x');  plt.ylabel('u')
                plt.axis([0, L, umin, umax])
                plt.legend(['t=%f' % t[n]], loc='lower left')
            else:
                self.lines[0].set_ydata(u)
                plt.legend(['t=%f' % t[n]], loc='lower left')
                plt.draw()
            if t[n] == 0:
                time.sleep(2) 
            else:
                time.sleep(0.2)
            plt.savefig('frame_%04d.png' % n)  # для генерации видео



    # Удаляем старые кадры
    for filename in glob.glob('frame_*.png'):
       os.remove(filename)

   #  Вызываем солвер и выполняем расчет
    plot = PlotMatplotlib()
    if animate:
        Action = plot
    else: 
        Action = None
    u, x, t, E = solveRS(
        dx,dt,L,a,f,psi,fi,mu1,mu2,T, Action,solution)

    # Генерируем видео файлы
    
    cmd = 'ffmpeg -r 24 -i frame_%04d.png -c:v libx264 movie_5.mp4'
    os.system(cmd)


    return E

def test1():
    
    gamma = 0.25
    
    L = math.pi
    a = 1
    N = 40  # Используем грубую сетку
    dx = L/N
    dt = dx*gamma/a
    T = 7

    def u_exact(x, t):
        return 1.5*np.sin(x+t)-0.5*np.cos(x+t)-1.5*np.sin(x-t)-0.5*np.cos(x-t)+np.cos(x)*((np.e)**(-3*t))

    def psi(x):
        return 0

    def fi(x):
        return 0
    
    def mu1(t):
        return 0
    
    def mu2(t):
        return 0

    def f(x, t):
        return 10*math.e**(-3*t)*math.cos(x)
    
    umax = 5
    umin = -umax
    E = viz(dx,dt,L,a,f,psi,fi,mu1,mu2,T, umin, umax,
              animate=True, solution = u_exact)
    return E, dt, dx

def test2():
    L = math.pi
    gamma = 1
    a = 1
    N = 160 
    dx = L/N
    dt = dx*gamma/a
    T = 4

    def u_exact(x, t):
        return -1/3*np.cos(4*(x+t))-1/3*np.cos(4*(x-t))+2/3*np.cos(2*t)*np.cos(4*x)

    def psi(x):
        return 0

    def fi(x):
        return 0
    
    def mu1(x):
        return 0
    
    def mu2(x):
        return 0

    def f(x, t):
        return 8*math.cos(2*t)*math.cos(4*x)
    
    umax = 5
    umin = -umax
    E = viz(dx,dt,L,a,f,psi,fi,mu1,mu2,T, umin, umax,
              animate=True, solution = u_exact)
    return E, dt, dx

def text_impuls_1():
    L = 2
    gamma = 1
    a = 1
    N = 100 
    dx = L/N
    dt = dx*gamma/a
    T = 10


    def psi(x):
        return 0

    def fi(x):
        return 0
    
    def mu1(x):
        return 0
    
    def mu2(x):
        return 0
    
    eps = 0.2
    ksi = 1/2
    def f(x, t):
        return g(x,ksi,eps)
    def g(x,ksi,eps):
        if ksi - eps <= x and x <= ksi + eps:
            return 2*eps
        else:
            return 0
        
    umax = 5
    umin = -umax
    E = viz(dx,dt,L,a,f,psi,fi,mu1,mu2,T, umin, umax,
              animate=True)
    return E, dt, dx

def text_impuls_2():
    L = 2
    gamma = 1
    a = 1
    N = 50 
    dx = L/N
    dt = dx*gamma/a
    T = 7


    def psi(x):
        return 0

    def fi(x):
        return 0
    
    def mu1(x):
        return 0
    
    def mu2(x):
        return 0
    
    eps = 0.2
    ksi = 1
    def f(x, t):
        return g(x,ksi,eps)
    def g(x,ksi,eps):
        if ksi - eps < x and x < ksi + eps:
            return 2*eps
        else:
            return 0
        
    umax = 5
    umin = -umax
    E = viz(dx,dt,L,a,f,psi,fi,mu1,mu2,T, umin, umax,
              animate=True)
    return E, dt, dx

def text_impuls_3():
    L = 2
    gamma = 1
    a = 1
    N = 50 
    dx = L/N
    dt = dx*gamma/a
    T = 5


    def psi(x):
        return 0

    def fi(x):
        return 0
    
    def mu1(x):
        return 0
    
    def mu2(x):
        return 0
    
    eps = 0.5
    ksi = 5/3
    def f(x, t):
        return g(x,ksi,eps)
    def g(x,ksi,eps):
        if ksi - eps < x and x < ksi + eps:
            return 2*eps
        else:
            return 0
        
    umax = 8
    umin = -1
    E = viz(dx,dt,L,a,f,psi,fi,mu1,mu2,T, umin, umax,
              animate=True)
    return E, dt, dx

def text_impuls_4():
    L = 2
    gamma = 0.5
    a = 1
    N = 50 
    dx = L/N
    dt = dx*gamma/a
    T = 10


    def psi(x):
        return 0

    def fi(x):
        return 0
    
    def mu1(x):
        return 0
    
    def mu2(x):
        return 0
    
    eps = 0.02
    ksi = 2/3
    def f(x, t):
        return g(x,ksi,eps)
    def g(x,ksi,eps):
        if ksi - eps < x and x < ksi + eps:
            return 2*eps
        else:
            return 0
        
    umax = 5
    umin = -umax
    E = viz(dx,dt,L,a,f,psi,fi,mu1,mu2,T, umin, umax,
              animate=True)
    return E, dt, dx
'''
def test_3():

    L = math.pi
    a = 1
    N = 40  # Используем грубую сетку
    dx = L/N
    dt = (dx)/(4*a)
    T = 3

    def u_exact(x, t):
        u_e = 1.5*np.sin(x+t)-0.5*np.cos(x+t)-1.5*np.sin(x-t)-0.5*np.cos(x-t)+np.cos(x)*((np.e)**(-3*t))

        return u_e

    def psi(x):
        return 0

    def fi(x):
        return 0
    
    def mu1(t):
        return 0
    
    def mu2(t):
        return 0

    def f(x, t):
        return 10*(math.e**(-3*t))*math.cos(x)


    
    def find_pogr(u, x, t, n):
        
        u_e = u_exact(x, t[n])
        diff = 0
        for i in range(len(u)):
            diff += (u[i] - u_e[i])*(u[i] - u_e[i])
            
        u_e = u_exact(x, t[n])
        diff = np.abs(u - u_e).max()
        Eps = 1E-13
        #print(diff)
       
        
    
        

    y,x,t,sq = solveRS(dx,dt,L,a,f,psi,fi,mu1,mu2,T, Action = None, solution = u_exact)
    #print(y)
    #u_e = u_exact(x, t[int(T/dt)])
    #print(np.abs(y-u_e).max())
    return sq,dt,dx
'''
E, dt, dx = text_impuls_3()
#print('Погрешность = ',E)
print('Шаг по времени = ', dt)
print('Шаг по пространству = ', dx)
