import numpy as np
import matplotlib.pyplot as plt
import string
x = np.linspace(0.1,2,100000)

#y=0.00005*((1+(x/(-2)))**(2))*((1-(x/5))**(-5))
#y=((1+(x**2)/4)**(-2))*(np.exp(0*np.arctan(x/2)))
#y=15*((x-1)**(-2))*(x**(-3))
#y=((1+(x**2)/4)**(3))
#y=np.exp(-(x**2)/2)
#y=10*(x**(-2))*np.exp(2/x)
#y=(1+x/(-4))**(5)
#y=(10)*np.exp(-x/2)
y=0.002*x**(-5)
plt.plot(x,y,'r')
#plt.text(1.7,180,u'y\u2080=0.002, m=5 ', ha='center', va='center')
plt.show()
