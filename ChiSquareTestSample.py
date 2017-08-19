import scipy.stats
import matplotlib.pyplot as plt
fe=[6.97,17.75,21.04,18.41,13.82,9.45,6.23,3.51,1.97,0.85]
fo1=[7.30,20.45,23.11,18.40,13.05,8.44,5.07,2.58,1.20,0.40]
fo2=[7.0,18,21.11,18.40,13.25,9.44,6.07,3.58,1.8,0.70]
a=scipy.stats.chisquare(fo1, fe, ddof=0, axis=0)
b=scipy.stats.chisquare(fo2, fe, ddof=0, axis=0)
print a[1]
print b[1]
if a[1]>b[1]:
    print "Curve one fits better"
else:
    print "Curve two fits better"
plt.plot(fe,'ro')
plt.plot(fo1,'go')
plt.plot(fo2,'bo')
plt.show()
