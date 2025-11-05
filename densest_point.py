IC = open("nfw-haloIC.txt", "r").readlines()[2:]
import numpy as np
import matplotlib.pyplot as plt

mass = np.array([float(line.split()[1]) for line in IC])*1e10
x = np.array([float(line.split()[2]) for line in IC])
y = np.array([float(line.split()[3]) for line in IC])
z = np.array([float(line.split()[4]) for line in IC])

CMx = np.sum(mass * x) / np.sum(mass)
CMy = np.sum(mass * y) / np.sum(mass)
CMz = np.sum(mass * z) / np.sum(mass)

r = np.sqrt(x**2 + y**2 + z**2)
radius = np.max(r)
while radius > 0.001:
    x = x - CMx
    y = y - CMy
    z = z - CMz
    r = np.sqrt(x**2 + y**2 + z**2)
    
    keep = r < radius
    mass = mass[keep]
    x = x[keep]
    y = y[keep]
    z = z[keep]
    r = r[keep]

    CMx = np.sum(mass * x) / np.sum(mass)
    CMy = np.sum(mass * y) / np.sum(mass)
    CMz = np.sum(mass * z) / np.sum(mass)
    if len(r) <10:
        break
    radius = radius*0.9

print (CMx, CMy, CMz)

mass0 = np.array([float(line.split()[1]) for line in IC])*1e10
x0 = np.array([float(line.split()[2]) for line in IC])
y0 = np.array([float(line.split()[3]) for line in IC])
z0 = np.array([float(line.split()[4]) for line in IC])

x1 = x0 - CMx
y1 = y0 - CMy
z1 = z0 - CMz
r = np.sqrt(x1**2 + y1**2 + z1**2)
max = np.max(r)
min = np.min(r)
bins = 10
den = []
r_bins = np.logspace(np.log10(min), np.log10(max), bins)
for j, _ in enumerate(r_bins):
    if j == 0:
        keep = (r < r_bins[j]) & (r > 0)
        volume = (4/3) * np.pi * (r_bins[j])**3
    else:
        keep = (r < r_bins[j]) & (r > r_bins[j-1])
        volume1 = (4/3) * np.pi * (r_bins[j-1])**3
        volume2 = (4/3) * np.pi * (r_bins[j])**3
        volume = volume2 - volume1
    mass_in_r = np.sum(mass0[keep])
    print(keep)
    print(mass_in_r)
    den.append(mass_in_r / volume)

density = np.array(den)
print(density)
print(r_bins)
plt.scatter(r_bins, density, s=1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Distance (kpc)')
plt.ylabel('Density (kg/kpc^3)') 
plt.show()
 


