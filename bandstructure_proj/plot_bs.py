import numpy as np
import matplotlib.pyplot as plt

DOSpath = "./DOSCAR"
DOSfile = open(DOSpath, "r")
nline = 0
for line in DOSfile:
    words = line.split()
    nline += 1
    if nline == 6:
        Ef = float(words[3])
        break

Evalpath = "./EIGENVAL"
Efile = open(Evalpath, "r")
nline = 0
iread = False
eigens = []
for line in Efile:
    words = line.split()
    nline += 1
    if nline == 6:
        nbands = int(words[2])
        nkpoints = int(words[1])
    if nline == 7:
        iread = True
    if iread:
        if len(words) == 3:
            eigens.append(float(words[1]) - Ef)

eigens = np.array(eigens).reshape([nkpoints, nbands])
#for i in range(nkpoints):
#    print("  ".join([str(w) for w in eigens[i]]))

nk = eigens.shape[0]
nb = eigens.shape[1]
PROpath = "./PROCAR"
PROfile = open(PROpath,  "r")
nline = 0
iread = False
data = []
for line in PROfile:
    words = line.split()
    nline += 1
    if nline == 2:
        nion = int(words[11])
        norbital = 9
        data = np.zeros([nk, nb, nion, norbital])
        #print(data.shape)
        ik = 0
        ib = 0
        ii = 0
    if not words: continue
    if iread:
        if words[0] == "tot":
            iread = False
            ib += 1
            continue
        data[ik, ib, ii, :] = words[1:10] 
        ii += 1
    if  words[0] == "ion":
        iread = True
    if words[0] == "band":
        ii = 0
    if  words[0] == "k-point" and nline > 5:
        ib = 0
        ik += 1




#bands to show : 0-8(Cu_dx2-y2); 0-4(Cu_dxy); 0-5+0-7(Cu_xz+yz); 0-6(Cu_z2); ion-3-6-total 
#print(eigens.shape)
#print(data.shape)
k_points = np.linspace(0, 2.889, 120)
colors = {
        'blue': np.array([0, 0, 1]),
        'green': np.array([0, 1, 0]),
        'yellow': np.array([1, 1, 0]),
        'red': np.array([1, 0, 0]),
        'pink': np.array([1, 0.75, 0.8])
    }
colors = list(colors.values())

proj = [data[:, :, 0, 8], data[:, :, 0, 4], data[:, :, 0, 5], data[:, :, 0, 6], np.sum(np.sum(data[:, :, 3:7, :], axis=-1), axis=-1)]

proj = np.array(proj).reshape([5, nk, nb])
for ik in range(nk):
    for ib in range(nb):
        if np.linalg.norm(proj[:, ik, ib]) < 0.0001:
            continue
        proj[:, ik, ib] /= np.linalg.norm(proj[:, ik, ib])

fig, ax = plt.subplots()
for ista in range(proj.shape[0]):
    for ib in range(nb):
        band = eigens[:, ib]
        #print(colors[ista].repeat(nk).reshape(3, nk).T)
        color = np.column_stack((colors[ista].repeat(nk).reshape(3, nk).T, proj[ista, :, ib]))
        points = ax.scatter(k_points, band, c=color) 
# Adding colorbar
cbar = plt.colorbar(points, ax=ax)
cbar.set_label('Property value')

# Labels and title
ax.set_ylim(-3, 3)
ax.set_xlabel('Wave Vector k')
ax.set_ylabel('Energy (eV)')
ax.set_title('Band Structure')

plt.show()
    
