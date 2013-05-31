import numpy as np

f = file('tscores.npz','rb')

T = np.load(f)

print T
nrows, nlv = np.shape(T)


###################################33


cent = np.empty(nlv,dtype=np.float64)
for a in range (nlv):
    cent[a]=np.mean(T[:,a])
    
# compute distances to X centroid and percentil 95

dcent = np.empty(nrows, dtype=np.float64)
for i in range(nrows):
    dcent[i] = np.sqrt(np.sum(np.square(cent-T[i,:])))

dcent = np.sort(dcent)
p95dcent = dcent[np.ceil((nrows*95)/100)-1]

print dcent, p95dcent

# compute mutual distances in X and percentil 95
dmut = np.empty((nrows*nrows-nrows)/2, dtype=np.float64)
k=0
for i in range (nrows):
    for j in range (nrows):
        if j>i:
            dmut[k]= np.sqrt(np.sum(np.square(T[j,:]-T[i,:])))
            k+=1
            
dmut = np.sort(dmut)
p95dmut = dmut[np.ceil((nrows*95)/100)-1]

print dmut, p95dmut

# compute percentil 96 dmodx

# compute distance to Y centroid and percentil 95

# compute mutual distance in Y and percentil 95

# write in a file, Am -> critical distances -> centroid -> t

