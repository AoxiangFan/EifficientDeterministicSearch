import numpy as np

def normalize2dpts(X):
	finiteind = np.where(np.abs(X[2,:]) > 1e-10)[0]
	X[0,finiteind] = X[0,finiteind]/X[2,finiteind]
	X[1,finiteind] = X[1,finiteind]/X[2,finiteind]
	X[2,finiteind] = 1
	c = np.mean(X[0:2,finiteind], axis=1).reshape((2,1))
	Xt = X[0:2,:] - c
	meandist = np.mean(np.sqrt(np.sum(Xt[0:2,:]**2,axis=0)))
	scale = np.sqrt(2)/meandist
	T = np.array([[scale, 0, -scale*c[0][0]],[0, scale, -scale*c[1][0]],[0, 0, 1]])
	Xn = np.dot(T,X)

	return Xn, T

def norm4Point(X, Y):
	N = X.shape[1]
	X, T1 = normalize2dpts(X)
	Y, T2 = normalize2dpts(Y)

	X[np.isnan(X)] = 1
	Y[np.isnan(Y)] = 1

	M = np.zeros((2*N,9))
	ooo = np.zeros((3))
	for i in range(N):
		p1 = X[:,i]
		p2 = Y[:,i]
		M[2*i,:] = np.hstack((p1.T*p2[2],ooo,-p1.T*p2[0]))
		M[2*i+1,:] = np.hstack((ooo,p1.T*p2[2],-p1.T*p2[1]))

	_, _, Vm = np.linalg.svd(M)
	Vm = Vm.T
	H = np.reshape(Vm[:,-1], (3,3))

	# evals, evecs = np.linalg.eig(np.dot(M.T,M))
	# sorted_indices = np.argsort(evals)
	# v = evecs[:, sorted_indices[0]]
	# H = np.reshape(v, (3,3))

	H = np.linalg.inv(T2).dot(H).dot(T1)

	return H

def weightedNorm4Point(X, Y, w):
	N = X.shape[1]
	X, T1 = normalize2dpts(X)
	Y, T2 = normalize2dpts(Y)

	X[np.isnan(X)] = 1
	Y[np.isnan(Y)] = 1

	M = np.zeros((2*N,9))
	ooo = np.zeros((3))
	for i in range(N):
		p1 = X[:,i]
		p2 = Y[:,i]
		M[2*i,:] = np.hstack((p1.T*p2[2],ooo,-p1.T*p2[0]))
		M[2*i+1,:] = np.hstack((ooo,p1.T*p2[2],-p1.T*p2[1]))

	wT = np.kron(np.diag(np.sqrt(w)),[[1,0],[0,1]])
	M = np.dot(wT, M)

	_, _, Vm = np.linalg.svd(M)
	Vm = Vm.T
	H = np.reshape(Vm[:,-1], (3,3))

	# evals, evecs = np.linalg.eig(np.dot(M.T,M))
	# sorted_indices = np.argsort(evals)
	# v = evecs[:, sorted_indices[0]]
	# H = np.reshape(v, (3,3))

	H = np.linalg.inv(T2).dot(H).dot(T1)

	return H

def norm7Point(X, Y):
	N = X.shape[1]

	X, T1 = normalize2dpts(X)
	Y, T2 = normalize2dpts(Y)

	X[np.isnan(X)] = 1
	Y[np.isnan(Y)] = 1

	M = np.zeros((N,9))
	for i in range(N):
		M[i,:] = [X[0,i]*Y[0,i], X[1,i]*Y[0,i], Y[0,i], X[0,i]*Y[1,i], X[1,i]*Y[1,i], Y[1,i], X[0,i], X[1,i], 1]

	_, S, V = np.linalg.svd(M)
	V = V.T
	F = np.reshape(V[:,-1], (3,3))

	f11=F[0,0]; f12=F[0,1]; f13=F[0,2]; f21=F[1,0]; f22=F[1,1]; f23=F[1,2]; f31=F[2,0]; f32=F[2,1]; f33=F[2,2]
	G = np.reshape(V[:,-2], (3,3))

	# evals, evecs = np.linalg.eig(np.dot(M.T,M))
	# sorted_indices = np.argsort(evals)
	# V = evecs[:, sorted_indices[0:2]]

	# F = np.reshape(V[:,0], (3,3))
	# f11=F[0,0]; f12=F[0,1]; f13=F[0,2]; f21=F[1,0]; f22=F[1,1]; f23=F[1,2]; f31=F[2,0]; f32=F[2,1]; f33=F[2,2]
	# G = np.reshape(V[:,1], (3,3))

	g11=G[0,0]; g12=G[0,1]; g13=G[0,2]; g21=G[1,0]; g22=G[1,1]; g23=G[1,2]; g31=G[2,0]; g32=G[2,1]; g33=G[2,2]
	a = (g11*g22*g33 - g11*g23*g32 - g12*g21*g33 + g12*g23*g31 + g13*g21*g32 - g13*g22*g31)
	b = (f11*g22*g33 - f11*g23*g32 - f12*g21*g33 + f12*g23*g31 + f13*g21*g32 - f13*g22*g31 - f21*g12*g33 + f21*g13*g32 + f22*g11*g33 - f22*g13*g31
    - f23*g11*g32 + f23*g12*g31 + f31*g12*g23 - f31*g13*g22 - f32*g11*g23 + f32*g13*g21 + f33*g11*g22 - f33*g12*g21)
	c = (f11*f22*g33 - f11*f23*g32 - f11*f32*g23 + f11*f33*g22 - f12*f21*g33 + f12*f23*g31 + f12*f31*g23 - f12*f33*g21 + f13*f21*g32 - f13*f22*g31
    - f13*f31*g22 + f13*f32*g21 + f21*f32*g13 - f21*f33*g12 - f22*f31*g13 + f22*f33*g11 + f23*f31*g12 - f23*f32*g11)
	d = f11*f22*f33 - f11*f23*f32 - f12*f21*f33 + f12*f23*f31 + f13*f21*f32 - f13*f22*f31

	s0 = np.roots([a, b, c, d])
	s = s0[np.isreal(s0)]
	s = s.real

	Fm = np.zeros((3,3,len(s)))
	for i in range(len(s)):
		Fm[:,:,i] = T2.T.dot(F+s[i]*G).dot(T1)

	return Fm

def norm8Point(X, Y):
	N = X.shape[1]
	X, T1 = normalize2dpts(X)
	Y, T2 = normalize2dpts(Y)

	X[np.isnan(X)] = 1
	Y[np.isnan(Y)] = 1

	M = np.zeros((N,9))
	for i in range(N):
		M[i,:] = [X[0,i]*Y[0,i], X[1,i]*Y[0,i], Y[0,i], X[0,i]*Y[1,i], X[1,i]*Y[1,i], Y[1,i], X[0,i], X[1,i], 1]

	_, _, Vm = np.linalg.svd(M)
	Vm = Vm.T
	F = np.reshape(Vm[:,-1], (3,3))

	# evals, evecs = np.linalg.eig(np.dot(M.T,M))
	# sorted_indices = np.argsort(evals)
	# v = evecs[:, sorted_indices[0]]
	# F = np.reshape(v, (3,3))

	U, S, V = np.linalg.svd(F)
	S[-1] = 0
	F = U.dot(np.diag(S)).dot(V)

	F = F/np.linalg.norm(F.flatten())
	F = T2.T.dot(F).dot(T1)

	return F


def weightedNorm8Point(X, Y, w):
	N = X.shape[1]
	X, T1 = normalize2dpts(X)
	Y, T2 = normalize2dpts(Y)

	X[np.isnan(X)] = 1
	Y[np.isnan(Y)] = 1

	M = np.zeros((N,9))
	for i in range(N):
		M[i,:] = [X[0,i]*Y[0,i], X[1,i]*Y[0,i], Y[0,i], X[0,i]*Y[1,i], X[1,i]*Y[1,i], Y[1,i], X[0,i], X[1,i], 1]

	M = np.dot(np.diag(np.sqrt(w)), M)

	_, _, Vm = np.linalg.svd(M)
	Vm = Vm.T
	F = np.reshape(Vm[:,-1], (3,3))

	# evals, evecs = np.linalg.eig(np.dot(M.T,M))
	# sorted_indices = np.argsort(evals)
	# v = evecs[:, sorted_indices[0]]
	# F = np.reshape(v, (3,3))

	U, S, V = np.linalg.svd(F)
	S[-1] = 0
	F = U.dot(np.diag(S)).dot(V)

	F = F/np.linalg.norm(F.flatten())
	F = T2.T.dot(F).dot(T1)

	return F


def sampsonDistanceF(X, Y, F):
	pfp = np.dot(Y.T, F).T
	pfp = pfp*X
	d = sum(pfp,0)**2
	epl1 = np.dot(F, X)
	epl2 = np.dot(F.T, Y)

	d = d/(epl1[0,:]**2 + epl1[1,:]**2 + epl2[0,:]**2 + epl2[1,:]**2)
	d = np.sqrt(d)
	
	return d



def sampsonDistanceH(X, Y, H):
	N = X.shape[1]
	D1 = np.zeros((N,9))
	D2 = np.zeros((N,9))
	ooo = np.zeros((3))
	for i in range(N):
		p1 = X[:,i]
		p2 = Y[:,i]
		D1[i,:] = np.hstack((p1.T*p2[2],ooo,-p1.T*p2[0]))
		D2[i,:] = np.hstack((ooo,p1.T*p2[2],-p1.T*p2[1]))	
	h = np.array([[H[0,0],H[0,1],H[0,2],H[1,0],H[1,1],H[1,2],H[2,0],H[2,1],H[2,2]]]).T
	r1 = np.dot(D1, h)
	r2 = np.dot(D2, h)

	algError = np.vstack((r1.T, r2.T))

	G1 = np.vstack((H[0,0] - Y[0,:]*H[2,0], H[0,1] - Y[0,:]*H[2,1], -X[0,:]*H[2,0] - X[1,:]*H[2,1] - H[2,2], np.zeros((1,N))))
	G2 = np.vstack((H[1,0] - Y[1,:]*H[2,0], H[1,1] - Y[1,:]*H[2,1], np.zeros((1,N)), -X[0,:]*H[2,0] - X[1,:]*H[2,1] - H[2,2]))

	magG1 = np.sqrt(np.sum(G1**2,axis=0))
	magG2 = np.sqrt(np.sum(G2**2,axis=0))
	magG1G2 = np.sum(G1*G2,axis=0)

	alpha = np.arccos(magG1G2/(magG1*magG2))
	D1 = algError[0,:]/magG1
	D2 = algError[1,:]/magG2

	d = (D1*D1 + D2*D2 - 2*D1*D2*np.cos(alpha))/np.sin(alpha)
	d = np.sqrt(d)
	
	return d

def preSampsonDistanceH_all(X, Y):
	global D1, D2
	N = X.shape[1]
	D1 = np.zeros((N,9))
	D2 = np.zeros((N,9))
	ooo = np.zeros((3))
	for i in range(N):
		p1 = X[:,i]
		p2 = Y[:,i]
		D1[i,:] = np.hstack((p1.T*p2[2],ooo,-p1.T*p2[0]))
		D2[i,:] = np.hstack((ooo,p1.T*p2[2],-p1.T*p2[1]))
	return D1, D2

def sampsonDistanceH_all(X, Y, H):
	global D1, D2
	N = X.shape[1]
	h = np.array([[H[0,0],H[0,1],H[0,2],H[1,0],H[1,1],H[1,2],H[2,0],H[2,1],H[2,2]]]).T

	r1 = np.dot(D1, h)
	r2 = np.dot(D2, h)

	algError = np.vstack((r1.T, r2.T))

	G1 = np.vstack((H[0,0] - Y[0,:]*H[2,0], H[0,1] - Y[0,:]*H[2,1], -X[0,:]*H[2,0] - X[1,:]*H[2,1] - H[2,2], np.zeros((1,N))))
	G2 = np.vstack((H[1,0] - Y[1,:]*H[2,0], H[1,1] - Y[1,:]*H[2,1], np.zeros((1,N)), -X[0,:]*H[2,0] - X[1,:]*H[2,1] - H[2,2]))

	magG1 = np.sqrt(np.sum(G1**2,axis=0))
	magG2 = np.sqrt(np.sum(G2**2,axis=0))
	magG1G2 = np.sum(G1*G2,axis=0)

	alpha = np.arccos(magG1G2/(magG1*magG2))
	D1t = algError[0,:]/magG1
	D2t = algError[1,:]/magG2

	d = (D1t*D1t + D2t*D2t - 2*D1t*D2t*np.cos(alpha))/np.sin(alpha)
	d = np.sqrt(d)
	
	return d


def minimalSampleF(X, Y, threshold):

	indices = np.random.choice(range(X.shape[1]), size=7, replace=False)
	F = norm7Point(X[:,indices],Y[:,indices])
	bestF = []
	bestS = -1
	bestInliers = []
	for i in range(F.shape[2]):
		d = sampsonDistanceF(X, Y, F[:,:,i])
		inliers = np.where(d<=threshold)[0]
		if len(inliers) > bestS:
			bestS = len(inliers)
			bestF = F[:,:,i]
			bestInliers = inliers
	return bestF, indices, bestInliers


# def minimalSampleF(X, Y):
# 	indices = np.random.choice(range(X.shape[1]), size=7, replace=False)
# 	F = norm7Point(X[:,indices],Y[:,indices])
	
# 	return F, indices

def minimalSampleH(X, Y):
	indices = np.random.choice(range(X.shape[1]), size=4, replace=False)
	H = norm4Point(X[:,indices],Y[:,indices])
	
	return H, indices


def updateMaxTrials(numCurInliers, maxTrials, N, confidence, order):

    ratioOfInliers = numCurInliers/N
    if ratioOfInliers > 1 - 1e-16:
        newNum = 1
    else:
        ratio_t = ratioOfInliers**order
        if ratio_t > 1e-16:
            logOneMinusRatio_t = np.log(1 - ratio_t)
            newNum = np.ceil(np.log(1-confidence) / logOneMinusRatio_t)
        else:
            newNum = 1e16

    if maxTrials > newNum:
        maxTrials = newNum
    return maxTrials
