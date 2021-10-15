import numpy as np
import BASE_utils
import LO_utils
import DEGEN_utils

def findFundamentalMatrix(X, Y, embeddingType, th, maxIter):

    label = np.zeros((X.shape[0],1))
    _, u1 = np.unique(Y, axis=0, return_index=True)
    X = X[u1,:]
    Y = Y[u1,:]
    _, u2 = np.unique(X, axis=0, return_index=True)
    X = X[u2,:]
    Y = Y[u2,:]

    N = X.shape[0]

    if N < 8:
        F = np.ones((3,3))
        return F, label.astype(bool).flatten()

    if embeddingType == 'A':
        D = generateEmbeddings(X, Y, 'A')
        if np.isnan(D).any():
            idx = np.ones((N), dtype=bool)
            idx = np.where(idx)[0]
        else:
            d1 = DPCP_solver(D,2)
            idxa = (d1 <= 0.25)
            idxb = (idxa == False)
            if any(idxb):
                d2 = DPCP_solver(D[:,idxb],2)
                idxc = (d2 <= 0.15)
                idxb[idxb] = idxc
                idx = np.zeros((N), dtype=bool)
                idx[idxa] = True
                idx[idxb] = True
            else:
                idx = np.zeros((N), dtype=bool)
                idx[idxa] = True
        idx = np.where(idx)[0]
    elif embeddingType == 'H':
        D = generateEmbeddings(X, Y, 'H')
        if np.isnan(D).any():
            idx = np.ones((N), dtype=bool)
            idx = np.where(idx)[0]
        else:
            d1 = DPCP_solver(D,1)
            d1 = (d1[range(0,2*N,2)]+d1[range(1,2*N,2)])/2
            idxa = (d1 <= 0.25)
            idxb = (idxa == False)
            tmp = np.zeros((N))
            tmp[idxb] = 1
            idxb2 = np.kron(tmp,np.array([1,1])).astype('bool')
            if any(idxb):
                d2 = DPCP_solver(D[:,idxb2],2)
                d2 = (d2[range(0,2*sum(tmp).astype('int'),2)]+d2[range(1,2*sum(tmp).astype('int'),2)])/2
                idxc = (d2 <= 0.15)
                idxb[idxb] = idxc
                idx = np.zeros((N), dtype=bool)
                idx[idxa] = True
                idx[idxb] = True
            else:
                idx = np.zeros((N), dtype=bool)
                idx[idxa] = True
        idx = np.where(idx)[0]
    elif embeddingType == 'F':
        D = generateEmbeddings(X, Y, 'F')
        if np.isnan(D).any():
            idx = np.ones((N), dtype=bool)
            idx = np.where(idx)[0]
        else:
            d1 = DPCP_solver(D,1)
            idx = (d1 <= 0.15)
            idx = np.where(idx)[0]
    
    if len(idx) < 8:
        idx = np.array(range(N))
    
    DEGEN = True
    LO = True

    F, _ = postFundamentalMatrix(X[idx,:], Y[idx,:], maxIter, th, DEGEN, LO)

    N = X.shape[0]
    Xt = np.hstack((X,np.ones((N,1)))).T
    Yt = np.hstack((Y,np.ones((N,1)))).T
    d = BASE_utils.sampsonDistanceF(Xt, Yt, F)
    inliers = np.where(d<=th)[0]

    label[u1[u2[inliers]]] = 1

    return F, label.astype(bool).flatten()

def findHomography(X, Y, embeddingType, th, maxIter):

    label = np.zeros((X.shape[0],1))
    _, u1 = np.unique(Y, axis=0, return_index=True)
    X = X[u1,:]
    Y = Y[u1,:]
    _, u2 = np.unique(X, axis=0, return_index=True)
    X = X[u2,:]
    Y = Y[u2,:]

    N = X.shape[0]

    if N < 4:
        H = np.ones((3,3))
        return H, label.astype(bool).flatten()

    if embeddingType == 'A':
        D = generateEmbeddings(X, Y, 'A')
        if np.isnan(D).any():
            idx = np.ones((N), dtype=bool)
            idx = np.where(idx)[0]
        else:
            d1 = DPCP_solver(D,2)
            idx = (d1 <= 0.15)
        idx = np.where(idx)[0]
    elif embeddingType == 'H':
        D = generateEmbeddings(X, Y, 'H')
        if np.isnan(D).any():
            idx = np.ones((N), dtype=bool)
            idx = np.where(idx)[0]
        else:
            d1 = DPCP_solver(D,1)
            d1 = (d1[range(0,2*N,2)]+d1[range(1,2*N,2)])/2
            idx = (d1 <= 0.15)
        idx = np.where(idx)[0]
    elif embeddingType == 'F':
        D = generateEmbeddings(X, Y, 'F')
        if np.isnan(D).any():
            idx = np.ones((N), dtype=bool)
            idx = np.where(idx)[0]
        else:
            d1 = DPCP_solver(D,3)
            idx = (d1 <= 0.15)
            idx = np.where(idx)[0]
    
    if len(idx) < 4:
        idx = np.array(range(N))
    
    LO = True
    H, _ = postHomography(X[idx,:], Y[idx,:], maxIter, th, LO)
    N = X.shape[0]
    Xt = np.hstack((X,np.ones((N,1)))).T
    Yt = np.hstack((Y,np.ones((N,1)))).T
    d = BASE_utils.sampsonDistanceH(Xt, Yt, H)
    inliers = np.where(d<=th)[0]

    label[u1[u2[inliers]]] = 1

    return H, label.astype(bool).flatten()


def generateEmbeddings(X, Y, embeddingType):
	X0 = X
	Y0 = Y
	N = X.shape[0]
	Xn = dataNorm(X0)
	Yn = dataNorm(Y0)
	Xt = np.vstack((Xn.T, np.ones((1,N))))
	Yt = np.vstack((Yn.T, np.ones((1,N))))
	if embeddingType == 'F':
		D = np.zeros((9,N))
		for i in range(N):
			D[:,i] = np.kron(Yt[:,i],Xt[:,i])
		D = D/np.sqrt(np.sum(D**2,axis=0))
	elif embeddingType == 'H':
		D = np.zeros((2*N,9))
		ooo = np.zeros((3))
		for i in range(N):
			p1 = Xt[:,i]
			p2 = Yt[:,i]
			D[2*i,:] = np.hstack((p1.T*p2[2],ooo,-p1.T*p2[0]))
			D[2*i+1,:] = np.hstack((ooo,p1.T*p2[2],-p1.T*p2[1]))
		D = D.T
		D = D/np.sqrt(np.sum(D**2,axis=0))
	elif embeddingType == 'At':
		D = np.zeros((2*N,9))
		ooo = np.zeros((3))
		for i in range(N):
			p1 = Xt[:,i]
			p2 = Yt[:,i]
			D[2*i,:] = np.hstack((p1.T,ooo,-p2[0]))
			D[2*i+1,:] = np.hstack((ooo,p1.T,-p2[1]))
		D = D.T
		D = D/np.sqrt(np.sum(D**2,axis=0))
	elif embeddingType == 'A':
		D = np.hstack((Xn,Yn,np.ones((N,1)))).T
		D = D/np.sqrt(np.sum(D**2,axis=0))

	return D


def DPCP_solver(Xtilde, c):
    mu_min = 1e-15 
    maxiter = 200

    mu_0 = 1e-2
    alpha = 1e-3
    beta = 1/2

    D, N = Xtilde.shape
	
    evals, evecs = np.linalg.eig(np.dot(Xtilde,Xtilde.T))
    sorted_indices = np.argsort(evals)
	
    B_0 = evecs[:, sorted_indices[0:c]]
	
    for j in range(0, c):
	    i = 0
	    b = B_0[:,j]
	    mu = mu_0
	    if j == 0:
		    obj_old = np.sum(np.abs(np.dot(Xtilde.T,b)))
		    while mu > mu_min and i <= maxiter:
			    i = i + 1
			    grad = np.sum(np.sign(np.dot(b.T,Xtilde))*Xtilde,axis=1)
			    grad_norm = np.square(np.linalg.norm(grad))
			    bk = b - mu*grad
			    while (np.sum(np.abs(np.dot(Xtilde.T,bk/np.linalg.norm(bk)))) > obj_old - alpha*mu*grad_norm) and mu > mu_min:
				    mu = mu*beta
				    bk = b - mu*grad
			    b = bk/np.linalg.norm(bk)
			    obj_old = np.sum(np.abs(np.dot(Xtilde.T,b)))
		    B = b.reshape((D,1))
	    else:
		    b = b - B.dot(B.T.dot(b))
		    b = b/np.linalg.norm(b)
		    obj_old = np.sum(np.abs(np.dot(Xtilde.T,b)))
		    while mu > mu_min and i < maxiter:
			    i = i + 1
			    grad = np.sum(np.sign(np.dot(b.T,Xtilde))*Xtilde,axis=1)
			    grad_norm = np.square(np.linalg.norm(grad))
			    bk = b - mu*grad
			    while (np.sum(np.abs(np.dot(Xtilde.T,bk/np.linalg.norm(bk)))) > obj_old - alpha*mu*grad_norm) and mu > mu_min:
				    mu = mu*beta
				    bk = b - mu*grad
			    bk = bk - B.dot(B.T.dot(bk))
			    b = bk/np.linalg.norm(bk)
			    obj_old = np.sum(np.abs(np.dot(Xtilde.T,b)))
		    B = np.concatenate((B,b.reshape((D,1))),axis=1)
	
    distances = np.zeros((N))
    for j in range(0,N):
	    distances[j] = np.linalg.norm(B.T.dot(Xtilde[:,j]))
	
    return distances

def dataNorm(X):
    Xn = (X - np.mean(X,axis=0))/np.std(X,axis=0)
    return Xn

def postFundamentalMatrix(X, Y, maxTrials, threshold, DEGEN, LO):
	N = X.shape[0]
	X = np.hstack((X,np.ones((N,1)))).T
	Y = np.hstack((Y,np.ones((N,1)))).T
	curTrials = 0
	bestInliers = []
	numBestInliers = 0

	numGoodInliers = 0

	BASE_utils.preSampsonDistanceH_all(X, Y)

	while curTrials <= maxTrials:
		F, indices, curInliers = BASE_utils.minimalSampleF(X, Y, threshold)
		numCurInliers = len(curInliers)
		if numGoodInliers < numCurInliers:
			numGoodInliers = numCurInliers
			# maxTrials = max(BASE_utils.updateMaxTrials(numGoodInliers, maxTrials, N, 0.999999, 7), 100)
			maxTrials = BASE_utils.updateMaxTrials(numGoodInliers, maxTrials, N, 0.999999, 7)
			if numBestInliers < numCurInliers:
				bestInliers = curInliers
				numBestInliers = numCurInliers
			if LO:
				if numCurInliers >= 16:
					F, curInliers = LO_utils.LO_F(X, Y, F, threshold, 3, 50)
					numCurInliers = len(curInliers)
					if numBestInliers < numCurInliers:
						bestInliers = curInliers
						numBestInliers = numCurInliers
			if DEGEN:
				flag, H = DEGEN_utils.checkSample(X[:,indices], Y[:,indices], F, 2*threshold)
				if flag:
					dH = BASE_utils.sampsonDistanceH_all(X, Y, H)
					inliersH = np.where(dH<=3*threshold)[0]
					if len(inliersH) >= 8:
						H, inliersH = LO_utils.LO_H(X, Y, H, threshold, 3, 50)
						if len(inliersH) >= 6:
							F, inliersF = DEGEN_utils.H2F(X, Y, H, threshold)
							if numBestInliers < len(inliersF):
								bestInliers = inliersF
								numBestInliers = len(inliersF)				
		curTrials = curTrials + 1
	if len(bestInliers) >= 8:
		F = BASE_utils.norm8Point(X[:, bestInliers], Y[:, bestInliers])
		d = BASE_utils.sampsonDistanceF(X, Y, F)
		bestInliers = np.where(d<=threshold)[0]
	else:
		F = np.ones((3,3))

	return F, bestInliers

def postHomography(X, Y, maxTrials, threshold, LO):
	N = X.shape[0]
	X = np.hstack((X,np.ones((N,1)))).T
	Y = np.hstack((Y,np.ones((N,1)))).T
	curTrials = 0
	bestInliers = []
	numBestInliers = 0

	numGoodInliers = 0

	BASE_utils.preSampsonDistanceH_all(X, Y)
	
	while curTrials <= maxTrials:
		H, _ = BASE_utils.minimalSampleH(X, Y)
		d = BASE_utils.sampsonDistanceH_all(X, Y, H)
		curInliers = np.where(d<=threshold)[0]
		numCurInliers = len(curInliers)
		if numGoodInliers < numCurInliers:
			numGoodInliers = numCurInliers
			# maxTrials = max(BASE_utils.updateMaxTrials(numGoodInliers, maxTrials, N, 0.999999, 7), 100)
			maxTrials = BASE_utils.updateMaxTrials(numGoodInliers, maxTrials, N, 0.999999, 4)
			if numBestInliers < numCurInliers:
				bestInliers = curInliers
				numBestInliers = numCurInliers
			if LO:
				if numCurInliers >= 8:
					H, curInliers = LO_utils.LO_H(X, Y, H, threshold, 3, 50)
					numCurInliers = len(curInliers)
					if numBestInliers < numCurInliers:
						bestInliers = curInliers
						numBestInliers = numCurInliers
		curTrials = curTrials + 1
		
	if len(bestInliers) >= 4:
		H = BASE_utils.norm4Point(X[:, bestInliers], Y[:, bestInliers])
		d = BASE_utils.sampsonDistanceH_all(X, Y, H)
		bestInliers = np.where(d<=threshold)[0]
	else:
		H = np.ones((3,3))

	return H, bestInliers