import numpy as np
import BASE_utils
import LO_utils

def skew_sym(a):
    ax = np.array([ [0,-a[2],a[1]], [a[2],0,-a[0]], [-a[1],a[0],0] ])
    return ax

def Hdetect(X7, Y7, F, IDXS):

    ec = np.cross(F[:,0],F[:,1])
    if (np.abs(ec[0]) < 1e-10) and (np.abs(ec[1]) < 1e-10) and (np.abs(ec[2]) < 1e-10):
        ec = np.cross(F[:,1],F[:,2])
    ec = ec/ec[2]
    Ex = skew_sym(ec)

    A = np.dot(np.linalg.inv(Ex), F)
    # A = np.dot(Ex, F)

    X3 = X7[:, IDXS]
    Y3 = Y7[:, IDXS]

    b = np.zeros((3,1))
    for i in range(3):
        p1 = np.cross(Y3[i,:], np.dot(A, X3[i,:]))
        p1 = p1/p1[2]
        p2 = np.cross(Y3[i,:], ec)
        p2 = p2/p2[2]
        p2 = p2/(np.linalg.norm(p2)**2)
        b[i] = np.dot(p1.T, p2)

    M = X3.T
    H = A - ec.reshape((3,1)).dot(np.dot(np.linalg.inv(M), b).T)

    if np.isnan(H[0,0]) or np.isinf(H[0,0]):
        H = np.eye(3)
    
    return H

def checkSample(X7, Y7, F, th):
    IDXS = np.array([[0,1,2],[3,4,5],[0,1,6],[3,4,6],[2,5,6]])
    for i in range(5):
        H = Hdetect(X7, Y7, F, IDXS[i,:])
        d = BASE_utils.sampsonDistanceH(X7, Y7, H)
        idx = np.argsort(d)
        H = BASE_utils.norm4Point(X7[:,idx[0:5]], Y7[:,idx[0:5]])
        d2 = BASE_utils.sampsonDistanceH(X7, Y7, H)
        inliers = np.where(d2<=th)[0]
        if len(inliers) > 4:
            return True, H
    
    return False, H

def innerFH(Xh, Yh, Xo, Yo, X, Y, th, num, sam_sizH, sam_sizO):
    lenH = Xh.shape[1]
    lenO = Xo.shape[1]
    max_i = 0
    max_s = 0

    bestF = []
    bestInliers = []

    for rep in range(num):
        hsam = np.random.choice(range(lenH), size=sam_sizH, replace=False)
        osam = np.random.choice(range(lenO), size=sam_sizO, replace=False)

        Xh_sam = Xh[:, hsam]
        Yh_sam = Yh[:, hsam]
        Xo_sam = Xo[:, osam]
        Yo_sam = Yo[:, osam]

        X_sam = np.hstack((Xh_sam,Xo_sam))
        Y_sam = np.hstack((Yh_sam,Yo_sam))

        aF = BASE_utils.norm8Point(X_sam, Y_sam)
        d = BASE_utils.sampsonDistanceF(X, Y, aF)
        inliers = np.where(d<=th)[0]

        nI = len(inliers)

        if max_i < nI:
            bestInliers = inliers
            bestF = aF
            max_i = nI

        if nI > max_s:
            max_s = nI
            aF, inliers = LO_utils.iterF(X, Y, aF, th, 2, 50)
            nI = len(inliers)

        if max_i < nI:
            bestInliers = inliers
            bestF = aF
            max_i = nI
    
    return bestF, bestInliers

def H2F(X, Y, H, th):
    MAX_SAM = 100
    sam_sizH = 6
    sam_sizO = 4

    hinl = np.where(BASE_utils.sampsonDistanceH_all(X, Y, H) < 4*th)[0]
    nhinl = np.where(BASE_utils.sampsonDistanceH_all(X, Y, H) > 5*th)[0]

    if (len(hinl) < sam_sizH) or (len(nhinl) < sam_sizO):
        F = np.ones((3,3))
        inliers = []
        return F, inliers

    Xo = X[:, nhinl]
    Yo = Y[:, nhinl]

    Xs = np.dot(H, Xo)
    Ys = Yo

    Xh = X[:, hinl]
    Yh = Y[:, hinl]

    len_xs = Xs.shape[1]
    m_i = sam_sizO
    max_sam = MAX_SAM
    no_sam = 0

    numBestInliers = 0
    bestInliers = []
    bestF = []

    while no_sam <= max_sam:
        sample = np.random.choice(range(len_xs), size=2, replace=False)
        no_sam = no_sam + 1
        ec = np.cross(np.cross(Xs[:,sample[0]],Ys[:,sample[0]]), np.cross(Xs[:,sample[1]],Ys[:,sample[1]]))
        ec = ec/ec[2]
        aFt = skew_sym(ec/np.linalg.norm(ec))
        Ftmp = np.dot(aFt,H)
        d = BASE_utils.sampsonDistanceF(Xo, Yo, Ftmp)
        inliersO = np.where(d<=2*th)[0]
        no_i = len(inliersO)
        if no_i > m_i:
            m_i = no_i
            max_sam = BASE_utils.updateMaxTrials(no_i, max_sam, len(nhinl), 0.999999, 2)
            F, inliers = innerFH(Xh, Yh, Xo[:, inliersO], Yo[:,inliersO], X, Y, th, 15, sam_sizH, sam_sizO)
            if len(inliers) > numBestInliers:
                numBestInliers = len(inliers)
                bestInliers = inliers
                bestF = F
        

    return bestF, bestInliers





