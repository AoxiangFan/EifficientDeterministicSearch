import numpy as np
import BASE_utils

def iterH(X, Y, H, th, m, inlLimit):
    ILSQ_ITERS = 4
    ths = m*th
    dth = (ths - th) / ILSQ_ITERS

    error = BASE_utils.sampsonDistanceH_all(X, Y, H)
    inliers = np.where(error<=th)[0]

    numBestInliers = len(inliers)
    bestInliers = inliers
    bestH = H

    if len(inliers) <= inlLimit:
        if len(inliers) >= 4:
            H = BASE_utils.norm4Point(X[:,inliers], Y[:,inliers])
        else:
            return bestH, bestInliers
    else:
        inliersSubset = np.random.choice(list(inliers), size=inlLimit, replace=False)
        H = BASE_utils.norm4Point(X[:,inliersSubset], Y[:,inliersSubset])
    
    for i in range(ILSQ_ITERS):
        error = BASE_utils.sampsonDistanceH_all(X, Y, H)

        inliersA = np.where(error<=th)[0]
        inliersB = np.where(error<=ths)[0]

        if len(inliersA) > numBestInliers:
            numBestInliers = len(inliersA)
            bestH = H
            bestInliers = inliersA

        if len(inliersB) <= inlLimit:
            if len(inliersB) >= 4:
                H = BASE_utils.norm4Point(X[:,inliersB], Y[:,inliersB])
            else:
                return bestH, bestInliers
        else:
            inliersSubset = np.random.choice(list(inliersB), size=inlLimit, replace=False)
            H = BASE_utils.norm4Point(X[:,inliersSubset], Y[:,inliersSubset])
        
        ths = ths - dth

        

    return bestH, bestInliers


def LO_H(X, Y, H, th, m, inlLimit):

    RAN_REP = 10

    error = BASE_utils.sampsonDistanceH_all(X, Y, H)
    inliers = np.where(error<=th)[0]
    bestH = H
    bestInliers = inliers
    numBestInliers = len(bestInliers)
    
    iniltialInliers = np.where(error<=m*th)[0]


    # initialH = BASE_utils.norm4Point(X[:,iniltialInliers], Y[:,iniltialInliers])

    if len(iniltialInliers) <= inlLimit:
        if len(iniltialInliers) >=8:
            initialH = BASE_utils.norm4Point(X[:,iniltialInliers], Y[:,iniltialInliers])
        else:
            return bestH, bestInliers
    else:
        inliersSubset = np.random.choice(list(iniltialInliers), size=inlLimit, replace=False)
        initialH = BASE_utils.norm4Point(X[:,inliersSubset], Y[:,inliersSubset])


    error = BASE_utils.sampsonDistanceH_all(X, Y, initialH)
    baseInliers = np.where(error<=th)[0]

    if len(baseInliers) < 16:
        return bestH, bestInliers

    siz = min(int(np.ceil(len(baseInliers)/2)), 14)

    for i in range(RAN_REP):
        inliersSubset = np.random.choice(list(baseInliers), size=siz, replace=False)
        H = BASE_utils.norm4Point(X[:,inliersSubset], Y[:,inliersSubset])
        H, inliers = iterH(X, Y, H, th, m, inlLimit)

        if len(inliers) > numBestInliers:
            bestH = H
            bestInliers = inliers
            numBestInliers = len(inliers)

    return bestH, bestInliers


def iterF(X, Y, F, th, m, inlLimit):
    ILSQ_ITERS = 4

    ths = m*th
    dth = (ths - th) / ILSQ_ITERS

    error = BASE_utils.sampsonDistanceF(X, Y, F)
    inliers = np.where(error<=th)[0]

    numBestInliers = len(inliers)
    bestInliers = inliers
    bestF = F

    if len(inliers) <= inlLimit:
        if len(inliers) >=8:
            F = BASE_utils.norm8Point(X[:,inliers], Y[:,inliers])
        else:
            return bestF, bestInliers
    else:
        inliersSubset = np.random.choice(list(inliers), size=inlLimit, replace=False)
        F = BASE_utils.norm8Point(X[:,inliersSubset], Y[:,inliersSubset])
    
    for i in range(ILSQ_ITERS):
        error = BASE_utils.sampsonDistanceF(X, Y, F)

        inliersA = np.where(error<=th)[0]
        inliersB = np.where(error<=ths)[0]

        if len(inliersA) > numBestInliers:
            numBestInliers = len(inliersA)
            bestF = F
            bestInliers = inliersA

        if len(inliersB) <= inlLimit:
            if len(inliersB) >= 8:
                F = BASE_utils.norm8Point(X[:,inliersB], Y[:,inliersB])
            else:
                return bestF, bestInliers
        else:
            inliersSubset = np.random.choice(list(inliersB), size=inlLimit, replace=False)
            F = BASE_utils.norm8Point(X[:,inliersSubset], Y[:,inliersSubset])
        
        ths = ths - dth

        

    return bestF, bestInliers


def LO_F(X, Y, F, th, m, inlLimit):
    
    RAN_REP = 10

    error = BASE_utils.sampsonDistanceF(X, Y, F)
    inliers = np.where(error<=th)[0]
    bestF = F
    bestInliers = inliers
    numBestInliers = len(bestInliers)
    
    iniltialInliers = np.where(error<=m*th)[0]


    # initialF = BASE_utils.norm8Point(X[:,iniltialInliers], Y[:,iniltialInliers])

    if len(iniltialInliers) <= inlLimit:
        if len(iniltialInliers) >=8:
            initialF = BASE_utils.norm8Point(X[:,iniltialInliers], Y[:,iniltialInliers])
        else:
            return bestF, bestInliers
    else:
        inliersSubset = np.random.choice(list(iniltialInliers), size=inlLimit, replace=False)
        initialF = BASE_utils.norm8Point(X[:,inliersSubset], Y[:,inliersSubset])

    error = BASE_utils.sampsonDistanceF(X, Y, initialF)
    baseInliers = np.where(error<=th)[0]

    if len(baseInliers) < 16:
        return bestF, bestInliers

    siz = min(int(np.ceil(len(baseInliers)/2)), 14)


    for i in range(RAN_REP):
        inliersSubset = np.random.choice(list(baseInliers), size=siz, replace=False)
        F = BASE_utils.norm8Point(X[:,inliersSubset], Y[:,inliersSubset])
        F, inliers = iterF(X, Y, F, th, m, inlLimit)

        if len(inliers) > numBestInliers:
            bestF = F
            bestInliers = inliers
            numBestInliers = len(inliers)

    return bestF, bestInliers