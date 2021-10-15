import scipy.io
import cv2
import numpy as np

import sys
sys.path.append('./src')
import EAS
import BASE_utils


def draw_match(img1, img2, corr1, corr2):

    corr1 = [cv2.KeyPoint(corr1[i, 0], corr1[i, 1], 1) for i in range(corr1.shape[0])]
    corr2 = [cv2.KeyPoint(corr2[i, 0], corr2[i, 1], 1) for i in range(corr2.shape[0])]

    assert len(corr1) == len(corr2)

    draw_matches = [cv2.DMatch(i, i, 0) for i in range(len(corr1))]

    display = cv2.drawMatches(img1, corr1, img2, corr2, draw_matches, None,
                              matchColor=(0, 255, 0),
                              singlePointColor=(0, 0, 255),
                              flags=0
                              )
    return display

if __name__ == "__main__":

    data = scipy.io.loadmat('data/box.mat')
    X = data['X0']
    Y = data['Y0']
    I1 = data['I1']
    I2 = data['I2']
    vpts = data['vpts']

    embeddingtype = 'A'
    th = 2.0
    maxIter = 500
    F, mask = EAS.findFundamentalMatrix(X, Y, embeddingtype, th, maxIter)
    print("Mean Geometric Error: {} pixels".format(np.mean(BASE_utils.sampsonDistanceF(vpts[0:3,:],vpts[3:6,:],F))))

    display = draw_match(I1, I2, X[mask,:], Y[mask,:])
    cv2.imshow("fundamental matrix estimation visualization", display)
    print('please press any key to terminate window')
    k = cv2.waitKey(0)
    cv2.destroyAllWindows()

    data = scipy.io.loadmat('data/graf.mat')
    X = data['X0']
    Y = data['Y0']
    I1 = data['I1']
    I2 = data['I2']
    vpts = data['vpts']

    embeddingtype = 'A'
    th = 1.0
    maxIter = 500
    H, mask = EAS.findHomography(X, Y, 'A', th, maxIter)
    print("Mean Geometric Error: {} pixels".format(np.mean(BASE_utils.sampsonDistanceH(vpts[0:3,:],vpts[3:6,:],H))))

    display = draw_match(I1, I2, X[mask,:], Y[mask,:])
    cv2.imshow("homography estimation visualization", display)
    print('please press any key to terminate window')
    k = cv2.waitKey(0)
    cv2.destroyAllWindows()

