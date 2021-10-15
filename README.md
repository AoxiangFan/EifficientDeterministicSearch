# Description

Implementation of ECCV 2020 paper "Geometric Estimation via Robust Subspace Recovery" and IEEE TPAMI 2021 paper "Efficient Deterministic Search with Robust Loss
Functions for Geometric Model Fitting". 

EES_MATLAB contains the code of the proposed Efficient Exact Search method for homography, fundamental matrix and linear estimation.

EAS_MATLAB contains the code of the proposed Efficient Approximate Search method for homography and fundamental matrix estimation. The EAS method is recommended for geometric estimation task which has better performance in our evaluation.

EAS_Python contains the Python code of the proposed Efficient Approximate Search method.

If you find this code useful for your research, please cite our paper:

```
@inproceedings{fan2020geometric,
  title={Geometric Estimation via Robust Subspace Recovery},
  author={Fan, Aoxiang and Jiang, Xingyu and Wang, Yang and Jiang, Junjun and Ma, Jiayi},
  booktitle={Computer Vision--ECCV 2020: 16th European Conference, Glasgow, UK, August 23--28, 2020, Proceedings, Part XXII 16},
  pages={462--478},
  year={2020},
  organization={Springer}
}
```
and
```
@article{fan2021efficient,
  title={Efficient Deterministic Search with Robust Loss Functions for Geometric Model Fitting},
  author={Fan, Aoxiang and Ma, Jiayi and Jiang, Xingyu and Ling, Haibin},
  journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
  year={2021},
  publisher={IEEE}
}
```

# Usage of MATLAB code

To use the MATLAB code, run initialization.m first and run demo for simple examples.

# Usage of Python code
To use the Python code, since the core functions are dependent solely on numpy, one only requires to add src files into system path to invoke EAS.

The file demo.py is given to show a simple example of fundamental matrix and homography estimation.

Note that scipy and opencv-python packages are additionally required to run the demo. 

# Performance
Since the code is implemented solely using numpy package, it can be a bit slow, taking a few seconds for challenging data.

We have tested the performance of EAS-A in the renowned image matching benchmark (https://github.com/ubc-vision/image-matching-benchmark). Below is the comparison result of EAS-A and several state-of-the-art robust estimators from the benchmark.

The test sequences in the figure are reichstag, sacre coeur and st peters square from Phototourism dataset. The competitors follow the recommended parameter settings as in https://www.cs.ubc.ca/research/image-matching-challenge/2021/submit/.
![ransac-rt-example](https://user-images.githubusercontent.com/59504874/137460726-6756bcfc-4a20-4d98-878f-7e7a17a8c654.png)



