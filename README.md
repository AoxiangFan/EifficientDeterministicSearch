# EifficientDeterministicSearch

Implementation of ECCV 2020 paper "Geometric Estimation via Robust Subspace Recovery" and IEEE TPAMI 2021 paper "Efficient Deterministic Search with Robust Loss
Functions for Geometric Model Fitting". 

This repo is now a rough version which only implements MATLAB code. A python implementation is coming soon. 

EES_MATLAB contains the code of the proposed Efficient Exact Search method for homography, fundamental matrix and linear estimation.
EAS_MATLAB contains the code of the proposed Efficient Approximate Search method for homography and fundamental matrix estimation. The EAS method is recommended for geometric estimation task which has better performance in our evaluation.

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

# usage

To use the code, run initialization.m first and run demo for simple examples.


