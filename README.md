## LDS-toolbox: a matlab toolbox for linear dynamical systems (LDSs) modeling

#### Authors: [Wenbing Huang](https://sites.google.com/site/wenbinghuangshomepage/)

### Overview

Linear Dynamical Systems (LDSs) are fundamental tools for modeling spatiotemporal data in various disciplines. 
Though rich in modeling, analyzing LDSs is not free of difficulty, mainly because LDSs do not comply with Euclidean geometry
and hence conventional learning techniques can not be applied directly. Specifically, LDSs apply parametric equations to model
the spatio-temporal data. The optimal system parameters (i.e. the system tuple (A,C)) learned from the input are employed as
the descriptor of each spatio-temporal sequence. 

With this toolbox, you can: 1) learn the stable LDS system tuple for any given sequence via various methods; 2) perform clustering or sparse coding on the space of LDSs; 3) conduct classification on spatiotemporal data, e.g. videos and tactile sequences.   

If you make use of this toolbox in your work, please cite the following papers:

   @inproceedings{huang2017efficient, 
         title={Efficient Optimization for Linear Dynamical Systems with Applications to Clustering and Sparse Coding.}, 
         author={Huang, Wenbing and Mehrtash, Harandi and Tong, Zhang and Lijie, Fan and Fuchun, Sun and Junzhou, Huang}, 
         booktitle={Advances in Neural Information Processing Systems (NIPS)}, 
         pages={3446--3456}, 
         year={2017} 
     } 

   @inproceedings{wenbing2016sparse, 
         title={Sparse coding and dictionary learning with linear dynamical systems}, 
         author={Huang, Wenbing and Sun, Fuchun and Cao, Lele and Zhao, Deli and Liu, Huaping and Harandi, Mehrtash}, 
         booktitle={IEEE Conference on Computer Vision and Pattern Recognition (CVPR)}, 
         pages={3938--3947}, 
         year={2016}, 
         organization={IEEE} 
     } 

   @inproceedings{huang2015scalable, 
         title={ Learning Stable Linear Dynamical Systems with the Weighted Least Square Method.}, 
         author={Huang, Wenbing and Cao, Lele and Sun, Fuchun and Zhao, Deli and Liu, Huaping and Yu, Shanshan}, 
         booktitle={Proceedings of the International Joint Conference on Artificial Intelligence (IJCAI)}, 
         pages={1599--1605}, 
         year={2016} 
     }

### Requirements
This package is almost self-contained, unless you want to recomplie the svm source codes in ./util/liblinear or ./util/libsvm.

### Running the code
1. To learn stable LDS tuples (A,C) from any input sequence, please run demo_dyntex;
2. For clustering of LDS tuples, please run demo_clustering;
3. For sparse coding, please run demo_dictionarylearning.

#### Input format
See the example in ./dataset/dyntex_plus/dyntex++lbp8-6.mat. The data are recorded in the matlab cell format. Different cells represent different sequences. Each sequence is restored by the SxT matrix, where S and T are the spatial and temporal dimensions, respectively.
