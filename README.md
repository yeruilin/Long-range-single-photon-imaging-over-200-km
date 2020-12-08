# Single-photon imaging over 200 km
 
This demo includes data and MATLAB codes used in the paper "Single-photon imaging over 200 km"
by Zheng-Ping Li, Jun-Tian Ye, Xin Huang, Peng-Yu Jiang, Yuan Cao, Yu Hong, Chao Yu, Jun Zhang, Qiang Zhang, Cheng-Zhi Peng, Feihu Xu, and Jian-Wei Pan
Under review

Corresponding author: feihuxu@ustc.edu.cn.

To try the codes, just download the zip and run the "demo_12_08_2020.m" in MATLAB. 
Warning: the code was tested using MATLAB 2017b and it might be incompatible with older versions.  

## Attribution

### Data
The experimental lidar data is obtained by our new single-photon lidar system, 
including a building 'tower.mat' over 9.8 km, a mountain over  124.1 km (waiting), and another mountain over 201.5 km (waiting).

Instructions for the experimental data:
- 'data_processed': (type: cell), it contains the TOF information of the arrival photons (unit: ps).
- 'mm' and 'MM' record the start time and stop time of the data.

### Code
We use the following three photon-efficient algorithms for reconstruction:
- Shin Dongeek, et al., "Photon-efficient computational 3-D and reflectivity imaging with single-photon detectors.", IEEE Transactions on Computational Imaging 1, 112 (2015).
--The related code can be found at https://github.com/photon-efficient-imaging.

- Rapp Joshua, et al., "A few photons among many: Unmixing signal and noise for photon-efficient active imaging.", IEEE Transactions on Computational Imaging 3, 445 (2017).
--The related code can be found at https://codeocean.com/capsule/6269424/tree/v1.

- Li, Zheng-Ping, et al., "Single-photon computational 3D imaging at 45 km.", Photonics Research 8, 1532 (2020).
--The related code can be found at https://github.com/quantum-inspired-lidar/long-range-photon-efficient-imaging




