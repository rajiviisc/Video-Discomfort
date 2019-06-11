IISc Video Discomfort Dataset

This dataset consists of:
1. Videos that can be generated using Blender files containing the mapped camera trajectories
2. Discomfort scores provided by subjects after watching the videos

References

If you use our dataset or code, please cite our paper:

@inproceedings{vdvr,
  title={Prediction of Discomfort due to Egomotion in Immersive Videos for Virtual Reality},
  author={Balasubramanian, Suprith and Soundararajan, Rajiv},
  booktitle={IEEE International Symposium on Mixed and Augmented Reality (ISMAR)},
  year={2019}
}

Videos

1. The videos can be generated using Blender and the .blend files in the ‘BlenderFiles’ folder to be downloaded from http://ece.iisc.ac.in/~rajivs/databases/BlenderFiles.rar 

2. Several of the videos can be directly generated using the .blend files. Some videos require dependencies. Different subfolders under BlenderFiles contain the necessary dependencies. The .blend files are appropriately located so that the dependencies are taken care of. Note that the .blend files for 100 videos are located across different folders and subfolders in BlenderFiles. 

3. The .blend files were generated after mapping specific trajectories from different datasets to individual scenes in the files. Camera trajectories were then modified using Blender parameters to generate the final sequences contained in the files.
 
Discomfort scores

The ‘StudyData’ folder contains 43 csv files corresponding to 43 subjects who took the study, with each csvfile containing the sequence of videos watched by the subject as well as the corresponding discomfort scores in the format [Video ID, Discomfort Score].

For password details to open the StudyData folder, please email rajivs@iisc.ac.in

Videos with ID = {5, 15, 25, 35, 45, 60, 70, 80, 90, 100} are Anchor videos, which are natural 360-degree videos.

