The file Sequence_map.csv contains the mapping of Video IDs to the trajectories in the sequence.
The raw discomfort scores and Z-scores are stored in the files scores.mat and zscores.mat respectively. The code for generating this data from the files in StudyData is in the file 'full_ratings.m'.
The discomfort prediction for both the average user as well as the user-specific case is computed using the regmodel.m file by setting flag values to 1 and 2 respectively.
For the user-specific case, the number of training videos can be specified by changing the value of the variable 'tvNum'. 
