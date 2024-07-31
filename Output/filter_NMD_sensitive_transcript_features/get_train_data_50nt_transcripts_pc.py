########################################################################################
#get_train_data_50nt_transcripts.py extracts NMD negative and positive sets that follow#
#the 50nt rule                                                                         #
#The user needs to provide the calculated NMD features for all transcripts, as et of   #
#differentially upregulated, downregulated and insignificant transcripts               #
########################################################################################

#usage: python get_train_data_50nt_transcripts.py features.csv up.csv down.csv insign.csv outname 0/1

import sys
import pandas as pd

feature_frame = pd.read_csv(sys.argv[1], header = 0, index_col = 0)
up_tids = pd.read_csv(sys.argv[2], header = 0, index_col = 0)
insign_tids = pd.read_csv(sys.argv[3], header = 0, index_col = 0)

outname = sys.argv[4]
num = float(sys.argv[5])

positive_features = feature_frame[feature_frame.index.isin(up_tids['x'])]
positive_features_50nt = positive_features[positive_features['50_nt'] == num]
#positive_features_50nt.to_csv(outname + 'positive_set.csv', index = True, header = True)
print(positive_features_50nt['50_nt'].head())

insign_features = feature_frame[feature_frame.index.isin(insign_tids.index)]
insign_features_50nt = insign_features[insign_features['50_nt'] == num]
print(len(insign_features_50nt['50_nt']), len(positive_features_50nt['50_nt']))

if len(insign_features_50nt.index) <=  len(positive_features_50nt.index):
    positive_features_50nt = positive_features_50nt.sample(n = len(insign_features_50nt.index), random_state = 100)
    negative_set = insign_features_50nt 
else:
    negative_set = insign_features_50nt.sample(n = len(positive_features_50nt.index), random_state = 100)
    positive_set = positive_features_50nt

positive_features_50nt.to_csv(outname + 'positive_set.csv', index = True, header = True)
negative_set.to_csv(outname + 'negative_set.csv', index = True, header = True)


