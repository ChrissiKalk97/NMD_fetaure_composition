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
down_tids = pd.read_csv(sys.argv[3], header = 0, index_col = 0)
insign_tids = pd.read_csv(sys.argv[4], header = 0, index_col = 0)

outname = sys.argv[5]
num = float(sys.argv[6])

positive_features = feature_frame[feature_frame.index.isin(up_tids['x'])]
positive_features_50nt = positive_features[positive_features['50_nt'] == num]
positive_features_50nt.to_csv(outname + 'positive_set.csv', index = True, header = True)
print(positive_features_50nt['50_nt'].head())

negative_features = feature_frame[feature_frame.index.isin(down_tids['x'])]
negative_features_50nt = negative_features[negative_features['50_nt'] == num]
if len(negative_features_50nt.index) <=  len(positive_features_50nt.index):
    number_upsample = len(positive_features_50nt.index) - len(negative_features_50nt.index)
    print(number_upsample)
    features_insign = feature_frame[feature_frame.index.isin(insign_tids['TXNAME'])]
    features_insign_50nt = features_insign[features_insign['50_nt'] == num]
    print(features_insign_50nt['50_nt'].head())
    negative_set = features_insign_50nt.sample(n = number_upsample, random_state = 100)
    negative_set = pd.concat([negative_set, negative_features_50nt], axis=0, join='outer')
    negative_set.to_csv(sys.argv[5] + 'negative_set.csv', index = True, header = True)
    print(negative_set['50_nt'].head())
else:
    negative_set = negative_features_50nt.sample(n = len(positive_features_50nt.index), random_state = 100)
    print(len(positive_features_50nt.index))
    negative_set.to_csv(sys.argv[5] + 'negative_set.csv', index = True, header = True)
    print(negative_set['50_nt'].head())

