import sys
import pandas as pd

feature_frame = pd.read_csv(sys.argv[1], header = 0, index_col = 0)
up_tids = pd.read_csv(sys.argv[2], header = 0, index_col = 0)
down_tids = pd.read_csv(sys.argv[3], header = 0, index_col = 0)

number_transcripts = sum(feature_frame.index.isin(up_tids.index)) - len(down_tids.index)
print(number_transcripts)

insign_tids = pd.read_csv(sys.argv[4], header = 0, index_col = 0)
print(insign_tids.head())
try:
	features_selected = feature_frame[feature_frame.index.isin(insign_tids['TXNAME'])]
except:
	features_selected = feature_frame[feature_frame.index.isin(insign_tids['x'])]
print(len(features_selected.index))
features_down = feature_frame[feature_frame.index.isin(down_tids.index)]
negative_set = features_selected.sample(n = number_transcripts, random_state = 100)
negative_set = pd.concat([negative_set, features_down], axis=0, join='outer')

negative_set.to_csv(sys.argv[5], index = True, header = True)
print(len(negative_set.index))