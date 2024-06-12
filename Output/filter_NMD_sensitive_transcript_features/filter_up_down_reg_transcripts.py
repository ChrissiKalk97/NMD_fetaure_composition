import sys
import pandas as pd

feature_frame = pd.read_csv(sys.argv[1], header = 0, index_col = 0)
tids = pd.read_csv(sys.argv[2], header = 0, index_col = 0)
print(len(tids.index))

features_selected = feature_frame[feature_frame.index.isin(tids['x'])]
features_selected.to_csv(sys.argv[3], index = True, header = True)
print(len(features_selected.index))