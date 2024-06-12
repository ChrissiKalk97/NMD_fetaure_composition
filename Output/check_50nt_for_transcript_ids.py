import sys
import pandas as pd 

def main():
    pipeline_output = pd.read_csv(sys.argv[1], header = 0, index_col = 0)
    #transcripts_of_interest = pd.read_csv(sys.argv[2], header = 0)
    transcripts_of_interest = pd.read_csv(sys.argv[2], header = 0, index_col = 0)

    print(transcripts_of_interest.head())
    print(pipeline_output.head())
    pipeline_output_of_interest = pipeline_output[pipeline_output.index.isin(list(transcripts_of_interest.iloc[:,0]))]
    print(len(pipeline_output_of_interest.index))
    pipeline_output_of_interest = pipeline_output_of_interest[pipeline_output_of_interest['50_nt'].notna()]
    pipeline_output_of_interest['50_nt'] = pipeline_output_of_interest['50_nt'].astype(int)
    print("nr of transcripts", len(pipeline_output_of_interest.index))

    print('of which 50 nt trnascripts', sum(pipeline_output_of_interest['50_nt']))
    




if __name__ == "__main__":
    main()