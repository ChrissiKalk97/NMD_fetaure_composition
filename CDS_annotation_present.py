import pandas as pd
from helper_functions_CDS_present import compose_transcript, find_termination_codon


def handle_cds_transcripts(transcript_gtftk_object, transcript_ids, NMD_features_df):
    transcript_ids_wo_cds = []
    counter = 0
    for transcript_id in transcript_ids:
        transcript_info = transcript_gtftk_object[transcript_id]
        #build sublist for each feature of the transcript
        transcript_info = [transcript_info[x:x+8] for x in range(0, len(transcript_info), 8)]
        #extract cds features
        cds = [sub_list for sub_list in transcript_info if 'CDS' in sub_list]
        #if there are CDS features: CDS is defined
        if len(cds) > 0:
            NMD_features_df.loc[transcript_id,'has_cds'] = 1
            stop_pos_genome, cds_length = find_termination_codon(transcript_info, cds)
            stop_pos_transcript, last_ejc = compose_transcript(transcript_info, stop_pos_genome)
            if stop_pos_transcript['stop_position'] == 0 or None:
                counter += 1
                NMD_features_df.loc[transcript_id, '50_nt'] = pd.NA
            elif (last_ejc - stop_pos_transcript['stop_position']) >= 50:
                NMD_features_df.loc[transcript_id, '50_nt'] = 1
                
            else:
                NMD_features_df.loc[transcript_id, '50_nt'] = 0
            NMD_features_df.loc[transcript_id, 'last_exon_length'] = last_ejc
            NMD_features_df.loc[transcript_id, 'end_ORF'] = stop_pos_transcript['stop_position']
            NMD_features_df.loc[transcript_id, 't_length'] = stop_pos_transcript['length']
            NMD_features_df.loc[transcript_id, 'start_ORF'] = stop_pos_transcript['stop_position'] - cds_length
        else: 
            NMD_features_df.loc[transcript_id,'has_cds'] = 0
            transcript_ids_wo_cds.append(transcript_id)
    NMD_features_df.drop(transcript_ids_wo_cds, inplace = True)
    return transcript_ids_wo_cds, NMD_features_df