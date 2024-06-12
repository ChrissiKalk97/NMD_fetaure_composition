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
            stop_pos_genome, cds_length = find_termination_codon(transcript_info, cds)
            if stop_pos_genome is None:
                print(transcript_info)
                continue
            stop_pos_transcript, last_ejc, exon_containing_stop_length, dist_stop_next_EJC\
                  = compose_transcript(transcript_info, stop_pos_genome)
            
            if stop_pos_transcript['stop_position'] == 0 or None:
                counter += 1
                NMD_features_df.loc[transcript_id, '50_nt'] = pd.NA
            elif (last_ejc - stop_pos_transcript['stop_position']) >= 50:
                NMD_features_df.loc[transcript_id, '50_nt'] = 1
                
            else:
                NMD_features_df.loc[transcript_id, '50_nt'] = 0
            NMD_features_df.loc[transcript_id, 'last_exon_length'] = stop_pos_transcript['length'] - last_ejc
            NMD_features_df.loc[transcript_id, 'end_ORF'] = stop_pos_transcript['stop_position']
            NMD_features_df.loc[transcript_id, 'distance_stop_EJC'] = stop_pos_transcript['length'] -\
                NMD_features_df.loc[transcript_id, 'last_exon_length'] -\
                NMD_features_df.loc[transcript_id, 'end_ORF']
            NMD_features_df.loc[transcript_id, 't_length'] = stop_pos_transcript['length']
            NMD_features_df.loc[transcript_id, 'start_ORF'] = stop_pos_transcript['stop_position'] - cds_length
            NMD_features_df.loc[transcript_id, 'exon_with_stop_length'] = exon_containing_stop_length
            NMD_features_df.loc[transcript_id, 'distance_stop_next_EJC'] = dist_stop_next_EJC
        
        else: 
            transcript_ids_wo_cds.append(transcript_id)
    NMD_features_df.drop(transcript_ids_wo_cds, inplace = True)
    return transcript_ids_wo_cds, NMD_features_df