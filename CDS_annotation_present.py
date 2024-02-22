from helper_functions_CDS_present import compose_transcript, find_termination_codon


def handle_cds_transcripts(transcript_gtftk_object, transcript_ids, NMD_features_df):
    transcript_ids_wo_cds = []
    counter = 0
    for transcript_id in transcript_ids:
        transcript_info = transcript_gtftk_object[transcript_id]
        #build sublist for each feature of the transcript
        transcript_info = [transcript_info[x:x+8] for x in range(0, len(transcript_info), 8)]
        #extract cds features
        cds = [sub_list for sub_list in transcript_info if "CDS" in sub_list]
        #if there are CDS features: CDS is defined
        if len(cds) > 0:
            NMD_features_df.loc[transcript_id,"has_cds"] = 1
            stop_pos_genome = find_termination_codon(transcript_info, cds)
            stop_pos_transcript, last_ejc = compose_transcript(transcript_info, stop_pos_genome)
            if stop_pos_transcript == 0:
                counter += 1
            if (last_ejc - stop_pos_transcript) >= 50:
                NMD_features_df.loc[transcript_id, "50_nt"] = 1
            else:
                NMD_features_df.loc[transcript_id, "50_nt"] = 0
        else: 
            NMD_features_df.loc[transcript_id,"has_cds"] = 0
            transcript_ids_wo_cds.append(transcript_id)
    print("counter", counter)
    return transcript_ids_wo_cds, NMD_features_df