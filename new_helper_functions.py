import regex as re
from typing import List

def get_transcript_string(transcript_ids: List[str]) -> str:
    transcript_string = ""
    for transcript_id in transcript_ids[:-1]:
        transcript_string += transcript_id+","
    transcript_string += transcript_ids[-1]
    return transcript_string

def get_cds_start(cds : List[List[str]], strand : str):
    start_position_plus = float("inf")
    start_position_minus = 0
    for partial_cds in cds:
        if strand == "+":
            five_prime = int(partial_cds[1])
            if start_position_plus < five_prime:
                start_position_plus = five_prime
                partial_cds[2] = str(int(partial_cds[1])+2)
        else:
            five_prime = int(partial_cds[2])
            if start_position_minus > five_prime:
                start_position_minus = five_prime
                partial_cds[1] = str(int(partial_cds[2])-2)
    return partial_cds





    
   