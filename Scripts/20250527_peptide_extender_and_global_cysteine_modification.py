#Created Date: 20220805
#Modified Date: 20250527
#Author: Cynthia M. Galicia-Medina

#Import Peptide List CSV file to Python
import pandas as pd
peptide_list = pd.read_csv (r'C:\Users\Yourname\Documents\peptidelist.csv')
print (peptide_list)

#Import Mouse Proteome Database from "UNITPROT" to Python
MouseProteomeDatabase = pd.read_csv(r'C:\Users\Yourname\Documents\Mouse_Proteome2023.csv')


#Merge Peptide List and Mouse Proteome Database to obtain the protein full sequences.
fullpetideseq = pd.merge(peptide_list, MouseProteomeDatabase, how='inner')
fullpetideseq.to_csv(r'C:\Users\Yourname\Documents\Merge_table_full_protein_sequence.csv', index=False)
print(fullpetideseq)

#Find the digested peptide sequences in the full protein sequence and extended desired peptide into a 15 length peptide sequence. 
alignpeptidesequence= []
for result in fullpetideseq.iloc:
    peptide_sequence = result ['Peptide Sequence']
    sequence = result ['Sequence']
    index = sequence.find(peptide_sequence)
    position_modification = int(result ['First Cysteine Extracted']) #Change depending on the position of the cysteine "First","Second", "Third"
    interest_character_index = (index + position_modification -1)
    alignpeptidesequence.append((sequence[interest_character_index - 7: interest_character_index + 8]))
print(alignpeptidesequence)
dict1 ={'Extended Peptide Sequence':alignpeptidesequence}
AlignsequenceDF=pd.DataFrame(dict1)
#Save the aligned peptide sequence into a csv file
AlignsequenceDF.to_csv(r'C:\Users\Yourname\Documents\alignpeptidesequences.csv', index=False)

#Find the global position in the protein of the cysteine of interested
globalcysteineposition= []
for result2 in fullpetideseq.iloc:
    peptide_sequence2 = result2 ['Peptide Sequence']
    sequence2 = result2 ['Sequence']
    index2 = sequence2.find(peptide_sequence2)
    position_modification2 = int(result2 ['First Cysteine Extracted']) #Change depending on the position of the cysteine "First","Second", "Third"
    interest_character_index2 = (index2 + position_modification2)
    globalcysteineposition.append(interest_character_index2 )
print(globalcysteineposition)
dict2 ={'Cysteine Modification Position':globalcysteineposition}
CysteineModificationPositionDF=pd.DataFrame(dict2)
#Save the aligned global cysteine position into a csv file
CysteineModificationPositionDF.to_csv(r'C:\Users\Yourname\Documents\globalcysteineposition.csv', index=False)
