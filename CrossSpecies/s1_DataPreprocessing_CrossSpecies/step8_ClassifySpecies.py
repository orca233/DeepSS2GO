import pandas as pd
from collections import Counter
from step0_DataPreprocessingSetting import *



print('\n################## a long, long time ago ... ##################\n')
print('# starting step8_ClassifySpecies #')

for aa_ss in ['aa', 'ss8', 'ss3']:
    for keep_species in ['HUMAN', 'ARATH', 'MOUSE', 'YEAST', 'SCHPO', 'RAT', 'ECOLI', 'DROME', 'CAEEL', 'MYCTU']:
        data_df = pd.read_pickle(path_pub_data + 'swissprot_clean_ALL00_' + aa_ss + '.pkl')

        prot_species_lst = []

        k = 0  #
        for index, row in data_df.iterrows():
            prot_name_species = row['proteins']
            prot_species = prot_name_species.split('_')
            prot_species_lst.append(prot_species[1])

            if prot_species[1] != keep_species:
                print('del')
                data_df = data_df.drop(index=index)
            else:
                print(keep_species + str(k))
                k += 1

        data_df.to_pickle(path_pub_data + 'swissprot_clean_' + keep_species + '_' + aa_ss + '.pkl')
        species = Counter(prot_species_lst)

        print('check specific proteins species::::::::::::' + keep_species)
        print(data_df['proteins'])
        print('\n')
        print('species count:::::::')
        print(species)

print('\n--------------------------------------------\n\n')

output_species_num = []
for keep_species in ['HUMAN', 'ARATH', 'MOUSE', 'YEAST', 'SCHPO', 'RAT', 'ECOLI', 'DROME', 'CAEEL', 'MYCTU']:
    fpath = '%sswissprot_clean_%s_aa.pkl' % (path_pub_data, keep_species)
    df = pd.read_pickle(fpath)
    species_num = '%s: %d' % (keep_species, len(df))
    output_species_num.append(species_num)
print(output_species_num)


print('done')
print('\n################## And they all lived happily ever after! ##################\n')


