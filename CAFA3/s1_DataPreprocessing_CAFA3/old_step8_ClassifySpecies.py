import pandas as pd
from collections import Counter
from step0_DataPreprocessingSetting import *

# swissprot_x_cullPDB_all_aa/ss8.pkl按物种分类： HUIMAN/MOUSE/RAT/ECOLI。。。
# swissprot 'protein' 名字里有分类  e.g. TEN1L_HUMAN
# 以及选中5个要保留的物种类别 keep_species_0/1/2/3/4/5 = xxxxxx， 执行 2 遍，分别执行aa_ss = 'aa' / 'ss8'


# 对于cullPDB：
# 'HUMAN': 198, 'ECOLI': 198, 'MOUSE': 113, 'YEAST': 79, 'RAT': 71, 'MYCTU': 33, 'PSEAE': 25, 'BOVIN': 21,
# 'DROME': 17, 'BACSU': 13, 'ARATH': 13, 'CHICK': 11, 'BPT4': 10, 'SCHPO': 9, 'SALTY': 8, 'PYRHO': 6, 'CAEEL': 6,



# 对于全集swissprot
# 这一行代码：species = Counter(prot_species_lst)  # 统计每种prot_species出现多少次 ！！！！！！！可以再输出最下面部分查看
# Counter({'HUMAN': 13238, 'MOUSE': 10152, 'ARATH': 10002, 'YEAST': 4964, 'SCHPO': 4296, 'RAT': 4217, 'ECOLI': 3288,
# 'DROME': 2472, 'CAEEL': 2042, 'MYCTU': 1456,
# 'DANRE': 1108, 'DICDI': 985, 'ORYSJ': 911, 'BOVIN': 669, 'CHICK': 511, 'XENLA': 474, 'PSEAE': 412, 'CANAL': 400,
# 'EMENI': 316, 'BACSU': 272, 'PIG': 216, 'MAIZE': 170, 'ASPFU': 152, 'ORYSI': 126, 'RABIT': 120, 'PLAF7': 119,
# 'SALTY': 115, 'CANGB': 106, 'CANLF': 104, 'SYNY3': 101,


# 对于经过 step2筛选后的 swissprot 物种个数species_num （除去过长的）
# 最后一段代码
# ['HUMAN: 13238', 'ARATH: 10002', 'MOUSE: 10152', 'YEAST: 4964', 'SCHPO: 4296', 'RAT: 4217', 'ECOLI: 3288',
# 'DROME: 2472', 'CAEEL: 2042', 'MYCTU: 1456']



### backup
# 'ARATH'  # Arabidopsis thaliana 拟南芥
# 'YEAST'  # 酵母
# 'SCHPO'  # 另一种酵母
# 'ECOLI'  # 大肠杆菌
# 'DROME'  # 黑腹果蝇
# 'CAEEL'  # Caenorhabditis elegans 秀丽隐杆线虫
# 'MYCTU'  # Mycobacterium tuberculosis 结核杆菌  原核？？？



for aa_ss in ['aa', 'ss8', 'ss3']:
    for keep_species in ['HUMAN', 'ARATH', 'MOUSE', 'YEAST', 'SCHPO', 'RAT', 'ECOLI', 'DROME', 'CAEEL', 'MYCTU']:  #
        data_df = pd.read_pickle(path_pub_data + 'swissprot_clean_ALL00_' + aa_ss + '.pkl')  # 每次循环要刷新，
        # data_df = data_df[0:100]  # 测试行，运行时注释掉  !!!!!!!!!!!!!!!!!!!!!

        prot_species_lst = []  # 统计每种prot_species出现多少次，每次循环要刷新，

        k = 0  # 计数 该物种，保留下来的
        for index, row in data_df.iterrows():
            prot_name_species = row['proteins']  # e.g. prot_name_species = TEN1L_HUMAN, prot_species[1]=HUMAN
            prot_species = prot_name_species.split('_')
            prot_species_lst.append(prot_species[1])  # 统计prot_species出现多少次

            if prot_species[1] != keep_species:  # # 不为“keep_species“，则删除改行
                print('del')
                data_df = data_df.drop(index=index)
            else:
                print(keep_species + str(k))  # 该物种保留下来了，统计个数
                k += 1
        # 保存pkl
        data_df.to_pickle(path_pub_data + 'swissprot_clean_' + keep_species + '_' + aa_ss + '.pkl')
        species = Counter(prot_species_lst)  # 统计每种prot_species出现多少次 ！！！！！！！，可以再输出最下面部分查看

        print('check specific proteins species::::::::::::' + keep_species)
        print(data_df['proteins'])  # 看看是否为同一种species
        print('\n')
        print('species count:::::::')
        print(species)


print('\n')
print('--------------------------------------------\n')
print('\n')

# 对于经过 step2筛选后的 swissprot 物种个数 species_num （除去过长的），这是我们已经知道了这十类最靠前
output_species_num = []
for keep_species in ['HUMAN', 'ARATH', 'MOUSE', 'YEAST', 'SCHPO', 'RAT', 'ECOLI', 'DROME', 'CAEEL', 'MYCTU']:
    fpath = '%sswissprot_clean_%s_aa.pkl' % (path_pub_data, keep_species)
    df = pd.read_pickle(fpath)
    species_num = '%s: %d' % (keep_species, len(df))
    output_species_num.append(species_num)
print(output_species_num)



print('done')


