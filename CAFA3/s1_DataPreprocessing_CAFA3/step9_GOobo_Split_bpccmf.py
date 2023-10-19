# 该py，把go.obo分成 go_bp/cc/mf.obo，保存在pub_data里
from step0_DataPreprocessingSetting import *

def split_obo_file(file_path):  # 把go.obo分成 go_bp/cc/mf.obo，保存在pub_data里
    with open(file_path, 'r') as file:
        data = file.read()

    go_all_nb = data.count('[Term]')
    print('go_all_nb = ', go_all_nb)  # len=47497

    sections = data.split('[Term]')
    mf_data = ''
    bp_data = ''
    cc_data = ''

    for section in sections:
        if 'namespace: biological_process' in section:
            bp_data += '[Term]' + section
        elif 'namespace: cellular_component' in section:
            cc_data += '[Term]' + section
        elif 'namespace: molecular_function' in section:
            mf_data += '[Term]' + section

    with open(path_pub_data + 'go_bp.obo', 'w') as file:
        file.write(bp_data)
        go_bp_nb = bp_data.count('[Term]')
        print('go_bp_nb = ', go_bp_nb)

    with open(path_pub_data + 'go_cc.obo', 'w') as file:
        file.write(cc_data)
        go_cc_nb = cc_data.count('[Term]')
        print('go_cc_nb = ', go_cc_nb)

    with open(path_pub_data + 'go_mf.obo', 'w') as file:
        file.write(mf_data)
        go_mf_nb = mf_data.count('[Term]')
        print('go_mf_nb = ', go_mf_nb)


# 指定 GO 合集文件的路径
file_path = path_pub_data+'go.obo'

# 调用函数将 GO 合集文件分割成 MF、BP 和 CC 三个文件
split_obo_file(file_path)


# go_bp_nb =  30543
# go_cc_nb =  4474
# go_mf_nb =  12480
# all = 47497


def extract_go_terms(filename):  # 提取每一个go.obo, go_bp/cc/mf.obo，并存在set里。
    go_terms = set()
    with open(filename, 'r') as file:
        contents = file.read()
        term_sections = contents.split('[Term]')
        for term in term_sections:
            lines = term.strip().split('\n')
            term_id = None
            for line in lines:
                if line.startswith('id:'):
                    term_id = line.split('id:')[1].strip()
                    break
            if term_id:
                go_terms.add(term_id)
    return go_terms


print('\n--------- extract_go_terms from each go(bp/cc/mf).obo ----------')
go_file = path_pub_data + 'go.obo'  # 替换成你的GO文件名
go_terms_set = extract_go_terms(go_file)

go_bp_file = path_pub_data + 'go_bp.obo'
go_bp_terms_set = extract_go_terms(go_bp_file)

go_cc_file = path_pub_data + 'go_cc.obo'
go_cc_terms_set = extract_go_terms(go_cc_file)

go_mf_file = path_pub_data + 'go_mf.obo'
go_mf_terms_set = extract_go_terms(go_mf_file)


print('go_terms_set = ', len(go_terms_set))
print('go_bp_terms_set = ', len(go_bp_terms_set))
print('go_cc_terms_set = ', len(go_cc_terms_set))
print('go_mf_terms_set = ', len(go_mf_terms_set))


