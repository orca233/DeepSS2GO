import os

def rename_folders(root_folder):
    # 获取root文件夹下的所有子文件夹
    subfolders = [f.path for f in os.scandir(root_folder) if f.is_dir()]

    for folder in subfolders:
        # 用新名称替换子文件夹名称中的"NOF65536"
        new_folder_name = folder.replace("NOF65536", "onlyF65536")

        # 重命名文件夹
        os.rename(folder, new_folder_name)
        print(f"文件夹 {folder} 已重命名为 {new_folder_name}")

if __name__ == "__main__":
    root_folder = "/home/fsong/work/py_proj/prot_algo/DeepSS2GO/CrossSpecies/s2_HPC_ARATH_ECOLI_MYCTU_YEAST_onlyF65536"  # 请替换为你实际的root文件夹路径
    rename_folders(root_folder)
