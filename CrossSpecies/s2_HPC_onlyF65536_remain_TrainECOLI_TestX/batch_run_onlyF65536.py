import os

# 指定root文件夹
#root_folder = "/scem/work/songfu/py_proj/prot_algo/DeepSS2GO/CrossSpecies/s2_HPC_ARATH_ECOLI_MYCTU_YEAST_onlyF65536/"

root_folder = os.getcwd()

# 遍历子文件夹
for folder in os.listdir(root_folder):
    folder_path = os.path.join(root_folder, folder)

    if os.path.isdir(folder_path):
        # 进入子文件夹
        os.chdir(folder_path)

        # 查找并运行run_KernelX_FilterY.py脚本
        script_path = "run_KernelX_FilterY.py"
        if os.path.exists(script_path):
            print(f"Running {script_path} in {folder_path}")
            os.system("python " + script_path)
        else:
            print(f"No {script_path} found in {folder_path}")

        # 返回上级文件夹
        os.chdir("..")
