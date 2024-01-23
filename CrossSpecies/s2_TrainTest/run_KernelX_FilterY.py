from step0_TrainTestSetting_global import *

current_dir = os.getcwd()
print('current_dir = ', current_dir)

os.chdir(path_base)
os.system('mkdir -p output')
os.chdir('output')

os.system('mkdir -p %s' % params_global_constant['dir0'])
dir0_path = '%soutput/%s' % (path_base, params_global_constant['dir0'])


######################### TrainX_TestY ########################
for kernel_temp in params_global_dynamic['kernels']:
    for filter_temp in params_global_dynamic['filters']:
        for ont_temp in params_global_dynamic['onts']:
            def process_input(x):
                if isinstance(x, list):
                    y = '_'.join(str(num) for num in x)
                else:
                    y = x
                return y

            dir1 = f'DeepSS2GO_Kernel{process_input(kernel_temp)}_Filter{process_input(filter_temp)}_Onts{ont_temp}/'
            print('dir1 = ', dir1)
            dir1_path = os.path.join(dir0_path, dir1)
            os.makedirs(dir1_path, exist_ok=True)

            os.system(f'cp {current_dir}/* {dir1_path}')

            os.chdir(dir1_path)
            os.system('mkdir -p data')
            os.system('cp %s/pub_data/go.obo data/' % path_base)
            os.system('mkdir -p results')
            os.system('mkdir -p results_alpha')
            os.system('mkdir -p results_alphabeta')


            step0_TrainTestSetting_local_fpath = os.path.join(dir1_path, "step0_TrainTestSetting_local.py")
            with open(step0_TrainTestSetting_local_fpath, "w") as local_file:
                local_file.write("params_local = {\n")
                local_file.write("    # from params_global_constant\n")

                # global_constant
                for key, value in params_global_constant.items():
                    local_file.write(
                        f"    '{key}': {repr(value)},\n")

                # global_dynamic
                local_file.write("\n    # split: params_global_dynamic\n")
                if isinstance(kernel_temp, int):
                    local_file.write(f"    'kernels': [{kernel_temp}],\n")
                else:
                    local_file.write(f"    'kernels': {kernel_temp},\n")
                if isinstance(filter_temp, int):
                    local_file.write(f"    'filters': [{filter_temp}],\n")
                else:
                    local_file.write(f"    'filters': {filter_temp},\n")
                local_file.write(f"    'onts': {repr(ont_temp)},\n")
                local_file.write("}\n")



            os.chdir(dir1_path)
            if params_global_constant['run_step1_SplitTrainTest_Terms'] == 'T':
                os.system('python step1_SplitTrainTest_Terms.py')
            if params_global_constant['run_step2_Train'] == 'T':
                os.system('python step2_Train.py')
            if params_global_constant['run_step3_Test'] == 'T':
                os.system('python step3_Test.py')
            if params_global_constant['run_step4_pkl2fa'] == 'T':
                os.system('python step4_pkl2fa.py &')

            if params_global_constant['run_step7.1_EvaluateWithoutAlpha'] == 'T':
                os.system('bash step7.1_EvaluateWithoutAlpha.sh &')
            if params_global_constant['run_step8.1_PredictWithoutAlpha'] == 'T':
                os.system('bash step8.1_PredictWithoutAlpha.sh')

            if params_global_constant['run_step5_Diamond4CrossSpecies'] == 'T':
                os.system('bash step5_Diamond4CrossSpecies.sh')
            if params_global_constant['run_step6_FindAlpha'] == 'T':
                os.system('bash step6_FindAlpha.sh')
            if params_global_constant['run_step7_EvaluateAlpha'] == 'T':
                os.system('bash step7_EvaluateAlpha.sh')
            if params_global_constant['run_step8_PredictAlpha'] == 'T':
                os.system('bash step8_PredictAlpha.sh')

            os.chdir(dir0_path)


print('run_TrainTest done')

