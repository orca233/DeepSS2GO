import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from CrossSpecies_db import *



# 下面6个选一个激活

##### Fmax #####
# data, vmin, vmax, cmap, top_title = data_bp_Fmax_aa, 0.25, 0.51, 'Greens', 'data_bp_Fmax_aa'
# data, vmin, vmax, cmap, top_title = data_cc_Fmax_aa, 0.45, 0.84, 'Blues', 'data_cc_Fmax_aa'
# data, vmin, vmax, cmap, top_title = data_mf_Fmax_aa, 0.27, 0.67, 'Oranges', 'data_mf_Fmax_aa'
# data, vmin, vmax, cmap, top_title = data_bp_Fmax_ss8, 0.25, 0.51, 'Greens', 'data_bp_Fmax_ss8'
# data, vmin, vmax, cmap, top_title = data_cc_Fmax_ss8, 0.45, 0.84, 'Blues', 'data_cc_Fmax_ss8'
# data, vmin, vmax, cmap, top_title = data_mf_Fmax_ss8, 0.27, 0.67, 'Oranges', 'data_mf_Fmax_ss8'

##### AUPR #####
# data, vmin, vmax, cmap, top_title = data_bp_AUPR_aa, 0.14, 0.45, 'Greys', 'data_bp_AUPR_aa'
# data, vmin, vmax, cmap, top_title = data_cc_AUPR_aa, 0.31, 0.85, 'copper_r', 'data_cc_AUPR_aa'
# data, vmin, vmax, cmap, top_title = data_mf_AUPR_aa, 0.2, 0.56, 'Reds', 'data_mf_AUPR_aa'
# data, vmin, vmax, cmap, top_title = data_bp_AUPR_ss8, 0.14, 0.45, 'Greys', 'data_bp_AUPR_ss8'
# data, vmin, vmax, cmap, top_title = data_cc_AUPR_ss8, 0.31, 0.85, 'copper_r', 'data_cc_AUPR_ss8'
data, vmin, vmax, cmap, top_title = data_mf_AUPR_ss8, 0.2, 0.56, 'Reds', 'data_mf_AUPR_ss8'



df = pd.DataFrame(data)
df.set_index('Test', inplace=True)

# 使用Seaborn创建热图
plt.figure(figsize=(10, 8))
# 设置vmin和vmax来调整颜色映射范围
heatmap = sns.heatmap(df, annot=True, cmap=cmap, fmt=".3f", linewidths=.5, cbar=True, vmin=vmin, vmax=vmax)  # cmap = 'coolwarm'

# 获取颜色条对象
cbar = heatmap.collections[0].colorbar

# 设置颜色条的刻度和标签
ticks = np.arange(vmin, vmax, 0.02)  # 设置刻度间隔为0.02
cbar.set_ticks(ticks)
cbar.set_ticklabels([f'{tick:.2f}' for tick in ticks])

# 在上方创建另一个坐标轴
ax2 = heatmap.twiny()
ax2.set_xticks(heatmap.get_xticks())
ax2.set_xticklabels(heatmap.get_xticklabels(), rotation=0)

# 设置上方和下方的坐标轴标签
heatmap.set_xlabel('X-axis (Bottom)', labelpad=15)
ax2.set_xlabel('X-axis (Top)', labelpad=15)

plt.title('Heatmap of Data = ' + top_title)
plt.show()
