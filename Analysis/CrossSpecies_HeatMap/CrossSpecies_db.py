
# 创建DataFrame

################################ bp ####################################

data_bp_Fmax_aa = {
    'Test': ['Test_ARATH', 'Test_ECOLI', 'Test_HUMAN', 'Test_MOUSE', 'Test_MYCTU', 'Test_YEAST'],
    'Train_ARATH': [0.413, 0.356, 0.337, 0.315, 0.325, 0.351],
    'Train_ECOLI': [0.302, 0.392, 0.321, 0.271, 0.352, 0.364],
    'Train_HUMAN': [0.319, 0.357, 0.419, 0.401, 0.319, 0.390],
    'Train_MOUSE': [0.318, 0.329, 0.410, 0.385, 0.298, 0.350],
    'Train_MYCTU': [0.258, 0.345, 0.289, 0.252, 0.472, 0.308],
    'Train_YEAST': [0.290, 0.350, 0.343, 0.286, 0.323, 0.449]
}

data_bp_Fmax_ss8 = {
    'Test': ['Test_ARATH', 'Test_ECOLI', 'Test_HUMAN', 'Test_MOUSE', 'Test_MYCTU', 'Test_YEAST'],
    'Train_ARATH': [0.452, 0.348, 0.342, 0.318, 0.334, 0.357],
    'Train_ECOLI': [0.295, 0.420, 0.327, 0.276, 0.366, 0.361],
    'Train_HUMAN': [0.315, 0.361, 0.440, 0.412, 0.334, 0.400],
    'Train_MOUSE': [0.324, 0.337, 0.434, 0.420, 0.328, 0.364],
    'Train_MYCTU': [0.282, 0.356, 0.300, 0.257, 0.502, 0.317],
    'Train_YEAST': [0.299, 0.357, 0.358, 0.293, 0.333, 0.474]
}

data_bp_AUPR_aa = {
    'Test': ['Test_ARATH', 'Test_ECOLI', 'Test_HUMAN', 'Test_MOUSE', 'Test_MYCTU', 'Test_YEAST'],
    'Train_ARATH': [0.351, 0.271, 0.259, 0.222, 0.248, 0.283],
    'Train_ECOLI': [0.190, 0.333, 0.229, 0.160, 0.257, 0.289],
    'Train_HUMAN': [0.228, 0.292, 0.374, 0.336, 0.238, 0.336],
    'Train_MOUSE': [0.232, 0.262, 0.361, 0.328, 0.215, 0.288],
    'Train_MYCTU': [0.165, 0.253, 0.195, 0.142, 0.386, 0.235],
    'Train_YEAST': [0.195, 0.296, 0.263, 0.191, 0.236, 0.416]
}

data_bp_AUPR_ss8 = {
    'Test': ['Test_ARATH', 'Test_ECOLI', 'Test_HUMAN', 'Test_MOUSE', 'Test_MYCTU', 'Test_YEAST'],
    'Train_ARATH': [0.390, 0.273, 0.264, 0.229, 0.258, 0.290],
    'Train_ECOLI': [0.197, 0.370, 0.237, 0.176, 0.281, 0.300],
    'Train_HUMAN': [0.230, 0.295, 0.401, 0.346, 0.252, 0.348],
    'Train_MOUSE': [0.238, 0.276, 0.391, 0.366, 0.246, 0.305],
    'Train_MYCTU': [0.173, 0.270, 0.195, 0.141, 0.409, 0.236],
    'Train_YEAST': [0.206, 0.302, 0.287, 0.208, 0.245, 0.444]
}

################################ cc ####################################

data_cc_Fmax_aa = {
    'Test': ['Test_ARATH', 'Test_ECOLI', 'Test_HUMAN', 'Test_MOUSE', 'Test_MYCTU', 'Test_YEAST'],
    'Train_ARATH': [0.751, 0.523, 0.607, 0.611, 0.492, 0.681],
    'Train_ECOLI': [0.511, 0.835, 0.480, 0.488, 0.689, 0.483],
    'Train_HUMAN': [0.706, 0.604, 0.674, 0.632, 0.477, 0.676],
    'Train_MOUSE': [0.707, 0.608, 0.622, 0.644, 0.497, 0.673],
    'Train_MYCTU': [0.538, 0.577, 0.494, 0.468, 0.759, 0.482],
    'Train_YEAST': [0.710, 0.554, 0.603, 0.607, 0.452, 0.744]
}


data_cc_Fmax_ss8 = {
    'Test': ['Test_ARATH', 'Test_ECOLI', 'Test_HUMAN', 'Test_MOUSE', 'Test_MYCTU', 'Test_YEAST'],
    'Train_ARATH': [0.754, 0.523, 0.599, 0.602, 0.500, 0.676],
    'Train_ECOLI': [0.510, 0.817, 0.475, 0.488, 0.700, 0.473],
    'Train_HUMAN': [0.698, 0.589, 0.672, 0.627, 0.500, 0.669],
    'Train_MOUSE': [0.699, 0.594, 0.621, 0.658, 0.521, 0.671],
    'Train_MYCTU': [0.501, 0.626, 0.496, 0.469, 0.761, 0.488],
    'Train_YEAST': [0.702, 0.523, 0.595, 0.601, 0.444, 0.732]
}

data_cc_AUPR_aa = {
    'Test': ['Test_ARATH', 'Test_ECOLI', 'Test_HUMAN', 'Test_MOUSE', 'Test_MYCTU', 'Test_YEAST'],
    'Train_ARATH': [0.791, 0.527, 0.594, 0.578, 0.495, 0.667],
    'Train_ECOLI': [0.448, 0.844, 0.400, 0.418, 0.670, 0.467],
    'Train_HUMAN': [0.695, 0.620, 0.721, 0.640, 0.509, 0.683],
    'Train_MOUSE': [0.701, 0.603, 0.644, 0.670, 0.522, 0.672],
    'Train_MYCTU': [0.412, 0.524, 0.341, 0.351, 0.850, 0.311],
    'Train_YEAST': [0.689, 0.518, 0.591, 0.564, 0.426, 0.771]
}

data_cc_AUPR_ss8 = {
    'Test': ['Test_ARATH', 'Test_ECOLI', 'Test_HUMAN', 'Test_MOUSE', 'Test_MYCTU', 'Test_YEAST'],
    'Train_ARATH': [0.789, 0.541, 0.588, 0.561, 0.499, 0.668],
    'Train_ECOLI': [0.431, 0.827, 0.389, 0.411, 0.668, 0.400],
    'Train_HUMAN': [0.690, 0.575, 0.715, 0.644, 0.504, 0.678],
    'Train_MOUSE': [0.691, 0.580, 0.640, 0.678, 0.536, 0.670],
    'Train_MYCTU': [0.365, 0.560, 0.353, 0.364, 0.847, 0.309],
    'Train_YEAST': [0.685, 0.528, 0.588, 0.563, 0.424, 0.767]
}



################################ mf ####################################

data_mf_Fmax_aa = {
    'Test': ['Test_ARATH', 'Test_ECOLI', 'Test_HUMAN', 'Test_MOUSE', 'Test_MYCTU', 'Test_YEAST'],
    'Train_ARATH': [0.534, 0.338, 0.392, 0.368, 0.322, 0.356],
    'Train_ECOLI': [0.342, 0.367, 0.331, 0.292, 0.348, 0.289],
    'Train_HUMAN': [0.405, 0.342, 0.471, 0.479, 0.325, 0.352],
    'Train_MOUSE': [0.402, 0.322, 0.468, 0.438, 0.307, 0.331],
    'Train_MYCTU': [0.297, 0.319, 0.314, 0.277, 0.369, 0.276],
    'Train_YEAST': [0.398, 0.337, 0.363, 0.345, 0.321, 0.396]
}


data_mf_Fmax_ss8 = {
    'Test': ['Test_ARATH', 'Test_ECOLI', 'Test_HUMAN', 'Test_MOUSE', 'Test_MYCTU', 'Test_YEAST'],
    'Train_ARATH': [0.602, 0.348, 0.409, 0.393, 0.362, 0.399],
    'Train_ECOLI': [0.333, 0.469, 0.335, 0.307, 0.394, 0.302],
    'Train_HUMAN': [0.446, 0.350, 0.515, 0.537, 0.363, 0.394],
    'Train_MOUSE': [0.438, 0.329, 0.504, 0.520, 0.350, 0.376],
    'Train_MYCTU': [0.352, 0.343, 0.338, 0.307, 0.400, 0.287],
    'Train_YEAST': [0.449, 0.333, 0.406, 0.390, 0.345, 0.423]
}

data_mf_AUPR_aa = {
    'Test': ['Test_ARATH', 'Test_ECOLI', 'Test_HUMAN', 'Test_MOUSE', 'Test_MYCTU', 'Test_YEAST'],
    'Train_ARATH': [0.491, 0.258, 0.291, 0.267, 0.257, 0.263],
    'Train_ECOLI': [0.243, 0.305, 0.218, 0.179, 0.280, 0.195],
    'Train_HUMAN': [0.340, 0.273, 0.438, 0.420, 0.266, 0.284],
    'Train_MOUSE': [0.328, 0.249, 0.414, 0.379, 0.243, 0.256],
    'Train_MYCTU': [0.204, 0.233, 0.196, 0.164, 0.272, 0.163],
    'Train_YEAST': [0.322, 0.268, 0.266, 0.248, 0.246, 0.335]
}

data_mf_AUPR_ss8 = {
    'Test': ['Test_ARATH', 'Test_ECOLI', 'Test_HUMAN', 'Test_MOUSE', 'Test_MYCTU', 'Test_YEAST'],
    'Train_ARATH': [0.558, 0.274, 0.322, 0.302, 0.287, 0.315],
    'Train_ECOLI': [0.252, 0.396, 0.243, 0.211, 0.336, 0.220],
    'Train_HUMAN': [0.404, 0.287, 0.478, 0.485, 0.309, 0.338],
    'Train_MOUSE': [0.378, 0.255, 0.452, 0.451, 0.274, 0.305],
    'Train_MYCTU': [0.256, 0.270, 0.224, 0.203, 0.298, 0.189],
    'Train_YEAST': [0.388, 0.268, 0.330, 0.310, 0.277, 0.372]
}




