import os
import math
import numpy as np


NBLADE = 3
RPM = 7000
BPF = (NBLADE * RPM) / 60.0

current_path = os.path.abspath(__file__)
dir_current = os.path.dirname(current_path)

EXP_dir = os.path.join(dir_current, "EXP_ISAE")
SONAR_dir = os.path.join(dir_current, "SONAR_ISAE")
BPF_dir = os.path.join(dir_current, "1stBPF")

EXP_fname = os.path.join(EXP_dir, 'exp_%sb_%sRPM_narrow.dat' % (
    NBLADE, RPM))
SONAR_fname = os.path.join(SONAR_dir, '%sb_%sRPM_tonal.dat' % (
    NBLADE, RPM))
out_fname = os.path.join(BPF_dir, '%sb_%sRPM_1stBPF.dat' % (
    NBLADE, RPM))


exp_freq_sets = []
exp_data_sets = []

with open(EXP_fname, 'r') as exp_file:
    lines = exp_file.readlines()
    lines = lines[1:]

    current_exp_freq = []
    current_exp_data = []

    for line in lines:
        line = line.strip()
        if line.startswith('zone t='):
            if current_exp_freq:
                exp_freq_sets.append(np.array(current_exp_freq))
                exp_data_sets.append(np.array(current_exp_data))
                current_exp_freq = []
                current_exp_data = []

        else:
            parts = line.split()
            if len(parts) >= 2:
                freq = float(parts[0])
                data = float(parts[1])
                current_exp_freq.append(freq)
                current_exp_data.append(data)

    if current_exp_freq:
        exp_freq_sets.append(np.array(current_exp_freq))
        exp_data_sets.append(np.array(current_exp_data))

sonar_freq_sets = []
sonar_data_sets = []

with open(SONAR_fname, 'r') as sonar_file:
    lines = sonar_file.readlines()
    lines = lines[1:]

    current_sonar_freq = []
    current_sonar_data = []

    for line in lines:
        line = line.strip()
        if line.startswith('zone t='):
            if current_sonar_freq:
                sonar_freq_sets.append(np.array(current_sonar_freq))
                sonar_data_sets.append(np.array(current_sonar_data))
                current_sonar_freq = []
                current_sonar_data = []
        else:
            parts = line.split()
            if len(parts) >= 4:
                freq = float(parts[0])
                data = float(parts[3])
                current_sonar_freq.append(freq)
                current_sonar_data.append(data)

    if current_sonar_freq:
        sonar_freq_sets.append(np.array(current_sonar_freq))
        sonar_data_sets.append(np.array(current_sonar_data))

exp_values_at_BPF = []
sonar_values_at_BPF = []

for exp_freq, exp_data in zip(exp_freq_sets, exp_data_sets):
    idx = (np.abs(exp_freq - BPF)).argmin()
    exp_values_at_BPF.append(exp_data[idx])

for sonar_freq, sonar_data in zip(sonar_freq_sets, sonar_data_sets):
    idx = (np.abs(sonar_freq - BPF)).argmin()
    # 인덱스 범위 설정 (-5 ~ +5)
    start_idx = max(idx - 5, 0)
    # +1은 슬라이싱 범위에 포함시키기 위함
    end_idx = min(idx + 5 + 1, len(sonar_data))
    # 해당 범위 내의 데이터 추출
    data_window = sonar_data[start_idx:end_idx]

    max_value = np.max(data_window)
    sonar_values_at_BPF.append(max_value)


print("EXP 파일에서 BPF에 가장 가까운 값:")
print(exp_values_at_BPF)

print("SONAR 파일에서 BPF에 가장 가까운 값:")
print(sonar_values_at_BPF)

map_fname = os.path.join(dir_current, 'map.dat')
map_data = np.genfromtxt(map_fname, delimiter='', skip_header=1)
map_x = map_data[:, 0]
map_y = map_data[:, 1]
map_z = map_data[:, 2]

hxy = []
r = []
el = []
az = []
for i in range(len(map_x)):
    hxy_temp = math.hypot(map_x[i], map_y[i])
    r_temp = math.hypot(hxy_temp, map_z[i])
    el_temp = math.atan2(map_z[i], hxy_temp)
    az_temp = math.atan2(map_y[i], map_x[i])

    hxy.append(hxy_temp)
    r.append(r_temp)
    el.append(el_temp)
    az.append(az_temp)

el = np.rad2deg(el)

with open(out_fname, 'w') as fo:
    fo.write('variables="Directivity[deg]",'
             '"Exp. SPL at BPF<sup>1st</sup> [dB]",'
             '"Present SPL at BPF<sup>1st</sup> [dB]"\n')
    fo.write(f'zone t="{NBLADE}b_{RPM}RPM"\n')
    for i in range(len(el)):
        fo.write(f'{el[i]}\t'
                 f'{exp_values_at_BPF[i]}\t'
                 f'{sonar_values_at_BPF[i]}\n')
