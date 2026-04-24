---
{"dg-publish":true,"permalink":"/计算相关/脚本的创建与使用（示例）/过渡态路径整合(pathway)/","tags":["服务器","脚本","DFT计算","代码"],"noteIcon":"","dg-note-properties":{"share_link":"https://share.note.sx/lhq6a348#fJwYzRCv/IdktT7mTdfpo2ePKpQBQCieasZ4swrKJW4","share_updated":"2026-04-25T04:05:16+08:00","tags":["服务器","脚本","DFT计算","代码"]}}
---

知识库入口：[[第一性原理概念/李院士知识库\|李院士知识库]]

### 一、 过渡态路径整合脚本：`pathway`

> [!note] 
> <font color="#f79646">使用该脚本可将过渡态计算（NEB）中 00 到 05 文件夹中发生解离的 Li 离子坐标自动识别并整合到一个新文件 Li-POSCAR 中，以便在 VESTA 中直观观察完整的解离轨迹。</font> 

#### 第一步：创建脚本

在 `~/my_scripts/` 目录下创建该文件（注意这里直接创建无后缀的文件，实现“文件名=命令名”）：

Bash

```
nano ~/my_scripts/pathway
```

> [!attention] 
>   若没有该文件夹则先创建一个：
>   

Bash

```
mkdir -p ~/my_scripts
```

#### 第二步：粘贴以下代码

```
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np

def read_poscar(filepath):
    """解析 POSCAR/CONTCAR 文件并提取晶格、元素及坐标信息"""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # 提取缩放系数和晶格矢量 (转换为真实笛卡尔矩阵)
    scale = float(lines[1].strip())
    lattice = np.array([list(map(float, line.split())) for line in lines[2:5]]) * scale

    # 提取元素种类和对应的原子数量 (VASP 5+ 标准格式)
    elements = lines[5].split()
    counts = list(map(int, lines[6].split()))

    # 判断是否存在 Selective dynamics (选择性动力学) 标识
    coord_start = 7
    has_selective = False
    if lines[7].strip()[0].lower() == 's':
        coord_start = 8
        has_selective = True

    # 读取坐标类型 (Direct or Cartesian)
    coord_type = lines[coord_start].strip()
    coord_start += 1 # 坐标数据正式开始的行号

    # 提取纯坐标数据 (只取前三列 x, y, z)
    coords = []
    for i in range(sum(counts)):
        parts = lines[coord_start + i].split()
        coords.append([float(x) for x in parts[:3]])

    return {
        'lines': lines,
        'lattice': lattice,
        'elements': elements,
        'counts': counts,
        'coord_start': coord_start,
        'has_selective': has_selective,
        'coords': np.array(coords)
    }

def main():
    print("正在解析初末态结构...")
    # 1. 读取初态和末态结构
    s00 = read_poscar(os.path.join('00', 'POSCAR'))
    s05 = read_poscar(os.path.join('05', 'POSCAR'))

    # 2. 定位 Li 原子的坐标索引范围
    if 'Li' not in s00['elements']:
        print("错误：在 00/POSCAR 中未找到 'Li' 元素，请检查文件格式是否为 VASP5 标准！")
        return
        
    li_index = s00['elements'].index('Li')
    li_start = sum(s00['counts'][:li_index])
    li_end = li_start + s00['counts'][li_index]

    # 3. 自动识别发生迁移的 Li 离子 (引入 PBC 校正)
    max_dist = 0.0
    target_atom_idx = -1

    for i in range(li_start, li_end):
        d_direct = s05['coords'][i] - s00['coords'][i]
        d_direct -= np.round(d_direct) 
        d_cart = np.dot(d_direct, s00['lattice'])
        dist = np.linalg.norm(d_cart)

        if dist > max_dist:
            max_dist = dist
            target_atom_idx = i

    print(f"成功锁定！发生迁移的 Li 离子序号为: 第 {target_atom_idx + 1} 号原子 (相对位移量: {max_dist:.4f} 埃)")

    # 4. 提取中间态和末态的轨迹坐标
    print("正在提取迁移轨迹...")
    trajectory_coords = []
    folders_files = [('01', 'CONTCAR'), ('02', 'CONTCAR'), ('03', 'CONTCAR'), ('04', 'CONTCAR'), ('05', 'POSCAR')]

    for folder, file in folders_files:
        path = os.path.join(folder, file)
        s = read_poscar(path)
        trajectory_coords.append(s['coords'][target_atom_idx])

    # 5. 重组并写入新的 Li-POSCAR
    print("正在生成包含轨迹的 Li-POSCAR...")
    with open('Li-POSCAR', 'w') as f:
        for i in range(6): 
            f.write(s00['lines'][i])
            
        new_counts = s00['counts'].copy()
        new_counts[li_index] += 5
        f.write("  " + "  ".join(map(str, new_counts)) + "\n")

        if s00['has_selective']:
            f.write(s00['lines'][s00['coord_start'] - 2])
        f.write(s00['lines'][s00['coord_start'] - 1])

        for i in range(len(s00['coords'])):
            f.write(s00['lines'][s00['coord_start'] + i])
            if i == li_end - 1:
                for idx, traj_coord in enumerate(trajectory_coords):
                    flag = "  F  F  F" if s00['has_selective'] else ""
                    f.write(f"  {traj_coord[0]:.6f}  {traj_coord[1]:.6f}  {traj_coord[2]:.6f}{flag}  # Traj_Image_{idx+1:02d}\n")

    print("\n任务完成！已生成文件：Li-POSCAR。您可以将其下载并在 VESTA 中打开查看。")

if __name__ == '__main__':
    main()
```

#### 第三步：给文件添加可执行权限（关键步骤，缺了必报错）

Linux 系统默认不给文本文件执行权限，必须手动添加：

Bash

```
chmod +x ~/my_scripts/pathway
```

#### 第四步：把脚本库加入系统全局环境变量 PATH（一劳永逸）

让系统在任何目录都能找到你的 `pathway` 命令，执行以下命令写入配置：

Bash (若该目录之前已经加入过[PATH](改变脚本名字.md#^39e1b2)则跳过该步骤)

```
echo 'export PATH="$HOME/my_scripts:$PATH"' >> ~/.bashrc
```

#### 第五步：刷新配置生效

Bash

```
source ~/.bashrc
```

 > [!caution] 
>  脚本文件格式应为 `UTF-8`，若在 Windows 下编辑过，可能变为 GBK 编码导致乱码报错，须输入以下命令转换格式：
 
Bash

``` 
iconv -f GBK -t UTF-8 ~/my_scripts/pathway -o ~/my_scripts/pathway
```
 _(这行代码的意思是：把该文件从 GBK 编码读取，转换为 UTF-8，然后重新输出并覆盖原文件)_

#服务器 #脚本 #DFT计算 #代码