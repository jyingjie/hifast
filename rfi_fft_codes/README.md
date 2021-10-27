# parameters meaning

Author: Xu Chen, NAOC, 2021.Jun ~ Oct

1. 先用cli_baseline多项式去基线
2. cli_mark_tRFI找时域RFI
3. cli_baseline_fft去驻波
4. cli_baseline再去一次基线
5. cli_multi找偏振很大的RFI、坐标系修正
6. （可以合并通道）
然后成图

有bug请联系stellarxu@qq.com或者xuchen21@mails.ucas.ac.cn

to do list:
0.92 MHz 周期脉冲RFI, FAST基地表示很难消除

参数具体设置参见notebook
(推荐)表示我推荐使用的方法。
(测试)表示最好用一个波束在notebook里测试，可以先试默认值
(默认)一般可以使用默认值



## cli_mark_tRFI标记rfi

### 时域

* --time_rfi：标记时域上突然出现的RFI，根据占据的频率长度，分为short-freq/long-freq。
    对该频率区间内的所有谱线做频率方向的平均后，画出一维的图。如果有--lf或--sf，可以省略--time_rfi

* --sf
* --sf_frange: 一个典型的1～2MHz宽时域RFI，经常出现在1380～1382MHz。
* --sf_file: 或者可以指定一个二维数组的npy文件路径，在里面循环找。
* --sf_frange_step(推荐):如果sf_frange和sf_file不指定，将以这个step(MHz)在频率上循环，对每个区间找一次时域RFI
```
sepname='xxxxxxx_path'
python cli_mark_tRFI.py $subname --outdir ./data \
        --mw_frange 1420.2 1420.45 --lf_sepname 'input_subname' \
        --lf --lf_frange 1385 1460 --lf_times 2 --lf_thr 0 --lf_rfi_last 50 --lf_ext 50 \
        --sf --sf_frange_step 30 --sf_times 3 --sf_thr 1 --sf_rfi_last 10 --sf_T_thr_times 2.5 \
        --plot --vmin_max -.05 .05 -f  --outdir $outdir || exit 1 
```

* --sf_times(测试): 大于各谱线之中，中值的*倍，初步认为异常
* --sf_thr(测试): 边缘处与附近值的差是中值的倍数，找出时间方向的边缘处陡峭的
* --sf_rfi_last(测试): 持续出现的谱线数
* --sf_T_thr_times(测试): 选定区域后，大于温度阈值的被标记
* --sf_ext(测试): 选定区域后，扩展边缘，一般为0

* --lf
* --lf_beams: 已知出现大卫星干扰的波束，其他不搜寻。如果不设置就对所有波束搜寻。
* --lf_frange: 大范围频率，[1400,1460]或者避开银河系[1421,1460]，RFI集中的频率
* --lf_times(测试): 大于各谱线之中，中值*倍，初步认为异常.
* --lf_thr(默认): 不需要边缘陡峭，所以设为0
* --lf_rfi_last(测试): 持续出现的谱线数
* --lf_ext(测试): 选定区域后，扩展边缘
* --lf_sepfile: 如果输入的是'input_subname'，None则从sub文件找lf；否则从sep文件找。
    自动搜索的路径顺序是上一级是/sep/、sub文件所在路径或者指定的路径

输出文件名含tr


输出的文件中```dict_out['is_rfi']```为二维的mask

### 其他参数
* --mw_frange(测试): 银河系(或者M31，M33都)存在的区域
* --ylim: 画谱线的上下限 
* --flux: 流量定标
* -c, --cali_fname: 定标源
* --keep_polar: 保留偏振
* --keep_rfi: 保留RFI
* --plot: 画图
* --vmin_max: 瀑布图的上下限

2021年7月28号前的数据有**8.1MHz周期RFI**，可以用以下方法mask为nan
## cli_mark_pdRFI标记rfi

### 频域周期RFI
* --period_rfi: 标记频率上周期变化的RFI

* --s_method_freq,--s_method_t:平滑方法，'gaussian', 'boxcar', 'median', 'None'
 
  --s_sigma_freq,--s_sigma_t: sigma同前面定义

* --rfi_thr: 大于rms rfi_thr倍的会被记为rfi
* --rms_frange: 用于计算rms的干净频率区间
* --rms_sigma: 如果区间有干扰，还是用该频率减去一个高斯平滑，得到准确的rms。此为高斯滤波的sigma。
* --mw_frange: 一个大致可能存在信号的区间，先保护起来，后面不会用它预估rfi位置
* --rfi_groups: rfi的形态呈现为8.1MHz间隔的几组，分类进行更精确的拟合。'two groups','three groups','all'
* --rfi_width_lim: 粗找rfi区域时，其宽度应该大于的通道数。太细的可能是驻波，可能是点源，排除
* --ext_sec: 粗找到rfi后，扩充边界的通道数
* --freq_thr: 对粗测的rfi中心分类时，可以偏离预计中心一段频率。
* --freq_step: 对粗测的rfi中心分类时，预计中心的间隔大约为8.1MHz。

分类后进行直线拟合，得到rfi的精确分布。如果大于rfi_thr的区域落在rfi精确理论频率上，则被标记为真rfi，可以防止高流量的源被标记。

输出文件名含pdr

### 标记

* --ext_edge: 结果扩展边缘的通道数
* --mask_thr: 只mask超过rms的几倍。很小的rfi会保留。
* --mask_RFI_method: 'fixed freq' 对于小rfi，固定mask宽度
* --mask_all_theory: 是否理论的全mask，这样会非常干净(输出文件名含strict)
      --freq_from_theory: 很小的去掉多少频率(固定宽度)
  
* --mask_RFI_method: '2 sides' 从中心向两边按step循环，确定mask边界，比较慢
      --small_rfi_times: 小于RMS这个倍数的不标记
      --chan_step: step通道数
      --freq_from_theory: mask rfi最大的宽度
例如
```
python cli_markRFI.py $subname --outdir ./data \
        --period_rfi \
        --s_method_freq gaussian --s_sigma_freq 3 --s_method_t boxcar --s_sigma_t 7 \
        --rfi_thr 14 --rms_frange 1390 1400 --mw_frange 1419 1425 --rfi_groups 'two groups' \
        --rfi_width_lim 20 --ext_sec 20 --freq_thr .7 --freq_step 8.1 \
        --mask_RFI_method 'fixed freq' --ext_edge 0 --mask_thr 16 --freq_from_theory .5 \
        --plot --ylim -1 5 --vmin_max -.05 .05 -f --keep_polar --time_coherent_per .90 || exit 1 
```
  
* --time_coherent_per: 如果该频率，例如 90% 都被标记了，那么整条都被标记为RFI。不怎么用。
* --save_rfi_list 是否保存各条谱线周期rfi理论值，```dict_out['rfi_list']```

输出的文件中```dict_out['is_rfi']```为二维的mask

需要注意，很难完全mask干净的，所以**仅供参考**

### 其他参数

* --ylim: 画谱线的上下限 
* --flux: 流量定标
* -c, --cali_fname: 定标源
* --keep_polar: 保留偏振
* --keep_rfi: 保留RFI
* --plot: 画图
* --vmin_max: 瀑布图的上下限



## cli_baseline_fft去驻波

### 对RFI
* --rfi_method: 防止大rfi/强源干扰fft，简化版只提供'near ripple'方法

* --rms_sigma,--rms_frange(测试): 如果区间有干扰，还是用该频率减去一个高斯平滑，得到准确的rms。此为高斯滤波的sigma。
* -rfi, --rfi_fname: 指定标记rfi的文件。无则自动到输入文件的上一级名为/rfi/的文件夹里找含-tr的文件，找不到则在当前路径找。

* --mw_frange(测试): 银河系(或者M31，M33都)存在的区域.如果两个大星系间隔较远，可以不给定区间或者只给定MW
* --frange: 希望得到结果的频率区间。无则默认全长
* --fft_frange: 用来fft的频率区间，一定>=frange。会根据数据的边缘自动决定向哪个方向做频率方向的周期延拓，以保证目标频率的效果较好

'near ripple'使用高流量处附近的驻波替代。例如
    ```
    fname = subnames[1]
    outdir = './data'
    sepname='xxxxxxx_path'
    python cli_baseline_fft.py $fname \
        -sep $sepname  --outdir $outdir  \
        --mw_frange  1420.2 1420.45 --rms_sigma 6 --rms_frange 1390 1400 \
        --rfi_method 'near ripple' --times_lower_thr 4.5 \
        --rfi_width_lim 10 --ext_sec 5 --ext_freq 1.3 \
        --fft_method rfft --chan_wide 7 --chan_narr 3 \
        --amp_thr_mean_factor 1.1 --amp_thr_factor 1.4 --choose_method 'all' \
        --rip_base --rip_1mhz --rip_2mhz --rip_0_04mhz --fft_ylim -5 160 \
        --plot --ylim -1 .5 --vmin_max -.05 .05 -f --fill_rfi rfi --keep_polar || exit 1 
    ```
    * 只需要已知时域rfi的文件，即*-tr*。    
      如果没有时域RFI，设置-rfi 'none'即可。
    * --times_lower_thr(测试):大于RMS这么多倍的会被替代
    * --rfi_width_lim(测试): 要替换的rfi/源，频率上的宽度应该大于一阈值(通道数)
    * --ext_sec(测试): (接上面的)并向两边扩展通道数
    * --ext_freq(默认): 先扩展边缘(MHz)，再寻找强流量两边的最低点，先使用左/右边的一段更干净一点的代替强流量处
    
    * 输出文件名含fft_blde.


### fft去驻波

* --fft_method(默认): 目前只提供rfft，mean和median的方法见‘cli_baseline_mean去基线‘

* --amp_thr_mean_factor(测试): fourier空间的平均振幅大于阈值的模，会被识别为已知的几种驻波。这里输入的是中值的倍数
* --amp_thr_factor(测试): fourier空间的每一谱线的振幅大于阈值，会被选中为驻波一部分被去除。这里输入的是中值的倍数
* --chan_wide(测试): 距离驻波的模的中心左右各(2 * chann - 1)个通道数的模，作为驻波的一部分被选中，这里为由于fourier空间中峰比较宽，多选一些
  --chan_narr(测试): 同理，窄一些的
* --choose_method(测试): 选择fourier模的方法，'all'或者'interpolate'
  
* --rip_base: 去除基频(实空间的常数)
* --rip_1mhz: fft去除1.08mhz 驻波，单镜都有的
* --rip_2mhz: fft去除1.92mhz 驻波，只有M06 yy有
* --rip_0_04mhz: fft去除0.039mhz 驻波，很弱的偏振驻波，由于光纤反射

* --fft_ylim(测试): 画fft模的上下限 

* --fill_rfi: 对输入的rfi区域，可以填nan或者保持原样。对应'nan','rfi'。需要注意只有输入对应RFI文件才能填nan。

### 其他参数
* smooth 参数同前面cli_baseline。一般不用平滑
* --ylim: 画谱线的上下限 
* --flux: 流量定标
* -c, --cali_fname: 定标源
* --keep_polar: 保留偏振
* --plot: 画图
* --vmin_max: 瀑布图的上下限
* --no_radec: 没有radec文件的话
* -sep, --sep_fname: 只分离光谱温度定标的文件。无则自动到输入文件的上一级名为/sep/的文件夹里找含-tr的文件，找不到则在当前路径找。

输出的hdf5包含去除驻波的数组(T或者flux)(单独的驻波```dict_out['ripple']```注释掉了)


## cli_baseline_mean去基线

* --sub_method: mean或median，取平均或中值前后若干条谱线作为基线
* --nspec: 平均谱线条数，默认200条。过少容易损失信号，过多容易驻波变化了。
* --func: 自己写的'iter'或卷积的'smooth'。默认'iter'，但是多了会比较慢。

```
python cli_baseline_mean.py $fname --outdir $outdir \
        --sub_method mean --nspec 200 --frange 1330 1360 \
        --plot --ylim -1 .5 --vmin_max -.05 .05 -f --keep_polar --fill_rfi nan || exit 1 
```

输出文件名含mean_bld/med_bld

需要注意，由于不同谱线基线不同，还需要使用cli_baseline.py多项式去一次基线。


## cli_RFIshape利用同一谱线不同频率的8.MHz周期RFI波形

利用同一谱线不同频率的RFI波形，分离污染了M31的RFI。新数据没有此RFI。

eg.
```
python cli_RFIshape.py $subname --outdir ./data \
        --m31_frange 1421 1425 --m31_times 1.5 --m31_thr 0 --m31_rfi_last 50 --m31_ext 3 \
        --no_m31_frange 1420 1423 \
        --rfi_thr 2 --rms_sigma 6 --rms_frange 1390 1400 --mw_frange 1419 1425 --rfi_groups 'three groups' \
        --rfi_width_lim 3 --ext_sec 4 --freq_thr .7 --freq_step 8.1 \
        --mask_RFI_method 'fixed freq' --ext_edge 0 --mask_thr 3 --freq_from_theory .5 \
        --plot --ylim -1 5 -f || exit 1 
```

* 其中``` --m31_frange 1421 1425 --m31_times 1.5 --m31_thr 0 --m31_rfi_last 50 --m31_ext 3 ```与cli_markRFI的lf含义相同，用来确定哪些谱线有M31。

* 下面这些与cli_markRFI的找周期RFI并mask的含义相同，用来收集其他频率的RFI波形
```   --rfi_thr 2 --rms_sigma 6 --rms_frange 1390 1400 --mw_frange 1419 1425 --rfi_groups 'three groups' \
        --rfi_width_lim 3 --ext_sec 4 --freq_thr .7 --freq_step 8.1 \
        --mask_RFI_method 'fixed freq' --ext_edge 0 --mask_thr 3 --freq_from_theory .5 \
```
* 用这些波形对齐峰值频率，插值成相同长度后取中值，用以代替不同频率的RFI
* --no_m31_frange: 周期RFI的区间里包含了M31，这里为了选出RFI的峰值用于对齐，需要把不使用的M31、MW设为False。
* --only_M31: 只把M31邻近处减去插值波形。指no_m31_frange及其左右3MHz
* 需要注意，很难完全减干净的，所以**仅供参考**
