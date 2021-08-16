# parameters meaning

Author: Xu Chen, NAOC, 2021.Jun ~ Aug

先用cli_baseline多项式去基线，然后的顺序看g15_pipe.sh

有bug请联系stellarxu@qq.com

to do list:
0.92 MHz 周期脉冲RFI


## cli_markRFI标记rfi

eg:
```
python cli_markRFI.py $subname --outdir ./data \
        --lf --sf --lf_beams ['05','06','13'] \
        --sf_frange 1380 1382 --sf_times 10 --sf_thr 10 --sf_rfi_last 20 --sf_T_thr .3 \
        --lf_frange 1400 1450 --lf_times 1.5 --lf_thr 0 --lf_rfi_last 50 --lf_ext 10 \
        --period_rfi \
        --s_method_freq gaussian --s_sigma_freq 3 --s_method_t boxcar --s_sigma_t 7 \
        --rfi_thr 3 --rms_frange 1400 1403 --mw_frange 1419 1425 --rfi_groups 'two groups' \
        --rfi_width_lim 20 --ext_sec 20 --freq_thr .3 --freq_step 8.1 \
        --mask_RFI_method 'fixed freq' --ext_edge 3 --mask_thr 2 --mask_all_theory --freq_from_theory .5 \
        --plot -f --save_sf || exit 1 
```
### 时域

* --time_rfi：标记时域上突然出现的RFI，根据占据的频率长度，分为short-freq/long-freq。对该频率区间内的所有谱线做频率方向的平均后，画出一维的图。

* --sf
* --sf_frange: 一个典型的1～2MHz宽时域RFI，经常出现在1380～1382MHz。
* --sf_file: 或者可以指定一个二维数组的npy文件路径，在里面循环找，例如以下是只找时域RFI:
```
time_rfi_file = '/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/time_rfi_list.npy'
python cli_markRFI.py $subname --outdir ./data \
        --sf \
        --sf_file $time_rfi_file --sf_times 1.8 --sf_thr 0 --sf_rfi_last 20 --sf_T_thr .3 \
        --plot --ylim -1 5 -f  || exit 1 
        
INFO: Looking for short-freq time RFI in [1378 1384] ... [markRFI]
Found :D
INFO: Looking for short-freq time RFI in [1310 1315] ... [markRFI]
Found :D
INFO: Looking for short-freq time RFI in [1320 1325] ... [markRFI]
Found :D
Finish
```

* --sf_times: 大于各谱线之中，中值的10倍，初步认为异常
* --sf_thr: 边缘处与附近值的差是中值的倍数，找出时间方向的边缘处陡峭的
* --sf_rfi_last: 持续出现的谱线数
* --sf_T_thr: 选定区域后，大于温度阈值的被标记
* --sf_ext: 选定区域后，扩展边缘

* --lf
* --lf_beams: 已知出现大卫星干扰的波束，其他不搜寻。如果不设置就对所有波束搜寻。
* --lf_frange: 大范围频率
* --lf_times: 大于各谱线之中，中值1.5倍，初步认为异常.
* --lf_thr: 不需要边缘陡峭，所以设为0
* --lf_rfi_last: 持续出现的谱线数
* --lf_ext: 选定区域后，扩展边缘

* --save_sf: 将sf单独保存为```dict_out['short_rfi']```

输出文件名含tr

### 频域
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
        --sf --lf --lf_beams ['05','06','13'] \
        --sf_frange 1380 1382 --sf_times 10 --sf_thr 10 --sf_rfi_last 20 --sf_T_thr .3 \
        --lf_frange 1400 1450 --lf_times 1.5 --lf_thr 0 --lf_rfi_last 50 --lf_ext 10 \
        --period_rfi \
        --s_method_freq gaussian --s_sigma_freq 3 --s_method_t boxcar --s_sigma_t 7 \
        --rfi_thr 3 --rms_frange 1400 1403 --mw_frange 1419 1425 --rfi_groups 'two groups' \
        --rfi_width_lim 20 --ext_sec 20 --freq_thr .3 --freq_step 8.1 \
        --mask_RFI_method '2 sides' --ext_edge 10 --mask_thr 2 --freq_from_theory 3\
        --small_rfi_times 2 --chan_step 5 \
        --plot -f  --save_tf|| exit 1 
  ```
  
* --time_coherent_per: 如果该频率，例如90%都被标记了，那么整条都被标记为RFI。不怎么用。
* --save_rfi_list: 是否保存各条谱线周期rfi理论值，```dict_out['rfi_list']```
* --save_tf:是否保存 time rfi，```dict_out['time_rfi']```

输出的文件中```dict_out['is_rfi']```为二维的mask

### 其他参数

* --ylim: 画谱线的上下限 
* --flux: 流量定标
* -c, --cali_fname: 定标源
* --keep_polar: 保留偏振
* --keep_rfi: 保留RFI
* --plot: 画图



## cli_baseline_fft去驻波

### 对RFI
* --rfi_method: 防止大rfi/强源干扰fft，提供7种方法：'subtract trpdr',
                      'subtract tr','subtract rfi', 'lower','set zeros','set noise','near ripple'
  推荐'near ripple'

* --rms_sigma: 如果区间有干扰，还是用该频率减去一个高斯平滑，得到准确的rms。此为高斯滤波的sigma。
* -rfi, --rfi_fname: 指定标记rfi的文件。无则自动到输入文件的上一级名为/rfi/的文件夹里找含-tr的文件，找不到则在当前路径找。

* --mw_frange: 银河系(或者M31，M33都)存在的区域

1. 'near ripple'使用高流量处附近的驻波替代。例如
    ```
    python cli_baseline_fft.py $subname --outdir $outdir \
        -sep $sepname   -rfi 'none' \
        --mw_frange 1420.2 1420.55 --rms_sigma 6 --rms_frange 1390 1400 \
        --rfi_method 'near ripple' --times_lower_thr 7 \
        --rfi_width_lim 15 --ext_sec 20 --ext_freq 1.3 \
        \
        --fft_method rfft --chan_wide 5 --chan_narr 3 \
        --amp_thr_mean_factor 1.05 --amp_thr_factor 1.4 \
        --rip_base --rip_1mhz --rip_2mhz --rip_0_04mhz --fft_ylim -5 130 \
        --plot --ylim -1 .5 --vmin_max -.15 .15 -f --fill_rfi nan --keep_polar || exit 1 
    ```
    * 只需要已知时域rfi的文件，即*-tr*。    
      如果没有时域RFI，设置-rfi 'none'即可。
    * --times_lower_thr:大于RMS这么多倍的会被替代
    * --rfi_width_lim,--ext_sec: 以上处理的范围，宽度阈值和扩展通道数(同周期rfi里的含义)
    * --ext_freq: 扩展边缘(MHz)，再寻找强流量两边的最低点，先使用左/右边的一段代替强流量处
    
    * 输出文件名含fft_blde.

2. 'subtract trpdr'用原始值减去滤波值。例如
    ```
    python cli_baseline_fft.py $subname --outdir ./data \
        --mw_frange 1420 1421 --rms_sigma 6 --rms_frange 1390 1400 \
        --rfi_method 'subtract trpdr' --sg_window 1.0 --sg_polyorder 7 --gauss_sigma 1 \
        # 这里略去相同的fft参数
    ```
    * 需要已知周期rfi和时域rfi的文件，即*-tr_pdr*
    * --sg_window --sg_polyorder: savgol_filter的窗口大小(MHz)和阶数. 存在周期rfi区域减去sg滤波值替代。存在时域RFI区域除以sg滤波值替代。
    * --gauss_sigma: 银河系存在的地方临时用一个高斯平滑减下去
    * 最后会除去所有超过3$\sigma$的强源或者异常点。
    
    * 输出文件名含fft_bldp.
    
3. 'subtract tr'用原始值减去滤波值。例如
    ```
    python cli_baseline_fft.py $subname --outdir ./data \
        --mw_frange 1420 1421 --rms_sigma 6 --rms_frange 1390 1400 \
        --rfi_method 'subtract tr' --sg_window 1.0 --sg_polyorder 7 --gauss_sigma 1 \
        --times_lower_thr 1.7   \
        # 这里略去相同的fft参数
    ```
    * 只需要已知时域rfi的文件，即*-tr*
    * --times_lower_thr: 大于RMS这么多倍的会被处理，低一些会干净。所以不需要已知pdr。
    * --sg_window --sg_polyorder: savgol_filter的窗口大小(MHz)和阶数.存在减去sg滤波值替代
    * --gauss_sigma: 银河系存在的地方临时用一个高斯平滑减下去

    * 输出文件名含fft_bldt.
    
 4. 'subtract rfi'用原始值减去滤波值。例如
    ```
    python cli_baseline_fft.py $subname --outdir ./data \
        --mw_frange 1420 1421 --rms_sigma 6 --rms_frange 1390 1400 \
        --rfi_method 'subtract rfi' --sg_window 1.0 --sg_polyorder 7 \
        --times_lower 1.0e4 --times_lower_thr 2   \
        # 这里略去相同的fft参数
    ```
    * 需要已知周期rfi的文件，即*pdr*
    * --sg_window --sg_polyorder: savgol_filter的窗口大小(MHz)和阶数.存在减去sg滤波值替代
    * --times_lower_thr: 减去滤波值后，大于RMS这么多倍的会被降低
    * --times_lower: 把高流量压低，比如原来的1/1e4
    
    * 输出文件名含fft_bldr.
 
5. 'lower'降低高流量。不需要rfi文件。例如
    ```
    python cli_baseline_fft.py $subname  --rms_frange 1390 1400 \
        --rfi_method 'lower' --times_lower 1.0e4 --times_lower_thr 2 \
        -# 这里略去相同的fft参数
    ```
    * --times_lower_thr:大于RMS这么多倍的会被降低
    * --times_lower: 把高流量压低，比如原来的1/1e4
    
    * 输出文件名含fft_bldl.
    
6.7. 'set zeros','set noise'高流量全部置零/标准差为RMS的正态分布噪声。不需要rfi文件。例如
    ```
    python cli_baseline_fft.py $subname --outdir ./data --rms_frange 1390 1400 \
        --rfi_method 'set zeros'  --times_lower_thr 2 \
        # 这里略去相同的fft参数
    ```
    * --times_lower_thr:大于RMS这么多倍的会被降低
    * 输出文件名含fft_bldz/fft_bldn.
    


### fft去驻波

* --fft_method: 目前只提供rfft，mean的方法见‘cli_baseline_mean去基线‘

* --amp_thr_mean_factor: fourier空间的平均振幅大于阈值的模，会被识别为已知的几种驻波。这里输入的是中值的倍数
* --amp_thr_factor: fourier空间的每一谱线的振幅大于阈值，会被选中为驻波一部分被去除。这里输入的是中值的倍数
* --chan_wide: 距离驻波的模的中心左右各(2 * chann - 1)个通道数的模，作为驻波的一部分被选中，这里为由于fourier空间中峰比较宽，多选一些
  --chan_narr: 同理，窄一些的
  
* --rip_base: 去除基频(实空间的常数)
* --rip_1mhz: fft去除1.08mhz 驻波，单镜都有的
* --rip_2mhz: fft去除1.92mhz 驻波，只有M06 yy有
* --rip_0_04mhz: fft去除0.039mhz 驻波，很弱的偏振驻波，由于光纤反射

* --rfi_8mhz: fft去除一部分8.1mhz rfi，2021.7.28后的新数据做好了压缩机电源屏蔽，已经没有了
* --rfi_8mhz_step: 残余的8.1MHzRFI，在fourier空间间隔为1/16.2$\mu$s左右

* --fft_ylim: 画fft模的上下限 

* --fill_rfi: 对输入的rfi区域，可以填nan或者保持原样。对应'nan','rfi'。需要注意只有输入对应RFI文件才能填nan。

### 其他参数
* smooth 参数同前面
* --ylim: 画谱线的上下限 
* --flux: 流量定标
* -c, --cali_fname: 定标源
* --keep_polar: 保留偏振
* --plot: 画图
* --no_radec: 没有radec文件的话
* -sep, --sep_fname: 只分离光谱温度定标的文件。无则自动到输入文件的上一级名为/sep/的文件夹里找含-tr的文件，找不到则在当前路径找。

输出的hdf5包含去除驻波的数组(T或者flux)(单独的驻波```dict_out['ripple']```注释掉了)


## cli_baseline_mean去基线

* --sub_method: mean，平均前后若干条谱线作为基线
* --nspec: 平均谱线条数，默认10条。过少容易损失信号，过多容易驻波变化了。

输出文件名含mean_bld.

需要注意，由于不同谱线基线不同，还需要使用cli_baseline.py多项式去一次基线。


## cli_RFIshape利用同一谱线不同频率的RFI波形

利用同一谱线不同频率的RFI波形，分离污染了M31的RFI

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
* 需要注意，很难完全减干净的
