# parameters meaning

先用cli_baseline多项式去基线，然后，以后再说

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
        --plot -f  || exit 1 
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
* --mask_all_theory: 是否理论的全mask，会非常干净(输出文件名含strict)
  --freq_from_theory: 很小的去掉多少频率(固定宽度)
  
* --mask_RFI_method: '2 sides' 从中心向两边按step循环，确定mask边界，比较慢
* --small_rfi_times: 小于RMS这个倍数的不标记
  --chan_step: step通道数
  --freq_from_theory: mask rfi最大的宽度
  例如
  ```
  python ../../cli_markRFI.py $subname --outdir ./data \
        --sf --lf --lf_beams ['05','06','13'] \
        --sf_frange 1380 1382 --sf_times 10 --sf_thr 10 --sf_rfi_last 20 --sf_T_thr .3 \
        --lf_frange 1400 1450 --lf_times 1.5 --lf_thr 0 --lf_rfi_last 50 --lf_ext 10 \
        --period_rfi \
        --s_method_freq gaussian --s_sigma_freq 3 --s_method_t boxcar --s_sigma_t 7 \
        --rfi_thr 3 --rms_frange 1400 1403 --mw_frange 1419 1425 --rfi_groups 'two groups' \
        --rfi_width_lim 20 --ext_sec 20 --freq_thr .3 --freq_step 8.1 \
        --mask_RFI_method '2 sides' --ext_edge 10 --mask_thr 2 --freq_from_theory 3\
        --small_rfi_times 2 --chan_step 5 \
        --plot -f  || exit 1 
  ```
* --time_coherent_per: 如果该频率，例如90%都被标记了，那么整条都被标记为RFI。不怎么用。
* --save_rfi_list: 是否保存各条谱线rfi理论值，```dict_out['rfi_list']```

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
* --rfi_method: 防止大rfi/强源干扰fft，提供三种方法：'lower','set zeros','subtract'
* --times_lower_thr: 大于RMS这么多倍的会被降低
* --times_lower: 把高流量压低，比如原来的1/1e4
* --rms_sigma: 如果区间有干扰，还是用该频率减去一个高斯平滑，得到准确的rms。此为高斯滤波的sigma。

1. 'lower'降低高流量。不需要rfi文件。例如
    ```
    python cli_baseline_fft.py $subname --outdir ./data \
        --rfi_method 'lower' --times_lower 1.0e4 --times_lower_thr 2 \
        --fft_method rfft --sw_freq 0.9254   --amp_thr 35  --sw_n 5 \
        --plot -f --keep_polar || exit 1 
    ```
    输出文件名含fft_bldl.
    
2. 'set zeros'高流量全部置零。不需要rfi文件。例如
    ```
    python ../../cli_baseline_fft.py $subname --outdir ./data \
        --rfi_method 'set zeros' --times_lower 1.0e4 --times_lower_thr 2 \
        --plot -f --keep_polar || exit 1 
    ```
    输出文件名含fft_bldz.
    
3. 'subtract'用原始值减去滤波值。例如
    ```
    python ../../cli_baseline_fft.py $subname --outdir ./data \
            -rfi $rfiname --mw_frange 1420.2 1420.55 \
            --rfi_method subtract --sg_window 1.0 --sg_polyorder 7 --times_lower 1.0e4 --times_lower_thr 2 \
            --fft_method rfft --sw_freq 0.9254   --amp_thr 35  --sw_n 5 \
            --plot -f --fill_rfi nan --keep_polar || exit 1 
    ```
    * -rfi, --rfi_fname: 指定标记rfi的文件。无则自动到输入文件的上一级名为/rfi/的文件夹里找含-tr的文件，找不到则在当前路径找。
    * --mw_frange: 银河系存在的区域
    * --sg_window --sg_polyorder: savgol_filter的窗口大小(MHz)和阶数
    输出文件名含fft_blds.

### fft去驻波

* --fft_method: 目前只提供rfft
* --sw_freq: 1mhz左右的驻波，在fourier空间为0.925$\mu$s左右
* --amp_thr: fourier空间的振幅大于阈值且处在1/16.2$\mu$s间隔RFI的模，被选中为驻波一部分
* --sw_n: 距离sw_freq中心左右各sw_n个通道数的模，作为0.925$\mu$s驻波的一部分被选中

* --rfi_8mhz: 是否fft去除8mhz rfi
* --rfi_freq_step: 残余的8.1MHzRFI，在fourier空间间隔为1/16.2$\mu$s左右

* --fill_rfi: 对输入的rfi区域，可以填nan或者保持原样。对应'nan','rfi'。需要注意只有输入对应RFI文件才能填nan。
* --keep_polar: 保留xx，yy，返回三维的结果
* -sep, --sep_fname: 只分离光谱温度定标的文件。无则自动到输入文件的上一级名为/sep/的文件夹里找含-tr的文件，找不到则在当前路径找。

输出的hdf5包含去除驻波的数组(T或者flux)(单独的驻波```dict_out['ripple']```注释掉了)


## cli_baseline_mean去基线

* --sub_method: mean，平均前后若干条谱线作为基线
* --nspec: 平均谱线条数，默认10条。过少容易损失信号，过多容易驻波变化了。

输出文件名含mean_bld.

需要注意，由于不同谱线基线不同，还需要使用cli_baseline.py多项式去一次基线。
