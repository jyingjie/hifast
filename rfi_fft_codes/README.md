# parameters meaning

Author: Xu Chen, NAOC, 2021.Jun ~ Aug

1. 先用cli_baseline多项式去基线
2. cli_mark_tRFI找时域RFI
3. cli_baseline_fft去驻波
4. cli_baseline再去一次基线
5. cli_multi找偏振很大的RFI、坐标系修正
6. （可以合并通道）
然后成图

有bug请联系stellarxu@qq.com

to do list:
0.92 MHz 周期脉冲RFI

参数具体设置参见notebook(需要配置可以可视化调参数jupyterlab环境)
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
python cli_mark_tRFI.py $subname \
        --mw_frange 1420.3 1423.3 --sf_frange_step 20 \
        --sf --sf_times 3 --sf_thr 0 --sf_rfi_last 10 --sf_T_thr_times 3 \
        --lf --lf_frange 1400 1450 --lf_times 1.5 --lf_thr 0 --lf_rfi_last 50 --lf_ext 10 \
        --plot --ylim -.5 .5 --vmin_max -.05 .05 -f  || exit 1 --outdir $outdir \
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


## cli_baseline_fft去驻波

### 对RFI
* --rfi_method: 防止大rfi/强源干扰fft，简化版只提供'near ripple'方法

* --rms_sigma,--rms_frange(测试): 如果区间有干扰，还是用该频率减去一个高斯平滑，得到准确的rms。此为高斯滤波的sigma。
* -rfi, --rfi_fname: 指定标记rfi的文件。无则自动到输入文件的上一级名为/rfi/的文件夹里找含-tr的文件，找不到则在当前路径找。

* --mw_frange(测试): 银河系(或者M31，M33都)存在的区域

 'near ripple'使用高流量处附近的驻波替代。例如
    ```
    fname = subnames[1]
    outdir = './data'
    sepname='/data/inspur_disk01/userdir/jyj/FAST/jingyj/M31_snapshot/data/M31_SnapShot_6_snapshot-M'
    python /data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/test_pipe/cli_baseline_fft.py $fname \
        --outdir $outdir -sep $sepname  \
        --mw_frange 1420.3 1423.3 --rms_sigma 6 --rms_frange 1390 1400 \
        --rfi_method 'near ripple' --times_lower_thr 4 \
        --rfi_width_lim 15 --ext_sec 20 --ext_freq 1.3 \
        --fft_method rfft --chan_wide 5 --chan_narr 3 \
        --amp_thr_mean_factor 1.05 --amp_thr_factor 1.4 --choose_method 'all' \
        --rip_base --rip_1mhz --rip_0_04mhz --fft_ylim -5 160 \
        --plot --ylim -1 .5 --vmin_max -.05 .05 -f --fill_rfi nan --keep_polar || exit 1 
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
* --nspec: 平均谱线条数，默认10条。过少容易损失信号，过多容易驻波变化了。
* --func: 自己写的'iter'或卷积的'smooth'。默认'iter'，但是多了会比较慢。

```
python cli_baseline_mean.py $fname --outdir $outdir \
        --sub_method mean --nspec 200 --frange 1330 1360 \
        --plot --ylim -1 .5 --vmin_max -.05 .05 -f --keep_polar --fill_rfi nan || exit 1 
```

输出文件名含mean_bld/med_bld

需要注意，由于不同谱线基线不同，还需要使用cli_baseline.py多项式去一次基线。


## cli_RFIshape利用同一谱线不同频率的8mhz 周期RFI波形

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
* 需要注意，很难完全减干净的
