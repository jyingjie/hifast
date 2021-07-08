# parameters meaning

先用cli_baseline多项式去基线，然后

## cli_markRFI标记rfi

G15,eg:
```
python cli_markRFI.py $subname --outdir ./data \
        --lf --sf --lf_beams ['05','06','13'] \
        --sf_frange 1380 1382 --sf_times 10 --sf_thr 10 --sf_rfi_last 20 --sf_T_thr .3 \
        --lf_frange 1400 1450 --lf_times 1.5 --lf_thr 0 --lf_rfi_last 50 --lf_ext 10 \
        \
        --s_method_freq gaussian --s_sigma_freq 3 --s_method_t boxcar --s_sigma_t 7 \
        --rfi_thr 3 --rms_frange 1400 1403 --mw_frange 1419 1425 --rfi_groups 'two groups' \
        --rfi_width_lim 20 --ext_sec 20 --freq_thr .3 --freq_step 8.1 \
        --mask_RFI_method 'fixed freq' --ext_edge 3 --mask_thr 2 --mask_all_theory --freq_from_theory .5 \
        --plot -f  || exit 1 
```
### 时域

* --time_rfi

标记时域上突然出现的RFI，根据占据的频率长度，分为short-freq/long-freq。对该频率区间内的所有谱线做频率方向的平均后，画出一维的图。

* --sf
* --sf_frange: 一个典型的1～2MHz宽时域RFI，经常出现在1380～1382MHz。
* --sf_times: 大于各谱线之中，中值10倍，初步认为异常
* --sf_thr: 边缘处与附近值的差是中值的倍数，找出时间方向的边缘处陡峭的
* --sf_rfi_last: 持续出现的谱线数
* --sf_T_thr: 选定区域后，大于温度阈值的被标记

* --lf
* --lf_beams: 已知出现大卫星干扰的波束，其他不搜寻。如果不设置就对所有波束搜寻。
* --lf_frange: 大范围频率
* --lf_times: 大于各谱线之中，中值1.5倍，初步认为异常.
* --lf_thr: 不需要边缘陡峭，所以设为0
* --lf_rfi_last: 持续出现的谱线数
* --lf_ext: 选定区域后，扩展边缘

### 频域

频率上周期变化的RFI

* --s_method_freq,--s_method_t:平滑方法，'gaussian', 'boxcar', 'median', 'None'
 
  --s_sigma_freq,--s_sigma_t: sigma同前面定义

* --rfi_thr: 大于rms rfi_thr倍的会被记为rfi
* --rms_frange: 用于计算rms的干净频率区间
* --mw_frange: 一个大致可能存在信号的区间，先保护起来，后面不会用它预估rfi位置
* --rfi_groups: rfi的形态呈现为8.1MHz间隔的几组，分类进行更精确的拟合。'two groups','three groups','all'
* --rfi_width_lim: 粗找rfi区域时，其宽度应该大于的通道数。太细的可能是驻波，可能是点源，排除
* --ext_sec: 粗找到rfi后，扩充边界的通道数
* --freq_thr: 对粗测的rfi中心分类时，可以偏离预计中心一段频率。
* --freq_step: 对粗测的rfi中心分类时，预计中心的间隔大约为8.1MHz。

分类后进行直线拟合，得到rfi的精确分布。如果大于rfi_thr的区域落在rfi精确理论频率上，则被标记为真rfi，可以防止高流量的源被标记。

### 标记

* --ext_edge: 结果扩展边缘的通道数
* --mask_thr: 只mask超过rms的几倍。很小的rfi会保留。
* --mask_RFI_method: 'fixed freq' 对于小rfi，固定mask宽度
* --mask_all_theory: 是否理论的全mask，会非常干净。
  --freq_from_theory: 很小的去掉多少频率(固定宽度)
  
* --mask_RFI_method: '2 sides' 从中心向两边按step循环，确定mask边界，比较慢
* --small_rfi_times: 小于RMS这个倍数的不标记
  --chan_step: step通道数
  --freq_from_theory: mask rfi最大的宽度
  ```
  python ../../cli_markRFI.py $subname --outdir ./data \
        --sf --lf --lf_beams ['05','06','13'] \
        --sf_frange 1380 1382 --sf_times 10 --sf_thr 10 --sf_rfi_last 20 --sf_T_thr .3 \
        --lf_frange 1400 1450 --lf_times 1.5 --lf_thr 0 --lf_rfi_last 50 --lf_ext 10 \
        \
        --s_method_freq gaussian --s_sigma_freq 3 --s_method_t boxcar --s_sigma_t 7 \
        --rfi_thr 3 --rms_frange 1400 1403 --mw_frange 1419 1425 --rfi_groups 'two groups' \
        --rfi_width_lim 20 --ext_sec 20 --freq_thr .3 --freq_step 8.1 \
        --mask_RFI_method '2 sides' --ext_edge 10 --mask_thr 2 --freq_from_theory 3\
        --small_rfi_times 2 --chan_step 5 \
        --plot -f  || exit 1 
  ```
  
* --save_rfi_list: 是否保存各条谱线rfi理论值，```dict_out['rfi_list']```





## cli_baseline_fft去驻波

G15,eg:
```
rfi_freq_step = 0.0617283950617284
python cli_baseline_fft.py $subname --outdir ./sub_once \
        -rfi './data/G15_drift_5_arcdrift-M15_W-20200606-specs_T-bld-pdrfi.hdf5' \
        --fft_method rfft --sw_freq 0.9254  --amp_thr 35  --sw_n 5 \
        --rfi_8mhz --rfi_freq_step $rfi_freq_step  \
        --rfi_method subtract --mw_frange 1420.2 1420.55 --sg_window 1.0 --sg_polyorder 7 --mw_lower 1.0e5 \
        --plot -f --fill_rfi nan --keep_polar|| exit 1 
```
### 对RFI

* -rfi, --rfi_fname: 指定标记rfi的文件。无则自动到上一级名为/cor_vel_once/的文件夹里找*pdrfi*
* --rfi_method: 防止大rfi干扰fft，用原始值减去滤波值。目前只提供减法
* --mw_frange: 保护银河系区域
* --sg_window --sg_polyorder: savgol_filter的窗口大小(MHz)和阶数
* --mw_lower: mw拟合不上，这里就把它的流量压低，比如原来的1/1e5

### fft去驻波

* --fft_method: 目前只提供rfft
* --sw_freq: 1mhz左右的驻波，在fourier空间为0.925$\mu$s左右
* --amp_thr: fourier空间的振幅大于阈值且处在1/16.2$\mu$s间隔RFI的模，被选中为驻波一部分
* --sw_n: 距离sw_freq中心左右各sw_n个通道数的模，作为0.925$\mu$s驻波的一部分被选中

* --rfi_8mhz: 是否fft去除8mhz rfi
* --rfi_freq_step: 残余的8.1MHzRFI，在fourier空间间隔为1/16.2$\mu$s左右

* --fill_rfi: 对输入的rfi区域，可以填nan，减去滤波后的噪声或者保持原样。'nan','noise','rfi'
* --keep_polar: 保留xx，yy，返回三维的结果

输出的hdf5包含去除驻波的数组(T或者flux)，和单独的驻波```dict_out['ripple']```


