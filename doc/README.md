#  安装
## 配置环境：
* 方法一：用conda配置环境

  如果没有conda，请安装[miniconda3](https://docs.conda.io/en/latest/miniconda.html)或[Miniforge3](https://github.com/conda-forge/miniforge)
  
  下载或者在解压后的hifast.xxx.tar.gz文件目录中找到环境配置文件[hifast_env.yml](../hifast_env.yml)，(ARM架构CPU使用[hifast_env.ARM64.yml](../hifast_env.ARM64.yml))
  
  * 用conda新建一个环境：
     ```
     $ conda env create -n myenv --file hifast_env.yml 
     ```
     可以修改myenv为其他字符。之后用
     ```
     $ conda activate myenv
     ```
     或
     ```
     $ source activate myenv
     ```
     来切换到对应的环境。 
  * 或 更新已存在的conda环境使其满足hifast的依赖库：
     ```
     $ conda env update --file hifast_env.yml -n ENV_NAME
     ```
     替换```ENV_NAME```为已存在的conda env名字，如果是主环境，则```ENV_NAME```为```base```

* 方法二：使用离线包配置环境（测试中，可能会有bug）
 下载 https://pan.cstcloud.cn/s/QmkBVYgyRO8 中的hifast_env.XXX.tar.gz
  
  首次配置：
  ```
  $ mkdir ~/hifast_env
  $ tar -zxvf hifast_env.XXX.tar.gz -C ~/hifast_env
  $ source ~/hifast_env/bin/activate
  (hifast_env) $ conda-unpack
  (hifast_env) $ source ~/hifast_env/bin/deactivate
  ```
  之后激活环境用：
  ```
   $ source ~/hifast_env/bin/activate
  ```
  移除环境用：
  ```
   (hifast_env) $ source ~/hifast_env/bin/deactivate
  ```
## 安装hifast

* 方法一： 用pip联网直接安装：

  见 [wikis/Install](../../wikis/Install)

* 方法二： 下载后安装

  下载安装包解压后**先 cd 切换到代码(setup.py)所在目录下**

  确认已激活之前配置好的环境

  * 安装
    ```
    $ python -m pip install . --upgrade 
    ```
  * 卸载
    ```
    $ pip uninstall hifast
    ```
  * 提示
    <br/>由于包文件里包括c代码，因此必须安装后才能正确 import。另外由于Python在 import 时是从当前路径开始搜索包，
    因此不要在安装包（setup.py）所在目录下执行 hifast相关命令。

# 示例
  
## 准备

  + 下载 Tcal文件夹(  https://pan.cstcloud.cn/s/AfnCB96cT2s 提取码：cqwy  ) 放到你的 home 目录 . (*2020-12-07后的代码需更新此文件*)

  **Note:**

   * 大部分hifast.xxx可以在终端中用 ```python -m hifast.xxx``` 来执行，后接输入的文件名和```-```加一个字母或```--```加多个字母的参数。
   *  ```python -m hifast.xxx -h``` 显示帮助。
   * 对于生成的hdf5文件，在终端执行
     * ```h5dump -g /Header XXX.hdf5``` 显示生成该文件时用的参数。
     * ```h5dump -n XXX.hdf5``` 显示文件中有什么内容
     
     也可以用[CARTA 2.0](https://carta.readthedocs.io/en/latest/index.html)打开。

## 1. 分离光谱并利用噪音管计算得出温度
 ```hifast.sep``` : 输入一个beam的谱线数据，得到 噪音管Cal on 和 off谱线的温度。

* 示例

  ```
  fname="/data/inspur_disk06/fast_data/3047/GAMA_G15/20191215/XXX_0001.fits"
  python -m hifast.sep $fname -d 0 -m 1 -n 120 --step 5 --frange 1329 1429 --smooth poly --s_deg 1 --outdir ./data
  # or
  python -m hifast.sep $fname -d 0 -m 1 -n 120 --step 5 --frange 1020 1445 --smooth gaussian --s_sigma 5 --outdir ./data
  # or 
  #对于F或N带，可以加 --dfactor W 降低采样到 W 带，减小文件大小，提高信噪比。
  python -m hifast.sep $fname -d 0 -m 1 -n 120 --step 5 --frange 1020 1445 --smooth gaussian --s_sigma 5 --dfactor W --outdir ./data
  ```
  程序先找到所有Cal on的谱线，然后减去相邻的Cal off的谱线得到Cal的Power。然后每条谱线使用时间上最近的Cal Power。
* 后面紧跟的是数据文件的绝对路径。
  <br/>FAST的数据里每个beam是存成了很多块文件，以0001.fits——9999.fits结尾。这个仅需要给任何一个块文件的文件名即可，程序会从第一个文件开始依次读取处理。
* 其他参数：
   * ```-d -m -n ```：分别为delay时间，Cal on时间和Cal off时间 除以谱线的采样时间。三个数都为整数。delay默认为0。
   * ```--noise_mode```: 噪音管强度。 high 或者 low，默认为high。
   * ```--noise_date```: 选用哪天的噪音温度文件。例如设为 20190115 或 20200531. 如果设为auto则选取与谱线观测时间最近的噪音管文件来定标。(20200531的噪音管文件Beam19 XX 偏振 在约1060MHz处有个大的gap。)
   * ```--step``` ：每次读入内存的块文件数量。
   * ```--frange```：程序提取和处理的频率范围：后接两个数，分别是频率的下限和上限。(需配合```--smooth```来设置)
   <br/> M31 HI 频率和速度（相对于LSR）的对应 见文件 freq_vs_vlsr_M31.txt。另外可以修改freq_vs_vlsr.py来估算你需要的 HI频率和速度的对应。不同的ra,dec和mjd，vlsr和freq的对应稍有差异。
   * ```--smooth```: 平滑方法。gaussian, poly 或者 mean（默认）。
     <br/>噪音管定标，在每一个Tcal on off的周期里, 温度在不同频率的值为
     <br/>   <img src="https://render.githubusercontent.com/render/math?math=T_{on}(\nu)= (\frac{P_{cal\_on}(\nu)}{(P_{cal\_on}(\nu)-P_{cal\_off}(\nu))_{smooth}}-1)*T_{cal}(\nu)_{smooth}">
     <br/>   <img src="https://render.githubusercontent.com/render/math?math=T_{off}(\nu)= \frac{P_{cal\_off}(\nu)}{(P_{cal\_on}(\nu)-P_{cal\_off}(\nu))_{smooth}}*T_{cal}(\nu)_{smooth}">
     * mean: 适用于freq区间（```--frange```）比较小（至少小于10?）。对freq内的所有流量值平均。
     * poly: 适用于freq区间比较小。多项式拟合，需要同时加```--s_deg```参数，```s_deg```为 1 时为线性拟合，为0时等同于mean，        默认值为1。
     * gaussian：适用于freq区间比较大，需要同时加```--s_sigma```参数，```s_sigma```需要小于freq区间的三分之一。平滑结果在        频率两端会不太准，因此推荐freq区间比需要的范围稍大一些。
   * ```--med_filter_size```：用于median filter 谱线的size，去除极窄的rfi。
   * ```--dfactor```：降低频率的采样率。```--dfactor 16```：每16个采样点平均以降低采样，```--dfactor W```：降低采样到W带。
   * ```--outdir``` ：输出文件存放的目录。
* 输出文件名以 ```specs_T.hdf5```结尾。可以用h5py来读取。例如：
  ```
   import h5py
   f= h5py.File('data/XXX_arcdrift-M02_F-specs_T.hdf5','r')
   S=f['S']
   print(S.keys())
   S['mjd'][:]
  ```
* 同时会输出pdf图片用来检查Cal on off的分离是否正确。

## 2. 转换馈源（KY）位置到RA DEC

```hifast.radec```: 转换馈源（KY）位置到RA DEC

* 示例

```
python -m hifast.radec data/XXX_arcdrift-M01_F-specs_T.hdf5
# or
python -m hifast.radec XXX_arcdrift_11_2020_XX_XX_22_54_09_000.xlsx
```

* 输入 ```hifast.sep```生成的hdf5文件，***仅需beam M01的即可***，其它beam的RA DEC会存在这同一个文件里。
  <br/>或者输入一个.xlsx结尾的馈源舱文件。
* 程序会先依次检测是否存在 *```your_HOME_dir/KY/```*, */data/inspur_disk06/fast_data/KY/*, */data31/KY/* (FAST服务器馈源文件所在目录)文件夹。 然后在最先检测到存在的文件夹里自动寻找对应的馈源舱文件。然后计算馈源舱文件文件里记录的RA DEC，由于谱线记录的时间采用点和馈源舱文件的一般不一致，因此会插值最后得到谱线的RA DEC。
* 其他参数：
  * ```--ky_files```：如果程序不能成功找到对应馈源舱文件，加此参数来手动指定馈源舱文件。参数后加馈源舱文件的路径。
  * ```--tol```: 正常情况下，馈源舱文件内记录的时间覆盖谱线记录的时间范围。加```--tol 5```可在馈源舱文件少记录```5```秒的情况下进行外插计算RA DEC。不过这部分谱线最终需扔掉。
  * ```--outdir``` ：指定输出文件存放的目录。如果输入hdf5文件默认与输入文件一致，如果输入.xlsx文件则默认为程序运行路径。
  * ```--ky_fixed```: 早期的一些Drift观测，馈源舱文件只记录开始的几分钟内的馈源舱位置，加此参数只利用这开始的几分钟来计算整个谱线的RA-DEC。由于Drift过程中馈源舱并不能完全稳定不动，这样计算出的RA DEC可能会不准确。(一般不建议使用)
  * ```--plot```：加此参数来画RA DEC分布图，保存为pdf图片。
* 存放radec的输出文件名是在输入文件名上加radec。可以h5py来读取，例如：
  ```
  import h5py
  f= h5py.File('data/XXX_arcdrift-M01_F-specs_T-radec.hdf5','r')
  S=f['S']
  print(S.keys())
  S['mjd'][:]
  ```

## 流量定标

```hifast.flux```: 流量定标Teff(K) &#8594; Flux(Jy/beam)

* 示例

```
python -m hifast.flux data/XXX_arcdrift-M01_F-specs_T.hdf5

```

* 输入未流量定标的hdf5文件。程序会整合RADEC，在输入的文件的所在目录下去读取对应的radec的文件。即```hifast.radec```输出的文件名中有'M01'的radec文件。
* 其他参数：
  * ```--cali_fname```：后接定标源文件路径，如果不加此参数，程序将使用https://arxiv.org/abs/2002.01786 中给出的Gain与天顶角的函数关系来流量定标。

## 扣除基线

```hifast.bld``` : 拟合并扣除基线(**b**ase**l**ine**d**)。

* 示例
```
python -m hifast.bld data/XXX_arcdrift-M01_F-specs_T-flux.hdf5 --method arPLS --nproc 5
```
* 输入 ```hifast.sep```产生的```specs_T.hdf5```文件或者```hifast.flux```产生的```specs_T-flux.hdf5```文件
    <br/>
    <br/> 也可以输入它自己生成的文件，进行二次去基线。
* 其他参数
  * ```--nproc```: 后接一个数，使用多少个进程来并行。
  * ```--frange```：
    <br/> 只用这个频率范围内谱线。后接两个数，空格隔开，下限在前。
    <br/> 范围越大，拟合用时越长，并且不是线性增长。
  * ```--method```: 拟合方法： 'arPLS', 'srPLS', 'Chebyshev', 'poly'。
  * ```--lam , --deg , --offset```：
    <br/> 去基线时的参数，--lam调整平滑度，越大越接近低阶多项式(poly)拟合。
  * 输出文件名会包含```bld```的hdf5文件。可以用h5py来读取。

## 拟合驻波

```hifast.sw``` : 拟合并扣除驻波(**s**tanding**w**ave)。

* 示例
```
python -m hifast.sw data/XXX_arcdrift-M01_F-specs_T-flux-bld.hdf5 --method fft --nproc 5
```
* 输入 ```hifast.bld```产生的```*-bld.hdf5```文件
* 其他参数
  * ```--nproc```: 后接一个数，使用多少个进程来并行。
  * --method {sin_poly,fft}

    sin_poly: 最小二乘法拟合一个正弦函数

    fft: 在傅里叶变换后的相空间里操作再反傅里叶变换得到驻波

## rfi标记
```hifast.rfi``` ：rfi标记

* 示例：
  ```
  fname=XXX-bld.hdf5
  python -m hifast.rfi $fname --pr True --pr_s_sigma 5 --pr_times 5 --pr_times_s 1.1
  #or
  python -m hifast.rfi $fname --tr True --tr_s_sigma 5 --tr_times 5 --pr_times_s 1.6 --tr_n_continue 50 --ext_add 0
  
  ```
* 输入去完基线得到的文件。
* 参数
  * ```--pr```, ```--tr``` 分别为 时域rfi和偏振rfi标记，会生成 is_rfi 这项存在输出文件里。

    ```--pr True```
    * ```--pr_s_sigma```: 沿时间维度高斯平滑（以pr_s_sigma为sigma）谱线以提高信噪比
    * ```--pr_times``` : 至少大于等于5
    * ```--pr_times_s```：大于1
    
         <br/>
    ```--tr True```
    * ```--tr_s_sigma```: 沿时间维度高斯平滑（以tr_s_sigma为sigma）谱线以提高信噪比
    * ```--tr_times```: 至少大于等于5
    * ```--tr_times_s```: 大于1.5
    * ```--tr_n_continue```: 一个“信号”持续多少条就标记为rfi

* 输出文件名中包含rfi。

## 静止坐标系（frame）修正等等
```hifast.multi``` ：静止坐标系（frame）修正等

* 示例：
  ```
   python -m hifast.multi XXX.hdf5 --fc True --frame LSRK
  ```
* 参数
  * ```--fc```从望远镜所在的地平参考系修正到太阳或者LSR为中心的参考系。同时会合并两个偏振；如果存在 is_rfi，则会把rfi替换为nan。

    ```--fc True```
    * ```--frame```: 参考系选择，HELIOCEN 或者 LSRK
  * ```--replace_rfi```: 如果设为True并且输入文件中存在is_rfi，则会把rfi的值替换为nan。默认为True
  * ```--merge_polar```: 如果设为True，合并两个偏振。默认为True

* 输出文件名根据输入参数改变，可能包含 fc

## 5. 生成fits cubes 文件

  ```
  python -m hifast.cube **/data/*-fc*.hdf5 --outname ./test_cubes.fits --bwidth 60 -p SIN
  ```
* 这里用 ```hifast.cube``` 来生成fits cubes文件，程序先生成ra、dec格点，然后找到距离格点中央为```--r_cut```范围内的谱线然后按```--method```处理谱线，最后保存在fits文件里。
* ```python -m hifast.cube```
  * 后面跟参考性修正后生成的hdf5文件，支持多个文件路径（空格隔开），支持通配符。程序运行后会首先输出要处理的文件路径，请检查无重复无错误。
  * ```--outname```: 输出的fits cubes文件路径，需要指定。
  * ```--bwidth```: ra dec 分格点时的间隔大小，单位为 角秒，默认为60。如果ra和dec采用一样间隔，参数后接一个数字即可，如果不一样，参数后接两个数，空格隔开。ra的间隔在前。
  * ```--r_cut```: 考虑距离格点中心r_cut范围内谱线。单位为 角秒，默认为90.
  * ```--method```: r_cut范围内谱线处理方法:
    * ```mean```: 对谱线求平均.
    * ```median```: 对谱线求median值.
    * ```reweight```: Barnes el. al. 2001, MNRAS 322, 486 https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..486B/abstract .
    * ```gaussian```: 待添加。
  * ```--proj```: 投影方式: SIN, AIT, TAN. 
  * ```--ra_range```: ra的范围，后接两个数，空格隔开，下限在前，单位为度。默认值为输入文件里ra的最小值和最大值。
  * ```--dec_range```: 类似```--ra_range```。
  * ```--range3```: 限制第三轴（速度）的范围，类似```--ra_range```。

# hifast.sh 组合脚本

#### 用hifast.sh来对多个文件执行相同的```python -m hifast.xxx```操作

  ```
  hifast.sh /data/inspur_disk06/fast_data/3047/G15_drift_5/20200606/*_W*0001.fits -n 10 -c "python -m hifast.sep | -d 0 -m 1 -n 1  --step 5  --frange 1369 1394 --smooth poly --s_deg 1 --outdir ./data"
  ```

* 这里用通配符来指定19个beam的文件.
* ```-c``` 参数指定用到的```python -m hifast.xxx```操作。
  <br/>```-c``` 后引号内的内容用 ```|``` 分割hifast命令和输入参数，省略输入的文件名。
* ```-n``` 指定多少个文件同时执行。

* 另外可以用.txt文件来输入文件列表，用.par文件来输入要执行的命令。
  
  例如:
  
  * *files.txt* 内容如下：
  
    ```
    /data/inspur_disk06/fast_data/3047/G15_drift_5/20200606/G15_drift_5_arcdrift-M01_W_0001.fits
    /data/inspur_disk06/fast_data/3047/G15_drift_5/20200606/G15_drift_5_arcdrift-M02_W_0001.fits
    /data/inspur_disk06/fast_data/3047/G15_drift_5/20200606/G15_drift_5_arcdrift-M03_W_0001.fits
    /data/inspur_disk06/fast_data/3047/G15_drift_5/20200606/G15_drift_5_arcdrift-M04_W_0001.fits
    /data/inspur_disk06/fast_data/3047/G15_drift_5/20200606/G15_drift_5_arcdrift-M05_W_0001.fits
    ```
  
  
  * *commands.par* 内容如下：
  
  ```
  python -m hifast.cli_sep | -d 0 -m 1 -n 1  --step 5  --frange 1369 1394 --smooth poly --s_deg 1 --outdir ./data

  ```
    
  * 可以执行：
  
  ```hifast.sh -i files.txt -n 10 -c commands.par```
  
##### hifast.sh计算radec
  
  ```
  hifast.sh data/*M01*specs_T.hdf5 -c "python -m hifast.radec |  "
  ```

##### hifast.sh 组合 ```hifast.bld``` 和 ```hifast.multi```

  * *commands.par* 内容如下：
  
  ```
  # 续行符 "\" 之后不能有空格
  python -m hifast.bld | --method arPLS --lam 1e7 --nproc 5 --outdir ./
  python -m hifast.multi  |  --tr --tr_method smooth --tr_s_sigma 5 --tr_n_continue 100 --fc \
        --keep_rfi --keep_polar

  ```
  * 执行：
  
  ```hifast.sh data/*M*specs_T.hdf5  -n 10 -c commands.par```
  <br/> 先运行```python -m hifast.bld```这行去基线，然后```hifast.sh```会自动把生成的文件输入到```python -m hifast.multi``` 这行。
  <br/> 注意 ```python -m hifast.bld```这行执行时会占用 ```3 * 5 = 15``` 个线程。
