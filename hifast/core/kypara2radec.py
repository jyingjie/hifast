#!/usr/bin/env python
# coding: utf-8

import numpy as np
import math
import erfa

#地理
hm = 1110.03;
elong = 1.8650006658875442;
phi = 0.44772847428643703;
#观测波长
wl = 300000.0;
gap = 0.27 #多波束中心间距 米

#属性-可更新的参数：
dUT1 =  0.1
phpa = 925. #气压
temperature =  15. #气温
humidity =  0.8 #相对湿度

feedLocalCoord = [  [ 0.0  , 0.0  , 0.0   ],
                    [ -0.27, 0.0, 0.0],
                    [ -0.135, 0.233826859, 0.0],
                    [ 0.135, 0.233826859, 0.0],
                    [ 0.27, 0.0, 0.0],
                    [ 0.135, -0.233826859, 0.0],
                    [ -0.135, -0.233826859, 0.0],
                    [ -0.54, 0.0, 0.0],
                    [ -0.405, 0.233826859, 0.0],
                    [ -0.27, 0.467653718, 0.0],
                    [ 0.0, 0.467653718, 0.0],
                    [ 0.27, 0.467653718, 0.0],
                    [ 0.405, 0.233826859, 0.0],
                    [ 0.54, 0.0, 0.0],
                    [ 0.405, -0.233826859, 0.0],
                    [ 0.27, -0.467653718, 0.0],
                    [ 0.0, -0.467653718, 0.0],
                    [ -0.27, -0.467653718, 0.0],
                    [ -0.405, -0.233826859, 0.0] ]


def kypara2radec(mjd, multibeamAngle, nB, globalCenterX,  globalCenterY, globalCenterZ, globalYaw, globalPitch, globalRoll):
    """
    Return: ra, dec in radians.
    """
    #1.多波束转角带来的旋转
    rotationMatrixMultiBeam= CalMultiBeamRotationMatrix(multibeamAngle)

    #2.Stewart下平台带来的旋转、
    rotationMatrixPlatform= CalPlatformRotationMatrix(globalYaw, globalPitch, globalRoll)

    #3.计算旋转后的位置
    useFeed = feedLocalCoord[nB-1]
    posRelative = VectorTransform(rotationMatrixPlatform, VectorTransform(rotationMatrixMultiBeam, useFeed))

    #4.计算全局坐标
    posAbsolute = posRelative + np.array([globalCenterX,  globalCenterY, globalCenterZ])

    #5.计算地平坐标
    aob = math.atan2(-posAbsolute[0], -posAbsolute[1])
    zob = math.atan2(math.sqrt(posAbsolute[0]**2 + posAbsolute[1] **2), -posAbsolute[2])

    #6.准备天文参数
    timeJD= mjd+ 2400000.5
    utc1 = math.floor(timeJD) + 0.5
    utc2 = timeJD - utc1
     
    xp = 0.
    yp = 0.
    return erfa.atoc13(b"A", aob, zob, utc1, utc2, dUT1, elong, phi, hm, xp, yp, phpa, temperature, humidity, wl)


#other functions
#从多波束转角计算旋转矩阵
def CalMultiBeamRotationMatrix(multibeamAngle):          
    """
    multibeamAngle: scalar  float64  Radian 
    """    
    ca= math.cos(multibeamAngle) 
    sa= math.sin(multibeamAngle) 
    rotationMatrixMultiBeam= np.zeros((3,3))
    
    rotationMatrixMultiBeam[0,0] = ca 
    rotationMatrixMultiBeam[0,1] = -sa 
    rotationMatrixMultiBeam[0,2] = 0.0 
    rotationMatrixMultiBeam[1,0] = sa 
    rotationMatrixMultiBeam[1,1] = ca 
    rotationMatrixMultiBeam[1,2] = 0.0 
    rotationMatrixMultiBeam[2,0] = 0.0 
    rotationMatrixMultiBeam[2,1] = 0.0 
    rotationMatrixMultiBeam[2,2] = 1.0 
    return rotationMatrixMultiBeam

#从偏航、俯仰、翻滚计算旋转矩阵。
def CalPlatformRotationMatrix(globalYaw, globalPitch, globalRoll):
    """
    globalYaw, globalPitch, globalRoll: scalar  float64
    """
    
    cy = math.cos(globalYaw) 
    sy = math.sin(globalYaw) 
    cp = math.cos(globalPitch) 
    sp = math.sin(globalPitch) 
    cr = math.cos(globalRoll) 
    sr = math.sin(globalRoll) 
    
    rotationMatrixPlatform = np.zeros((3,3))
    rotationMatrixPlatform[0,0] = cy * cp 
    rotationMatrixPlatform[0,1] = cy * sp * sr - sy * cr 
    rotationMatrixPlatform[0,2] = sy * sr + cy * sp * cr 
    rotationMatrixPlatform[1,0] = sy * cp 
    rotationMatrixPlatform[1,1] = cy * cr + sy * sp * sr 
    rotationMatrixPlatform[1,2] = sy * sp * cr - cy * sr 
    rotationMatrixPlatform[2,0] = -sp 
    rotationMatrixPlatform[2,1] = cp * sr 
    rotationMatrixPlatform[2,2] = cp * cr 
        
    return rotationMatrixPlatform

def VectorTransform(matrixA, vectorB):
    """
    """ 
    
    m11 = matrixA[0, 0] * vectorB[0] + matrixA[0, 1] * vectorB[1] + matrixA[0, 2] * vectorB[2]
    m21 = matrixA[1, 0] * vectorB[0] + matrixA[1, 1] * vectorB[1] + matrixA[1, 2] * vectorB[2]
    m31 = matrixA[2, 0] * vectorB[0] + matrixA[2, 1] * vectorB[1] + matrixA[2, 2] * vectorB[2]

   
    return  np.array([m11, m21, m31])
