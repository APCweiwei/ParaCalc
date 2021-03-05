# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 11:20:45 2018

@author: David

calculate paras of filters.

"""
import skrf as rf
import numpy as np
import os
import re
import csv

def select_s1p_filename(path,suffix_name): #提取所有.s1p的文件名但不带路径
    allfile=(os.listdir(path))
    namelist=[]
    for name in allfile:  
        if os.path.splitext(name)[1]==suffix_name:
            namelist.append(name)
    return(namelist)


def path_convert(path):
    path_temp=path+'\\'
    path_new=path_temp.replace('\\','/')
    return path_new

def filter_para_extract(Snn,IL_backoff1,IL_backoff2,filename,f_low,f_up):
    
    try:
        Row=eval(re.findall('-?\d+',filename[-14:-4])[-2])
        Column=eval(re.findall('-?\d+',filename[-14:-4])[-1])
    except IndexError:
        Row=0
        Column=0
    
    print(filename)
    
    #截取一段频率范围，在这个范围内计算fmax
    low_edge=np.where((Snn.f<f_low*1000000))[0][-1]
    up_edge=np.where((Snn.f>f_up*1000000))[0][0]
    
    Mag_s  =Snn.s21.s_db #提取s21的幅度
    Mag_s11=Snn.s11.s_db

    Mag_s_list=([])
    for jj in range(len(Mag_s)):
        Mag_s_list.append(Mag_s[jj][0][0])
        
    Mag_s11_list=([])
    for jj in range(len(Mag_s11)):
        Mag_s11_list.append(Mag_s11[jj][0][0])
        
    #list转成元祖进行最大值索引查询。
    Mag_s_array=np.array(Mag_s_list)
    
    #计算IL_max
    db_s21_max = np.max(Mag_s_array[low_edge:up_edge])
    db_s21_n1db=db_s21_max-IL_backoff1 #for example -20db
    
    #计算f_max
    f_max_index= np.where(Mag_s_array==db_s21_max)[0][0]
    fmax=Snn.f[f_max_index]

    try:
        #右侧30db滚降点计算。
        db_s21_30db=db_s21_max-30#回退30db
        db30_index_R = np.where((Mag_s_array<=db_s21_30db)&(Snn.f>fmax)&(Snn.f<f_up*1000000))[0][0]
        f_db30_R=round(Snn.f[db30_index_R]/1000000,2)
        
    except IndexError: 
        f_db30_R=0
    
    try:
        #左右两侧滚降20db计算
        f_left_ndb_index = np.where((Mag_s_array<=db_s21_n1db)&(Snn.f<fmax))[0][-1]
        f_right_ndb_index = np.where((Mag_s_array>=db_s21_n1db)&(Snn.f<f_up*1000000))[0][-1]
        
    except IndexError: 
        f_left_ndb_index=f_max_index-1
        f_right_ndb_index=f_max_index+1
        
    f_L_20db=round((Snn.f[f_left_ndb_index])/1000000,3)
    f_R_20db=(Snn.f[f_right_ndb_index])/1000000
       

    #fc计算
    index_fc=(f_left_ndb_index+f_right_ndb_index)/2
    if ((index_fc%2)!=0): #预防不能整除
        index_fc=(f_left_ndb_index+f_right_ndb_index+1)/2

    IL_fc=round(Mag_s_list[int(index_fc)] ,2) 
    f_c=(f_L_20db+f_R_20db)/2
    
    #others计算
    fb_fc_indexlist=np.where((Snn.f<f_c*1000000))[0][-1]
    fc_ft_indexlist=np.where((Snn.f>f_c*1000000))[0][0]
    
    IL_Max_left=round(np.max(Mag_s_list[0:fb_fc_indexlist]),2)
    IL_Max_right=round(np.max(Mag_s_list[fc_ft_indexlist:-1]),2)
    
    S11_fc=round(Mag_s11_list[int(index_fc)] ,2)

    #计算f1,f2处的抑制。
    f1_rejection=1500000000
    f1_rejection_index=np.where((Snn.f<f1_rejection))[0][-1]
    IL_f1_rejection=round(Mag_s_list[f1_rejection_index] ,2)
    
    f2_rejection=4500000000
    f2_rejection_index=np.where((Snn.f<f2_rejection))[0][-1]
    IL_f2_rejection=round(Mag_s_list[f2_rejection_index] ,2)
    
    BW_1db=0
    try:
        #计算1db BW。
        IL_backoff=IL_fc-IL_backoff2 
        BW_f_L_index = np.where((Mag_s_array<=IL_backoff)&(Snn.f<(f_c*1000000)))[0][-1]
        BW_f_R_index = np.where((Mag_s_array<=IL_backoff)&(Snn.f>(f_c*1000000)))[0][0]
    
        BW_f_L=(Snn.f[BW_f_L_index])/1000000
        BW_f_R=(Snn.f[BW_f_R_index])/1000000
        BW_1db=round(BW_f_R-BW_f_L,2)
    except IndexError: 
        BW_1db=0
    
    BW_3db=0
    try:
        #计算3db BW。
        IL_backoff=IL_fc-3 
        BW_f_L_index = np.where((Mag_s_array<=IL_backoff)&(Snn.f<(f_c*1000000)))[0][-1]
        BW_f_R_index = np.where((Mag_s_array<=IL_backoff)&(Snn.f>(f_c*1000000)))[0][0]
    
        BW_f_L=(Snn.f[BW_f_L_index])/1000000
        BW_f_R=(Snn.f[BW_f_R_index])/1000000
        BW_3db=round(BW_f_R-BW_f_L,2)  
        
    except IndexError: 
        BW_3db=0
        BW_f_L=f_c
        BW_f_R=f_c
    
    #右侧notch点计算
    #200MHz平分20分作为计算最小值网格。
    devide=20
    f_grid_step=(200)//devide
    
    min_list=[]
    for ii in range(devide):
        grid_low_edge=np.where((Snn.f<(f_R_20db+f_grid_step*ii)*1000000))[0][-1]
        grid_up_edge=np.where((Snn.f>(f_R_20db+f_grid_step*(ii+1))*1000000))[0][0]
        min_s21_grid=np.min(Mag_s_array[grid_low_edge:grid_up_edge])
        min_list.append(min_s21_grid)

    notch_peak_values=[]

    temp_grid_min=min_list[0]
    flag_notch1=1
    for ii in range(len(min_list)):

        if flag_notch1==1:
            notch_peak_freqs=[0]
            if min_list[ii]>temp_grid_min:
                notch_peak_freqs=[]
                notch_peak_values.append(min_list[ii-1])
                f_grid_index= np.where((Mag_s_array==min_list[ii-1]) )[0][0]
                f_grid=Snn.f[f_grid_index]/1000000
                notch_peak_freqs.append(f_grid)#这里用的是=号不是append
                flag_notch1=0

    
    print(notch_peak_freqs)
    
    #左侧notch点计算
    #notch点是3db带宽处往外扩50MHz，这段频率内的S21最小值。
    bound_start_L = np.where((Snn.f>((BW_f_L-50)*1000000)))[0][0]
    bound_stop_L = np.where((Snn.f>((BW_f_L)*1000000)))[0][0] 
    db_notch_left = np.min(Mag_s_array[bound_start_L:bound_stop_L])        
    db_notch_left_index=np.where(Mag_s_array==db_notch_left)[0][-1]
    
    
    f_notch_left=round(Snn.f[db_notch_left_index]/1000000,2)
    f_notch_right=round(notch_peak_freqs[0],2)


    delta_left=f_c-f_notch_left
    delta_right=f_notch_right-f_c
    db_s21_max=round(db_s21_max,2)

    paras=[filename,Row,Column,f_L_20db,f_R_20db,f_db30_R,f_c,S11_fc,IL_fc,IL_Max_left,IL_Max_right,\
                                  db_s21_max,IL_f1_rejection,IL_f2_rejection,BW_1db,\
                            BW_3db,f_notch_left,f_notch_right,delta_left,delta_right]    

    return paras


def save_filters_paras(path,name,data):   
    
    with open(path+name+'.csv', 'w', newline='') as f:

        csv_writer=csv.writer(f)
        csv_writer.writerow(['filename','Row','Column','f_L_20db','f_R_20db','f_R_30db','f_c','S11_fc',\
                             'IL_fc','IL_Max_left','IL_Max_right','IL_max','1500MHz_Rejection',\
                             '4500MHz_Rejection','BW_1db','BW_3db','f_notch_left','f_notch_right',\
                             'fc-notchL','notchR-fc'])          
        csv_writer.writerows(data)#写入多行时不要用for循环writerow，不然每次csv都被占用。  


if __name__ =='__main__':
    #r是告诉python转义字符无效。
    path = path_convert(r'C:\20200910_CBAW2A\2_1st_Test_before_trimming\Big_designs\WF01 apart result')
    namelist=select_s1p_filename(path,'.s2p')
    
    csv_name1='HLWF12-mid_filter_paras'
    #csv_name2='HH2_WF16_filter_paras_after_trimming_1st'
    
    paras=([])
    for ii in range(len(namelist)):
        Snn=rf.Network(path+namelist[ii])
        para=filter_para_extract(Snn,20,1,namelist[ii],4200,5200)  #10db是确定fc，1db是确定BW
        paras.append(para)
    print(paras)
    save_filters_paras(path,csv_name1,paras)        
      





    
    
    
    
    
    
    
    
    
    
    
    
    
