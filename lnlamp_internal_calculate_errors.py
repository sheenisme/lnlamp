#!/usr/bin/env python3

import sys
from pathlib import Path
from decimal import *
import pandas as pd
from functools import *



def ReadValues(filename):
  with open(filename, 'r') as f:
    l = f.readline()
    while l != '':
      for v in l.strip().split():
        if v != '':
          yield v
      l = f.readline()



def ComputeDifference(float_data, double_data):
  n = 0
  acc_err = Decimal(0)
  acc_val = Decimal(0)
  max_abs_err=Decimal('-Infinity')
  max_prc_err=Decimal('-Infinity')
  double_no_fl = 0
  float_no_fl = 0

  thres_ofl_cp = Decimal('0.01')

  for sv_float, sv_double in zip(float_data, double_data):
    v_double, v_float = Decimal(sv_double), Decimal(sv_float)

    if v_double.is_nan() or v_double.is_infinite() :
      double_no_fl += 1
    elif v_float.is_nan() or v_float.is_infinite():
      float_no_fl += 1
      double_no_fl += 1
    elif ((v_float + v_double).copy_abs() - (v_float.copy_abs() + v_double.copy_abs())) > thres_ofl_cp:
      double_no_fl += 1
    else:
      n += 1
      abs_err_temp = (v_double-v_float).copy_abs()
      acc_err += abs_err_temp
      # 更新最大绝对误差
      if(abs_err_temp > max_abs_err) :
        max_abs_err = abs_err_temp
      
      # 如果真值是零的话，则继续
      if v_double.copy_abs().is_zero():
        continue;
      else:
        prc_err_temp=abs_err_temp / v_double.copy_abs() * 100
        # 更新最大相对误差
        if(prc_err_temp > max_prc_err) :
          max_prc_err = prc_err_temp
        acc_val += v_double.copy_abs()
      
  e_prc = (acc_err / acc_val * 100) if acc_val > 0 and n > 0 else -1
  e_abs = (acc_err / n) if n > 0 else -1
      
  return {'double_no_fl': double_no_fl, \
          'float_no_fl': float_no_fl, \
          'mean_pec': e_prc,\
          'mean_abs': e_abs,\
          'max_pec': max_prc_err,\
          'max_abs': max_abs_err}               



def PrettyPrint(table):
  df = pd.DataFrame.from_dict(table,orient='index').T
  df = df.transpose()
  return df.to_string()



if __name__ == "__main__":
  if len(sys.argv) != 3:
    print('参数个数为:', len(sys.argv), '个参数。')
    print('参数列表:', str(sys.argv))
    print('脚本名为：', sys.argv[0])
    for i in range(1, len(sys.argv)):
        print('参数 %s 为：%s' % (i, sys.argv[i]))
    print("lnlamp: error: lnlamp_internal_calculate_errors.py 接收的参数存在错误!!!")
    exit()

  float_data_file = sys.argv[1]
  float_data = ReadValues(str(float_data_file))
  double_data_file = sys.argv[2]
  double_data = ReadValues(str(double_data_file))
  try:
    res = ComputeDifference(float_data, double_data)
  except Exception as inst:
    print("Problem were encountered: ", inst)
    
  # print(PrettyPrint(res))
  print(res['double_no_fl'])
  print(res['float_no_fl'])
  print(res['mean_pec'])
  print(res['mean_abs'])
  print(res['max_pec'])
  print(res['max_abs'])
