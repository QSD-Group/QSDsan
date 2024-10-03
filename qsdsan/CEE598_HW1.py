# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 21:03:25 2024

@author: Junhyung Park
"""

import pandas as pd

# 파일 경로 설정 (파일 경로는 자신의 경로로 설정해야 함)
faostat_file_path = 'FAOSTAT1.xlsx'
vwc_file_path = 'vwc.xlsx'

# FAOSTAT 데이터 불러오기 (CSV 파일)
faostat_data = pd.read_csv(faostat_file_path)

# VWC 데이터 불러오기 (Excel 파일)
vwc_data = pd.read_excel(vwc_file_path)

# 필요한 열만 사용하여 VWC 데이터 정리 (필요한 경우 열 이름 변경)
vwc_data_cleaned = vwc_data.rename(columns={'Product description (HS)': 'Item', 'vwc': 'VWC'})

# FAOSTAT와 VWC 데이터 매칭
merged_data = pd.merge(faostat_data, vwc_data_cleaned[['Item', 'VWC']], on='Item', how='inner')

# 가상 물 무역(VWT) 계산 (수출/수입량 * VWC)
merged_data['VWT'] = merged_data['Value'] * merged_data['VWC']

# 최종 데이터 출력 (또는 파일로 저장)
print(merged_data.head())

# 최종 데이터 CSV 파일로 저장하기 (원하는 경우)
merged_data.to_csv('output_vwt_data.csv', index=False)